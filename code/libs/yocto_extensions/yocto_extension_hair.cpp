//
// Implementation for Yocto/Extension.
//

//
// LICENSE:
//
// Copyright (c) 2020 -- 2020 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#include "yocto_extension_hair.h"

#include <iostream>
#include <chrono>

namespace yocto::extension{
    using math::abs;
    using math::asin;
    using math::clamp;
    using math::exp;
    using math::fresnel_dielectric;
    using math::log;
    using math::max;
    using math::min;
    using math::normalize;
    using math::sqrt;
    using math::zero3f;
}

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR EXTENSION
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// HAIR SHADER
// -----------------------------------------------------------------------------
namespace yocto::extension{

    // BSDF inline functions
    inline float AbsCosTheta(const vec3f& w) { return abs(w.z); }
    inline bool isNan(const vec3f& w) { return std::isnan(w.x) || std::isnan(w.y) || std::isnan(w.z);}

    // -------------------------
    // GENERAL UTILITY FUNCTIONS
    // -------------------------

        static uint32_t Compact1By1(uint32_t x) {
        // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
        x &= 0x55555555;
        // x = --fe --dc --ba --98 --76 --54 --32 --10
        x = (x ^ (x >> 1)) & 0x33333333;
        // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
        x = (x ^ (x >> 2)) & 0x0f0f0f0f;
        // x = ---- ---- fedc ba98 ---- ---- 7654 3210
        x = (x ^ (x >> 4)) & 0x00ff00ff;
        // x = ---- ---- ---- ---- fedc ba98 7654 3210
        x = (x ^ (x >> 8)) & 0x0000ffff;
        return x;
    }

    static vec2f Demuxfloat(float f) {
        uint64_t v = f * (1ull << 32);
        uint32_t bits[2] = {Compact1By1(v), Compact1By1(v >> 1)};
        return {bits[0] / float(1 << 16), bits[1] / float(1 << 16)};
    }

    inline float Sqr(float v) { return v * v; }
    inline vec3f Sqr(vec3f v) { return v * v; }

    template <int n>
    static float Pow(float v) {
        static_assert(n > 0, "Power can’t be negative");
        float n2 = Pow<n / 2>(v);
        return n2 * n2 * Pow<n & 1>(v);
    }
    template <> float Pow<1>(float v) { return v; }
    template <> float Pow<0>(float v) { return 1; }

    inline float SafeASin(float x) {
        return asin(clamp(x, -1.0f, 1.0f));
    }

    inline float SafeSqrt(float x) {
        return sqrt(max(float(0), x));
    }

    // -----------------------
    // HAIR LOCAL DECLARATIONS
    // -----------------------

    // compute the values of the modified Bessel function of the first kind
    inline float I0(float x) {
        float val = 0;
        float x2i = 1;
        int64_t ifact = 1;
        int i4 = 1;
        for (int i = 0; i < 10; ++i) {
            if (i > 1) ifact *= i;
            val += x2i / (i4 * Sqr(ifact));
            x2i *= x * x;
            i4 *= 4;
        }
        return val;
    }

    inline float LogI0(float x) {
        if (x > 12)
            return x + 0.5 * (-log(2 * pi) + log(1 / x) + 1 / (8 * x));
        else
            return log(I0(x));
    }


    // ---------------------
    // HAIR LOCAL FUNCTIONS
    // ---------------------
    
    // Longitudinal Scattering Function
    static float Mp(float cosThetaI, float cosThetaO, float sinThetaI, float sinThetaO, float v) {
        float a = cosThetaI * cosThetaO / v;
        float b = sinThetaI * sinThetaO / v;
        float mp = (v <= .1) ?
            (exp(LogI0(a) - b - 1/v + 0.6931f + log(1 / (2*v)))) :
            (exp(-b) * I0(a)) / (std::sinh(1 / v) * 2 * v);
        if(std::isinf(mp) || std::isnan(mp)) std::cerr<<"CHECK FAILED : MP"<<std::endl;
        return mp;  
    }

    // Attenuation Function
    static std::array<vec3f, pMax + 1> Ap(float cos_theta_o, float eta, float h, const vec3f &T) {
        std::array<vec3f, pMax + 1> ap;
        //Compute p = 0 attenuation at initial cylinder intersection
        auto cosGamma_o = SafeSqrt(1 - h * h);
        auto cosTheta = cos_theta_o * cosGamma_o;
        auto f = fresnel_dielectric(eta, {0, 0, 1}, {0, 0, cosTheta});
        if (f == 1.f)  {
            std::cout<<"cosTheta "<<cosTheta<<std::endl;
            std::cout<<"cos_theta_o "<<cos_theta_o<<std::endl;
            std::cout<<"cosGamma_o "<<cosGamma_o<<std::endl;
            std::cout<<"h "<<h<<std::endl;}
        ap[0] = vec3f(f);
        //Compute p = 1 attenuation term
        ap[1] = Sqr(1 - f) * T;
        //Compute attenuation terms up to p = pMax
        for (int p = 2; p < pMax; ++p)
            ap[p] = ap[p - 1] * T * f;
        //Compute attenuation term accounting for remaining orders of scattering
        ap[pMax] = ap[pMax - 1] * f * T / (vec3f(1.f) - T * f);

        return ap;
    }

    // Azimuthal Scattering Function
    inline float Np(float phi, int p, float s, float gammaO, float gammaT) {
        float dphi = phi - Phi(p, gammaO, gammaT);
        //Remap dphi to [−π, π]
        while (dphi > pi) dphi -= 2 * pi;
        while (dphi < -pi) dphi += 2 * pi;
        return TrimmedLogistic(dphi, s, -pi, pi);
    }
    // PHI : gives the net change in azimuthal direction
    inline float Phi(int p, float gammaO, float gammaT) {
        return 2 * p * gammaT - 2 * gammaO + p * pi;
    }
    //  way to represent surface roughness
    inline float Logistic(float x, float s) {
        x = abs(x);
        return exp(-x / s) / (s * Sqr(1 + exp(-x / s)));
    }

    inline float LogisticCDF(float x, float s) { return 1 / (1 + exp(-x / s)); }

    inline float TrimmedLogistic(float x, float s, float a, float b) { 
        if(a >= b) std::cerr<<"CHECK FAILED : TRIMMED LOGISTIC"<<std::endl;
        return Logistic(x, s) / (LogisticCDF(b, s) - LogisticCDF(a, s)); }

    static float SampleTrimmedLogistic(float u, float s, float a, float b) {
        if(a >= b) std::cerr<<"CHECK FAILED :  SAMPLE TRIMMED LOGISTIC"<<std::endl;
        float k = LogisticCDF(b, s) - LogisticCDF(a, s);
        float x = -s * std::log(1 / (u * k + LogisticCDF(a, s)) - 1);
        if(std::isnan(x)) std::cerr<<"CHECK FAILED : SAMPLE TRIMMED LOGISTIC"<<std::endl;
        return clamp(x, a, b);
    }

    // ---------------------
    // HAIR BSDF CONSTRUCTOR
    // ---------------------

    HAIR_BSDF::HAIR_BSDF(float h_, float eta_, vec3f sigma_a_, float beta_m_, float beta_n_, float alpha_, 
        vec3f normal_, vec3f tangent_, vec3f color_, float eumelanin_, float pheomelanin_){

        material.beta_m = beta_m_;
        material.beta_n = beta_n_;
        material.alpha = alpha_;
        material.sigma_a = sigma_a_;
        material.color = color_;
        material.eumelanin = eumelanin_;
        material.pheomelanin = pheomelanin_;
        material.eta = eta_;

        #ifdef YOCTO_EMBREE
            h = h_;
        #else
            h = -1 + 2. * h_;
        #endif

        if(h < -1 || h > 1) std::cerr<<"CHECK FAILED : H "<<std::endl;
        if(beta_m_ < 0 || beta_m_ > 1) std::cerr<<"CHECK FAILED : BETA_M "<<std::endl;
        if(beta_n_ < 0 || beta_n_ > 1) std::cerr<<"CHECK FAILED : BETA_N "<<std::endl;

        gammaO = SafeASin(h);
        
        normal = normal_;
        tangent = tangent_;

        if (material.sigma_a != zero3f){
            sigma_a = material.sigma_a;
        }
        else if(material.color != zero3f){
            sigma_a = sigma_a_from_reflectance(material.color, material.beta_n);
        }
        else if(material.eumelanin != 0 ||material.pheomelanin != 0){
            sigma_a = sigma_a_from_concentration(material.eumelanin, material.pheomelanin);
        }

        // Compute longitudinal variance from βm
        v[0] = Sqr(0.726f * material.beta_m + 0.812f * Sqr(material.beta_m) + 3.7f * Pow<20>(material.beta_m));
        v[1] = .25 * v[0]; 
        v[2] = 4 * v[0];
        for (int p = 3; p <= pMax; ++p)
            v[p] = v[2];

        //Compute azimuthal logistic scale factor from βn
        s = SqrtPiOver8 * (0.265f * material.beta_n + 1.194f * Sqr(material.beta_n) + 5.372f * Pow<22>(material.beta_n));
        if(std::isnan(s)) std::cerr<<"CHECK FAILED : S "<<std::endl;

        //Compute α terms for hair scales
        sin2kAlpha[0] = sin(material.alpha);
        cos2kAlpha[0] = SafeSqrt(1 - Sqr(sin2kAlpha[0]));
        for (int i = 1; i < 3; ++i) {
            sin2kAlpha[i] = 2 * cos2kAlpha[i - 1] * sin2kAlpha[i - 1];
            cos2kAlpha[i] = Sqr(cos2kAlpha[i - 1]) - Sqr(sin2kAlpha[i - 1]);
            }
    }

    // ---------------------
    // HAIR BSDF METHODS
    // ---------------------

    // GLOBAL

    vec3f HAIR_BSDF::f(vec3f& wo, vec3f& wi, const bool convert) {

        // get variables from bsdf
        if(convert){
            yocto_to_hair(&wo, normal, tangent);
            yocto_to_hair(&wi, normal, tangent);
        }

        //Compute hair coordinate system terms related to wo
        float sin_theta_o = wo.x;
        float cos_theta_o = SafeSqrt(1 - Sqr(sin_theta_o));
        float phiO = atan2(wo.z, wo.y);

        //Compute hair coordinate system terms related to wi
        float sin_theta_i = wi.x;
        float cos_theta_i = SafeSqrt(1-Sqr(sin_theta_i));
        float phiI = atan2(wi.z, wi.y);

        //Compute cos θt for refracted ray
        float sinThetaT = sin_theta_o / material.eta;
        float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

        //Compute γt for refracted ray
        float etap = sqrt(material.eta * material.eta - Sqr(sin_theta_o)) / cos_theta_o;
        float sinGammaT = h / etap;
        float cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));
        float gammaT = SafeASin(sinGammaT);

        //Compute the transmittance T of a single path through the cylinder 
        vec3f T = exp(-sigma_a * (2 * cosGammaT / cosThetaT));

        //Evaluate hair BSDF
        float phi = phiI - phiO;
        auto ap = Ap(cos_theta_o, material.eta, h, T);
        vec3f fsum(0.);
        for (int p = 0; p < pMax; ++p) {
            // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
            float sin_theta_op, cos_theta_op;
            if (p == 0) {
                sin_theta_op = sin_theta_o * cos2kAlpha[1] - cos_theta_o * sin2kAlpha[1];
                cos_theta_op = cos_theta_o * cos2kAlpha[1] + sin_theta_o * sin2kAlpha[1];
            }

            // Handle remainder of $p$ values for hair scale tilt
            else if (p == 1) {
                sin_theta_op = sin_theta_o * cos2kAlpha[0] + cos_theta_o * sin2kAlpha[0];
                cos_theta_op = cos_theta_o * cos2kAlpha[0] - sin_theta_o * sin2kAlpha[0];
            } else if (p == 2) {
                sin_theta_op = sin_theta_o * cos2kAlpha[2] + cos_theta_o * sin2kAlpha[2];
                cos_theta_op = cos_theta_o * cos2kAlpha[2] - sin_theta_o * sin2kAlpha[2];
            } else {
                sin_theta_op = sin_theta_o;
                cos_theta_op = cos_theta_o;
            }

            // Handle out-of-range $\cos \thetao$ from scale adjustment
            cos_theta_op = abs(cos_theta_op);
            fsum += Mp(cos_theta_i, cos_theta_op, sin_theta_i, sin_theta_op, v[p]) * ap[p] *
                    Np(phi, p, s, gammaO, gammaT);
        }

        // Compute contribution of remaining terms after _pMax_
        fsum += Mp(cos_theta_i, cos_theta_o, sin_theta_i, sin_theta_o, v[pMax]) * ap[pMax] / (2.f * pi);

        if(std::isinf(luminance(fsum)) || std::isnan(luminance(fsum))) std::cerr<<"CHECK FAILED : HAIR::F"<<std::endl;

        return fsum;
    }

    vec3f HAIR_BSDF::sample_f(vec3f &WO, vec3f *WI, const vec2f &u2, float *pdf_, const bool convert) {

        
        if(convert) yocto_to_hair(&WO, normal, tangent);

        if (*WI == zero3f){

            // Compute hair coordinate system terms related to wo
            float sin_theta_o = WO.x;
            float cos_theta_o = SafeSqrt(1 - Sqr(sin_theta_o));
            float phiO = atan2(WO.z, WO.y);

            // Derive four random samples from u2
            vec2f u[2] = { Demuxfloat(u2[0]), Demuxfloat(u2[1]) };

            // Determine which term p to sample for hair scattering
            std::array<float, pMax + 1> apPdf = compute_ap_pdf(cos_theta_o);
            int p;
            for (p = 0; p < pMax; ++p) {
                if (u[0][0] < apPdf[p]) break;
                u[0][0] -= apPdf[p];
            }
                // Rotate $\sin \thetao$ and $\cos \thetao$ to account for hair scale tilt
            float sin_theta_op, cos_theta_op;
            if (p == 0) {
                    sin_theta_op = sin_theta_o * cos2kAlpha[1] - cos_theta_o * sin2kAlpha[1];
                    cos_theta_op = cos_theta_o * cos2kAlpha[1] + sin_theta_o * sin2kAlpha[1];
            }
            else if (p == 1) {
                    sin_theta_op = sin_theta_o * cos2kAlpha[0] + cos_theta_o * sin2kAlpha[0];
                    cos_theta_op = cos_theta_o * cos2kAlpha[0] - sin_theta_o * sin2kAlpha[0];
            } else if (p == 2) {
                    sin_theta_op = sin_theta_o * cos2kAlpha[2] + cos_theta_o * sin2kAlpha[2];
                    cos_theta_op = cos_theta_o * cos2kAlpha[2] - sin_theta_o * sin2kAlpha[2];
            } else {
                    sin_theta_op = sin_theta_o;
                    cos_theta_op = cos_theta_o;
            }

            // Sample Mp to compute θi
            u[1][0] = max(u[1][0], float(1e-5));
            float cosTheta = 1 + v[p] * log(u[1][0] + (1 - u[1][0]) * exp(-2 / v[p]));
            float sinTheta = SafeSqrt(1 - Sqr(cosTheta));
            float cosPhi = cos(2 * pi * u[1][1]);
            float sin_theta_i = -cosTheta * sin_theta_op + sinTheta * cosPhi * cos_theta_op;
            float cos_theta_i = SafeSqrt(1 - Sqr(sin_theta_i));
            // Sample Np to compute ∆φ
                // Compute γt for refracted ray
                float etap = sqrt(material.eta * material.eta - Sqr(sin_theta_o)) / cos_theta_o;
                float sinGammaT = h / etap;
                float gammaT = SafeASin(sinGammaT);
                float dphi;
                if (p < pMax) dphi = Phi(p, gammaO, gammaT) + SampleTrimmedLogistic(u[0][1], s, -pi, pi);
                else dphi = 2 * pi * u[0][1];
            // Compute wi from sampled hair scattering angles
            float phiI = phiO + dphi;
            // auto wi = vec3f(sin_theta_i, cos_theta_i * cos(phiI), cos_theta_i * sin(phiI));
            *WI = vec3f(sin_theta_i, cos_theta_i * cos(phiI), cos_theta_i * sin(phiI));}

        else yocto_to_hair(WI, normal, tangent);

        *pdf_ = pdf(WO, *WI, false);
        auto f_ = f(WO, *WI, false);

        if(!convert && AbsCosTheta(*WI) > 0) f_ /= AbsCosTheta(*WI);

        if(convert) {
            hair_to_yocto(WI, normal, tangent);
            hair_to_yocto(&WO, normal, tangent);
            }

        return f_;
    }

    float HAIR_BSDF::pdf(vec3f &wo, vec3f &wi, const bool convert) {

        if(convert) {
            yocto_to_hair(&wo, normal, tangent);
            yocto_to_hair(&wi, normal, tangent);
        }
        
        // Compute hair coordinate system terms related to _wo_
        float sin_theta_o = wo.x;
        float cos_theta_o = SafeSqrt(1 - Sqr(sin_theta_o));
        float phiO = atan2(wo.z, wo.y);

        // Compute hair coordinate system terms related to _wi_
        float sin_theta_i = wi.x;
        float cos_theta_i = SafeSqrt(1 - Sqr(sin_theta_i));
        float phiI = std::atan2(wi.z, wi.y);

        // Compute $\gammat$ for refracted ray
        float etap = sqrt(material.eta * material.eta - Sqr(sin_theta_o)) / cos_theta_o;
        float sinGammaT = h / etap;
        float gammaT = SafeASin(sinGammaT);

        // Compute PDF for $A_p$ terms
        std::array<float, pMax + 1> apPdf = compute_ap_pdf(cos_theta_o);

        // Compute PDF sum for hair scattering events
        float phi = phiI - phiO;
        float pdf = 0;
        for (int p = 0; p < pMax; ++p) {
            // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
            float sin_theta_op, cos_theta_op;
            if (p == 0) {
                sin_theta_op = sin_theta_o * cos2kAlpha[1] - cos_theta_o * sin2kAlpha[1];
                cos_theta_op = cos_theta_o * cos2kAlpha[1] + sin_theta_o * sin2kAlpha[1];
            }

            // Handle remainder of $p$ values for hair scale tilt
            else if (p == 1) {
                sin_theta_op = sin_theta_o * cos2kAlpha[0] + cos_theta_o * sin2kAlpha[0];
                cos_theta_op = cos_theta_o * cos2kAlpha[0] - sin_theta_o * sin2kAlpha[0];
            } else if (p == 2) {
                sin_theta_op = sin_theta_o * cos2kAlpha[2] + cos_theta_o * sin2kAlpha[2];
                cos_theta_op = cos_theta_o * cos2kAlpha[2] - sin_theta_o * sin2kAlpha[2];
            } else {
                sin_theta_op = sin_theta_o;
                cos_theta_op = cos_theta_o;
            }

            // Handle out-of-range $\cos \thetao$ from scale adjustment
            cos_theta_op = std::abs(cos_theta_op);
            pdf += Mp(cos_theta_i, cos_theta_op, sin_theta_i, sin_theta_op, v[p]) *
                apPdf[p] * Np(phi, p, s, gammaO, gammaT);
        }
        pdf += Mp(cos_theta_i, cos_theta_o, sin_theta_i, sin_theta_o, v[pMax]) *
            apPdf[pMax] * (1 / (2 * pi));
        return pdf;
    }

    // LOCAL

    std::array<float, pMax + 1> HAIR_BSDF::compute_ap_pdf(float cos_theta_o){

        // Compute array of Ap values for cos_theta_o
            float sin_theta_o = SafeSqrt(1 - cos_theta_o * cos_theta_o);
            //Compute cos θt for refracted ray
            float sinThetaT = sin_theta_o / material.eta;
            float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));
            //Compute γt for refracted ray
            float etap = std::sqrt(material.eta * material.eta - Sqr(sin_theta_o)) / cos_theta_o;
            float sinGammaT = h / etap;
            float cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));
            float gammaT = SafeASin(sinGammaT);
            // Compute the transmittance T of a single path through the cylinder 
            vec3f T = exp(-sigma_a * (2 * cosGammaT / cosThetaT));
            auto ap = Ap(cos_theta_o, material.eta, h, T);
        // Compute Ap PDF from individual Ap terms
            std::array<float, pMax + 1> apPdf;
            float sumY = 0.0f; 
            //for (int i = 0; i <= pMax; ++i) sumY += ap[i].y;
            for (int i = 0; i <= pMax; ++i) sumY += luminance(ap[i]);
            //for (int i = 0; i <= pMax; ++i) apPdf[i] = ap[i].y / sumY;
            for (int i = 0; i <= pMax; ++i) apPdf[i] = luminance(ap[i]) / sumY;
            return apPdf;
    }

    vec3f HAIR_BSDF::sigma_a_from_concentration(float ce, float cp) {
        auto eumelaninSigmaA = vec3f{0.419f, 0.697f, 1.37f};
        auto pheomelaninSigmaA = vec3f{0.187f, 0.4f, 1.05f};
        return ce * eumelaninSigmaA + cp * pheomelaninSigmaA;
    }

    vec3f HAIR_BSDF::sigma_a_from_reflectance(const vec3f& c, float beta_n) {
        vec3f sigma_a;
        sigma_a = Sqr(log(c) /
        (5.969f - 0.215f * beta_n + 2.532f * Sqr(beta_n) -
        10.73f * Pow<3>(beta_n) + 5.574f * Pow<4>(beta_n) +
        0.245f * Pow<5>(beta_n)));
        return sigma_a;
    }
    

    // -------------------------
    // AUXILIARY FUNCTIONS
    // -------------------------

    // yocto brdf uses a different coordinate system from the one used in the hair shader
    // the following functions convert the two coordinate systems to make them work together

    void yocto_to_hair(vec3f *ray, vec3f z, vec3f x){
        auto frame = inverse(frame_fromzx(zero3f, z, x));
        *ray = transform_direction(frame, *ray);
    }

    void hair_to_yocto(vec3f *ray, vec3f z, vec3f x){
        auto frame = frame_fromzx(zero3f, z, x);
        *ray = transform_direction(frame, *ray);
    }
}

// -----------------------------------------------------------------------------
// TESTS 
// -----------------------------------------------------------------------------

namespace yocto::extension{

    // -----------------------------------------------------------------------------
    // HAIR SHADER
    // -----------------------------------------------------------------------------
    // Almost the same as math::sample_sphere(), only using the double form of \pi, used only to ensure 
    // the success of the test functions
    vec3f uniform_sample_sphere(const vec2f &u) {
    float z = 1 - 2 * u[0];
    float r = sqrt(max((float)0, (float)1 - z * z));
    float phi = 2 * pi * u[1];
    return vec3f(r * cos(phi), r * sin(phi), z);
    }

    void white_furnace(){
        rng_state rng;
        auto wo = uniform_sample_sphere(rand2f(rng));
        for (float beta_m = .1; beta_m < 1; beta_m += .2) {
            for (float beta_n = .1; beta_n < 1; beta_n += .2) {

                //Estimate reflected uniform incident radiance from hair 
                vec3f sum = zero3f;
                int count = 300000;
                std::vector<float> test;
                for (int i = 0; i < count; ++i) {

                    #ifdef YOCTO_EMBREE
                            auto h = -1 + 2 * rand1f(rng);
                    #else
                            auto h = rand1f(rng);
                    #endif

                    if (h == 0) h += math::flt_eps;

                    auto hair = HAIR_BSDF(h, 1.55, zero3f, beta_m, beta_n, 0.f, {0, 0, 1}, {1, 0, 0});
                    auto wi = uniform_sample_sphere(rand2f(rng));

                    if(AbsCosTheta(wi) <= 0) sum += hair.f(wo, wi, false)  * AbsCosTheta(wi);
                    else sum += hair.f(wo, wi, false);
                    }

                    float avg = luminance(sum) / (count * sample_sphere_pdf(wo));

                    if(!(avg >= .95 && avg <= 1.05)) throw std::runtime_error("TEST WHITE FURNACE FAILED!");
                }
            }
        std::cout << "TEST WHITE FURNACE SUCCEDED" << std::endl;
    }

    void white_furnace_sampled(){
        rng_state rng;
        vec3f wo = uniform_sample_sphere(math::rand2f(rng));
        for (float beta_m = .1; beta_m < 1; beta_m += .2) {
            for (float beta_n = .1; beta_n < 1; beta_n += .2) {
                //Estimate reflected uniform incident radiance from hair 
                vec3f sum = zero3f;
                int count = 300000;
                for (int i = 0; i < count; ++i) {
                    
                    #ifdef YOCTO_EMBREE
                            auto h = -1 + 2 * rand1f(rng);
                    #else
                            auto h = rand1f(rng);
                    #endif

                    if (h == 0) h += math::flt_eps;

                    auto hair = HAIR_BSDF(h, 1.55, zero3f, beta_m, beta_n, 0.f, {0, 0, 1}, {1, 0, 0});

                    vec3f wi = zero3f;
                    float pdf;
                    vec3f f = hair.sample_f(wo,&wi, math::rand2f(rng), &pdf, false);
                    if (pdf > 0) sum += f * AbsCosTheta(wi) / pdf;
                    }

                    float avg = luminance(sum) / count;
                    if (!(avg >= 0.99f && avg <= 1.01f)) throw std::runtime_error("TEST WHITE FURNACE SAMPLED FAILED!");
                }
            }
        std::cout<<"TEST WHITE FURNACE SAMPLED SUCCEDED\n";
    }

    void sampling_weights() {
        rng_state rng;
        for (float beta_m = .1; beta_m < 1; beta_m += .2) {
            for (float beta_n = .4; beta_n < 1; beta_n += .2) {

                //Check sample_f() sample weight
                int count = 300000;
                for (int i = 0; i < count; ++i) {

                    #ifdef YOCTO_EMBREE
                            auto h = -1 + 2 * rand1f(rng);
                    #else
                            auto h = rand1f(rng);
                    #endif

                    if (h == 0) h += math::flt_eps;

                    auto hair = HAIR_BSDF(h, 1.55, zero3f, beta_m, beta_n, 0.f, {1.f, 0.f, 0.f}, {0.f, 0.f, 1.f});

                    vec3f wo = uniform_sample_sphere(math::rand2f(rng));
                    vec2f u = math::rand2f(rng);

                    vec3f wi = zero3f;
                    float pdf;
                    vec3f f  = hair.sample_f(wo, &wi, u, &pdf, false);
                    if (pdf > 0) {
                        if (luminance(f) * AbsCosTheta(wi) / pdf < 0.999 || luminance(f) * AbsCosTheta(wi) / pdf> 1.001)
                            throw std::runtime_error("TEST SAMPLING WEIGHTS FAILED!");
                        }
                    }
                }
            }
        std::cout << "TEST SAMPLING WEIGHTS SUCCEDED" << std::endl;
    }

    void sampling_consistency() {
        rng_state rng;
        for (float beta_m = .2; beta_m < 1; beta_m += .2) {
            for (float beta_n = .4; beta_n < 1; beta_n += .2) {
                    // Declare variables for hair sampling test
                    const int count = 64*1024;
                    vec3f sigma_a = vec3f(.25f);
                    vec3f wo = uniform_sample_sphere(math::rand2f(rng));
                    auto Li = [](const vec3f &w) -> vec3f {return vec3f(w.z * w.z);};
                    vec3f fImportance = zero3f, fUniform = zero3f;
                    for (int i = 0; i < count; ++i) {
                        // Compute estimates of scattered radiance for hair sampling test
                        #ifdef YOCTO_EMBREE
                                auto h = -1 + 2 * rand1f(rng);
                        #else
                                auto h = rand1f(rng);
                        #endif

                        auto hair = HAIR_BSDF(h, 1.55, zero3f, beta_m, beta_n, 0.f, {0, 0, 1}, {1, 0, 0});
                        auto u = math::rand2f(rng);
                        vec3f wi = zero3f;
                        float pdf;
                        auto f = hair.sample_f(wo, &wi, u, &pdf, false);
                        if (pdf > 0) fImportance += f * Li(wi) * AbsCosTheta(wi) / (count * pdf);
                        wi = uniform_sample_sphere(u);
                        if (AbsCosTheta(wi) <= 0) fUniform += hair.f(wo, wi, false) * Li(wi) * AbsCosTheta(wi) / (count * sample_sphere_pdf(wo));
                        else fUniform += hair.f(wo, wi, false) * Li(wi) / (count * sample_sphere_pdf(wo));
                        
                    }
                // Verify consistency of estimated hair reflected radiance values
                float err = std::abs(luminance(fImportance)- luminance(fUniform)) / luminance(fUniform);
                if (err >= 0.05) throw std::runtime_error("TEST SAMPLING CONSISTENCY FAILED!");
                }
        }
        std::cout << "TEST SAMPLING CONSISTENCY SUCCEDED" << std::endl;
    }   
}
