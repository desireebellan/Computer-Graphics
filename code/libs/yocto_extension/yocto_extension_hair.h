//
// # Yocto/Extension: Tiny Yocto/GL extension
//
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
//

#ifndef _YOCTO_EXTENSION_HAIR_H_
#define _YOCTO_EXTENSION_HAIR_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <yocto/yocto_math.h>
#include <yocto/yocto_shape.h>

// -----------------------------------------------------------------------------
// ALIASES
// -----------------------------------------------------------------------------

namespace yocto::extension{

    // math
    using math::luminance;
    using math::make_rng;
    using math::identity3x3f;
    using math::mat3f;
    using math::vec2f;
    using math::vec3f;
    using math::vec4f;
    using math::vec2i;
    using math::vec3i;
    using math::vec4i;
    using math::zero3f;
    using math::pi;
    using math::pif;
    using math::ray3f;
    using math::rand1f;
    using math::rng_state;
    using math::sample_hemisphere;

    //shape
    using shape::bvh_tree;
    using shape::bvh_node;

}

// -----------------------------------------------------------------------------
// HIGH LEVEL API - HAIR SHADER
// -----------------------------------------------------------------------------

namespace yocto::extension{

    // HairBSDF Constants
    inline const float SqrtPiOver8 = 0.626657069f;
    inline const int pMax = 3;

    // Hair BSDF definition
    class HAIR_BSDF{
        public:
            // constructor
            HAIR_BSDF(float h, float eta, vec3f sigma_a, float beta_m, float beta_n, float alpha, vec3f normal, vec3f tangent, 
                vec3f color_ = zero3f, float eumelainin = 0, float pheomelanin = 0);

            // method

            static math::vec3f sigma_a_from_concentration(float ce, float cp);
            static math::vec3f sigma_a_from_reflectance(const math::vec3f &c, float beta_n);
            math::vec3f f(math::vec3f& wo, math::vec3f& wi, const bool convert = true);
            math::vec3f sample_f(math::vec3f& wo, math::vec3f *wi, const math::vec2f &u, float *pdf, const bool convert = true);
            float pdf(math::vec3f& wo, math::vec3f& wi, const bool convert = true); 
            std::array<float, pMax + 1> compute_ap_pdf(float cos_theta_o);

        private:
            // Parameters
            float h = 0;
            float gammaO = 0;
            vec3f sigma_a = zero3f;

            math::vec3f normal, tangent;
            std::array<float, pMax + 1> v;
            float s;
            float sin2kAlpha[3], cos2kAlpha[3];

            struct hair_material {
                vec3f sigma_a = zero3f;
                vec3f color = zero3f;
                float eumelanin = 0;
                float pheomelanin = 0;
                float eta = 1.55;
                float beta_m = 0.3;
                float beta_n = 0.3;
                float alpha = 2;
            };

            hair_material material;
    };

    // Hair Local Functions

    inline float Phi(int p, float gammaO, float gammaT);
    inline float LogI0(float x);
    inline float I0(float x);
    inline float LogisticCDF(float x, float s);
    inline float TrimmedLogistic(float x, float s, float a, float b);

    // YOCTO-HAIR conversion functions

    void yocto_to_hair(vec3f *ray, vec3f z, vec3f x);
    void hair_to_yocto(vec3f *ray, vec3f z, vec3f x);
}

// -----------------------------------------------------------------------------
// HIGH LEVEL API - TESTS
// -----------------------------------------------------------------------------

namespace yocto::extension{
    void white_furnace();
    void white_furnace_sampled();
    void sampling_weights();
    void sampling_consistency();
}

#endif


