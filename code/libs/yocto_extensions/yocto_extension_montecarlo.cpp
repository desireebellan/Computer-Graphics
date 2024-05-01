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

#include "yocto_extension_montecarlo.h"


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
// MONTECARLO GEOMETRY PROCESSING 
// -----------------------------------------------------------------------------

namespace yocto::extension{

    static inline void print(vec3f x, std::string message){std::cout<<message<<x.x<<" "<<x.y<<" "<<x.z<<std::endl;}
     static inline void print(vec2f x, std::string message){std::cout<<message<<x.x<<" "<<x.y<<" "<<std::endl;}

    // Performs a linear interpolation between x and y using a weight h between them
    inline float mix(float x, float y, float h){return x*(1-h)+y*h;}
    inline vec3f mix(vec3f x, vec3f y, float h){return x*(1-h)+y*h;}
    // Computes fractional of a float value
    float fractional(float x){return x - floor(x);}
    vec3f fractional(vec3f x){return x - vec3f{(float)floor(x.x), (float)floor(x.y), (float)floor(x.z)};}

    //Plane
    float sdPlane(vec3f point, vec3f normal, float r ){return dot(point,normal) - r;}
    // Sphere
    float sdSphere(vec3f point, vec3f center, float radius) { return length(point - center) - radius;}
    // Torus of major radius R and minor radius r; centered in C and spun about the z-axis
    float sdTorus(vec3f point, vec3f C, float R, float r){
        auto xy = vec2f(point.x, point.y);
        auto Cxy = vec2f(C.x, C.y);
        return length(vec2f(length(xy - Cxy) - R, point.z - C.z)) - r;
    }

    // IMPLICIT SURFACES OPERATIONS
    // Union
    float opUnion( float d1, float d2 ) { return min(d1,d2); }
    //Subtraction
    float opSubtraction( float d1, float d2 ) {return max(-d1,d2); }
    // Intersection
    float opIntersection( float d1, float d2 ) {return max(d1,d2); }
    // Smooth Union
    float opSmoothUnion(float d1, float d2, float k) {
        float h = clamp( 0.5f + 0.5f*(d2-d1)/k, 0.0f, 1.0f);
        return mix(d2, d1, h) - k*h*(1.0-h); }
    // Smooth Subtraction
    float opSmoothSubtraction( float d1, float d2, float k ) {
        float h = clamp( 0.5f - 0.5f*(d2+d1)/k, 0.0f, 1.0f );
        return mix( d2, -d1, h ) + k*h*(1.0-h); }
    // Smooth Intersection
    float opSmoothIntersection( float d1, float d2, float k ) {
        float h = clamp( 0.5f - 0.5f*(d2-d1)/k, 0.0f, 1.0f );
        return mix( d2, d1, h ) + k*h*(1.0-h); }

    template <int n>
    static float Pow(float v) {
        static_assert(n > 0, "Power can’t be negative");
        float n2 = Pow<n / 2>(v);
        return n2 * n2 * Pow<n & 1>(v);
    }
    template <> float Pow<1>(float v) { return v; }
    template <> float Pow<0>(float v) { return 1; }

    // ---------------
    // WALK ON CIRCLES
    // ---------------
    // Estimate PDE in a 2D point using Montecarlo Walk On Circles for a single walk
    // 
    template <typename Func>
    vec3f wos_pde(std::function<float(vec2f)> sdf, vec2f x0, Func g, Func f, Func h, wos_params* params, rng_state rng){
        
        vec3f sum = vec3f(0.f); 
        float R;
        auto step = 0;
        vec3f f_;
        vec2f x_prec = x0;
        vec2f x = x0;
        vec2f x1;
        // russian roulette
        if(rand1f(rng) < params->q) return vec3f(params->b);
        if (sdf(x0) > params->tollerance){
            do{
                    // solve exterior problem using russian roulette
                    R = sdf(x_prec);
                    x = x_prec + R * sample_disk(rand2f(rng));
                    x_prec = x;
                    step ++;}
                while(R > params->tollerance && step < params->max_steps);}
            // x = closest point on surface?
            sum += g(x);
            // roussian soulette, account for skipped samples
            sum = (sum - params->q*params->b)/(1-params->q);
            return sum;
    }

    // ----------------
    // WALK ON SPHERE
    // ----------------

    // Estimate PDE in a 3D point using Montecarlo Walk On Spheres for a single walk
    // 
    template <typename SDF = std::function<float(vec3f)>, typename Func>
    float wos_pde(SDF sdf, vec3f x0, Func g, Func f, Func h, wos_params* params, rng_state rng){
        
        float sum = 0.f; 
        float R;
        auto step = 0;
        float f_;
        vec3f x_prec = x0;
        vec3f x = x0;
        vec3f x1;

        // russian roulette
        if(rand1f(rng) < params->q) return params->b;
            do{
                    // solve exterior problem using russian roulette
                    R = -sdf(x_prec);                          
                    if (R <= params->tollerance) break; 
                    x = x_prec + R * sample_sphere(rand2f(rng));
                    if (step == 0) x1 = x;
                    // Poisson and Biharmonic Solution
                    // Sample a point y inside the B(x) ball
                    auto y = importance_sample(f_, f, R, x_prec, *params, rng);
                    auto volume = 4/3 * pi * Pow<3>(R);
                    // compute the laplace estimator over the point y using the boundary function h
                    // This Nested approach is slow and increase variance
                    // f_ = wos_pde(sdf, y, h, [](vec3f point){return std;;nanf;}, [](vec3f point){return 0.f;}, params, rng);
                    // A more efficient approach is tree walking
                    if(std::isnan(f_)) f_ = h(x); 
                    if (params->c <= 0) sum += volume * f_ * Green_harmonic(x_prec, y, R);               
                    else {
                        auto screening = R * sqrt(params->c) / sinh(R * sqrt(params->c));
                        sum = sum * screening + volume * f_ * Yukawa_potential(x_prec, y, R, screening);}
                    x_prec = x;
                    step ++;}
                while(R > params->tollerance && step < params->max_steps);
            // x = closest point on surface?
            sum += g(x);
            // roussian soulette, account for skipped samples
            sum = (sum - params->q*params->b)/(1-params->q);
            if (params->variance)
                return sum - dot(wos_gradient(sdf, x0, g, f, h, *params, rng), x1 - x0);
            return sum;
    }

    // version of ptr::trace_samples() to test montecarlo for implicit surfaces
    // 3D solver
    template <typename SDF_intersect = std::function<float(ray3f)>, typename SDF = std::function<float(vec3f)>, typename Func>
    void solve_wos_samples(ptr::state* state, SDF_intersect sdf_intersect, SDF sdf, Func g, Func f, Func h, 
        const ptr::camera* camera, wos_params* params){

        if (params->noparallel) {
            for (auto j = 0; j < state->render.size().y; j++) {
                for (auto i = 0; i < state->render.size().x; i++) {
                    state->render[{i, j}] = solve_wos_sample(state, sdf_intersect, sdf, g, f, h, camera, {i, j}, params);
                }
            }
        } else {
            parallel_for(
                state->render.size(), [state, camera, params, sdf_intersect,  sdf, g, f, h](const vec2i& ij) {
                    state->render[ij] = solve_wos_sample(state, sdf_intersect, sdf, g, f, h, camera, ij, params);
                });
        }
    }

    // 2D solver
    template <typename Func>
    void solve_wos_samples(ptr::state* state, std::function<float(vec2f)> bound_sdf, Func g, Func f, Func h, wos_params* params){
        if (params->noparallel) {
            for (auto j = 0; j < state->render.size().y; j++) {
                for (auto i = 0; i < state->render.size().x; i++) {               
                    state->render[{i, j}] = solve_wos_sample(state, bound_sdf, g, f, h, {i, j}, params);
                }
            }
        } else {
            parallel_for(
                state->render.size(), [state, params, bound_sdf, g, f, h](const vec2i& ij) {
                    state->render[ij] = solve_wos_sample(state, bound_sdf, g, f, h, ij, params);
                });
        }
    }
    // 3D solver
    template <typename SDF_intersect = std::function<float(ray3f)>, typename SDF = std::function<float(vec3f)>, typename Func>
    static vec4f solve_wos_sample(ptr::state* state,  SDF_intersect sdf_intersect, SDF sdf, Func g, Func f, Func h, const ptr::camera* camera, 
        const vec2i& ij, wos_params* params){
        auto& pixel  = state->pixels[ij];
        ray3f ray = sample_camera(camera, ij, state->pixels.size(), rand2f(pixel.rng), rand2f(pixel.rng)); 
        vec4f solution;
        vec3f res;
        bool hit = false;
        auto wos_solver = [h, f, g, sdf, params, pixel](vec3f x0)->vec3f{
            auto res = wos_pde(sdf, x0, g, f, h, params, pixel.rng);
            // auto normal = sdNormal(sdf, x0, params->h);
            //return colormap(res, *params) * normal;
            // return cross(colormap(res, *params),normal);
            // return normalize(colormap(res, *params) + normal);
            return colormap(res, *params);
        };
        float dist = sdf_intersect(ray);
        if(dist < math::flt_max){
            auto x0 = ray.o + ray.d*dist;
            res =  wos_solver(x0);
            /*if (sdPlane(x0, vec3f(0, 0, 1), 0.03) < 0.f && sdPlane(x0, vec3f(0, 0, 1), 0.03) > -0.0001f)
                res = {0.f, 0.f, 0.f};
            else res = {0.94,0.3,0.07};*/
            // res = colormap(g(x0), *params);
            hit = true;
            }
        else{
            // outside the surface
            // implementare environment lights
            res = vec3f{1.0f, 1.0f, 1.0f};
            }   
        solution = {res, hit? 1.f : 0.f};
        if (!isfinite(xyz(solution))) xyz(solution) = zero3f;
        if (max(xyz(solution)) > params->clamp)
            xyz(solution) = xyz(solution) * (params->clamp / max(xyz(solution)));
        pixel.accumulated += solution;
        pixel.samples += 1;
        return pixel.accumulated / pixel.samples;
    }


    // 2D solver
    template <typename Func>
    static vec4f solve_wos_sample(ptr::state* state, std::function<float(vec2f)> bound_sdf, Func g, Func f, Func h, const vec2i& ij, wos_params* params){
        auto& pixel  = state->pixels[ij];
        auto x0 = vec2f((float)(2.f*ij.x/state->render.size().x - 1.f), (float)(2.f*(state->render.size().y - ij.y)/state->render.size().y - 1.f));
        float dist;
        bool hit = false;
        auto wos_solver = [h, f, g, bound_sdf, params, pixel](vec2f x0)->vec3f{
            return wos_pde(bound_sdf, x0, g, f, h, params, pixel.rng);
        };
        vec4f solution =  {wos_solver(x0), 1.f};
        if (!isfinite(xyz(solution))) xyz(solution) = zero3f;
        if (max(xyz(solution)) > params->clamp)
            xyz(solution) = xyz(solution) * (params->clamp / max(xyz(solution)));
        pixel.accumulated += solution;
        pixel.samples += 1;
        return pixel.accumulated / pixel.samples;
    }

    // Estimate gradient of PDE in a point using Montecarlo Walk On Sphere
    template <typename SDF, typename Func>
    static vec3f wos_gradient(SDF sdf, vec3f x0, Func g, Func f, Func h, wos_params params, rng_state rng){
        auto x = x0;
        auto sum = vec3f(0.f); 
        int step = 0;           
        // first point
        auto R = -sdf(x0);
        auto x1 = x0 + R * sample_sphere(rand2f(rng));
        vec3f normal = (x1 - x0)/R;
        float coefficient = 3/R;
        x0 = x1;
        float f_;

        do{
            // radius of the smallest ball around the point inside the volume  
            R = -sdf(x0);
            // Sample a point y inside the B(x) ball
            auto y = importance_sample(f_, f, R, x0, params, rng);
            x1 = x0 + R * sample_sphere(rand2f(rng));
            auto volume = 4/3 * pi * Pow<3>(R);
            // compute the laplace estimator over the point y using the boundary function h
            // This Nested approach is slow and increase variance
            // f_ = wos_pde(sdf, y, h, [](vec3f point){return std::nanf;}, [](vec3f point){return 0.f;}, params, rng);
            // A more efficient approach is tree walking
            if(std::isnan(f_)) f_ = h(x1);
            auto screening = R * sqrt(params.c) / sinh(R * sqrt(params.c));
            if (params.c <= 0) sum += volume * f_ * Green_harmonic_gradient(x0, y, R);
            else sum = sum * params.c + volume * f_ * Yukawa_potential_gradient(x0, y, R, screening);   
            x0 = x1;        
            step ++;
        } 
        while( R > params.tollerance && step < params.max_steps);
        // x = closest point on the surface to x0 ?
        return  coefficient*g(x0)*normal + sum ;
    }

    // Estimate the hessian using the simple wos algorithm
    template <typename Func>
    void wos_hessian(FWN_tree tree_node, vec3f x0, Func g, Func f, float eps, rng_state rng){
        auto R = tree_node.closest_triangle(x0);
        auto x1 = x0 + R * sample_sphere(rand2f(rng));
        auto y0 = R * rand1f(rng) * sample_sphere(rand2f(rng));
        // TEST
        float c=1000, q=0., b=0.;
        auto uf = wos_poisson(tree_node, x1, g, f, c, q, b, eps, rng);
        auto volume = 4/3 * pi * Pow<3>(R);
        mat3f x1x0 = {x1-x0, zero3f, zero3f};
        return uf * (transpose(x1x0)*9/Pow<4>(R)*x1x0-identity3x3f*3/Pow<2>(R)) + volume*Green_harmonic_hessian(f, x0, y0, R);
    }

    // -------------------
    // POTENTIAL FUNCTIONS
    // -------------------
    // HARMONIC GREEN FUNCTION
    float Green_harmonic(vec3f x, vec3f y, float R){
        auto r = math::distance(x,y);
        return clamp(1.f/(4*pif) * (R - r)/(r*R), math::flt_min, math::flt_max);
    }
    // DERIVATIVE GREEN FUNCTIONS
    vec3f Green_harmonic_gradient(vec3f x, vec3f y, float R){
        auto r = math::distance(x,y);
        return (y-x)/(4*pif)*(1/pow(r, 3) - 1/pow(R, 3));
    }

    // GREEN IMPORTANCE SAMPLING
    vec3f sample_green(rng_state rng){
        // sample uniformly a direction 
        auto d = sample_sphere(rand2f(rng));
        //polar angle
        auto theta = atan2(d.y, d.x);
        // Use Ulrich polar method to get r = |y-x|
        float S, U, V;
        do{
            U = rand1f(rng);
            V = 2*rand1f(rng) - 1;
            S = pow(U, 2) + pow(V, 2);
        }
        while(S > 1);
        float sinTheta = sin(theta);

        auto r = sinTheta/(4*pif)*(0.5 + U*V/S * sqrt(1.f-pow(S, 2.f/3.f)));
        return d * r ;
    }

    // YUKAWA POTENTIAL FOR SCREENED POISSON
    float Yukawa_potential(vec3f x, vec3f y, float R, float C){
        auto r = math::distance(x,y);
        return 1.f/(4*pi) * (sinh((R-r))/r) * C/R;
    }  
    // DERIVATIVE YUKAWA POTENTIAL
    vec3f Yukawa_potential_gradient(vec3f x, vec3f y, float R, float c){
        auto r = math::distance(x,y);
        return (y-x)/(4*pi)*(sqrt(c)*cosh((R - r)*sqrt(c))/(r*sinh(R*sqrt(c)))*(1/r - 1/R)
                            + sinh((R - r)*sqrt(c))/(r*sinh(R*sqrt(c)))*(1/pow(r,2) + cosh((R - r)*sqrt(c))/(r*sinh(R*sqrt(c)))));
    }


    // ----------------
    // GLOBAL FUNCTIONS 
    // ----------------

    // Sample random rotation matrix using PVR 
    mat3f sample_rotation(vec3f normal, rng_state rng){
        auto v = sample_hemisphere(normal, rand2f(rng));
        float phi = rand1f(rng) * 2 * pi;
        float cphi = cos(phi);
        float sphi = sin(phi);
        float sigma = 1 - cphi;
        vec3f x = {(float)pow(v[0], 2) * sigma + cphi, v[0] * v[1] * sigma + v[2] * sphi, v[0] * v[2] * sigma - v[1] * sphi};
        vec3f y = {v[1] * v[0] * sigma - v[2] * sphi, (float)pow(v[1], 2) * sigma + cphi, v[1] * v[2] * sigma + v[0] * sphi};
        vec3f z = {v[2] * v[0] * sigma + v[1] * sphi, v[2] * v[1] *  sigma - v[0] * sphi, (float)pow(v[2], 2) * sigma + cphi};
        return mat3f(x, y, z);
    }

    // Get importance sampling 3D
    template <typename Func>
    vec3f importance_sample(float &f_, Func f, float R, vec3f x, wos_params params, rng_state rng){
        vec3f y;
        auto importance = params.importance;
        if (!importance.sample_y) y = x + R * rand1f(rng) * sample_sphere(rand2f(rng));
        else if (params.c <= 0) y = x + R * sample_green(rng);
        else {
            std::cout<<"Screened Green sampling is not implemented.\n";
            y = x + R * rand1f(rng) * sample_sphere(rand2f(rng));}     
        if(importance.sample_f == "point"s) {
            if (length(x - params.z) < R){
                // f_ = params.cz;
                y = params.z;
            }
        }
        else if (importance.sample_f == "curve"s){
            std::cout<<"Curve source sampling is not implemented.\n";
        } 
        f_ = f(y);
        return y;
    }


    // Extract point from  set that are placed inside a volume
    // Uses General Winding Number
    void points_in_volume(std::vector<vec3f>& sample_points, std::vector<vec3f> points, std::vector<vec3f> positions, std::vector<vec3i> triangles){           
        for (auto& point: points)
            if (in_volume(point, triangles, positions)) 
                    sample_points.push_back(point);
    }

    // Check if a point p is inside a volume
    // Generalized Winding Number, correct for clean meshes
    bool in_volume(vec3f point, std::vector<vec3i> triangles, std::vector<vec3f> points){
        float w = 0.f;
        for (auto& triangle: triangles)
            w += 1.f/(4*pif)*solid_angle(points[triangle.x], points[triangle.y], points[triangle.z], point);
        if (w > 0.5) return true;
        else return false;
    }

    // Compute solid angle projection of a triangle over a point
    // Used to check if a point p is inside a volume

    float solid_angle(vec3f vi, vec3f vj, vec3f vk, vec3f p){
        //distances
        auto a = vi - p;
        auto b = vj - p;
        auto c = vk - p;
        //determinant
        auto det = dot(a,cross(b,c));
        auto norm_a = math::length(a);
        auto norm_b = math::length(b);
        auto norm_c = math::length(c);
        auto denom = norm_a*norm_b*norm_c + dot(a,b)*norm_c + dot(b,c)*norm_a + dot(a,c)*norm_b;
        return atan2(det,denom)*2;
    }

    // sample n points given a vector of points
    std::vector<vec3f> sample_box(std::vector<vec3f> vertices, vec3f& center, float& radius, int num, math::rng_state rng){
        std::vector<vec3f> min_max = get_min_max(vertices);
        center = (min_max[0] + min_max[1])/2;
        radius = math::distance(center, min_max[1])/2;
        // slice = (min_max[1].x - min_max[0].x)*slice + min_max[0].x;
        std::vector<vec3f> samples = {};
        for (int n=0; n<num; n++){
            auto rand_point = math::rand3f(rng)*(min_max[1]-min_max[0])+min_max[0];
            samples.push_back(rand_point);
        }
        return samples;
    }

    // Returns the minimum and maximum vertices of a mesh
    std::vector<vec3f> get_min_max(std::vector<vec3f> vertices){
        auto point_max = vec3f(-math::flt_max);
        auto point_min = vec3f(math::flt_max);
        for (auto point: vertices){
            point_max = max(point_max, point);
            point_min = min(point_min, point);
        }
        std::vector<vec3f> min_max = {point_min, point_max};
        return min_max;
    }

    // Check if a ray intersect an implicit surface
    template <typename SDF>
    float sphereTrace(ray3f ray, SDF sdf, float epsilon){
        auto t = ray.tmin;
        while(t < ray.tmax) {
        auto p = eval_ray(ray, t);
        auto d = abs(sdf(p));
        if(d < epsilon) {return t; }
        t += d;
        }
        return math::flt_max;
    }

    // estimate the normal to a SDF in a point given an increment h
    template <typename SDF>
    vec3f sdNormal(SDF sdf, vec3f x, float h){
        return normalize(vec3f(sdf(x + vec3f{h, 0, 0}) - sdf(x -vec3f{h, 0, 0}),
                sdf(x + vec3f{0, h, 0}) - sdf(x -vec3f{0, h, 0}),
                sdf(x + vec3f{0, 0, h}) - sdf(x -vec3f{0, 0, h})));}

    // SHAPES
    // IMPLICIT SURFACES
    // Plane given the normal of the plane and distance from the origin h
    // Diffusion curve
    void sdLine(vec2f point, vec3f& col,  float& d, vec2f a, vec2f b, vec3f cu, vec3f cv)
    {
        auto pa = point - a;
        auto ba = b - a;
        float h = clamp( dot(pa,ba)/dot(ba,ba),0.0,1.0);
        float distance = length(pa-h*ba);
        if(distance < d ){
            d = distance;
            if (distance < 0.01)
            {float s = pa.x*ba.y-pa.y*ba.x;
            col = s<0.0 ? cu : cv;}}
    }
    void sdLine(vec2f point, vec3f& color, float& d, vec2f a, vec2f b, vec3f cu0,vec3f cv0,vec3f cu1, vec3f cv1)
        {
            auto pa = point - a;
            auto ba = b - a;
            float h = clamp( dot(pa,ba)/dot(ba,ba),0.0,1.0);
            float distance = length(pa-h*ba);
            if(distance < d)
            {   if (distance < 0.01){
                float s = pa.x*ba.y-pa.y*ba.x;
                h = 1.0-h;
                color = s<0.0 ? mix(cu0,cu1,h) : mix(cv0,cv1,h);}
            }
        }


    void create_point_cloud(std::vector<vec4i>& quads, std::vector<vec3f>& positions, std::vector<vec3i>& triangles, std::vector<vec2i>& lines){
        rng_state rng;
        if(!lines.empty()){
            for (auto line : lines){
                auto x = positions[line[0]];
                auto y = positions[line[1]];
                auto barycenter = 1/2 * (x + y);
                positions.push_back(barycenter);
            }
        }
        if(!quads.empty()){
            for (auto quad : quads){
                auto x = positions[quad[0]];
                auto y = positions[quad[1]];
                auto z = positions[quad[2]];
                auto q = positions[quad[3]];
                auto bc0 = 1/3 * (x + y + z);
                auto bc1 = 1/3 * (y + z + q);
                auto bc2 = 1/3 * (z + q + x);
                auto bc3 = 1/3 * (q + x + y);
                auto v1 = 1/3 * (bc0 - bc2);
                auto v2 = 1/3 * (bc1 - bc3);
                auto w = 1/3 * (bc1 - bc0);
                auto s = (w.y*v1.x - w.x*v1.y) / (v2.x*v1.y - v2.y*v1.x);
                auto barycenter = bc1 + s * v2;
                positions.push_back(barycenter);  
            }       
        }
        if(!triangles.empty()){
            for (auto triangle : triangles){
                auto x = positions[triangle[0]];
                auto y = positions[triangle[1]];
                auto z = positions[triangle[2]];
                auto barycenter = 1/3 * (x + y + z);
                positions.push_back(barycenter);                
            }

        }
    }

    void create_polygon_soup(std::vector<vec4i>& quads, std::vector<vec3f>& positions, std::vector<vec3i>& triangles, std::vector<vec2i>& lines){
        // create vector with new points
        auto positions_o = positions;
        auto triangles_o = triangles;
        auto quads_o = quads;
        positions.clear();
        triangles.clear();
        quads.clear();
        // randomly rotate elements around their barycenter
        rng_state rng;
        if(!lines.empty()){
            // NOT IMPLEMENTED
        }
        if(!quads_o.empty()){
            for (auto quad : quads_o){
                auto x = positions_o[quad[0]];
                auto y = positions_o[quad[1]];
                auto z = positions_o[quad[2]];
                auto q = positions_o[quad[3]];
                auto bc0 = 1/3 * (x + y + z);
                auto bc1 = 1/3 * (y + z + q);
                auto bc2 = 1/3 * (z + q + x);
                auto bc3 = 1/3 * (q + x + y);
                auto v1 = 1/3 * (bc0 - bc2);
                auto v2 = 1/3 * (bc1 - bc3);
                auto w = 1/3 * (bc1 - bc0);
                auto s = (w.y*v1.x - w.x*v1.y) / (v2.x*v1.y - v2.y*v1.x);
                auto barycenter = bc1 + s * v2;
                vec3f normal = math::quad_normal(x, y, z, q);
                vec4i quad_new;
                auto R = sample_rotation(normal, rng);
                for (auto i=0; i<4; i++){
                    positions.push_back((barycenter - positions_o[quad[i]]) * R + barycenter);
                    quad_new[i] = positions.size() - 1;
                }
                quads.push_back(quad_new);

            }
        }
        if(!triangles_o.empty()){
            for(auto triangle : triangles_o){
                auto x = positions_o[triangle[0]];
                auto y = positions_o[triangle[1]];
                auto z = positions_o[triangle[2]];
                auto barycenter = 1/3 * (x + y + z);
                auto normal = math::triangle_normal(x, y, z);
                auto R = sample_rotation(normal, rng);
                vec3i triangle_new;
                for (auto i=0; i<3; i++){
                    positions.push_back((barycenter - positions_o[triangle[i]]) * R + barycenter);
                    triangle_new[i] = positions.size() - 1;
                }
                triangles.push_back(triangle_new);
            }
            
        }
    }

    // shape deformation (biharmonic)
    // surface reconstructions (poisson)
    // helmotz decomposition (poisson 2d e 3d)
    // symmetric direction fields (laplace 2d)
    // rendering

    template <typename Func>
    inline void parallel_for(const int size, Func&& func) {
        auto             futures  = std::vector<std::future<void>>{};
        auto             nthreads = std::thread::hardware_concurrency();
        std::atomic<int> next_idx(0);
        auto progress = vec2i{0, (int)size};
        for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
            futures.emplace_back(
                std::async(std::launch::async, [&func, &next_idx, &progress, size]() {
                    while(true){
                        auto i = next_idx.fetch_add(1);
                        if (i >= size) break;
                        func(i);
                    }
                }));
            }
            for (auto& f : futures) f.get();
    }
    template <typename Func>
    inline void parallel_for(const vec2i ij, Func&& func) {
        auto             futures  = std::vector<std::future<void>>{};
        auto             nthreads = std::thread::hardware_concurrency();
        std::atomic<int> next_idx(0);
        for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
            futures.emplace_back(
                std::async(std::launch::async, [&func, &next_idx, ij]() {
                while (true) {
                    auto i = next_idx.fetch_add(1);
                    if (i >= ij.x) break;
                    for (auto j = 0; j < ij.y; j++) {
                        func(vec2i{i, j});}
                }
                }));
        }
        for (auto& f : futures) f.get();
    }
}

// -----------------------------------------------------------------------------
// FAST WINDING NUMBERS
// -----------------------------------------------------------------------------

namespace yocto::extension{

    // Class Methods

    FWN_tree::FWN_tree(std::vector<vec3f> positions_, std::vector<vec3i> triangles_){
        
        triangles = triangles_;
        positions = positions_;

        make_triangles_bvh(bvh, triangles, positions, std::vector<float> (positions.size(), 1e-5f));

        nodes.resize(bvh.nodes.size());
        auto curr_id = 0;

        auto Ct = [](vec3f x, vec3f y, vec3f z, vec3f p) -> mat3f {
            auto xij = 1.f/2.f*(x + y) - p;
            auto xjk = 1.f/2.f*(y + z) - p;
            auto xik = 1.f/2.f*(x + z) - p;
            return (outer(xij, xij)* 1.f/3.f + outer(xjk, xjk) * 1.f/3.f + outer(xik, xik) * 1.f/3.f);
        };

        std::cout<<"Building taylor expansion for the bvh...\n";

        for (auto node : bvh.nodes){      
              
            nodes[curr_id].p = math::center(node.bbox); 
            nodes[curr_id].id = curr_id; 
            nodes[curr_id].r = length(bvh.nodes[curr_id].bbox.max - bvh.nodes[curr_id].bbox.min)/2.f;
            if(node.internal){
                nodes[curr_id].children.resize(node.num, nullptr);
                float distance;
                for(auto child = 0; child < node.num; child ++){
                    nodes[curr_id].children[child] = &nodes[node.start + child];
                    nodes[node.start + child].parent = &nodes[curr_id];
                    //distance = math::distance(nodes[curr_id].p, math::center(bvh.nodes[node.start + child].bbox)); 
                    //if (distance > nodes[curr_id].r) nodes[curr_id].r = distance;

                }     
            }
            else{
                // Compute the taylor terms for the node
                // For each leaf
                vec3f t1 = zero3f;
                mat3f t2 = {zero3f, zero3f, zero3f};
                auto t3 = mat9x3f();
                for (auto i=0; i < node.num; i++){
                    int idx = bvh.primitives[node.start + i];
                    auto x = positions[triangles[idx].x];
                    auto y = positions[triangles[idx].y];
                    auto z = positions[triangles[idx].z];                  
                    auto bbox = math::triangle_bounds(x, y, z);
                    auto r = math::distance(nodes[curr_id].p, math::center(bbox));
                    if (r > nodes[curr_id].r) nodes[curr_id].r = r;                           
                    auto at = math::triangle_area(x, y, z);
                    auto nt = math::triangle_normal(x, y, z);
                    // First order term of Taylor series
                    t1 += at * nt;
                    t2 += outer((1/3*(x + y + z) - nodes[curr_id].p), nt) * at;
                    auto xij = 1.f/2.f*(x + y) - nodes[curr_id].p;
                    auto xjk = 1.f/2.f*(y + z) - nodes[curr_id].p;
                    auto xik = 1.f/2.f*((x + z) - nodes[curr_id].p);
                    // auto Ct = (outer(xij, xij) + outer(xjk, xjk) + outer(xik, xik))*1.f/3.f;
                    // Third order term of Taylor expansion
                    t3 = t3 + at * outer(Ct(x, y, z, nodes[curr_id].p), nt);
                } 

                auto parent = &nodes[curr_id];
                while(parent != nullptr){
                    auto child = parent;
                    child -> t1 += t1;
                    child -> t2 += t2;
                    child -> t3 += t3;
                    parent = child -> parent;
                } 
            }
            curr_id ++;
        }

        std::cout<<"Taylor expansion finished.\n";
    
    }

    float FWN_tree::intersect(ray3f &ray){
        vec3i triangle;
        vec2f uv;
        shape::bvh_intersection intersection = shape::intersect_triangles_bvh(bvh, triangles, positions, ray, true);
        if (intersection.hit) return intersection.distance;
        return math::flt_max;
    }

    // Fast Winding number implementation for clean meshes, triangle soups and point clouds
    // Check if a point is inside a volume 
    // Implementation inspired by "Fast Winding Numbers for Soups and Clouds" by G. Barill et al.
    void FWN_tree::fast_points_in_volume(std::vector<vec3f>& sample_points, std::vector<vec3f> points, int approximation, float beta){
        
        std::vector<float> w(points.size());      
        std::cout<<"Checking points in volume ...\n";       
        auto in_volume = [&w, this, points, approximation, beta](int i){
            w[i] = fast_point_in_volume(points[i], beta, approximation);
        };

        for (auto i=0; i<points.size(); i++) in_volume(i);
        for (auto i = 0; i<w.size(); i++) {
            if (w[i] > 0.5) sample_points.push_back(points[i]);}
    }

    float FWN_tree::fast_point_in_volume(vec3f point, float beta, int approximation){
        float w = 0.;
        auto node = &nodes[0];
        auto checked = std::vector<int> (nodes.size(), 0);

        auto identity_product = [](vec3f r) -> mat9x3f{
            auto sum = mat9x3f();
            for (auto i = 0; i<3; i++) 
            sum += (outer(r, outer(identity3x3f[i], identity3x3f[i])) 
                    + outer(identity3x3f[i], outer(r, identity3x3f[i])) 
                    + outer(identity3x3f[i], outer(identity3x3f[i], r)));
            return sum;
        };

        while(true){
            if(bvh.nodes[node->id].internal){
                auto condition = math::distance(point, node->p) > beta*node->r;
                if(checked[node->id] == 0 && condition){                 
                            auto r = node->p - point;
                            auto length = math::length(r);
                            auto r3 = pow(length,3);
                            // first order derivative of Green's function
                            auto dG = r/r3*1.f/(4*pif);
                            w += math::dot(node -> t1,dG);
                            if (approximation > 1){
                                auto r5 = pow(length, 5);
                                auto crossr = outer(r, r);
                                // second order derivative of Green's function
                                auto d2G = (identity3x3f/r3 - crossr*2.f/r5)*1.f/(4*pif);
                                w += inner(node -> t2,d2G);
                                if (approximation == 3){
                                // third order derivative of Green's function
                                    auto d3G = -identity_product(r)/r5 + 15*outer(r, crossr)/pow(length, 7)*1.f/(4*pif);
                                    w += 1./2. * inner(node -> t3, d3G);
                                }
                            }
                             
                            if (node->parent != nullptr) 
                                node = node->parent;
                            else break;

                            continue;
                }
                if (checked[node->id] < bvh.nodes[node->id].num){
                    node = node->children[checked[node->id]];
                    if (node->parent != nullptr) checked[node->parent->id] += 1;
                    continue;
                }
                if (node->parent != nullptr) node = node->parent;
                else break;
            }
            else{
                for (auto i=0; i<bvh.nodes[node->id].num; i++) {
                    auto id = bvh.primitives[bvh.nodes[node->id].start + i];
                    //solid angle
                    w += solid_angle(positions[triangles[id].x], positions[triangles[id].y], positions[triangles[id].z], point);
                    // area-weighted dipole
                    // w +=  math::dot(math::triangle_area(positions[triangles[id].x], positions[triangles[id].y], positions[triangles[id].z]) * 
                    //    math::triangle_normal(positions[triangles[id].x], positions[triangles[id].y], positions[triangles[id].z]),
                    //    (fnode.p - point)/(4*pif*pow(math::length(fnode.p - point),3)));
                } 
                checked[node->parent->id] += 1;
                node = node->parent;
            }
        
        
        }
        return w;
    }

    float FWN_tree::closest_triangle(vec3f point){

        auto node = &nodes[0];
        float R, sdf;
        vec3f p;
        auto triangle_distance = [point](const vec3f p0, const vec3f p1, const vec3f p2, vec3f& pt) -> float{
            auto cuv = closestuv_triangle(point, p0, p1, p2);
            pt = interpolate_triangle(p0, p1, p2, cuv);
            return length(point - pt);
            /*auto normal = triangle_normal(p0, p1, p2);
            auto p   = p0 * (1 - cuv.x - cuv.y) + p1 * cuv.x + p2 * cuv.y;
            uv = cuv;
            //return dot(p-point, normal);
            auto dd  = dot(p - point, p - point);
            return sqrt(dd);*/
        };
        while(true){
            R = math::flt_max;
            if(bvh.nodes[node->id].internal){
                fwn_node *new_node; 
                for(auto child = 0; child < bvh.nodes[node->id].num; child ++){  
                    auto distance = math::distance(point, node->children[child]->p);
                    if(std::isinf(distance)){return R;}
                    if (distance < R) {
                        R = distance;
                        new_node = node->children[child];
                    }
                }
                node = new_node;
            }
            else{
                vec3f pt;
                for(auto child = 0; child < bvh.nodes[node->id].num; child ++){
                    auto idx = bvh.primitives[bvh.nodes[node->id].start + child];
                    auto x = positions[triangles[idx].x];
                    auto y = positions[triangles[idx].y];
                    auto z = positions[triangles[idx].z]; 
                    auto distance = triangle_distance(x, y, z, pt);
                    if (distance < R){
                        R = distance;
                        sdf = dot(triangle_normal(x, y, z), point - pt) >= 0 ? R : -R;
                    }
                }
                break;
            }
        }
        return - abs(sdf);
    }

    vec3f colormap(float x, wos_params params){
        switch(params.colormap_id){
            case 1:{
                // yellow-purple stripes
                auto color2 = vec3f{0.94,0.74,0.09};
                // auto color1 = vec3f{0.03, 0.02, 0.2};
                auto color1 = vec3f{0.04,0.08,0.47};
                // auto color2 = vec3f{1.0, 1.0, 0.28};
                return mix(color1, color2, x);
            }

           case 2:{
            // red-blue
            auto color1 = vec3f{1.0, 0.16, 0.17};
            auto color2 = vec3f{0.0, 0.3, 0.8};
            return mix(color1, color2, x);               
            }

            case 3:{
            // yellow-purple
            auto color1 = vec3f{1.0, 1.0, 0.28};
            auto color2 = vec3f{0.03, 0.02, 0.2};
            return mix(color1, color2, x);
            }}
        return {0, 0, 0};
    }

}

// -------------------
// IMAGES AND SCENES
// -------------------
namespace yocto::extension{
    float sdLine(vec2f point, vec3f& color){
        float d = math::flt_max;
        sdLine(point, color, d, vec2f( 0.032566, 0.636419),vec2f(-0.140558, 0.554919),vec3f(0.802,0.563,0.443),vec3f(0.865,0.792,0.690));
        sdLine(point, color, d, vec2f( 0.157917, 0.667506),vec2f( 0.032586, 0.636408),vec3f(0.824,0.575,0.453),vec3f(0.869,0.804,0.702));
        sdLine(point, color, d, vec2f( 0.292774, 0.683634),vec2f( 0.157637, 0.667307),vec3f(0.849,0.600,0.467),vec3f(0.871,0.804,0.702));
        sdLine(point, color, d, vec2f( 0.462459, 0.681116),vec2f( 0.293026, 0.683576),vec3f(0.890,0.647,0.488),vec3f(0.871,0.804,0.700));
        sdLine(point, color, d, vec2f( 0.532097, 0.663221),vec2f( 0.462281, 0.681257),vec3f(0.931,0.692,0.510),vec3f(0.871,0.804,0.698));
        sdLine(point, color, d, vec2f( 0.679623, 0.589789),vec2f( 0.532148, 0.663353),vec3f(0.951,0.712,0.518),vec3f(0.871,0.804,0.698));
        sdLine(point, color, d, vec2f( 0.737901, 0.539194),vec2f( 0.680024, 0.589419),vec3f(0.951,0.712,0.512),vec3f(0.871,0.804,0.698));
        sdLine(point, color, d, vec2f( 0.816252, 0.438381),vec2f( 0.737331, 0.539616),vec3f(0.933,0.700,0.502),vec3f(0.871,0.804,0.696));
        sdLine(point, color, d, vec2f( 0.855568, 0.362220),vec2f( 0.816666, 0.437574),vec3f(0.920,0.684,0.486),vec3f(0.871,0.802,0.694));
        sdLine(point, color, d, vec2f( 0.914063, 0.195570),vec2f( 0.855412, 0.362903),vec3f(0.904,0.667,0.467),vec3f(0.871,0.800,0.694));
        sdLine(point, color, d, vec2f( 0.906394,-0.050806),vec2f( 0.914184, 0.195399),vec3f(0.876,0.637,0.437),vec3f(0.871,0.800,0.692));
        sdLine(point, color, d, vec2f( 0.847684,-0.254043),vec2f( 0.906304,-0.051022),vec3f(0.841,0.600,0.402),vec3f(0.871,0.800,0.688));
        sdLine(point, color, d, vec2f( 0.753134,-0.419097),vec2f( 0.847548,-0.253827),vec3f(0.808,0.571,0.373),vec3f(0.871,0.798,0.686));
        sdLine(point, color, d, vec2f( 0.655326,-0.537575),vec2f( 0.752991,-0.419390),vec3f(0.796,0.557,0.357),vec3f(0.871,0.796,0.686),vec3f(0.698,0.439,0.271),vec3f(0.871,0.796,0.686));
        sdLine(point, color, d, vec2f( 0.589686,-0.601779),vec2f( 0.655427,-0.537635),vec3f(0.698,0.439,0.271),vec3f(0.871,0.796,0.686),vec3f(0.592,0.298,0.149),vec3f(0.871,0.796,0.686));
        sdLine(point, color, d, vec2f( 0.485769,-0.673529),vec2f( 0.589993,-0.601259),vec3f(0.553,0.269,0.131),vec3f(0.871,0.796,0.684));
        sdLine(point, color, d, vec2f( 0.420683,-0.695361),vec2f( 0.485736,-0.673833),vec3f(0.506,0.247,0.124),vec3f(0.876,0.806,0.692));
        sdLine(point, color, d, vec2f( 0.289074,-0.705976),vec2f( 0.420341,-0.695230),vec3f(0.482,0.275,0.155),vec3f(0.876,0.806,0.692));
        sdLine(point, color, d, vec2f( 0.261794,-0.694878),vec2f( 0.289218,-0.706199),vec3f(0.459,0.284,0.173),vec3f(0.871,0.796,0.680));
        sdLine(point, color, d, vec2f( 0.234336,-0.652438),vec2f( 0.261805,-0.694874),vec3f(0.443,0.249,0.153),vec3f(0.871,0.796,0.678));
        sdLine(point, color, d, vec2f( 0.229172,-0.613117),vec2f( 0.234298,-0.652477),vec3f(0.420,0.184,0.116),vec3f(0.871,0.796,0.678));
        sdLine(point, color, d, vec2f( 0.211059,-0.573956),vec2f( 0.229373,-0.612721),vec3f(0.404,0.145,0.094),vec3f(0.871,0.796,0.678),vec3f(0.384,0.090,0.059),vec3f(0.631,0.506,0.447));
        sdLine(point, color, d, vec2f( 0.794264, 0.332819),vec2f( 0.820637, 0.237829),vec3f(0.961,0.525,0.318),vec3f(0.941,0.394,0.175));
        sdLine(point, color, d, vec2f( 0.791799, 0.406427),vec2f( 0.794099, 0.333239),vec3f(0.967,0.565,0.365),vec3f(0.941,0.388,0.163));
        sdLine(point, color, d, vec2f( 0.773116, 0.456365),vec2f( 0.791677, 0.406208),vec3f(0.973,0.596,0.406),vec3f(0.941,0.382,0.153));
        sdLine(point, color, d, vec2f( 0.691214, 0.542967),vec2f( 0.772875, 0.456604),vec3f(0.980,0.631,0.455),vec3f(0.941,0.375,0.139));
        sdLine(point, color, d, vec2f( 0.584748, 0.612130),vec2f( 0.691856, 0.543066),vec3f(0.988,0.675,0.510),vec3f(0.941,0.367,0.124));
        sdLine(point, color, d, vec2f( 0.523270, 0.633009),vec2f( 0.583850, 0.612122),vec3f(0.994,0.702,0.543),vec3f(0.941,0.361,0.114));
        sdLine(point, color, d, vec2f( 0.416081, 0.649827),vec2f( 0.523622, 0.632948),vec3f(0.996,0.710,0.553),vec3f(0.941,0.357,0.110),vec3f(0.965,0.596,0.420),vec3f(0.941,0.353,0.094));
        sdLine(point, color, d, vec2f( 0.328754, 0.652317),vec2f( 0.416523, 0.649811),vec3f(0.955,0.573,0.390),vec3f(0.941,0.349,0.090));
        sdLine(point, color, d, vec2f( 0.152437, 0.628754),vec2f( 0.328509, 0.652168),vec3f(0.945,0.549,0.361),vec3f(0.941,0.345,0.086),vec3f(0.918,0.671,0.486),vec3f(0.929,0.404,0.192));
        sdLine(point, color, d, vec2f(-0.058755, 0.562422),vec2f( 0.152170, 0.628817),vec3f(0.918,0.671,0.486),vec3f(0.929,0.404,0.192),vec3f(0.847,0.643,0.471),vec3f(0.910,0.478,0.325));
        sdLine(point, color, d, vec2f(-0.152258, 0.413501),vec2f(-0.073990, 0.418056),vec3f(0.322,0.165,0.318),vec3f(0.839,0.357,0.294),vec3f(0.325,0.161,0.341),vec3f(0.875,0.427,0.400));
        sdLine(point, color, d, vec2f(-0.166656, 0.426903),vec2f(-0.152256, 0.413325),vec3f(0.325,0.161,0.341),vec3f(0.875,0.427,0.400),vec3f(0.329,0.157,0.361),vec3f(0.929,0.510,0.518));
        sdLine(point, color, d, vec2f(-0.150301, 0.459101),vec2f(-0.167019, 0.427044),vec3f(0.327,0.155,0.355),vec3f(0.902,0.504,0.522));
        sdLine(point, color, d, vec2f(-0.118651, 0.485569),vec2f(-0.150151, 0.459146),vec3f(0.322,0.151,0.329),vec3f(0.892,0.510,0.522));
        sdLine(point, color, d, vec2f(-0.002352, 0.524810),vec2f(-0.118560, 0.485419),vec3f(0.302,0.141,0.253),vec3f(0.918,0.525,0.484));
        sdLine(point, color, d, vec2f( 0.031055, 0.523623),vec2f(-0.002648, 0.525035),vec3f(0.329,0.161,0.216),vec3f(0.927,0.531,0.443));
        sdLine(point, color, d, vec2f( 0.698487, 0.039536),vec2f( 0.691611,-0.011391),vec3f(0.241,0.084,0.069),vec3f(0.894,0.271,0.043));
        sdLine(point, color, d, vec2f( 0.719461, 0.087611),vec2f( 0.698093, 0.039144),vec3f(0.208,0.076,0.112),vec3f(0.892,0.271,0.041));
        sdLine(point, color, d, vec2f( 0.756735, 0.118077),vec2f( 0.719793, 0.087918),vec3f(0.180,0.078,0.149),vec3f(0.890,0.271,0.037));
        sdLine(point, color, d, vec2f( 0.781480, 0.121070),vec2f( 0.756822, 0.117813),vec3f(0.208,0.100,0.151),vec3f(0.890,0.271,0.035));
        sdLine(point, color, d, vec2f( 0.801014, 0.097611),vec2f( 0.781005, 0.120993),vec3f(0.220,0.114,0.149),vec3f(0.888,0.269,0.033));
        sdLine(point, color, d, vec2f( 0.823724, 0.024839),vec2f( 0.801409, 0.097493),vec3f(0.210,0.129,0.149),vec3f(0.886,0.267,0.029));
        sdLine(point, color, d, vec2f( 0.823815,-0.104922),vec2f( 0.823541, 0.024940),vec3f(0.247,0.180,0.135),vec3f(0.882,0.265,0.024));
        sdLine(point, color, d, vec2f( 0.812986,-0.134971),vec2f( 0.823897,-0.104904),vec3f(0.275,0.206,0.120),vec3f(0.878,0.263,0.018));
        sdLine(point, color, d, vec2f( 0.781477,-0.159956),vec2f( 0.813047,-0.134843),vec3f(0.275,0.182,0.102),vec3f(0.876,0.263,0.014));
        sdLine(point, color, d, vec2f( 0.742899,-0.136473),vec2f( 0.781414,-0.160083),vec3f(0.273,0.161,0.084),vec3f(0.873,0.261,0.010));
        sdLine(point, color, d, vec2f( 0.707845,-0.099381),vec2f( 0.742671,-0.136282),vec3f(0.278,0.145,0.071),vec3f(0.871,0.259,0.006));
        sdLine(point, color, d, vec2f( 0.692791,-0.057184),vec2f( 0.707941,-0.099276),vec3f(0.278,0.127,0.057),vec3f(0.869,0.259,0.002));
        sdLine(point, color, d, vec2f( 0.691497,-0.011596),vec2f( 0.692977,-0.057240),vec3f(0.269,0.114,0.049),vec3f(0.859,0.257,0.000));
        sdLine(point, color, d, vec2f( 0.649558,-0.424193),vec2f( 0.570243,-0.515577),vec3f(0.394,0.220,0.129),vec3f(0.763,0.245,0.125));
        sdLine(point, color, d, vec2f( 0.691564,-0.410217),vec2f( 0.649794,-0.424293),vec3f(0.408,0.216,0.118),vec3f(0.800,0.247,0.122),vec3f(0.635,0.341,0.208),vec3f(0.773,0.318,0.176));
        sdLine(point, color, d, vec2f( 0.696625,-0.423080),vec2f( 0.691171,-0.410166),vec3f(0.635,0.341,0.208),vec3f(0.773,0.318,0.176),vec3f(0.569,0.345,0.220),vec3f(0.749,0.416,0.251));
        sdLine(point, color, d, vec2f( 0.689415,-0.442915),vec2f( 0.697144,-0.422830),vec3f(0.529,0.343,0.227),vec3f(0.729,0.441,0.273));
        sdLine(point, color, d, vec2f( 0.593770,-0.551001),vec2f( 0.689147,-0.443398),vec3f(0.441,0.322,0.225),vec3f(0.688,0.429,0.271));
        sdLine(point, color, d, vec2f( 0.555892,-0.565232),vec2f( 0.593883,-0.550881),vec3f(0.431,0.282,0.198),vec3f(0.665,0.351,0.220));
        sdLine(point, color, d, vec2f( 0.551358,-0.559085),vec2f( 0.555731,-0.565064),vec3f(0.471,0.263,0.180),vec3f(0.663,0.310,0.192),vec3f(0.345,0.263,0.192),vec3f(0.659,0.271,0.169));
        sdLine(point, color, d, vec2f( 0.570349,-0.515707),vec2f( 0.551420,-0.559235),vec3f(0.343,0.257,0.188),vec3f(0.661,0.257,0.157));
        sdLine(point, color, d, vec2f( 0.139656, 0.051930),vec2f( 0.125174,-0.011939),vec3f(0.212,0.176,0.188),vec3f(0.773,0.251,0.153),vec3f(0.322,0.298,0.314),vec3f(0.831,0.376,0.290));
        sdLine(point, color, d, vec2f( 0.156524, 0.072534),vec2f( 0.139337, 0.052148),vec3f(0.322,0.298,0.314),vec3f(0.831,0.376,0.290),vec3f(0.455,0.388,0.420),vec3f(0.863,0.510,0.447));
        sdLine(point, color, d, vec2f( 0.207624, 0.092959),vec2f( 0.156443, 0.072549),vec3f(0.455,0.388,0.420),vec3f(0.863,0.510,0.447),vec3f(0.608,0.612,0.635),vec3f(0.902,0.635,0.588));
        sdLine(point, color, d, vec2f( 0.277523, 0.089856),vec2f( 0.208031, 0.092954),vec3f(0.608,0.612,0.635),vec3f(0.902,0.635,0.588),vec3f(0.792,0.820,0.851),vec3f(0.945,0.776,0.753));
        sdLine(point, color, d, vec2f( 0.342574, 0.060277),vec2f( 0.277521, 0.089925),vec3f(0.792,0.820,0.851),vec3f(0.945,0.776,0.753),vec3f(0.851,0.867,0.894),vec3f(0.988,0.914,0.910));
        sdLine(point, color, d, vec2f( 0.395896, 0.011487),vec2f( 0.342325, 0.060273),vec3f(0.851,0.867,0.894),vec3f(0.988,0.914,0.910),vec3f(0.816,0.769,0.788),vec3f(0.973,0.792,0.769));
        sdLine(point, color, d, vec2f( 0.426593,-0.056929),vec2f( 0.395983, 0.010998),vec3f(0.816,0.769,0.788),vec3f(0.973,0.792,0.769),vec3f(0.784,0.690,0.702),vec3f(0.953,0.675,0.624));
        sdLine(point, color, d, vec2f( 0.429756,-0.097479),vec2f( 0.426646,-0.056343),vec3f(0.784,0.690,0.702),vec3f(0.955,0.676,0.627));
        sdLine(point, color, d, vec2f( 0.211977, 0.039860),vec2f( 0.152329, 0.066394),vec3f(0.251,0.255,0.322),vec3f(0.573,0.522,0.557),vec3f(0.239,0.251,0.325),vec3f(0.612,0.620,0.655));
        sdLine(point, color, d, vec2f( 0.253875,-0.009396),vec2f( 0.212258, 0.039555),vec3f(0.224,0.235,0.298),vec3f(0.627,0.651,0.694));
        sdLine(point, color, d, vec2f( 0.304442,-0.112422),vec2f( 0.254010,-0.009088),vec3f(0.190,0.196,0.251),vec3f(0.629,0.653,0.704));
        sdLine(point, color, d, vec2f( 0.333171,-0.133931),vec2f( 0.303949,-0.112429),vec3f(0.178,0.169,0.220),vec3f(0.612,0.620,0.671));
        sdLine(point, color, d, vec2f( 0.374078,-0.127448),vec2f( 0.333105,-0.133826),vec3f(0.222,0.190,0.220),vec3f(0.629,0.629,0.675));
        sdLine(point, color, d, vec2f( 0.429594,-0.097757),vec2f( 0.374129,-0.127644),vec3f(0.278,0.229,0.237),vec3f(0.678,0.659,0.694));
        sdLine(point, color, d, vec2f( 0.128499,-0.076818),vec2f( 0.113227,-0.128679),vec3f(0.131,0.071,0.118),vec3f(0.986,0.927,0.737));
        sdLine(point, color, d, vec2f( 0.214776,-0.207363),vec2f( 0.168035,-0.203218),vec3f(0.982,0.796,0.496),vec3f(0.124,0.039,0.092));
        sdLine(point, color, d, vec2f( 0.238606,-0.218984),vec2f( 0.215027,-0.207500),vec3f(0.990,0.782,0.447),vec3f(0.153,0.033,0.057));
        sdLine(point, color, d, vec2f( 0.252587,-0.370868),vec2f( 0.218651,-0.371118),vec3f(0.441,0.100,0.059),vec3f(0.845,0.380,0.137));
        sdLine(point, color, d, vec2f( 0.276234,-0.360153),vec2f( 0.252950,-0.370523),vec3f(0.496,0.124,0.057),vec3f(0.847,0.412,0.163));
        sdLine(point, color, d, vec2f( 0.296842,-0.316514),vec2f( 0.275936,-0.360588),vec3f(0.533,0.141,0.055),vec3f(0.847,0.431,0.180),vec3f(0.643,0.192,0.051),vec3f(0.910,0.561,0.271));
        sdLine(point, color, d, vec2f( 0.296875,-0.238281),vec2f( 0.296875,-0.316406),vec3f(0.643,0.192,0.051),vec3f(0.910,0.561,0.271),vec3f(0.741,0.239,0.047),vec3f(0.953,0.678,0.357));
        sdLine(point, color, d, vec2f(-0.018745,-0.009492),vec2f(-0.023346,-0.023500),vec3f(0.984,0.665,0.482),vec3f(0.824,0.235,0.180));
        sdLine(point, color, d, vec2f( 0.004475, 0.004419),vec2f(-0.019006,-0.009652),vec3f(0.976,0.686,0.506),vec3f(0.798,0.210,0.153));
        sdLine(point, color, d, vec2f( 0.078024,-0.000098),vec2f( 0.004545, 0.004548),vec3f(0.935,0.696,0.480),vec3f(0.792,0.229,0.122));
        sdLine(point, color, d, vec2f( 0.129006,-0.035229),vec2f( 0.078225,-0.000073),vec3f(0.902,0.686,0.447),vec3f(0.788,0.243,0.102),vec3f(0.937,0.769,0.545),vec3f(0.769,0.255,0.059));
        sdLine(point, color, d, vec2f( 0.214680,-0.370938),vec2f( 0.167824,-0.202999),vec3f(0.263,0.220,0.243),vec3f(0.941,0.643,0.373));
        sdLine(point, color, d, vec2f( 0.117434,-0.128738),vec2f(-0.023211,-0.023293),vec3f(0.718,0.698,0.690),vec3f(0.984,0.867,0.718));
        sdLine(point, color, d, vec2f(-0.027606,-0.023401),vec2f(-0.250232, 0.058604),vec3f(0.369,0.376,0.439),vec3f(0.824,0.357,0.349));
        sdLine(point, color, d, vec2f( 0.224314,-0.490840),vec2f( 0.214794,-0.371444),vec3f(0.110,0.076,0.102),vec3f(0.553,0.102,0.061));
        sdLine(point, color, d, vec2f( 0.210828,-0.578506),vec2f( 0.224257,-0.490867),vec3f(0.118,0.094,0.114),vec3f(0.467,0.108,0.084));
        sdLine(point, color, d, vec2f(-0.738216, 0.577331),vec2f(-0.835883, 0.550764),vec3f(0.647,0.639,0.914),vec3f(0.859,0.780,0.682));
        sdLine(point, color, d, vec2f(-0.428590, 0.601696),vec2f(-0.738185, 0.577395),vec3f(0.698,0.680,0.884),vec3f(0.865,0.778,0.667));
        sdLine(point, color, d, vec2f(-0.257866, 0.593737),vec2f(-0.428754, 0.601618),vec3f(0.763,0.733,0.845),vec3f(0.876,0.775,0.645));
        sdLine(point, color, d, vec2f(-0.136723, 0.554703),vec2f(-0.257809, 0.593759),vec3f(0.808,0.769,0.820),vec3f(0.878,0.775,0.657));
        sdLine(point, color, d, vec2f(-0.910838, 0.948616),vec2f(-0.980469, 0.968751),vec3f(0.782,0.753,0.910),vec3f(0.486,0.422,0.769));
        sdLine(point, color, d, vec2f(-0.849120, 0.905236),vec2f(-0.910830, 0.948608),vec3f(0.771,0.741,0.908),vec3f(0.494,0.431,0.773));
        sdLine(point, color, d, vec2f(-0.805151, 0.844085),vec2f(-0.849092, 0.905222),vec3f(0.751,0.724,0.904),vec3f(0.494,0.431,0.773));
        sdLine(point, color, d, vec2f(-0.792960, 0.808567),vec2f(-0.805191, 0.844089),vec3f(0.741,0.714,0.902),vec3f(0.494,0.431,0.773));
        sdLine(point, color, d, vec2f(-0.660730, 0.350115),vec2f(-0.687380, 0.324080),vec3f(0.569,0.539,0.929),vec3f(0.625,0.612,0.963));
        sdLine(point, color, d, vec2f(-0.581342, 0.390737),vec2f(-0.660931, 0.350300),vec3f(0.610,0.565,0.941),vec3f(0.671,0.641,0.963));
        sdLine(point, color, d, vec2f(-0.519169, 0.405706),vec2f(-0.581260, 0.390541),vec3f(0.653,0.590,0.953),vec3f(0.722,0.671,0.963));
        sdLine(point, color, d, vec2f(-0.325640, 0.410142),vec2f(-0.519014, 0.405733),vec3f(0.708,0.624,0.969),vec3f(0.784,0.706,0.965));
        sdLine(point, color, d, vec2f(-0.265662, 0.386779),vec2f(-0.325886, 0.410329),vec3f(0.767,0.661,0.982),vec3f(0.847,0.741,0.965));
        sdLine(point, color, d, vec2f(-0.890704, 0.394153),vec2f(-0.996080, 0.374882),vec3f(0.576,0.569,0.918),vec3f(0.635,0.637,0.951));
        sdLine(point, color, d, vec2f(-0.808329, 0.393799),vec2f(-0.890754, 0.394278),vec3f(0.575,0.569,0.922),vec3f(0.635,0.635,0.955));
        sdLine(point, color, d, vec2f(-0.746093, 0.374902),vec2f(-0.808295, 0.393687),vec3f(0.573,0.569,0.924),vec3f(0.635,0.629,0.961));
        sdLine(point, color, d, vec2f(-0.691314, 0.328133),vec2f(-0.746006, 0.375002),vec3f(0.571,0.569,0.927),vec3f(0.635,0.624,0.969));
        sdLine(point, color, d, vec2f(-0.713022, 0.223982),vec2f(-0.699215, 0.078105),vec3f(0.443,0.404,0.831),vec3f(0.490,0.471,0.855),vec3f(0.486,0.455,0.871),vec3f(0.576,0.573,0.929));
        sdLine(point, color, d, vec2f(-0.687633, 0.328193),vec2f(-0.713148, 0.224068),vec3f(0.502,0.473,0.884),vec3f(0.590,0.586,0.943));
        sdLine(point, color, d, vec2f(-0.962077,-0.021061),vec2f(-0.984470,-0.062241),vec3f(0.567,0.553,0.914),vec3f(0.506,0.502,0.843));
        sdLine(point, color, d, vec2f(-0.929730, 0.012170),vec2f(-0.962095,-0.021151),vec3f(0.563,0.543,0.914),vec3f(0.508,0.514,0.829));
        sdLine(point, color, d, vec2f(-0.842855, 0.056023),vec2f(-0.929597, 0.011852),vec3f(0.555,0.527,0.914),vec3f(0.506,0.522,0.804));
        sdLine(point, color, d, vec2f(-0.740058, 0.072302),vec2f(-0.842912, 0.056291),vec3f(0.549,0.518,0.914),vec3f(0.502,0.522,0.788),vec3f(0.537,0.498,0.914),vec3f(0.431,0.412,0.788));
        sdLine(point, color, d, vec2f(-0.648617, 0.066357),vec2f(-0.739962, 0.072359),vec3f(0.531,0.486,0.912),vec3f(0.406,0.369,0.786));
        sdLine(point, color, d, vec2f(-0.599562, 0.052212),vec2f(-0.648693, 0.066421),vec3f(0.522,0.465,0.912),vec3f(0.390,0.339,0.776));
        sdLine(point, color, d, vec2f(-0.523674,-0.011767),vec2f(-0.599782, 0.051913),vec3f(0.510,0.443,0.912),vec3f(0.416,0.371,0.761));
        sdLine(point, color, d, vec2f(-0.216909, 0.467215),vec2f(-0.261703, 0.386747),vec3f(0.943,0.525,0.555),vec3f(0.851,0.796,0.924));
        sdLine(point, color, d, vec2f(-0.140684, 0.554868),vec2f(-0.216981, 0.467356),vec3f(0.898,0.506,0.527),vec3f(0.851,0.786,0.865));
        sdLine(point, color, d, vec2f(-0.294609, 0.328762),vec2f(-0.320297, 0.253955),vec3f(0.945,0.480,0.496),vec3f(0.655,0.545,0.908));
        sdLine(point, color, d, vec2f(-0.261515, 0.386675),vec2f(-0.294411, 0.328660),vec3f(0.953,0.482,0.498),vec3f(0.678,0.576,0.945),vec3f(0.953,0.482,0.494),vec3f(0.776,0.694,0.961));
        sdLine(point, color, d, vec2f(-0.320312, 0.152344),vec2f(-0.320312, 0.250000),vec3f(0.227,0.159,0.249),vec3f(0.892,0.455,0.476));
        sdLine(point, color, d, vec2f(-0.285102, 0.093592),vec2f(-0.320160, 0.152217),vec3f(0.194,0.129,0.220),vec3f(0.878,0.441,0.463));
        sdLine(point, color, d, vec2f(-0.249821, 0.058499),vec2f(-0.285083, 0.093620),vec3f(0.206,0.137,0.220),vec3f(0.861,0.422,0.445));
        sdLine(point, color, d, vec2f(-0.370886, 0.200884),vec2f(-0.296881, 0.058849),vec3f(0.210,0.208,0.296),vec3f(0.506,0.441,0.829));
        sdLine(point, color, d, vec2f(-0.371080, 0.218669),vec2f(-0.370863, 0.200571),vec3f(0.216,0.208,0.275),vec3f(0.539,0.473,0.853));
        sdLine(point, color, d, vec2f(-0.320338, 0.250030),vec2f(-0.371125, 0.218779),vec3f(0.214,0.182,0.245),vec3f(0.573,0.504,0.873));
        sdLine(point, color, d, vec2f(-0.318397, 0.038871),vec2f(-0.296987, 0.058686),vec3f(0.408,0.345,0.751),vec3f(0.235,0.243,0.371));
        sdLine(point, color, d, vec2f(-0.335799, 0.000126),vec2f(-0.318147, 0.038906),vec3f(0.394,0.337,0.686),vec3f(0.225,0.229,0.355));
        sdLine(point, color, d, vec2f(-0.457045, 0.015479),vec2f(-0.343771, 0.003763),vec3f(0.404,0.341,0.757),vec3f(0.714,0.702,0.792),vec3f(0.427,0.369,0.796),vec3f(0.478,0.443,0.510));
        sdLine(point, color, d, vec2f(-0.524601,-0.011832),vec2f(-0.457157, 0.015632),vec3f(0.427,0.369,0.796),vec3f(0.478,0.443,0.510),vec3f(0.475,0.424,0.792),vec3f(0.627,0.596,0.741));
        sdLine(point, color, d, vec2f(-0.567845,-0.048929),vec2f(-0.524475,-0.011939),vec3f(0.475,0.424,0.792),vec3f(0.627,0.596,0.741),vec3f(0.545,0.490,0.871),vec3f(0.522,0.490,0.647));
        sdLine(point, color, d, vec2f(-0.687512,-0.202819),vec2f(-0.567697,-0.048786),vec3f(0.545,0.490,0.871),vec3f(0.522,0.490,0.647),vec3f(0.749,0.706,0.980),vec3f(0.702,0.651,0.710));
        sdLine(point, color, d, vec2f(-0.738370,-0.246033),vec2f(-0.687782,-0.202852),vec3f(0.749,0.706,0.980),vec3f(0.702,0.651,0.710),vec3f(0.690,0.651,0.929),vec3f(0.259,0.235,0.337));
        sdLine(point, color, d, vec2f(-0.819820,-0.267858),vec2f(-0.738322,-0.246131),vec3f(0.690,0.651,0.929),vec3f(0.259,0.235,0.337),vec3f(0.631,0.604,0.933),vec3f(0.380,0.333,0.416));
        sdLine(point, color, d, vec2f(-0.925652,-0.360953),vec2f(-0.819489,-0.267857),vec3f(0.631,0.604,0.933),vec3f(0.380,0.333,0.416),vec3f(0.525,0.502,0.808),vec3f(0.243,0.188,0.333));
        sdLine(point, color, d, vec2f(-0.924077,-0.382938),vec2f(-0.926012,-0.360824),vec3f(0.525,0.502,0.808),vec3f(0.243,0.188,0.333),vec3f(0.506,0.478,0.733),vec3f(0.541,0.486,0.631));
        sdLine(point, color, d, vec2f(-0.908709,-0.396678),vec2f(-0.923914,-0.383129),vec3f(0.502,0.473,0.714),vec3f(0.524,0.476,0.637));
        sdLine(point, color, d, vec2f(-0.877801,-0.327818),vec2f(-0.908834,-0.396482),vec3f(0.498,0.467,0.694),vec3f(0.506,0.467,0.643),vec3f(0.584,0.565,0.894),vec3f(0.490,0.439,0.659));
        sdLine(point, color, d, vec2f(-0.850665,-0.306523),vec2f(-0.877737,-0.328046),vec3f(0.584,0.565,0.894),vec3f(0.490,0.439,0.659),vec3f(0.553,0.510,0.765),vec3f(0.329,0.278,0.439));
        sdLine(point, color, d, vec2f(-0.842943,-0.336355),vec2f(-0.850644,-0.306423),vec3f(0.553,0.510,0.765),vec3f(0.329,0.278,0.439),vec3f(0.518,0.478,0.761),vec3f(0.537,0.494,0.643));
        sdLine(point, color, d, vec2f(-0.646520,-0.253453),vec2f(-0.842985,-0.336477),vec3f(0.518,0.478,0.761),vec3f(0.537,0.494,0.643),vec3f(0.671,0.620,0.957),vec3f(0.467,0.392,0.565));
        sdLine(point, color, d, vec2f(-0.613162,-0.227037),vec2f(-0.646645,-0.253104),vec3f(0.651,0.602,0.945),vec3f(0.469,0.398,0.555));
        sdLine(point, color, d, vec2f(-0.561901,-0.162148),vec2f(-0.613081,-0.227201),vec3f(0.594,0.549,0.910),vec3f(0.476,0.429,0.539));
        sdLine(point, color, d, vec2f(-0.500240,-0.050707),vec2f(-0.562155,-0.162105),vec3f(0.557,0.514,0.886),vec3f(0.482,0.455,0.533),vec3f(0.443,0.408,0.812),vec3f(0.286,0.267,0.369));
        sdLine(point, color, d, vec2f(-0.464788,-0.030266),vec2f(-0.500093,-0.050719),vec3f(0.443,0.408,0.812),vec3f(0.286,0.267,0.369),vec3f(0.408,0.376,0.780),vec3f(0.157,0.114,0.267));
        sdLine(point, color, d, vec2f(-0.439427,-0.045861),vec2f(-0.464616,-0.030239),vec3f(0.408,0.376,0.780),vec3f(0.157,0.114,0.267),vec3f(0.455,0.427,0.784),vec3f(0.416,0.353,0.435));
        sdLine(point, color, d, vec2f(-0.421871,-0.081990),vec2f(-0.439501,-0.045912),vec3f(0.443,0.418,0.773),vec3f(0.404,0.369,0.437));
        sdLine(point, color, d, vec2f(-0.386711,-0.105557),vec2f(-0.421868,-0.082119),vec3f(0.431,0.408,0.761),vec3f(0.392,0.384,0.439),vec3f(0.357,0.349,0.655),vec3f(0.376,0.427,0.447));
        sdLine(point, color, d, vec2f(-0.339923, 0.003866),vec2f(-0.386802,-0.101587),vec3f(0.200,0.196,0.267),vec3f(0.345,0.333,0.365));
        sdLine(point, color, d, vec2f(-0.407061,-0.201258),vec2f(-0.386665,-0.105627),vec3f(0.361,0.341,0.624),vec3f(0.792,0.812,0.751));
        sdLine(point, color, d, vec2f(-0.403323,-0.268221),vec2f(-0.406911,-0.201492),vec3f(0.363,0.329,0.606),vec3f(0.771,0.769,0.724));
        sdLine(point, color, d, vec2f(-0.390770,-0.273672),vec2f(-0.403673,-0.268054),vec3f(0.384,0.333,0.643),vec3f(0.771,0.749,0.704));
        sdLine(point, color, d, vec2f(-0.345881,-0.060745),vec2f(-0.386656,-0.101656),vec3f(0.776,0.835,0.784),vec3f(0.361,0.400,0.427),vec3f(0.792,0.839,0.780),vec3f(0.145,0.176,0.227));
        sdLine(point, color, d, vec2f(-0.328983,-0.056141),vec2f(-0.345818,-0.060743),vec3f(0.792,0.839,0.780),vec3f(0.145,0.176,0.227),vec3f(0.808,0.847,0.780),vec3f(0.220,0.251,0.310));
        sdLine(point, color, d, vec2f(-0.252749,-0.126588),vec2f(-0.329397,-0.055845),vec3f(0.837,0.859,0.775),vec3f(0.227,0.247,0.306));
        sdLine(point, color, d, vec2f(-0.230481,-0.160125),vec2f(-0.252484,-0.126757),vec3f(0.867,0.871,0.769),vec3f(0.235,0.243,0.302),vec3f(0.859,0.851,0.749),vec3f(0.161,0.153,0.204));
        sdLine(point, color, d, vec2f(-0.226341,-0.219167),vec2f(-0.230372,-0.160222),vec3f(0.882,0.869,0.757),vec3f(0.151,0.133,0.180));
        sdLine(point, color, d, vec2f(-0.277439,-0.343654),vec2f(-0.226534,-0.219015),vec3f(0.906,0.886,0.765),vec3f(0.141,0.114,0.157),vec3f(0.929,0.898,0.761),vec3f(0.251,0.212,0.239));
        sdLine(point, color, d, vec2f(-0.331270,-0.274257),vec2f(-0.390414,-0.273335),vec3f(0.076,0.061,0.096),vec3f(0.808,0.780,0.714));
        sdLine(point, color, d, vec2f(-0.307108,-0.285553),vec2f(-0.331497,-0.274239),vec3f(0.075,0.059,0.094),vec3f(0.829,0.806,0.718));
        sdLine(point, color, d, vec2f(-0.285269,-0.316297),vec2f(-0.307202,-0.285572),vec3f(0.075,0.059,0.094),vec3f(0.835,0.820,0.722));
        sdLine(point, color, d, vec2f(-0.277101,-0.348165),vec2f(-0.284918,-0.316904),vec3f(0.075,0.059,0.094),vec3f(0.831,0.820,0.722),vec3f(0.086,0.063,0.102),vec3f(0.620,0.584,0.522));
        sdLine(point, color, d, vec2f( 0.191818,-0.655109),vec2f( 0.210666,-0.574643),vec3f(0.137,0.122,0.129),vec3f(0.729,0.608,0.541),vec3f(0.145,0.133,0.133),vec3f(0.804,0.737,0.616));
        sdLine(point, color, d, vec2f( 0.177802,-0.664657),vec2f( 0.191779,-0.655002),vec3f(0.145,0.133,0.133),vec3f(0.804,0.737,0.616),vec3f(0.063,0.059,0.075),vec3f(0.796,0.753,0.620));
        sdLine(point, color, d, vec2f( 0.168515,-0.738957),vec2f( 0.178126,-0.664256),vec3f(0.096,0.092,0.100),vec3f(0.812,0.759,0.625));
        sdLine(point, color, d, vec2f( 0.155592,-0.752345),vec2f( 0.168584,-0.738793),vec3f(0.104,0.100,0.102),vec3f(0.824,0.757,0.627));
        sdLine(point, color, d, vec2f( 0.132513,-0.757861),vec2f( 0.155207,-0.752647),vec3f(0.073,0.069,0.073),vec3f(0.820,0.747,0.624));
        sdLine(point, color, d, vec2f( 0.086261,-0.759725),vec2f(-0.023197,-0.695559),vec3f(0.157,0.141,0.184),vec3f(0.063,0.055,0.067),vec3f(0.286,0.275,0.325),vec3f(0.055,0.043,0.055));
        sdLine(point, color, d, vec2f( 0.132843,-0.757744),vec2f( 0.086059,-0.759407),vec3f(0.286,0.275,0.325),vec3f(0.055,0.043,0.055),vec3f(0.761,0.722,0.522),vec3f(0.059,0.043,0.059));
        sdLine(point, color, d, vec2f(-0.008536,-0.523957),vec2f(-0.113095,-0.546922),vec3f(0.761,0.722,0.522),vec3f(0.051,0.033,0.051));
        sdLine(point, color, d, vec2f( 0.028417,-0.535309),vec2f(-0.008810,-0.523762),vec3f(0.792,0.767,0.563),vec3f(0.055,0.041,0.057));
        sdLine(point, color, d, vec2f( 0.049655,-0.559962),vec2f( 0.028517,-0.535386),vec3f(0.798,0.782,0.580),vec3f(0.076,0.065,0.073));
        sdLine(point, color, d, vec2f( 0.071714,-0.659300),vec2f( 0.049649,-0.559687),vec3f(0.771,0.751,0.557),vec3f(0.082,0.073,0.080));
        sdLine(point, color, d, vec2f( 0.062107,-0.673453),vec2f( 0.072068,-0.659420),vec3f(0.731,0.704,0.520),vec3f(0.088,0.080,0.088));
        sdLine(point, color, d, vec2f(-0.019664,-0.695497),vec2f( 0.061596,-0.673865),vec3f(0.722,0.694,0.514),vec3f(0.090,0.082,0.090),vec3f(0.722,0.690,0.502),vec3f(0.220,0.216,0.180));
        sdLine(point, color, d, vec2f(-0.103143,-0.577233),vec2f(-0.113704,-0.546335),vec3f(0.078,0.063,0.094),vec3f(0.669,0.586,0.406));
        sdLine(point, color, d, vec2f(-0.128872,-0.656514),vec2f(-0.102680,-0.578043),vec3f(0.078,0.063,0.094),vec3f(0.643,0.573,0.396),vec3f(0.078,0.063,0.094),vec3f(0.545,0.506,0.349));
        sdLine(point, color, d, vec2f(-0.110341,-0.684318),vec2f(-0.132708,-0.656415),vec3f(0.090,0.076,0.129),vec3f(0.712,0.680,0.494));
        sdLine(point, color, d, vec2f(-0.065097,-0.698489),vec2f(-0.110606,-0.684138),vec3f(0.106,0.094,0.133),vec3f(0.712,0.680,0.494));
        sdLine(point, color, d, vec2f(-0.023516,-0.695404),vec2f(-0.065017,-0.698593),vec3f(0.124,0.112,0.137),vec3f(0.718,0.686,0.498));
        sdLine(point, color, d, vec2f(-0.356409,-0.396602),vec2f(-0.363344,-0.417840),vec3f(0.667,0.684,0.514),vec3f(0.086,0.073,0.102));
        sdLine(point, color, d, vec2f(-0.329514,-0.372734),vec2f(-0.356301,-0.396902),vec3f(0.678,0.694,0.525),vec3f(0.098,0.086,0.110),vec3f(0.796,0.792,0.624),vec3f(0.075,0.059,0.094));
        sdLine(point, color, d, vec2f(-0.304495,-0.384405),vec2f(-0.329372,-0.372760),vec3f(0.794,0.786,0.625),vec3f(0.073,0.055,0.094));
        sdLine(point, color, d, vec2f(-0.304686,-0.445287),vec2f(-0.304686,-0.384176),vec3f(0.792,0.780,0.627),vec3f(0.071,0.051,0.094),vec3f(0.678,0.667,0.541),vec3f(0.071,0.047,0.094));
        sdLine(point, color, d, vec2f(-0.329703,-0.471680),vec2f(-0.304762,-0.445327),vec3f(0.678,0.667,0.541),vec3f(0.071,0.047,0.094),vec3f(0.533,0.522,0.427),vec3f(0.071,0.051,0.094));
        sdLine(point, color, d, vec2f(-0.352235,-0.461623),vec2f(-0.329757,-0.471344),vec3f(0.504,0.494,0.406),vec3f(0.071,0.051,0.094));
        sdLine(point, color, d, vec2f(-0.363307,-0.418071),vec2f(-0.352136,-0.462043),vec3f(0.475,0.467,0.384),vec3f(0.071,0.051,0.094),vec3f(0.561,0.565,0.439),vec3f(0.129,0.118,0.137));
        sdLine(point, color, d, vec2f(-0.276110,-0.521710),vec2f(-0.289045,-0.550497),vec3f(0.641,0.637,0.478),vec3f(0.035,0.033,0.063));
        sdLine(point, color, d, vec2f(-0.258075,-0.516729),vec2f(-0.275890,-0.521970),vec3f(0.684,0.682,0.514),vec3f(0.047,0.045,0.075));
        sdLine(point, color, d, vec2f(-0.206378,-0.541564),vec2f(-0.258683,-0.516738),vec3f(0.706,0.706,0.529),vec3f(0.047,0.047,0.078),vec3f(0.792,0.788,0.600),vec3f(0.051,0.055,0.063));
        sdLine(point, color, d, vec2f(-0.187895,-0.563812),vec2f(-0.206382,-0.541493),vec3f(0.757,0.753,0.576),vec3f(0.063,0.057,0.078));
        sdLine(point, color, d, vec2f(-0.189739,-0.573083),vec2f(-0.187133,-0.563528),vec3f(0.735,0.731,0.561),vec3f(0.075,0.059,0.094));
        sdLine(point, color, d, vec2f(-0.251815,-0.589992),vec2f(-0.189995,-0.573544),vec3f(0.749,0.745,0.569),vec3f(0.075,0.059,0.094),vec3f(0.498,0.502,0.345),vec3f(0.075,0.063,0.098));
        sdLine(point, color, d, vec2f(-0.272446,-0.579505),vec2f(-0.251715,-0.590033),vec3f(0.514,0.518,0.363),vec3f(0.075,0.063,0.098));
        sdLine(point, color, d, vec2f(-0.289115,-0.550792),vec2f(-0.272717,-0.579388),vec3f(0.549,0.551,0.398),vec3f(0.086,0.075,0.106));
        sdLine(point, color, d, vec2f(-0.192736,-0.469520),vec2f(-0.277483,-0.347082),vec3f(0.108,0.076,0.120),vec3f(0.065,0.057,0.090));
        sdLine(point, color, d, vec2f(-0.113128,-0.547007),vec2f(-0.192449,-0.470244),vec3f(0.080,0.045,0.080),vec3f(0.063,0.047,0.073));
        sdLine(point, color, d, vec2f(-0.399689,-0.363777),vec2f(-0.390465,-0.273580),vec3f(0.355,0.292,0.675),vec3f(0.088,0.075,0.112));
        sdLine(point, color, d, vec2f(-0.414173,-0.390763),vec2f(-0.399958,-0.363772),vec3f(0.347,0.288,0.651),vec3f(0.090,0.075,0.110));
        sdLine(point, color, d, vec2f( 0.402104,-0.164371),vec2f( 0.429754,-0.097640),vec3f(0.237,0.182,0.180),vec3f(0.784,0.296,0.143));
        sdLine(point, color, d, vec2f( 0.366996,-0.206976),vec2f( 0.401848,-0.164337),vec3f(0.210,0.145,0.145),vec3f(0.778,0.304,0.147));
        sdLine(point, color, d, vec2f( 0.296913,-0.238165),vec2f( 0.367225,-0.206912),vec3f(0.194,0.118,0.120),vec3f(0.743,0.304,0.149));
        sdLine(point, color, d, vec2f( 0.167900,-0.203204),vec2f( 0.113213,-0.128975),vec3f(0.706,0.690,0.678),vec3f(0.145,0.110,0.157));
        sdLine(point, color, d, vec2f( 0.454448,-0.085567),vec2f( 0.429490,-0.097974),vec3f(0.898,0.478,0.357),vec3f(0.945,0.643,0.580));
        sdLine(point, color, d, vec2f( 0.501070,-0.040927),vec2f( 0.454701,-0.085293),vec3f(0.898,0.475,0.349),vec3f(0.941,0.631,0.561),vec3f(0.902,0.439,0.278),vec3f(0.910,0.525,0.416));
        sdLine(point, color, d, vec2f( 0.536751, 0.019080),vec2f( 0.501035,-0.040870),vec3f(0.904,0.425,0.247),vec3f(0.902,0.502,0.384));
        sdLine(point, color, d, vec2f( 0.554739, 0.081950),vec2f( 0.536789, 0.018991),vec3f(0.906,0.412,0.216),vec3f(0.900,0.508,0.390));
        sdLine(point, color, d, vec2f( 0.535861, 0.166934),vec2f( 0.554677, 0.082064),vec3f(0.906,0.412,0.216),vec3f(0.906,0.537,0.427),vec3f(0.906,0.416,0.224),vec3f(0.925,0.647,0.557));
        sdLine(point, color, d, vec2f( 0.419826, 0.338301),vec2f( 0.535615, 0.166647),vec3f(0.906,0.414,0.220),vec3f(0.920,0.614,0.524));
        sdLine(point, color, d, vec2f( 0.406078, 0.390921),vec2f( 0.419883, 0.338831),vec3f(0.906,0.412,0.216),vec3f(0.910,0.565,0.473));
        sdLine(point, color, d, vec2f( 0.146477, 0.567224),vec2f( 0.031189, 0.523314),vec3f(0.863,0.475,0.404),vec3f(0.963,0.447,0.278));
        sdLine(point, color, d, vec2f( 0.238636, 0.574102),vec2f( 0.146891, 0.567233),vec3f(0.869,0.484,0.425),vec3f(0.957,0.427,0.241));
        sdLine(point, color, d, vec2f( 0.276654, 0.562718),vec2f( 0.238421, 0.574198),vec3f(0.876,0.496,0.445),vec3f(0.949,0.400,0.192));
        sdLine(point, color, d, vec2f( 0.328054, 0.527362),vec2f( 0.276439, 0.562762),vec3f(0.880,0.508,0.463),vec3f(0.939,0.382,0.165));
        sdLine(point, color, d, vec2f( 0.407058, 0.213498),vec2f( 0.371072, 0.195345),vec3f(0.988,0.973,0.969),vec3f(0.992,0.855,0.831));
        sdLine(point, color, d, vec2f( 0.458303, 0.213407),vec2f( 0.407282, 0.213566),vec3f(0.988,0.973,0.969),vec3f(0.992,0.855,0.831));
        sdLine(point, color, d, vec2f( 0.481278, 0.196816),vec2f( 0.457990, 0.213665),vec3f(0.988,0.973,0.967),vec3f(0.992,0.855,0.829));
        sdLine(point, color, d, vec2f( 0.499009, 0.163011),vec2f( 0.481392, 0.196333),vec3f(0.988,0.973,0.965),vec3f(0.992,0.855,0.827));
        sdLine(point, color, d, vec2f( 0.503604, 0.090461),vec2f( 0.499073, 0.163192),vec3f(0.988,0.971,0.965),vec3f(0.990,0.855,0.827));
        sdLine(point, color, d, vec2f( 0.489446, 0.076448),vec2f( 0.503554, 0.090238),vec3f(0.988,0.969,0.963),vec3f(0.986,0.847,0.816));
        sdLine(point, color, d, vec2f( 0.443782, 0.062339),vec2f( 0.489220, 0.076528),vec3f(0.986,0.969,0.961),vec3f(0.986,0.845,0.814));
        sdLine(point, color, d, vec2f( 0.378683, 0.066466),vec2f( 0.443770, 0.062469),vec3f(0.984,0.969,0.959),vec3f(0.988,0.851,0.822));
        sdLine(point, color, d, vec2f( 0.339578, 0.097522),vec2f( 0.378706, 0.066500),vec3f(0.984,0.967,0.957),vec3f(0.986,0.851,0.820));
        sdLine(point, color, d, vec2f( 0.329906, 0.140698),vec2f( 0.339788, 0.097651),vec3f(0.984,0.965,0.955),vec3f(0.984,0.851,0.818));
        sdLine(point, color, d, vec2f( 0.349815, 0.180720),vec2f( 0.329907, 0.140407),vec3f(0.984,0.965,0.953),vec3f(0.982,0.849,0.814));
        sdLine(point, color, d, vec2f( 0.371027, 0.195457),vec2f( 0.349730, 0.180934),vec3f(0.984,0.965,0.953),vec3f(0.982,0.849,0.814));
        sdLine(point, color, d, vec2f(-0.117425, 0.341220),vec2f(-0.074257, 0.417986),vec3f(0.847,0.318,0.216),vec3f(0.794,0.425,0.357));
        sdLine(point, color, d, vec2f(-0.131058, 0.286069),vec2f(-0.117280, 0.341245),vec3f(0.847,0.318,0.216),vec3f(0.814,0.429,0.355));
        sdLine(point, color, d, vec2f(-0.120979, 0.239269),vec2f(-0.131229, 0.286301),vec3f(0.847,0.318,0.216),vec3f(0.831,0.435,0.355));
        sdLine(point, color, d, vec2f(-0.094915, 0.195662),vec2f(-0.120748, 0.238869),vec3f(0.853,0.308,0.194),vec3f(0.847,0.439,0.355));
        sdLine(point, color, d, vec2f(-0.058454, 0.163892),vec2f(-0.094941, 0.195621),vec3f(0.861,0.296,0.171),vec3f(0.861,0.439,0.351));
        sdLine(point, color, d, vec2f( 0.096154, 0.101755),vec2f(-0.058568, 0.164000),vec3f(0.855,0.304,0.188),vec3f(0.884,0.443,0.349));
        sdLine(point, color, d, vec2f( 0.148627, 0.066296),vec2f( 0.096309, 0.101713),vec3f(0.845,0.316,0.214),vec3f(0.906,0.447,0.349));
        sdLine(point, color, d, vec2f( 0.238337,-0.218583),vec2f( 0.296931,-0.238114),vec3f(0.188,0.059,0.055),vec3f(0.941,0.714,0.384));
        sdLine(point, color, d, vec2f( 0.356732, 0.495659),vec2f( 0.328217, 0.527093),vec3f(0.475,0.510,0.604),vec3f(0.418,0.445,0.533));
        sdLine(point, color, d, vec2f( 0.406416, 0.390930),vec2f( 0.356795, 0.496206),vec3f(0.475,0.510,0.604),vec3f(0.398,0.414,0.496));
        sdLine(point, color, d, vec2f( 0.351847, 0.563942),vec2f( 0.328272, 0.527247),vec3f(0.439,0.475,0.569),vec3f(0.725,0.459,0.408),vec3f(0.443,0.475,0.569),vec3f(0.882,0.420,0.259));
        sdLine(point, color, d, vec2f( 0.386954, 0.567493),vec2f( 0.352309, 0.564131),vec3f(0.431,0.453,0.537),vec3f(0.884,0.416,0.247));
        sdLine(point, color, d, vec2f( 0.434360, 0.555848),vec2f( 0.386390, 0.567397),vec3f(0.420,0.431,0.506),vec3f(0.886,0.412,0.235),vec3f(0.361,0.341,0.396),vec3f(0.882,0.416,0.247));
        sdLine(point, color, d, vec2f( 0.481493, 0.527439),vec2f( 0.434614, 0.555661),vec3f(0.376,0.353,0.406),vec3f(0.882,0.416,0.247));
        sdLine(point, color, d, vec2f( 0.556668, 0.460958),vec2f( 0.481472, 0.527483),vec3f(0.392,0.365,0.416),vec3f(0.882,0.416,0.247),vec3f(0.522,0.502,0.557),vec3f(0.882,0.412,0.247));
        sdLine(point, color, d, vec2f( 0.568077, 0.433950),vec2f( 0.556441, 0.461051),vec3f(0.524,0.506,0.563),vec3f(0.880,0.414,0.247));
        sdLine(point, color, d, vec2f( 0.550912, 0.394412),vec2f( 0.568131, 0.433870),vec3f(0.494,0.480,0.541),vec3f(0.880,0.416,0.247));
        sdLine(point, color, d, vec2f( 0.490768, 0.382274),vec2f( 0.550937, 0.394469),vec3f(0.422,0.425,0.496),vec3f(0.882,0.416,0.247));
        sdLine(point, color, d, vec2f( 0.406039, 0.390810),vec2f( 0.490402, 0.382519),vec3f(0.375,0.382,0.457),vec3f(0.882,0.418,0.253));
        sdLine(point, color, d, vec2f( 0.330744, 0.478846),vec2f( 0.328238, 0.527133),vec3f(0.908,0.586,0.498),vec3f(0.475,0.510,0.604));
        sdLine(point, color, d, vec2f( 0.365538, 0.422127),vec2f( 0.330688, 0.479010),vec3f(0.910,0.588,0.498),vec3f(0.475,0.510,0.604));
        sdLine(point, color, d, vec2f( 0.406309, 0.390576),vec2f( 0.365537, 0.422123),vec3f(0.910,0.588,0.498),vec3f(0.475,0.502,0.592));
        sdLine(point, color, d, vec2f( 0.037170, 0.502175),vec2f( 0.031325, 0.523294),vec3f(0.357,0.220,0.255),vec3f(0.725,0.404,0.349),vec3f(0.357,0.220,0.255),vec3f(0.914,0.506,0.408));
        sdLine(point, color, d, vec2f( 0.009112, 0.456598),vec2f( 0.036991, 0.502416),vec3f(0.357,0.220,0.255),vec3f(0.900,0.484,0.404));
        sdLine(point, color, d, vec2f(-0.014264, 0.437994),vec2f( 0.009264, 0.456398),vec3f(0.357,0.220,0.255),vec3f(0.882,0.457,0.398));
        sdLine(point, color, d, vec2f(-0.074384, 0.418060),vec2f(-0.014479, 0.438192),vec3f(0.357,0.220,0.255),vec3f(0.873,0.443,0.394));
        sdLine(point, color, d, vec2f(-0.034865, 0.469095),vec2f(-0.074326, 0.418051),vec3f(0.357,0.220,0.255),vec3f(0.325,0.153,0.310));
        sdLine(point, color, d, vec2f( 0.031160, 0.519567),vec2f(-0.034850, 0.469051),vec3f(0.357,0.220,0.255),vec3f(0.312,0.147,0.269));
        sdLine(point, color, d, vec2f(-0.476788,-0.270699),vec2f(-0.413803,-0.390958),vec3f(0.400,0.333,0.724),vec3f(0.529,0.406,0.447));
        sdLine(point, color, d, vec2f(-0.515771,-0.226315),vec2f(-0.477190,-0.270110),vec3f(0.439,0.367,0.784),vec3f(0.571,0.437,0.476));
        sdLine(point, color, d, vec2f(-0.548925,-0.216192),vec2f(-0.515512,-0.226748),vec3f(0.475,0.408,0.808),vec3f(0.588,0.451,0.490),vec3f(0.580,0.529,0.878),vec3f(0.404,0.282,0.416));
        sdLine(point, color, d, vec2f(-0.563920,-0.221379),vec2f(-0.548931,-0.216115),vec3f(0.580,0.529,0.878),vec3f(0.404,0.282,0.416),vec3f(0.671,0.635,0.937),vec3f(0.533,0.455,0.702));
        sdLine(point, color, d, vec2f(-0.487346,-0.300005),vec2f(-0.564121,-0.221089),vec3f(0.671,0.635,0.937),vec3f(0.533,0.455,0.702),vec3f(0.502,0.459,0.749),vec3f(0.522,0.392,0.561));
        sdLine(point, color, d, vec2f(-0.456709,-0.351242),vec2f(-0.487322,-0.300058),vec3f(0.502,0.459,0.749),vec3f(0.522,0.392,0.561),vec3f(0.396,0.349,0.639),vec3f(0.412,0.290,0.471));
        sdLine(point, color, d, vec2f(-0.437438,-0.426111),vec2f(-0.456570,-0.351699),vec3f(0.396,0.349,0.639),vec3f(0.412,0.290,0.471),vec3f(0.294,0.243,0.537),vec3f(0.553,0.475,0.596));
        sdLine(point, color, d, vec2f(-0.426103,-0.437341),vec2f(-0.437819,-0.425625),vec3f(0.294,0.243,0.537),vec3f(0.553,0.475,0.596),vec3f(0.380,0.345,0.600),vec3f(0.659,0.600,0.761));
        sdLine(point, color, d, vec2f(-0.137942,-0.735671),vec2f(-0.234693,-0.652174),vec3f(0.390,0.329,0.331),vec3f(0.078,0.075,0.145));
        sdLine(point, color, d, vec2f(-0.101453,-0.787596),vec2f(-0.137531,-0.736006),vec3f(0.376,0.318,0.302),vec3f(0.067,0.067,0.141),vec3f(0.729,0.698,0.671),vec3f(0.051,0.055,0.137));
        sdLine(point, color, d, vec2f(-0.093467,-0.816814),vec2f(-0.101262,-0.787835),vec3f(0.729,0.698,0.671),vec3f(0.051,0.055,0.137),vec3f(0.592,0.549,0.522),vec3f(0.059,0.059,0.141));
        sdLine(point, color, d, vec2f(-0.101072,-0.835520),vec2f(-0.094078,-0.816636),vec3f(0.592,0.549,0.522),vec3f(0.059,0.059,0.141),vec3f(0.427,0.369,0.341),vec3f(0.067,0.067,0.145));
        sdLine(point, color, d, vec2f(-0.120105,-0.842187),vec2f(-0.100815,-0.835333),vec3f(0.427,0.369,0.341),vec3f(0.067,0.067,0.145),vec3f(0.541,0.467,0.439),vec3f(0.129,0.129,0.192));
        sdLine(point, color, d, vec2f(-0.139116,-0.832593),vec2f(-0.120099,-0.842135),vec3f(0.541,0.467,0.439),vec3f(0.129,0.129,0.192),vec3f(0.651,0.561,0.541),vec3f(0.220,0.212,0.298));
        sdLine(point, color, d, vec2f(-0.182775,-0.748898),vec2f(-0.139277,-0.832422),vec3f(0.651,0.561,0.541),vec3f(0.220,0.212,0.298),vec3f(0.176,0.090,0.090),vec3f(0.329,0.290,0.529));
        sdLine(point, color, d, vec2f(-0.242482,-0.691315),vec2f(-0.182860,-0.748982),vec3f(0.176,0.090,0.090),vec3f(0.329,0.290,0.529),vec3f(0.376,0.247,0.227),vec3f(0.247,0.212,0.412));
        sdLine(point, color, d, vec2f(-0.376186,-0.507215),vec2f(-0.414107,-0.390583),vec3f(0.341,0.347,0.365),vec3f(0.071,0.061,0.108));
        sdLine(point, color, d, vec2f(-0.334872,-0.560660),vec2f(-0.375961,-0.507315),vec3f(0.347,0.353,0.337),vec3f(0.061,0.051,0.096));
        sdLine(point, color, d, vec2f(-0.234322,-0.652779),vec2f(-0.334981,-0.561046),vec3f(0.333,0.337,0.322),vec3f(0.059,0.047,0.090),vec3f(0.251,0.251,0.255),vec3f(0.067,0.063,0.098));
        sdLine(point, color, d, vec2f(-0.443668,-0.478495),vec2f(-0.425835,-0.437392),vec3f(0.365,0.322,0.592),vec3f(0.392,0.392,0.510),vec3f(0.353,0.310,0.580),vec3f(0.235,0.227,0.365));
        sdLine(point, color, d, vec2f(-0.439128,-0.529675),vec2f(-0.443433,-0.478735),vec3f(0.353,0.310,0.580),vec3f(0.235,0.227,0.365),vec3f(0.282,0.275,0.471),vec3f(0.110,0.106,0.200));
        sdLine(point, color, d, vec2f(-0.375449,-0.627070),vec2f(-0.439486,-0.529481),vec3f(0.282,0.275,0.471),vec3f(0.110,0.106,0.200),vec3f(0.361,0.314,0.592),vec3f(0.110,0.102,0.133));
        sdLine(point, color, d, vec2f(-0.332113,-0.645086),vec2f(-0.375323,-0.626910),vec3f(0.361,0.314,0.592),vec3f(0.131,0.122,0.145));
        sdLine(point, color, d, vec2f(-0.265559,-0.648265),vec2f(-0.332002,-0.645122),vec3f(0.361,0.314,0.592),vec3f(0.165,0.153,0.165));
        sdLine(point, color, d, vec2f( 0.205231,-0.817938),vec2f( 0.109033,-0.765411),vec3f(0.429,0.459,0.576),vec3f(0.820,0.737,0.624));
        sdLine(point, color, d, vec2f( 0.293818,-0.895384),vec2f( 0.205808,-0.818524),vec3f(0.447,0.478,0.592),vec3f(0.820,0.737,0.624),vec3f(0.529,0.557,0.659),vec3f(0.820,0.737,0.624));
        sdLine(point, color, d, vec2f( 0.354093,-0.996256),vec2f( 0.293952,-0.895611),vec3f(0.529,0.557,0.659),vec3f(0.820,0.737,0.624),vec3f(0.612,0.635,0.722),vec3f(0.820,0.737,0.624));
        sdLine(point, color, d, vec2f( 0.367427,-1.054936),vec2f( 0.353965,-0.995872),vec3f(0.612,0.635,0.722),vec3f(0.820,0.737,0.624),vec3f(0.725,0.690,0.663),vec3f(0.816,0.737,0.624));
        sdLine(point, color, d, vec2f(-0.954948,-0.473706),vec2f(-0.984357,-0.410454),vec3f(0.161,0.147,0.363),vec3f(0.398,0.363,0.714));
        sdLine(point, color, d, vec2f(-0.903811,-0.638807),vec2f(-0.954931,-0.473614),vec3f(0.153,0.141,0.349),vec3f(0.404,0.369,0.718),vec3f(0.067,0.063,0.188),vec3f(0.420,0.392,0.722));
        sdLine(point, color, d, vec2f(-0.875055,-0.699107),vec2f(-0.903928,-0.638463),vec3f(0.073,0.069,0.202),vec3f(0.408,0.380,0.702));
        sdLine(point, color, d, vec2f(-0.629447,-0.956558),vec2f(-0.687403,-0.870736),vec3f(0.029,0.029,0.135),vec3f(0.120,0.114,0.286));
        sdLine(point, color, d, vec2f(-0.605535,-1.023209),vec2f(-0.629620,-0.956694),vec3f(0.029,0.029,0.125),vec3f(0.096,0.096,0.247));
        sdLine(point, color, d, vec2f(-0.265320,-0.663013),vec2f(-0.242156,-0.691406),vec3f(0.216,0.129,0.125),vec3f(0.361,0.310,0.592),vec3f(0.325,0.267,0.251),vec3f(0.361,0.310,0.592));
        sdLine(point, color, d, vec2f(-0.265459,-0.648465),vec2f(-0.265186,-0.663040),vec3f(0.355,0.298,0.288),vec3f(0.355,0.306,0.580));
        sdLine(point, color, d, vec2f(-0.630864,-0.893412),vec2f(-0.687487,-0.870996),vec3f(0.125,0.112,0.275),vec3f(0.259,0.253,0.482));
        sdLine(point, color, d, vec2f(-0.586520,-0.885655),vec2f(-0.630843,-0.893310),vec3f(0.108,0.092,0.233),vec3f(0.208,0.206,0.418));
        sdLine(point, color, d, vec2f(-0.487897,-0.819075),vec2f(-0.586679,-0.885778),vec3f(0.104,0.088,0.218),vec3f(0.176,0.176,0.375));
        sdLine(point, color, d, vec2f(-0.300556,-0.770950),vec2f(-0.487735,-0.819175),vec3f(0.098,0.082,0.188),vec3f(0.188,0.182,0.373));
        sdLine(point, color, d, vec2f(-0.274553,-0.751048),vec2f(-0.300678,-0.771091),vec3f(0.094,0.078,0.163),vec3f(0.202,0.188,0.363));
        sdLine(point, color, d, vec2f(-0.242248,-0.691634),vec2f(-0.274548,-0.751130),vec3f(0.090,0.075,0.147),vec3f(0.216,0.196,0.359));
        sdLine(point, color, d, vec2f(-0.800847,-0.800581),vec2f(-0.875072,-0.699008),vec3f(0.069,0.057,0.225),vec3f(0.351,0.331,0.616));
        sdLine(point, color, d, vec2f(-0.687515,-0.871181),vec2f(-0.800804,-0.800859),vec3f(0.055,0.049,0.192),vec3f(0.331,0.316,0.584));
        sdLine(point, color, d, vec2f(-0.727078, 0.837400),vec2f(-0.792967, 0.812525),vec3f(0.851,0.773,0.682),vec3f(0.453,0.390,0.737));
        sdLine(point, color, d, vec2f(-0.657094, 0.890944),vec2f(-0.727080, 0.837365),vec3f(0.859,0.780,0.682),vec3f(0.475,0.412,0.759));
        sdLine(point, color, d, vec2f(-0.603768, 0.961718),vec2f(-0.657083, 0.890933),vec3f(0.857,0.780,0.682),vec3f(0.508,0.445,0.794));
        sdLine(point, color, d, vec2f(-0.589831, 0.999997),vec2f(-0.603747, 0.961733),vec3f(0.857,0.780,0.682),vec3f(0.527,0.467,0.816));
        sdLine(point, color, d, vec2f(-0.902586, 0.451620),vec2f(-0.836020, 0.550739),vec3f(0.569,0.553,0.867),vec3f(0.618,0.624,0.949));
        sdLine(point, color, d, vec2f(-0.958157, 0.411244),vec2f(-0.902331, 0.451519),vec3f(0.506,0.494,0.859),vec3f(0.637,0.647,0.957));
        sdLine(point, color, d, vec2f(-0.988291, 0.406328),vec2f(-0.958317, 0.411469),vec3f(0.465,0.453,0.851),vec3f(0.649,0.663,0.961));
        sdLine(point, color, d, vec2f(-0.790561, 0.737649),vec2f(-0.792976, 0.808576),vec3f(0.722,0.694,0.898),vec3f(0.859,0.788,0.661));
        sdLine(point, color, d, vec2f(-0.835941, 0.550741),vec2f(-0.790568, 0.737610),vec3f(0.708,0.680,0.892),vec3f(0.849,0.778,0.680));
        sdLine(point, color, d,  vec2f( 0.124987,-0.012010),vec2f( 0.128540,-0.077336),vec3f(0.153,0.094,0.125),vec3f(0.980,0.929,0.725),vec3f(0.082,0.039,0.075),vec3f(0.773,0.251,0.153));
        return d;
    }

}

// -----------------------------------------------------------------------------
// TESTS 
// -----------------------------------------------------------------------------

namespace yocto::extension{

    // -----------------------------------------------------------------------------
    // FAST WINDING NUMBER
    // -----------------------------------------------------------------------------

    void fwn_clean_mesh(std::string filename, int num_points){
        // Obtain positions and triangles from shape
        std::vector<vec4i> quads;
        std::vector<vec3f> positions;
        auto scale = 10.f;
        shape::make_monkey(quads, positions, scale);
        auto triangles = shape::quads_to_triangles(quads);

        // Create the Fast Winding Number tree
        auto fwn = FWN_tree(positions, triangles);

        // Create square box around the mesh
        rng_state rng;
        vec3f center;
        float radius;
        auto points = sample_box(positions, center, radius, num_points, rng);

        // Check if points are in the volume

        // Exact solution
        std::cout<<"Exact Solution\n";
        std::vector<vec3f> new_points0;
        auto tstart = std::chrono::high_resolution_clock::now();
        points_in_volume(new_points0, points, positions, triangles);

        // Check time of execution
        auto tend = std::chrono::high_resolution_clock::now();
        auto exe_time = std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count();
        std::cout<<"Solution found. Time : "<<exe_time*1e-3<<"sec \n";
        std::cout<<"Number of points : "<<new_points0.size()<<std::endl;

        // -------------------------
        // First order approximation
        // -------------------------

        std::cout<<"First order approximation\n";
        std::vector<vec3f> new_points1;
        tstart = std::chrono::high_resolution_clock::now();
        fwn.fast_points_in_volume(new_points1, points, 1);

        // Check time of execution
        tend = std::chrono::high_resolution_clock::now();
        exe_time = std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count();
        std::cout<<"Solution found. Time : "<<exe_time*1e-3<<"sec \n";
        std::cout<<"Number of points : "<<new_points1.size()<<std::endl;

        // Check correctness of the result
        auto sum1 = 0.f;
        for(auto point : new_points1)
            sum1 += in_volume(point, triangles, positions) ? 1 : 0;
        std::cout<<"Average Accuracy : "<<sum1/new_points1.size()<<std::endl;

        // --------------------------
        // Second order approximation
        // --------------------------
        std::cout<<"Second order approximation\n";
        std::vector<vec3f> new_points2;
        tstart = std::chrono::high_resolution_clock::now();
        fwn.fast_points_in_volume(new_points2, points, 2);

        // Check time of execution
        tend = std::chrono::high_resolution_clock::now();
        exe_time = std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count();
        std::cout<<"Solution found. Time : "<<exe_time*1e-3<<"sec \n";
        std::cout<<"Number of points : "<<new_points2.size()<<std::endl;

        // Check correctness of the result
        auto sum2 = 0.f;
        for(auto point : new_points2)
            sum2 += in_volume(point, triangles, positions) ? 1 : 0;
        std::cout<<"Average Accuracy : "<<sum2/new_points2.size()<<std::endl;

        // -------------------------
        // Third order approximation
        // -------------------------
        std::cout<<"Third order approximation\n";
        std::vector<vec3f> new_points3;
        tstart = std::chrono::high_resolution_clock::now();
        fwn.fast_points_in_volume(new_points3, points, 3);

        // Check time of execution
        tend = std::chrono::high_resolution_clock::now();
        exe_time = std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count();
        std::cout<<"Solution found. Time : "<<exe_time*1e-3<<"sec \n";
        std::cout<<"Number of points : "<<new_points3.size()<<std::endl;

        // Check correctness of the result
        auto sum3 = 0.f;
        for(auto point : new_points3)
            sum3 += in_volume(point, triangles, positions) ? 1 : 0;
        std::cout<<"Average Accuracy : "<<sum3/new_points3.size()<<std::endl;
    }

    // DA TESTARE
    void closest_point_clean_mesh (std::string filename, int num_points){
        // Obtain positions and triangles from shape
        std::vector<vec4i> quads;
        std::vector<vec3f> positions;
        std::vector<vec3i> triangles;
        auto scale = 10.f;

        if(filename == "monkey"){
            shape::make_monkey(quads, positions, scale);
            auto triangles = shape::quads_to_triangles(quads);
        }

        // Create the Fast Winding Number tree
        auto fwn = FWN_tree(positions, triangles);

        // Create square box around the mesh
        rng_state rng;
        vec3f center;
        float radius;
        auto points = sample_box(positions, center, radius, num_points, rng);

        // Compute closest triangle for each point 
        // Costume made function 
        float R;
        vec3i triangle;
        vec2f uv;
        vec3f sdf;
        auto tstart = std::chrono::high_resolution_clock::now();
        for (auto point: points){
            R = abs(fwn.closest_triangle(point));
        }

        // Check time of execution
        auto tend = std::chrono::high_resolution_clock::now();
        auto exe_time = std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count();
        std::cout<<"Solution found. Time : "<<exe_time*1e-3<<"sec \n";

        // Yocto overlap triangle function
        bvh_tree bvh;
        auto triangles_radius = std::vector<float> (positions.size(), 1e-5f);
        shape::make_triangles_bvh(bvh, triangles, positions, triangles_radius);
        tstart = std::chrono::high_resolution_clock::now();
        for (auto point: points){
            shape::overlap_triangles_bvh(bvh, triangles, positions, triangles_radius, point, radius, true);
        }    
        // Check time of execution
        tend = std::chrono::high_resolution_clock::now();
        exe_time = std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count();
        std::cout<<"Solution found. Time : "<<exe_time*1e-3<<"sec \n";    
    }

    void laplace_implicit(std::string imfilename, ptr::state* state, bool poisson, ptr::camera* camera, wos_params* params){
        // smooth blend of two tori
        std::function<float(vec3f)> SDF_torus =[](vec3f point)->float{
            auto R = 0.125;
            auto r = 0.06;
            vec3f center = {0, 0, -1};
            vec3f center1 = center - R*vec3f(1, 0, 0);
            vec3f center2 = center + R*vec3f(1, 0, 0);
            auto torus1 = sdTorus(vec3f(point.x, point.z, point.y), vec3f(center1.x, center1.z, center1.y), R, r);
            auto torus2 = sdTorus(point, center2, R, r);
            auto k = 0.1;
            return opSmoothUnion(torus1, torus2, k);
        };
        // yellow and blue stripes bounding function
        std::function<float(vec3f)> g = [](vec3f x){return fractional((x.y+0.008)*18.75 + 0.25)>0.5? 1.f : 0.f;};
        std::function<float(vec3f)> f;
        if(poisson)
            f = [params](vec3f x){return 1000*max(-sdSphere(x, {0, 0, -1}, 0.2), 0.f);};
        else if(params->importance.sample_f == "point"s) f = [params](vec3f x) {if (x == params->z) return params->cz; else return NAN;};
        else
            f = [](vec3f x){return NAN;};
        std::function<float(vec3f)> h = [](vec3f x){return 0.f;};
        std::function<float(ray3f)> sdfPlane = [params](ray3f ray)->float{
            return sphereTrace(ray, [params](vec3f point){return sdPlane(point, vec3f(0, 0, 1), -1);}, params->epsilon);};
        std::function<float(vec3f)> sdplane = [](vec3f point)->float{return sdPlane(point, vec3f(0, 0, 1), -1);};
        std::function<float(ray3f)> sdTorus_intersect = [params, SDF_torus](ray3f ray)->float{
            return sphereTrace(ray, SDF_torus, params->epsilon);};
        std::function<float(ray3f)> sdf = [sdplane, SDF_torus, params](ray3f ray){
            return sphereTrace(ray, [SDF_torus, sdplane](vec3f point){return opIntersection(SDF_torus(point), sdplane(point));}, params->epsilon);};
        std::function<float(vec3f)> bound_sdf = SDF_torus;
        // render
        cli::print_progress("render image", 0, params->num_walks);
        for(auto sample = 0; sample < params->num_walks; sample ++) {
            cli::print_progress("render image", sample, params->num_walks);

            solve_wos_samples(state, sdf, bound_sdf, g, f, h, camera, params);}
            cli::print_progress("render image", params->num_walks, params->num_walks);

            // save image
            cli::print_progress("save image", 0, 1);
            auto ioerror = ""s;
            if (!img::save_image(imfilename, state->render, ioerror)) cli::print_fatal(ioerror);
            cli::print_progress("save image", 1, 1);
        }

    void wos_diffusion_curves(std::string imfilename, ptr::state* state, wos_params* params){
        // boundary curves
        std::function<float(vec2f)> bound_sdf = [](vec2f point){vec3f color; return sdLine(point, color);};
        // boundary function
        std::function<vec3f(vec2f)> g = [](vec2f point){vec3f color; float d = sdLine(point, color); return color;};
        std::function<vec3f(vec2f)> f = [](vec2f point){return vec3f(NAN);};
        std::function<vec3f(vec2f)> h = [](vec2f point){return vec3f{0, 0, 0};};
        // render
        cli::print_progress("render image", 0, params->num_walks);
        for(auto sample = 0; sample < params->num_walks; sample ++) {
            cli::print_progress("render image", sample, params->num_walks);

            solve_wos_samples(state, bound_sdf, g, f, h, params);}
        cli::print_progress("render image", params->num_walks, params->num_walks);

        // save image
        cli::print_progress("save image", 0, 1);
        auto ioerror = ""s;
        if (!img::save_image(imfilename, state->render, ioerror)) cli::print_fatal(ioerror);
        cli::print_progress("save image", 1, 1);

    }

    void laplace_mesh_boolean(std::string imfilename, ptr::scene* scene, ptr::state* state, ptr::camera* camera, wos_params* params){
        std::function<float(ray3f)> bound_sdf = [](ray3f ray){return math::flt_min;};
        std::function<float(vec3f)> sdf = [](vec3f point){return math::flt_min;};
        std::vector<std::function<bool(vec3f)>> vinside = {[](vec3f){return true;}};
        auto fwn_trees = std::vector<FWN_tree> (scene->objects.size());
        auto id = 0;

        // from shapes to bvh
        for(auto object : scene->objects){
            auto shape = object->shape;
            fwn_trees[id] = FWN_tree(shape->positions, shape->triangles);
            std::function<float(vec3f)> sdf_mesh = [&fwn_trees, id](vec3f point)->float{
                vec3i triangle; vec2f uv; 
                auto fwn_tree = fwn_trees[id];
                return fwn_tree.closest_triangle(point);};
            std::function<float(ray3f)> bound_mesh = [&fwn_trees, id](ray3f ray)->float{return fwn_trees[id].intersect(ray);};
            sdf = [sdf, sdf_mesh](vec3f point){return opSmoothIntersection(sdf(point), sdf_mesh(point), 0.1);};
            bound_sdf = [bound_sdf, bound_mesh](ray3f ray){return opSmoothIntersection(bound_sdf(ray), bound_mesh(ray), 0.1);};
            vinside.push_back([&fwn_trees, id](vec3f point)->bool{return fwn_trees[id].fast_point_in_volume(point) >= 0.5;});
            id += 1;
        }
        // functions 
        std::function<float(vec3f)> g = [](vec3f x){return fractional((x.x+0.004)*18.75 + 0.25)>0.5? 1.f : 0.f;};
        std::function<float(vec3f)> f = [](vec3f x){return NAN;};
        std::function<float(vec3f)> h = [](vec3f x){return 0.f;};
        // external boundaries
        // implicit surfaces
        auto plane = [](vec3f point ){return sdPlane(point, vec3f(0, 0, 1), 0.03);};
        auto sphere = [](vec3f point){return sdSphere(point, vec3f(0, 0.04, 0.03), 0.03);};
        std::function<float(ray3f)> sdf_intersect = [sphere, plane, params](ray3f ray)->float{
            return sphereTrace(ray,[plane, sphere](vec3f point) {return opSubtraction(sphere(point), plane(point));}, params->epsilon);};
        std::function<float(vec3f)> sdf_implicit = [sphere, plane](vec3f point)->float{return opSubtraction(sphere(point), plane(point));};
        // implicit surface and mesh        
        std::function<float(ray3f)> sdf_tot = [bound_sdf, sdf_intersect, vinside, sdf_implicit, params](ray3f ray){
            auto d = bound_sdf(ray);
            if (!(d<math::flt_max)) return d;
            if(sdf_implicit(ray.o + ray.d * d) < 0.f) return d;
            d = sdf_intersect(ray);
            auto point = ray.o + ray.d * d;
            if (std::all_of(vinside.begin(), vinside.end(), [point](std::function<bool(yocto::math::vec3f)> v){return v(point);}))
                return d;
            else return math::flt_max;};
        // internal boundaries
        sdf = [sdf, sphere](vec3f point){return opSubtraction(sphere(point), sdf(point));};

        // render
        auto ioerror = ""s;
        cli::print_progress("render image", 0, params->num_walks);
        for(auto sample = 0; sample < params->num_walks; sample ++) {
            cli::print_progress("render image", sample, params->num_walks);

            solve_wos_samples(state, sdf_tot, sdf, g, f, h, camera, params);
            if (params->save_batch > 0) {
                auto output = imfilename + "." + std::to_string(sample) + ".png";
                if (!img::save_image(output, state->render, ioerror)) cli::print_fatal(ioerror);
            }}
            cli::print_progress("render image", params->num_walks, params->num_walks);

            // save image
            cli::print_progress("save image", 0, 1);
            
            if (!img::save_image(imfilename, state->render, ioerror)) cli::print_fatal(ioerror);
            cli::print_progress("save image", 1, 1);
    }
}