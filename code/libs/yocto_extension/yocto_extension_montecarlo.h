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

#ifndef _YOCTO_EXTENSION_MONTECARLO_H_
#define _YOCTO_EXTENSION_MONTECARLO_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <yocto/yocto_math.h>
#include <yocto/yocto_shape.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_sceneio.h>
#include "yocto_pathtrace/yocto_pathtrace.h"
#include "yocto/yocto_commonio.h"

#include <memory>
using namespace std::string_literals;

// -----------------------------------------------------------------------------
// ALIASES
// -----------------------------------------------------------------------------

namespace yocto::extension{
        // Namespace aliases
    namespace cli = yocto::commonio;
    namespace img = yocto::image;
    namespace ptr = yocto::pathtrace;
    namespace ext = yocto::extension;
    namespace sio = yocto::sceneio;

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
    using math::zero2f;
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
// MATH FUNCTIONS
// -----------------------------------------------------------------------------

namespace yocto::extension{

    // Extension to yocto::math for 6x6 matrices needed to compute the taylor series third order approximation for Fast Winding Numbers
    struct mat9x3f{
        mat3f m1 = mat3f(zero3f, zero3f, zero3f);
        mat3f m2 = mat3f(zero3f, zero3f, zero3f);
        mat3f m3 = mat3f(zero3f, zero3f, zero3f);

        mat9x3f();
        mat9x3f(mat3f m1, mat3f m2, mat3f m3);

    };

    inline float inner(mat3f a, mat3f b){return dot(a.x, b.x) + dot(a.y, b.y) + dot(a.z, b.z);}
    inline mat3f outer(vec3f u, vec3f v){
        return mat3f({u.x*v.x, u.x*v.y, u.x*v.z}, {u.y*v.x, u.y*v.y, u.y*v.z}, {u.z*v.x, u.z*v.y, u.z*v.z});}
    inline float inner(mat9x3f a, mat9x3f b){return inner(a.m1, b.m1) + inner(a.m2, b.m2) + inner(a.m3, b.m3);}   
    inline mat9x3f outer(vec3f u, mat3f V){return mat9x3f(V*u.x, V*u.y, V*u.z);}
    inline mat9x3f outer(mat3f V, vec3f u){return mat9x3f(V*u.x, V*u.y, V*u.z);}
    
    inline mat3f operator-(const mat3f& a){return{-a.x, -a.y, -a.z};}

    inline mat9x3f::mat9x3f() {}
    inline mat9x3f::mat9x3f(mat3f m1, mat3f m2, mat3f m3) : m1{m1}, m2{m2}, m3{m3} {} 
   
    inline mat9x3f operator+(const mat9x3f& a, const mat9x3f& b){return {a.m1 + b.m1, a.m2+b.m2, a.m3+b.m3};}
    inline mat9x3f operator-(const mat9x3f& a, const mat9x3f& b){return {a.m1 - b.m1, a.m2 - b.m2, a.m3 - b.m3};}
    inline mat9x3f operator-(const mat9x3f& a){return {-a.m1, -a.m2, -a.m3};}
    inline mat9x3f operator*(const float a, const mat9x3f& b){return {b.m1*a, b.m2*a, b.m3*a};}
    inline mat9x3f operator*(const mat9x3f& a, const float b){return b*a;}
    inline mat9x3f operator/(const mat9x3f& a, const float b){return {a.m1*1.f/b, a.m2/b, a.m3/b};}
    inline mat9x3f& operator+=(mat9x3f& a, const mat9x3f& b){return a = a + b;};

}

// -----------------------------------------------------------------------------
// HIGH LEVEL API - FAST WINDING NUMBERS
// -----------------------------------------------------------------------------

namespace yocto::extension{

    // Structure for Fast Winding Numbers
    class FWN_tree{
        private:          
            struct fwn_node;
            bvh_tree              bvh;
            std::vector<fwn_node> nodes      = {};
            std::vector<vec3f>    positions  = {};
            std::vector<vec3i>    triangles  = {};

        public:
            // Default Constructor
            FWN_tree(){};
            // Constructor
            FWN_tree(std::vector<vec3f> positions, std::vector<vec3i> triangles);

            // Methods
            void fast_points_in_volume(std::vector<vec3f> &sample_points, std::vector<vec3f> points, int approximation = 1, float beta = 2.);
            float fast_point_in_volume(vec3f point, float beta = 2., int approximation = 1);
            float closest_triangle(vec3f point);
            float intersect(ray3f &ray);

    };
    struct FWN_tree::fwn_node{         
        std::vector<fwn_node*>  children;  
        vec3f                   p;
        int id;

        fwn_node* parent =  nullptr                  ;
        float   r        =  0.f                      ;
        vec3f   t1       =  zero3f                   ;
        mat3f   t2       =  {zero3f, zero3f, zero3f} ; 
        mat9x3f t3       =  {{zero3f, zero3f, zero3f}, 
                             {zero3f, zero3f, zero3f}, 
                             {zero3f, zero3f, zero3f}};

        // Constructor
        fwn_node(){}
     
    };
}

// -----------------------------------------------------------------------------
// HIGH LEVEL API - MONTECARLO GEOMETRY PROCESSING
// -----------------------------------------------------------------------------

namespace yocto::extension{

    struct importance_sampling{
        bool        sample_y   = false;
        std::string sample_f   = ""s;
    };

    struct wos_params{

        float       clamp       = 100;
        int         save_batch  = -1;
        bool        noparallel  = false;      
        float       q           = 0  ;
        float       b           = 0  ;
        float       c           = -1.f; // screening parameter
        vec3f       z           = zero3f; // point source
        float       cz          = 0; // point source coefficient
        float       tollerance  = 1e-5;
        float       epsilon     = 1e-5;
        int         num_walks   = 128;
        float       max_steps   = 128;
        vec2f       bounds      = {0, 1};
        
        int         colormap_id = 1;
        // std::vector<vec3f> colormap = std::vector<vec3f>(256);
        
        bool variance          = false;
        importance_sampling importance;

        // increment used to compute the normal to the sdf 
        float       h          = 1e-3;
    };

    struct colors{
        std::vector<std::vector<int>> colormap = {};
        std::vector<vec3f>            colors   = {};
    };

    // GLOBAL FUNCTIONS
    // Compute solid angle projection of a triangle over a point
    // Used to check if a point p is inside a volume
    float solid_angle(vec3f vi, vec3f vj, vec3f vk, vec3f p);
    // sample n points inside a box given its minimum and maximum vertices
    std::vector<vec3f> sample_box(std::vector<vec3f> vertices, vec3f& center, float& radius, int num = 10000, math::rng_state rng = {});
    // Returns the minimum and maximum vertices of a mesh
    std::vector<vec3f> get_min_max(std::vector<vec3f> points);
    // Check if a point p is inside a volume
    // Generalized Winding Number, correct for clean meshes
    bool in_volume(vec3f point, std::vector<vec3i> triangles, std::vector<vec3f> points);
    template <typename Func>
    vec3f importance_sample(float &f_, Func f, float R, vec3f x, wos_params params, rng_state rng);
    // Check if a ray intersect an implicit surface
    template <typename SDF>
    bool sphereTrace(ray3f ray, SDF &sdf, float &dist, float epsilon);
    template <typename SDF>
    vec3f sdNormal(SDF sdf, vec3f x, float h);

    vec3f colormap(float x, wos_params params);

    // Compute ray point
    inline vec3f eval_ray(ray3f ray, float t){return ray.o + t*ray.d;}

    // MONTECARLO FUNCTIONS
    template <typename SDF_intersect = std::function<float(ray3f)>, typename SDF = std::function<float(vec3f)>, typename Func>
    void solve_wos_samples(ptr::state* state, SDF_intersect sdf_intersect, SDF sdf, Func g, Func f, Func h, 
        const ptr::camera* camera, wos_params* params);
    template <typename SDF_intersect = std::function<float(ray3f)>, typename SDF = std::function<float(vec3f)>, typename Func>
    static vec4f solve_wos_sample(ptr::state* state,  SDF_intersect sdf_intersect, SDF sdf, Func g, Func f, Func h, const ptr::camera* camera, 
        const vec2i& ij, wos_params* params);

    // Estimate gradient of PDE in a point using Montecarlo Walk On Sphere
    template <typename SDF, typename Func>
    static vec3f wos_gradient(SDF sdf, vec3f x0, Func g, Func f, Func h, wos_params params, rng_state rng);

    // -------------------
    // POTENTIAL FUNCTIONS
    // -------------------
    // HARMONIC GREEN FUNCTION
    float Green_harmonic(vec3f x, vec3f y, float R);
    // DERIVATIVE GREEN FUNCTIONS
    vec3f Green_harmonic_gradient(vec3f x, vec3f y, float R);
    // YUKAWA POTENTIAL FOR SCREENED POISSON
    float Yukawa_potential(vec3f x, vec3f y, float R, float C);
    // DERIVATIVE YUKAWA POTENTIAL
    vec3f Yukawa_potential_gradient(vec3f x, vec3f y, float R, float c);

    template <typename Func>
    inline void parallel_for(const int size, Func&& func);
    template <typename Func>
    inline void parallel_for(const vec2i ij, Func&& func);


    struct scene{
        std::vector<vec2f> positions;
        std::vector<vec2i> lines;
        std::vector<vec3f> colors;
        std::vector<std::vector<int>>  color_code;
    };

}

// ---------
// TESTS
// ---------
namespace yocto::extension{

    void fwn_clean_mesh(std::string filename = "monkey", int num_points = 10000);
    void closest_point_clean_mesh(std::string filename = "monkey", int num_points = 10000);
    void laplace_implicit(std::string imfilename, ptr::state* state, bool poisson, ptr::camera* camera, wos_params* params);
    void wos_diffusion_curves(std::string imfilename, ptr::state* state, wos_params* params);
    void laplace_mesh_boolean(std::string imfilename, ptr::scene* scene, ptr::state* state, ptr::camera* camera, wos_params* params);
}

#endif


