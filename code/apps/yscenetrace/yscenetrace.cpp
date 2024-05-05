//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
#include <iostream>

#include <yocto/yocto_commonio.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_shape.h>
#include <yocto/yocto_sceneio.h>
#include <yocto_pathtrace/yocto_pathtrace.h>
using namespace yocto::math;
namespace ptr = yocto::pathtrace;
namespace cli = yocto::commonio;
namespace sio = yocto::sceneio;
namespace shp = yocto::shape;
namespace ext = yocto::extension;

#include <map>
#include <memory>
using namespace std::string_literals;

#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;

// construct a scene from io
void init_scene(ptr::scene* scene, sio::model* ioscene, ptr::camera*& camera,
    sio::camera* iocamera, sio::progress_callback progress_cb = {}) {
  // handle progress
  auto progress = vec2i{
      0, (int)ioscene->cameras.size() + (int)ioscene->environments.size() +
             (int)ioscene->materials.size() + (int)ioscene->textures.size() +
             (int)ioscene->shapes.size() + (int)ioscene->subdivs.size() +
             (int)ioscene->objects.size()};

  auto camera_map     = std::unordered_map<sio::camera*, ptr::camera*>{};
  camera_map[nullptr] = nullptr;
  for (auto iocamera : ioscene->cameras) {
    if (progress_cb) progress_cb("convert camera", progress.x++, progress.y);
    auto camera = add_camera(scene);
    set_frame(camera, iocamera->frame);
    set_lens(camera, iocamera->lens, iocamera->aspect, iocamera->film);
    set_focus(camera, iocamera->aperture, iocamera->focus);
    camera_map[iocamera] = camera;
  }

  auto texture_map     = std::unordered_map<sio::texture*, ptr::texture*>{};
  texture_map[nullptr] = nullptr;
  for (auto iotexture : ioscene->textures) {
    if (progress_cb) progress_cb("convert texture", progress.x++, progress.y);
    auto texture = add_texture(scene);
    if (!iotexture->colorf.empty()) {
      set_texture(texture, iotexture->colorf);
    } else if (!iotexture->colorb.empty()) {
      set_texture(texture, iotexture->colorb);
    } else if (!iotexture->scalarf.empty()) {
      set_texture(texture, iotexture->scalarf);
    } else if (!iotexture->scalarb.empty()) {
      set_texture(texture, iotexture->scalarb);
    }
    texture_map[iotexture] = texture;
  }

  auto material_map = std::unordered_map<sio::material*, ptr::material*>{};
  material_map[nullptr] = nullptr;
  for (auto iomaterial : ioscene->materials) {
    if (progress_cb) progress_cb("convert material", progress.x++, progress.y);
    auto material = add_material(scene);
    set_emission(material, iomaterial->emission,
        texture_map.at(iomaterial->emission_tex));
    set_color(
        material, iomaterial->color, texture_map.at(iomaterial->color_tex));
    set_specular(material, iomaterial->specular,
        texture_map.at(iomaterial->specular_tex));
    set_ior(material, iomaterial->ior);
    set_metallic(material, iomaterial->metallic,
        texture_map.at(iomaterial->metallic_tex));
    set_transmission(material, iomaterial->transmission, iomaterial->thin,
        iomaterial->trdepth, texture_map.at(iomaterial->transmission_tex));
    set_roughness(material, iomaterial->roughness,
        texture_map.at(iomaterial->roughness_tex));
    set_opacity(
        material, iomaterial->opacity, texture_map.at(iomaterial->opacity_tex));
    set_thin(material, iomaterial->thin);
    set_scattering(material, iomaterial->scattering, iomaterial->scanisotropy,
        texture_map.at(iomaterial->scattering_tex));
    set_normalmap(material, texture_map.at(iomaterial->normal_tex));

    // Add hair properties
    set_hair(material, iomaterial->eta, iomaterial->alpha, iomaterial->sigma_a, iomaterial->beta_m, iomaterial->beta_n, 
      iomaterial->eumelanin, iomaterial->pheomelanin);

    material_map[iomaterial] = material;


  }

  auto subdiv_map = std::unordered_map<sio::subdiv*, ptr::shape*>{};
  for (auto iosubdiv : ioscene->subdivs) {
    if (progress_cb)
      progress_cb("convert subdiv", progress.x++, progress.y);
    auto subdiv = add_shape(scene);
    set_subdiv_quadspos(subdiv, iosubdiv->quadspos);
    set_subdiv_quadstexcoord(subdiv, iosubdiv->quadstexcoord);
    set_subdiv_positions(subdiv, iosubdiv->positions);
    set_subdiv_texcoords(subdiv, iosubdiv->texcoords);
    subdiv_map[iosubdiv] = subdiv;
  }

  auto shape_map     = std::unordered_map<sio::shape*, ptr::shape*>{};
  shape_map[nullptr] = nullptr;
  for (auto ioshape : ioscene->shapes) {
    if (progress_cb) progress_cb("convert shape", progress.x++, progress.y);
    auto shape = add_shape(scene);
    set_points(shape, ioshape->points);
    set_lines(shape, ioshape->lines);
    set_triangles(shape, ioshape->triangles);
    if(!ioshape->quads.empty())
      set_triangles(shape, shp::quads_to_triangles(ioshape->quads));
    set_positions(shape, ioshape->positions);
    set_normals(shape, ioshape->normals);
    set_texcoords(shape, ioshape->texcoords);
    set_radius(shape, ioshape->radius);
    shape_map[ioshape] = shape;
  }

  for (auto ioobject : ioscene->objects) {
    if (progress_cb) progress_cb("convert object", progress.x++, progress.y);
    if(ioobject->instance) {
      for(auto frame : ioobject->instance->frames) {
        auto object = add_object(scene);
        set_frame(object, frame * ioobject->frame);
        if (ioobject->shape) set_shape(object, shape_map.at(ioobject->shape));
        if (ioobject->subdiv) {
          set_shape(object, subdiv_map.at(ioobject->subdiv));
          set_subdiv_subdivision(subdiv_map.at(ioobject->subdiv),
              ioobject->material->subdivisions, true);
          set_subdiv_displacement(subdiv_map.at(ioobject->subdiv),
              ioobject->material->displacement,
              texture_map.at(ioobject->material->displacement_tex));
        }
        set_material(object, material_map.at(ioobject->material));
      }
    } else {
      auto object = add_object(scene);
      set_frame(object, ioobject->frame);
      if (ioobject->shape) set_shape(object, shape_map.at(ioobject->shape));
      if (ioobject->subdiv) {
        set_shape(object, subdiv_map.at(ioobject->subdiv));
        set_subdiv_subdivision(subdiv_map.at(ioobject->subdiv),
            ioobject->material->subdivisions, true);
        set_subdiv_displacement(subdiv_map.at(ioobject->subdiv),
            ioobject->material->displacement,
            texture_map.at(ioobject->material->displacement_tex));
      }
      set_material(object, material_map.at(ioobject->material));
    }
  }

  for (auto ioenvironment : ioscene->environments) {
    if (progress_cb)
      progress_cb("convert environment", progress.x++, progress.y);
    auto environment = add_environment(scene);
    set_frame(environment, ioenvironment->frame);
    set_emission(environment, ioenvironment->emission,
        texture_map.at(ioenvironment->emission_tex));
  }

  // done
  if (progress_cb) progress_cb("convert done", progress.x++, progress.y);

  // get camera
  camera = camera_map.at(iocamera);
}

int main(int argc, const char* argv[]) {
  // options
  auto params      = ptr::trace_params{};
  auto save_batch  = false;
  auto camera_name = ""s;
  auto imfilename  = "out.hdr"s;
  auto filename    = "scene.json"s;
  auto dimension   = 3;
  auto montecarlo  = false;
  int num_points   = 1e4;
  auto scale       = 1.f;
  auto bounding    = "fraction-norm"s;
  auto skeleton    = false;
  auto screening   = 0.f;
  auto solver      = "laplace"s;
  auto slice       = 1.f;

  // parse command line
  auto cli = cli::make_cli("yscntrace", "Offline path tracing");
  add_option(cli, "--camera", camera_name, "Camera name.");
  add_option(cli, "--resolution,-r", params.resolution, "Image resolution.");
  add_option(cli, "--samples,-s", params.samples, "Number of samples.");
  add_option(cli, "--shader,-t", params.shader, "Shader type.", ptr::shader_names);
  add_option(cli, "--bounces,-b", params.bounces, "Maximum number of bounces.");
  add_option(cli, "--clamp", params.clamp, "Final pixel clamping.");
  add_option(cli, "--noparallel", params.noparallel, "Don't use parallel for loops. ");
  add_option(cli, "--save-batch/--no-save-batch", save_batch, "Save images progressively");
  add_option(cli, "--output-image,-o", imfilename, "Image filename");

  // Hair
  // add_option(cli, "--hair/--no-hair", params.has_hair, "Use or not the hair shader.");

  // monte carlo
  // add_option(cli, "--num-points", num_points, "Number of points sampled for montecarlo");
  // add_option(cli, "--montecarlo/--no-montecarlo", montecarlo, "Show only the edges of the mesh");
  // add_option(cli, "--scale", scale, "Scaling value");
  // add_option(cli, "--num_walks", params.num_walks, "Number of walks for montecarlo");
  // add_option(cli, "--max_steps", params.max_steps, "Max number of steps during each montecarlo walk.");
  // add_option(cli, "--bound", bounding, "Bound function : [fraction-norm, norm, stripes]");
  // add_option(cli, "--skeleton", skeleton, "Visualize objects skeleton");
  // add_option(cli, "--screen", screening, "Screening value for screened Poisson solver");
  //add_option(cli, "--solver", solver, "Type of montecarlo solver [laplace, poisson, nested-biharmonic]");
  // add_option(cli, "--importance", ext::importance_sampling, "Use importance sampling");
  // add_option(cli, "--colormap, -c", ext::COLORMAP_ID, "Index of chosen colormap {1: 'yellow-purple', 2: 'red-blue'.");
  // add_option(cli, "--slice", slice, "float between [0,1] representing how much you want to slice from the meshn");
  
  add_option(cli, "scene", filename, "Scene filename", true);

  parse_cli(cli, argc, argv);

  // scene loading
  auto ioscene_guard = std::make_unique<sio::model>();
  auto ioscene       = ioscene_guard.get();
  auto ioerror       = ""s;
  if (!load_scene(filename, ioscene, ioerror, cli::print_progress))
    cli::print_fatal(ioerror);

  // get camera
  auto iocamera = get_camera(ioscene, camera_name);

  // convert scene
  auto scene_guard = std::make_unique<ptr::scene>();
  auto scene       = scene_guard.get();
  auto camera      = (ptr::camera*)nullptr;
  init_scene(scene, ioscene, camera, iocamera, cli::print_progress);

  // cleanup
  if (ioscene_guard) ioscene_guard.reset();

  /*if (montecarlo){
    auto size = scene->objects.size();
    for (auto i = 0; i<size; i++){
    //for (auto& object: scene->objects){
      auto scale_frame = scaling_frame({scale, scale, scale});
      auto& object = scene->objects[i];
      for (auto& p : object -> shape -> positions) p = transform_point(scale_frame, p);
      for (auto& n : object -> shape -> normals) n = transform_normal(scale_frame, n, false);

      // SAMPLE POINTS INSIDE THE OBJECT
      auto quads_points = std::vector<std::vector<vec4i>> {};
      auto positions_points = std::vector<std::vector<vec3f>>{};
      auto normals_points = std::vector<std::vector<vec3f>>{};
      auto texcoords_points = std::vector<std::vector<vec2f>>{};

      //auto vertices = object -> shape -> positions;
      std::vector<vec3f> U = {};
      ext::montecarlo_volume(object -> shape -> positions, object -> shape -> triangles, object -> shape -> normals, num_points,
                                quads_points, positions_points, normals_points, texcoords_points, screening, U, bounding, solver,
                                slice, params.num_walks, params.max_steps);

      // CREATE OBJECT SKELETON
      auto edges = shp::get_edges(object->shape->triangles);

      std::vector<vec4i> quads; 
      std::vector<vec3f> qpositions;
      std::vector<vec3f> qnormals;
      std::vector<vec2f> qtexcoords;

      if (skeleton){
        // VISUALIZE SKELETON
        ext::lines_to_cylinders(quads, qpositions, qnormals, qtexcoords, edges, object -> shape -> positions, 4, 0.0001f);

        object -> shape -> triangles = shp::quads_to_triangles(quads);
        object -> shape -> positions = qpositions;
        object -> shape -> normals = qnormals;
        object -> shape -> texcoords = qtexcoords;
        object -> material -> color = {0,0,0};
      }
      //else std::remove(scene->objects.begin(), scene->objects.end(), object);
      else scene->objects.erase(scene->objects.begin()+i);


      // VISUALIZE POINTS

      for (auto i = 0; i<quads_points.size();i++){
        auto point_obj = add_object(scene);
        auto point_shape = add_shape(scene);
        point_shape->triangles = shp::quads_to_triangles(quads_points[i]);
        point_shape->positions = positions_points[i];
        point_shape->normals = normals_points[i];
        point_shape->texcoords = texcoords_points[i];
        point_obj->shape = point_shape;
        auto point_material = add_material(scene);
        point_material->color = U[i];
        point_obj->material = point_material;
      }
    }
  }*/

  // init subdivs
  init_subdivs(scene, params, cli::print_progress);

  // build bvh
  init_bvh(scene, params, cli::print_progress);

  // build lights
  init_lights(scene, params, cli::print_progress);

  // init state
  auto state_guard = std::make_unique<ptr::state>();
  auto state = state_guard.get();
  init_state(state, scene, camera, params);

  // render
  cli::print_progress("render image", 0, params.samples);
  for(auto sample = 0; sample < params.samples; sample ++) {
    cli::print_progress("render image", sample, params.samples);
    trace_samples(state, scene, camera, params);
    if(save_batch) {
        auto ext = "-s" + std::to_string(sample) +
                   fs::path(imfilename).extension().string();
        auto outfilename = fs::path(imfilename).replace_extension(ext).string();
        auto ioerror     = ""s;
        cli::print_progress("save image", sample, params.samples);
        if (!save_image(outfilename, state->render, ioerror))
          cli::print_fatal(ioerror);
    }
  }
  cli::print_progress("render image", params.samples, params.samples);

  // save image
  cli::print_progress("save image", 0, 1);
  if (!save_image(imfilename, state->render, ioerror)) cli::print_fatal(ioerror);
  cli::print_progress("save image", 1, 1);

  // done
  return 0;
}
