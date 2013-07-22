//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <optix_world.h>
#include <optixu/optixu_math_namespace.h>
#include "OPTIX_COMMONSTRUCTS.h"
#include "OPTIX_HELPERS.h"
#include "OPTIX_RAY_STRUCTS.h"
using namespace optix;

rtTextureSampler<float, 3, cudaReadModeElementType> smoke_dencities_tex;
rtTextureSampler<float, 3, cudaReadModeElementType> photon_map_tex;

rtDeclareVariable(float3,       low_corner, , );
rtDeclareVariable(float3,       up_corner, , );
rtDeclareVariable(float,       scene_epsilon, , );
rtDeclareVariable(float,       step, , );
rtDeclareVariable(float,       density_multiplier, , );
rtDeclareVariable(float,       exp_multiplier, , );
rtDeclareVariable(float,       exp_multiplier_2, , );
rtDeclareVariable(float,       exp_tex_resolution, , );

rtDeclareVariable(rtObject,          top_object, , );

rtDeclareVariable(rtObject,          top_opaque_object, , );

rtDeclareVariable(PerRayData_radiance, prd, rtPayload, );
rtDeclareVariable(float, t_hit, rtIntersectionDistance, );
rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );

rtDeclareVariable(float, front_box_hit_point, attribute front_box_hit_point, );
rtDeclareVariable(float, back_box_hit_point, attribute back_box_hit_point, );

__device__ bool isInsideDomain(float3 point) {
    return !(point.x < low_corner.x || point.x > up_corner.x ||
        point.y < low_corner.y || point.y > up_corner.y ||
        point.z < low_corner.z || point.z > up_corner.z);
}

__device__ float3 operator*(float3 a, uint3 b) {
    return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}

__device__ float getIntensity(float3 point) {
    if (!isInsideDomain(point)) {
        return -1.f;
    }
    float3 tex_coord = (point - low_corner) / (up_corner - low_corner);
    return tex3D(smoke_dencities_tex, tex_coord.x, tex_coord.y, tex_coord.z);
}

RT_PROGRAM void any_hit_shadow() {
  // phongShadowed();
  rtIgnoreIntersection();
}

#if 0
RT_PROGRAM void closest_hit_radiance() {
  float3 hit_point = ray.origin + (t_hit + scene_epsilon) * ray.direction;

  PerRayData_radiance new_prd;
  new_prd.importance = prd.importance;
  new_prd.depth = prd.depth;
  optix::Ray next_ray = optix::make_Ray(hit_point, ray.direction, 0, scene_epsilon, RT_DEFAULT_MAX);
  rtTrace(top_opaque_object, next_ray, new_prd);

  float back_point = back_box_hit_point;//new_prd.distance > 0 ? (back_box_hit_point < new_prd.distance ? back_box_hit_point : new_prd.distance) : back_box_hit_point;
  /*
  if (new_prd.distance > 0) {
    back_box_hit_point = new_prd.distance;//back_box_hit_point < new_prd.distance ? back_box_hit_point : new_prd.distance;//fminf(back_box_hit_point, new_prd.distance);
  }
  */
  if (back_point - front_box_hit_point < scene_epsilon) {
      prd.result = new_prd.result;
      return;
  }

  float density = getIntensity(hit_point);
  float accumulative_density = 0;
  float accumulative_step = front_box_hit_point;
  if (density > 0) {
      do {
         accumulative_density += density * density_multiplier;
         accumulative_step += step;
         // hit_point += ray.direction * ;
         if (accumulative_density >= 1) {
            accumulative_density = 1;
            break;
         }
         density = getIntensity(hit_point + ray.direction * accumulative_step);
      } while (density >= 0.f && accumulative_step < front_box_hit_point + 2);

      // accumulative_density = back_point - front_box_hit_point;
      prd.result = make_float3(accumulative_density, accumulative_density, accumulative_density);// + new_prd.result;
  } else {
      prd.result = make_float3(1.f, 0.8f, 0.f);
  }
}

#endif

#if 0
RT_PROGRAM void closest_hit_radiance() {
  float3 hit_point = ray.origin + (t_hit + scene_epsilon) * ray.direction;

  if (!isInsideDomain(hit_point)) {
     prd.result = make_float3(1.f, 0.f, 0.f);
  }

  float density = getIntensity(hit_point);
  float accumulative_density;
  do {
     accumulative_density += density * density_multiplier;
     hit_point += ray.direction * step;
     if (accumulative_density >= 1) {
        accumulative_density = 1;
        break;
     }
     density = getIntensity(hit_point);
  } while (density >= 0.f);

  PerRayData_radiance new_prd;
  new_prd.importance = prd.importance * (1.f - accumulative_density);
  new_prd.depth = prd.depth;

  optix::Ray next_ray = optix::make_Ray(hit_point, ray.direction, 0, scene_epsilon, RT_DEFAULT_MAX);
  rtTrace(top_object, next_ray, new_prd);
  prd.result = make_float3(accumulative_density, accumulative_density, accumulative_density) + new_prd.result;
}
#endif

#if 0
RT_PROGRAM void closest_hit_radiance() {
  float3 hit_point = ray.origin + (t_hit + scene_epsilon) * ray.direction;


  float t = t_hit - front_box_hit_point;
  if (t > scene_epsilon || t < -scene_epsilon) {
    prd.result = make_float3(1.f, 0.f, 0.f);
    return;
  }

  PerRayData_radiance new_prd;
  new_prd.importance = prd.importance;
  new_prd.depth = prd.depth;

  optix::Ray next_ray = optix::make_Ray(ray.origin, ray.direction, 0, scene_epsilon, RT_DEFAULT_MAX);
  rtTrace(top_opaque_object, next_ray, new_prd);
  // prd.result = make_float3(accumulative_density, accumulative_density, accumulative_density) + new_prd.result;

  if (!isInsideDomain(hit_point)) {
     prd.result = new_prd.result;
     return;
  }

  float back_point = back_box_hit_point;
  if (new_prd.distance > 0) {
    back_point = fminf(back_box_hit_point, new_prd.distance);
  }

  float density = getIntensity(hit_point);
  float accumulative_density = 0;
  float accumulative_step = front_box_hit_point;
  do {
     accumulative_density += density * step;
     hit_point += ray.direction * step;
     if (accumulative_density >= 1) {
        accumulative_density = 1;
        break;
     }
     density = getIntensity(hit_point);
     accumulative_step += step;
  } while (density >= 0.f && accumulative_step < back_point);

  // PerRayData_radiance new_prd;
  /*new_prd.importance = prd.importance * (1.f - accumulative_density);
  new_prd.depth = prd.depth;

  optix::Ray next_ray = optix::make_Ray(hit_point, ray.direction, 0, scene_epsilon, RT_DEFAULT_MAX);
  rtTrace(top_object, next_ray, new_prd);
  */
  accumulative_density *= density_multiplier;
  prd.result = make_float3(accumulative_density, accumulative_density, accumulative_density) + new_prd.result;
}
#endif

RT_PROGRAM void closest_hit_radiance() {
  float3 hit_point = ray.origin + (t_hit + scene_epsilon) * ray.direction;

  float t = t_hit - front_box_hit_point;
  if (t > scene_epsilon || t < -scene_epsilon) {
    prd.result = make_float3(1.f, 0.f, 0.f);
    return;
  }

  PerRayData_radiance new_prd;
  new_prd.importance = prd.importance;
  new_prd.depth = prd.depth;

  optix::Ray next_ray = optix::make_Ray(ray.origin, ray.direction, 0, scene_epsilon, RT_DEFAULT_MAX);
  rtTrace(top_opaque_object, next_ray, new_prd);
  // prd.result = make_float3(accumulative_density, accumulative_density, accumulative_density) + new_prd.result;

  if (!isInsideDomain(hit_point)) {
     prd.result = new_prd.result;
     return;
  }

  float back_point = back_box_hit_point;
  if (new_prd.distance > 0) {
    back_point = fminf(back_box_hit_point, new_prd.distance);
  }

  float density = getIntensity(hit_point);
  float result_color_radiance = 0;
  float alpha = 1;

  // float accumulative_density = 0;
  float accumulative_step = front_box_hit_point;
  // float result_intensity = 0;
  float3 normalized_hit_point;

  do {
     normalized_hit_point = (hit_point - low_corner) / (up_corner - low_corner);
     result_color_radiance += alpha * density_multiplier * density * step * tex3D(photon_map_tex, normalized_hit_point.x, normalized_hit_point.y, normalized_hit_point.z);
     if (result_color_radiance > 1.0f) {
        result_color_radiance = 1.0f;
        break;
     }
     alpha *= (1 - density_multiplier * density * step);
     if (alpha < 0) {
        alpha = 0.f;
        break;
     }

     // result_intensity += tex3D(photon_map_tex, normalized_hit_point.x, normalized_hit_point.y, normalized_hit_point.z);// * tex1D(exp_tex, accumulative_density / exp_tex_resolution);

     // accumulative_density += density;
     hit_point += ray.direction * step;
     density = getIntensity(hit_point);
     accumulative_step += step;
  } while (density >= 0.f && accumulative_step < back_point);

  // result_intensity *= density_multiplier;
  // normalized_hit_point = (hit_point - low_corner) / (up_corner - low_corner);
  prd.result = make_float3(result_color_radiance, result_color_radiance, result_color_radiance) + alpha * new_prd.result;// + new_prd.result * exp_multiplier_2 * tex3D(photon_map_tex, normalized_hit_point.x, normalized_hit_point.y, normalized_hit_point.z) * tex1D(exp_tex, accumulative_density / exp_tex_resolution);
}
