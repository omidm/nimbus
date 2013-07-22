//#####################################################################
// Copyright 2011, Valeria Nikolaenko
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <optix_world.h>
#include <optixu/optixu_math_namespace.h>
#include "OPTIX_COMMONSTRUCTS.h"
#include "OPTIX_HELPERS.h"
#include "OPTIX_RAY_STRUCTS.h"
using namespace optix;

rtTextureSampler<float, 3, cudaReadModeElementType> smoke_dencities_tex;

rtDeclareVariable(float3,       low_corner, , );
rtDeclareVariable(float3,       up_corner, , );
rtDeclareVariable(float,       scene_epsilon, , );
rtDeclareVariable(float,       step, , );
rtDeclareVariable(float,       color_multiplier, , );
// rtDeclareVariable(float,       exp_multiplier, , );

// rtDeclareVariable(rtObject,          smoke_object, , );
rtDeclareVariable(rtObject,          top_opaque_object, , );

rtDeclareVariable(PerRayData_photon, prd, rtPayload, );
rtDeclareVariable(float, t_hit, rtIntersectionDistance, );
rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );

//rtDeclareVariable(float, distance, attribute distance, );
rtDeclareVariable(float, front_box_hit_point, attribute front_box_hit_point, );

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

RT_PROGRAM void closest_hit_radiance() {
  float3 hit_point = ray.origin + (t_hit + scene_epsilon) * ray.direction;

  /*
  float t = t_hit - front_box_hit_point;
  if (t > scene_epsilon || t < -scene_epsilon) {
    prd.result = make_float3(1.f, 0.f, 0.f);
    return;
  }*/

  PerRayData_radiance new_prd;
  new_prd.depth = 0;
  new_prd.importance = 1.f;

  optix::Ray next_ray = optix::make_Ray(ray.origin, ray.direction, 0, scene_epsilon, RT_DEFAULT_MAX);
  rtTrace(top_opaque_object, next_ray, new_prd);

  if (!isInsideDomain(hit_point)) {
     prd.result = 0.f;
     return;
  }

  if (new_prd.distance > 0 && new_prd.distance < prd.stop_distance) {
    prd.result = 0.f;
    return;
  }

  float density = getIntensity(hit_point);
  float accumulative_density = 0;
  float accumulative_step = front_box_hit_point;
  do {
     accumulative_density += density * step;
     // hit_point += ray.direction * step;
     density = getIntensity(hit_point);
     accumulative_step += step;
  } while (density >= 0.f && accumulative_step < prd.stop_distance);

  float res = color_multiplier * accumulative_density;//exp(-density_multiplier * accumulative_density);
  if (res > 1.0f) {
    res = 1.0f;
  } else if (res < 0.f) {
    res = 0.f;
  }
  prd.result = /*getIntensity(ray.origin + prd.stop_distance * ray.direction) * */(1.0f - res);
}
