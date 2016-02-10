//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <optix_world.h>
#include "OPTIX_COMMONSTRUCTS.h"
#include "OPTIX_RAY_STRUCTS.h"
#include "OPTIX_HELPERS.h"
using namespace optix;

rtTextureSampler<float, 3, cudaReadModeElementType> smoke_dencities_tex2;
rtBuffer<float, 3> photon_map;

rtBuffer<BasicLight>                 lights;
rtDeclareVariable(float3,       low_corner, , );
rtDeclareVariable(float3,       up_corner, , );
rtDeclareVariable(float,       scene_epsilon, , );
rtDeclareVariable(float3,       absorption, , );
rtDeclareVariable(float3,       scattering, , );
rtDeclareVariable(float,       step_size, , );
rtDeclareVariable(rtObject, smoke_object, , );
rtDeclareVariable(float,       step, , );
rtDeclareVariable(float,       color_multiplier, , );

rtDeclareVariable(uint3, launch_index, rtLaunchIndex, );
rtDeclareVariable(uint3, launch_dim, rtLaunchDim, );
rtDeclareVariable(rtObject,          top_opaque_object, , );

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
    return tex3D(smoke_dencities_tex2, tex_coord.x, tex_coord.y, tex_coord.z);
}

__device__ float intersect(optix::Ray ray) {
    float3 t1 = (up_corner - ray.origin) / ray.direction;
    float3 t2 = (low_corner - ray.origin) / ray.direction;

    float3 T1 = fminf(t1, t2);
    float3 T2 = fmaxf(t1, t2);

    float Tnear = fmaxf(T1);
    float Tfar = fminf(T2);

    if (Tnear < Tfar && Tfar > 0) {
        if (Tnear < 0) {
            Tnear = scene_epsilon * 10;
        }
        return Tnear;
    }
    return -1;
}

__device__ float attenuation(optix::Ray ray, float t_hit, float stop_distance) {
  float3 hit_point = ray.origin + (t_hit + scene_epsilon) * ray.direction;

  PerRayData_radiance new_prd;
  new_prd.depth = 0;
  new_prd.importance = 1.f;

  rtTrace(top_opaque_object, ray, new_prd);

  if (!isInsideDomain(hit_point)) {
     return 0.f;
  }

  float alpha = 1.0f;
  if (new_prd.distance > 0 && new_prd.distance < stop_distance) {
    alpha = 0.5f;
  }

  float density = getIntensity(hit_point);
  float accumulative_density = 0;
  float accumulative_step = t_hit;
  do {
     accumulative_density += density * step;
     hit_point += ray.direction * step;
     density = getIntensity(hit_point);
     accumulative_step += step;
  } while (density >= 0.f && accumulative_step < stop_distance);

  float res = color_multiplier * accumulative_density;//exp(-density_multiplier * accumulative_density);
  if (res > 1.0f) {
    res = 1.0f;
  } else if (res < 0.f) {
    res = 0.f;
  }
  return (1.0f - res) * alpha;
}

RT_PROGRAM void ray_gen() {
  unsigned int num_lights = lights.size();
  float3 offsets[4] = {make_float3(0.25f, 0.25f, 0.25f),
                       // make_float3(-0.25f, 0.25f, 0.25f),
                       make_float3(-0.25f, -0.25f, 0.25f),
                       // make_float3(0.25f, -0.25f, 0.25f),

                       // make_float3(0.25f, 0.25f, -0.25f),
                       make_float3(-0.25f, 0.25f, -0.25f),
                       // make_float3(-0.25f, -0.25f, -0.25f),
                       make_float3(0.25f, -0.25f, -0.25f)};

  photon_map[launch_index] = 0.f;
  for(int i = 0; i < num_lights; ++i) {
      // photon_map[launch_index] += 1.0f;

      for (int j = 0; j < 4; j++) {
          float3 grid_point = ((make_float3(launch_index) + make_float3(0.5f, 0.5f, 0.5f) + offsets[j]) / make_float3(launch_dim)) * (up_corner - low_corner) + low_corner;
          optix::Ray ray = optix::make_Ray(lights[i].pos, normalize(grid_point - lights[i].pos), 0, scene_epsilon, RT_DEFAULT_MAX);

          // rtTrace(smoke_object, ray, prd);
          float distance = intersect(ray);
          if (distance > 0) {
            photon_map[launch_index] += 0.25f * attenuation(ray, distance, length(lights[i].pos - grid_point));
          }
      }
      // accumulating photons from all light sources
      // prd.result = 1.0f;
      // photon_map[launch_index] += prd.result;
  }
}

RT_PROGRAM void exception() {
  const unsigned int code = rtGetExceptionCode();
  rtPrintf( "Caught exception during photon mapping 0x%X at launch index (%d,%d,%d)\n", code, launch_index.x, launch_index.y, launch_index.z);
  photon_map[launch_index] = 0.f;
}


