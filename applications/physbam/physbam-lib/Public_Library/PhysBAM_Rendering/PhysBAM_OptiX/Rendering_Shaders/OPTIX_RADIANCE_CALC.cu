//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth, Bo Zhu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// OptiX phong material program radiance_calc.cu
//#####################################################################

#include <optix_world.h>
#include <optixu/optixu_math_namespace.h>
#include <optix.h>
#include <optixu/optixpp_namespace.h>
#include <optix_world.h>
#include "OPTIX_COMMONSTRUCTS.h"
#include "OPTIX_HELPERS.h"
#include "OPTIX_RAY_STRUCTS.h"

using namespace optix;

rtDeclareVariable(float3,Ka,,);
rtDeclareVariable(float3,Kd,,);
rtDeclareVariable(float3,Ks,,);
rtDeclareVariable(float3,Kt,,);
rtDeclareVariable(float3,reflectivity,,);
rtDeclareVariable(float,phong_exp,,);
rtDeclareVariable(float,refract_coef,,);
rtDeclareVariable(float,transpar_atten_coef,,);
rtDeclareVariable(float3,cutoff_color,,);
rtDeclareVariable(float3,ambient_light_color,,);

rtDeclareVariable(float3,shading_normal,attribute shading_normal,);
rtDeclareVariable(float3,geometric_normal,attribute geometric_normal,);
rtDeclareVariable(int,max_depth,,);
rtDeclareVariable(float,scene_epsilon,,);

rtDeclareVariable(float,t_hit,rtIntersectionDistance,);
rtDeclareVariable(float3,extinction_constant,,);
rtDeclareVariable(rtObject,top_object,,);
rtDeclareVariable(rtObject,top_opaque_object,,);

rtDeclareVariable(optix::Ray,ray,rtCurrentRay,);
rtDeclareVariable(PerRayData_radiance,prd,rtPayload,);
rtDeclareVariable(PerRayData_shadow,prd_shadow,rtPayload,);

rtBuffer<BasicLight> lights;

RT_PROGRAM void any_hit_shadow_phong() 
{
    rtIgnoreIntersection();
    prd_shadow.attenuation*=0.5f;
    rtTerminateRay();
}

RT_PROGRAM void closest_hit_radiance_phong() {
  float3 world_shading_normal   = normalize( rtTransformNormal( RT_OBJECT_TO_WORLD, shading_normal ) );
  float3 world_geometric_normal = normalize( rtTransformNormal( RT_OBJECT_TO_WORLD, geometric_normal ) );

  prd.distance = t_hit;
  float3 hit_point = ray.origin + t_hit * ray.direction;
  float3 beer_attenuation = make_float3(1);

  float3 result = make_float3(0.f);
  result += Ka * ambient_light_color;

  // compute direct lighting
  uint num_lights = lights.size();
  for (int i = 0; i < num_lights; ++i) {
    BasicLight light = lights[i];
    float Ldist = optix::length(light.pos - hit_point);
    float3 L = optix::normalize(light.pos - hit_point);
    float nDl = optix::dot(world_shading_normal, L);

    // cast shadow ray
    float3 light_attenuation = make_float3(static_cast<float>( nDl > 0.0f ));
    /*
    // no shadows for speeding
    if (nDl > 0.0f && light.casts_shadow) {
      PerRayData_shadow shadow_prd;
      shadow_prd.attenuation = make_float3(1.0f);
      optix::Ray shadow_ray = optix::make_Ray(fhp, L, 1, scene_epsilon, Ldist);
      rtTrace(top_opaque_object, shadow_ray, shadow_prd);
      light_attenuation = shadow_prd.attenuation;
    }
    */

    // If not completely shadowed, light the hit point
    if (fmaxf(light_attenuation) > 0.0f) {
      float3 Lc = light.color * light_attenuation;
      result += Kd * nDl * Lc;

      float3 H = optix::normalize(L - ray.direction);
      float nDh = optix::dot(world_shading_normal, H);
      if(nDh > 0 && fmaxf(Ks) > 0) {
        float power = pow(nDh, phong_exp);
        result += Ks * power * Lc;
      }
    }
  }

  if (fmaxf(reflectivity) > 0 && prd.depth < max_depth) {
    // ray tree attenuation
    PerRayData_radiance new_prd;
    new_prd.importance = prd.importance * optix::luminance(reflectivity);
    new_prd.depth = prd.depth + 1;
    new_prd.result = make_float3(0.f, 0.f, 0.f);

    // reflection ray
    if (new_prd.importance > 0.01f) {
      float3 R = optix::reflect(ray.direction, world_shading_normal);
      optix::Ray refl_ray = optix::make_Ray(hit_point, R, 0, scene_epsilon, RT_DEFAULT_MAX);
      rtTrace(top_object, refl_ray, new_prd);
      result += reflectivity * new_prd.result * beer_attenuation;
    }
  }

  // Refraction
  if (fmaxf(Kt) > 0.0f && prd.depth < max_depth) {

    float3 t;
    if (refract(t, ray.direction, world_shading_normal, refract_coef)) {
      // check for external or internal reflection
      float cos_theta = dot(ray.direction, world_shading_normal);
      if (cos_theta < 0.0f)
        cos_theta = -cos_theta;
      else
        cos_theta = dot(t, world_shading_normal);

      // float fresnel_exponent = 3.0f, fresnel_minimum = 1.0f, fresnel_maximum = 1.0f

      float reflection = 0.f;//fresnel_schlick(cos_theta, 4.0f, 0.1f, 1.f);
      PerRayData_radiance new_prd;
      new_prd.importance = prd.importance * (1.0f-reflection) * optix::luminance(Kt * beer_attenuation);
      new_prd.depth = prd.depth + 1;

      if (new_prd.importance > 0.01f) {
        optix::Ray refl_ray = optix::make_Ray(hit_point, t, 0, scene_epsilon, RT_DEFAULT_MAX);
        rtTrace(top_object, refl_ray, new_prd);
        // attenuate by depth if we are coming into the surface
        if (new_prd.distance >= 0 && dot(world_shading_normal, ray.direction) < 0) {
            result += (1.0f-reflection) * Kt * beer_attenuation * new_prd.result * lerp(cutoff_color/*make_float3(0.3f, 0.7f, 0.9f)*/, make_float3(1.0f, 1.0f, 1.0f), exp(-transpar_atten_coef * new_prd.distance));
        } else {
            result += (1.0f-reflection) * Kt * beer_attenuation * new_prd.result;
        }
      }
    } else { // full inner reflection
        // ray tree attenuation
        PerRayData_radiance new_prd;
        new_prd.importance = prd.importance * optix::luminance(Kt);
        new_prd.depth = prd.depth + 1;
        new_prd.result = make_float3(0.f, 0.f, 0.f);

        // reflection ray
        if (new_prd.importance > 0.01f) {
          float3 R = optix::reflect(ray.direction, world_shading_normal);
          optix::Ray refl_ray = optix::make_Ray(hit_point, R, 0, scene_epsilon, RT_DEFAULT_MAX);
          rtTrace(top_object, refl_ray, new_prd);
          if (new_prd.distance >= 0 && dot(world_shading_normal, ray.direction) < 0) {
            result += Kt * new_prd.result * beer_attenuation * lerp(cutoff_color, make_float3(1.0f, 1.0f, 1.0f), exp(-transpar_atten_coef * new_prd.distance));
          } else {
            result += Kt * new_prd.result * beer_attenuation;
          }
        }
    }
  }

  // cutoff color
  if (prd.depth >= max_depth && fmaxf(result) == 0) {
      result = cutoff_color;
  }
  prd.result = result;
}

RT_PROGRAM void any_hit_shadow_lambertian()
{
  prd_shadow.attenuation=make_float3(0);
  rtTerminateRay();
}

RT_PROGRAM void closest_hit_radiance_floor_lambertian()
{
    float3 ffnormal=make_float3(0,1.0f,0);
    float3 color=Ka*ambient_light_color;

    float3 hit_point=ray.origin+t_hit*ray.direction;
    for(int i=0;i<lights.size();i++){
        BasicLight light=lights[i];
        float3 light_hit_dir=normalize(light.pos-hit_point);
        float lam=dot(ffnormal,light_hit_dir);
        if(lam>0)color+=Kd*lam*light.color;
    }
    prd.result=color;
}

RT_PROGRAM void closest_hit_radiance_floor_lambertian_with_shadow()
{
    float3 ffnormal=make_float3(0,1.0f,0);
    float3 color=Ka*ambient_light_color;

    float3 hit_point=ray.origin+t_hit*ray.direction;
    for(int i=0;i<lights.size();i++){
        BasicLight light=lights[i];
        float3 light_hit_dir=normalize(light.pos-hit_point);
        float light_hit_distance=length(light.pos-hit_point);
        float lam=dot(ffnormal,light_hit_dir);
        if(lam>0){
            optix::Ray ray_shadow=optix::make_Ray(light.pos,-light_hit_dir,RAY_TYPE_SHADOW,RAY_T_MIN,light_hit_distance-RAY_T_MIN);
			PerRayData_shadow prd_shadow;prd_shadow.attenuation=make_float3(1.0f);
			rtTrace(top_object,ray_shadow,prd_shadow);
			float3 light_attenuation=prd_shadow.attenuation;
			color+=Kd*lam*light.color*light_attenuation;
        }
    }
    prd.result=color;
}

RT_PROGRAM void closest_hit_radiance_lambertian()
{
    float3 world_shading_normal=normalize(rtTransformNormal(RT_OBJECT_TO_WORLD,shading_normal));
    float3 world_geometric_normal=normalize(rtTransformNormal(RT_OBJECT_TO_WORLD,geometric_normal));
    float3 ffnormal=faceforward(world_shading_normal,-ray.direction,world_geometric_normal);
    float3 color=Ka*ambient_light_color;

    float3 hit_point=ray.origin+t_hit*ray.direction;
    for(int i=0;i<lights.size();i++){
        BasicLight light=lights[i];
        float3 light_hit_dir=normalize(light.pos-hit_point);
        float lam=dot(ffnormal,light_hit_dir);
        if(lam>0)color+=Kd*lam*light.color;
    }
    prd.result=color;
}

RT_PROGRAM void closest_hit_radiance_lambertian_with_shadow()
{
    float3 world_shading_normal=normalize(rtTransformNormal(RT_OBJECT_TO_WORLD,shading_normal));
    float3 world_geometric_normal=normalize(rtTransformNormal(RT_OBJECT_TO_WORLD,geometric_normal));
    float3 ffnormal=faceforward(world_shading_normal,-ray.direction,world_geometric_normal);
    float3 color=Ka*ambient_light_color;

    float3 hit_point=ray.origin+t_hit*ray.direction;
    for(int i=0;i<lights.size();i++){
        BasicLight light=lights[i];
        float3 light_hit_dir=normalize(light.pos-hit_point);
        float light_hit_distance=length(light.pos-hit_point);
        float lam=dot(ffnormal,light_hit_dir);
        if(lam>0){
            optix::Ray ray_shadow=optix::make_Ray(light.pos,-light_hit_dir,RAY_TYPE_SHADOW,RAY_T_MIN,light_hit_distance-RAY_T_MIN);
			PerRayData_shadow prd_shadow;prd_shadow.attenuation=make_float3(1.0f);
			rtTrace(top_object,ray_shadow,prd_shadow);
			float3 light_attenuation=prd_shadow.attenuation;
			color+=Kd*lam*light.color*light_attenuation;
        }
    }
    prd.result=color;
}