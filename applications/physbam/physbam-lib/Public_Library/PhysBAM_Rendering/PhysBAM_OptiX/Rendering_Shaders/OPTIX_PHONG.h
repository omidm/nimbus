//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// OptiX phong material program radiance_calc.cu
//#####################################################################


#ifdef __CUDACC__

#include <optix.h>
#include <optixu/optixu_math_namespace.h>
#include <optixu/optixpp_namespace.h>
#include <optix_world.h>
#include "OPTIX_COMMONSTRUCTS.h"
#include "OPTIX_HELPERS.h"
#include "OPTIX_RAY_STRUCTS.h"
using namespace optix;

rtDeclareVariable(int,               max_depth, , );
rtBuffer<BasicLight>                 lights;
rtDeclareVariable(float3,            ambient_light_color, , );
// rtDeclareVariable(unsigned int,      radiance_ray_type, , );
// rtDeclareVariable(unsigned int,      shadow_ray_type, , );
// rtDeclareVariable(float,             scene_epsilon, , );
rtDeclareVariable(rtObject,          top_object, , );
rtDeclareVariable(float,          scene_epsilon, , );

// rtDeclareVariable(rtObject,          top_shadower, , );

rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );
rtDeclareVariable(float, t_hit, rtIntersectionDistance, );
rtDeclareVariable(PerRayData_radiance, prd, rtPayload, );
// rtDeclareVariable(float3, front_hit_point, attribute front_hit_point, );
// rtDeclareVariable(float3, back_hit_point, attribute back_hit_point, );
rtDeclareVariable(float3, extinction_constant, , );

rtDeclareVariable(PerRayData_shadow,   prd_shadow, rtPayload, );
rtDeclareVariable(rtObject,          top_opaque_object, , );

__device__ void phongShadowed(float3 p_normal) {
    rtIgnoreIntersection();
    prd_shadow.attenuation *= 0.5f;
    rtTerminateRay();
    /*
  float nDi = fabs(dot(p_normal, ray.direction));

  prd_shadow.attenuation *= 1 - fresnel_schlick(nDi, 5, make_float3(0.6f, 0.3f, 0.6f), make_float3(1));
  if(optix::luminance(prd_shadow.attenuation) < 0.01f)
    rtTerminateRay();
  else
    rtIgnoreIntersection();
    */
}

__device__ bool my_refract(float3& r, float3 i, float3 n, float ior)
{
  float3 nn = n;
  float negNdotV = dot(i,nn);
  float eta;

  if (negNdotV > 0.0f) {
    eta = ior;
    nn = -n;
    negNdotV = -negNdotV;
  } else {
    eta = 1.f / ior;
  }

  const float k = 1.f - eta*eta * (1.f - negNdotV * negNdotV);

  if (k < 0.0f) {
    return false;
  } else {
    r = normalize(eta*i - (eta*negNdotV + sqrtf(k)) * nn);
    return true;
  }
}

__device__ __inline__ float3 exp(const float3& x) {
  return make_float3(exp(x.x), exp(x.y), exp(x.z));
}

__device__ void phongShade( float3 p_Kd,
                            float3 p_Ka,
                            float3 p_Ks,
                            float3 p_Kt,
                            float3 p_normal,
                            float  p_phong_exp,
                            float3 p_reflectivity ) {
  prd.distance = t_hit;
  float3 hit_point = ray.origin + t_hit * ray.direction;
  // hit_point = rtTransformPoint(RT_OBJECT_TO_WORLD, hit_point);

  // ambient contribution

  // float3 result = make_float3(1.f, 1.f, 1.f);//p_Ka * ambient_light_color;
  float3 beer_attenuation = make_float3(1);
  /*
  if(dot(p_normal, ray.direction) > 0) {
    // Beer's law attenuation
    beer_attenuation = exp(extinction_constant * t_hit);
  } else {
    beer_attenuation = make_float3(1);
  }
  */


  float3 result = make_float3(0.f);
  result += p_Ka * ambient_light_color;

  // compute direct lighting
  unsigned int num_lights = lights.size();
  for(int i = 0; i < num_lights; ++i) {
    BasicLight light = lights[i];
    float Ldist = optix::length(light.pos - hit_point);
    float3 L = optix::normalize(light.pos - hit_point);
    float nDl = optix::dot( p_normal, L);

    // cast shadow ray
    float3 light_attenuation = make_float3(static_cast<float>( nDl > 0.0f ));
    if (nDl > 0.0f/* && light.casts_shadow*/ ) {
        /*
      PerRayData_shadow shadow_prd;
      shadow_prd.attenuation = make_float3(1.0f);
      optix::Ray shadow_ray = optix::make_Ray(fhp, L, 1, scene_epsilon, Ldist);
      rtTrace(top_opaque_object, shadow_ray, shadow_prd);
      light_attenuation = shadow_prd.attenuation;
      */
    }

    // If not completely shadowed, light the hit point
    if( fmaxf(light_attenuation) > 0.0f ) {
      float3 Lc = light.color * light_attenuation;
      result += p_Kd * nDl * Lc;

      float3 H = optix::normalize(L - ray.direction);
      float nDh = optix::dot( p_normal, H );
      if(nDh > 0 && fmaxf(p_Ks) > 0) {
        float power = pow(nDh, p_phong_exp);
        result += p_Ks * power * Lc;
      }
    }
  }

/*
if (prd.depth == 1) {
  prd.result = result;
  return;
}*/
  if( fmaxf( p_reflectivity ) > 0 && prd.depth < max_depth) {

    // ray tree attenuation
    PerRayData_radiance new_prd;
    new_prd.importance = prd.importance * optix::luminance( p_reflectivity );
    new_prd.depth = prd.depth + 1;
    new_prd.result = make_float3(0.f, 0.f, 0.f);

    // reflection ray
    if( new_prd.importance > 0.01f) {
      float3 R = optix::reflect( ray.direction, p_normal );
      optix::Ray refl_ray = optix::make_Ray(hit_point, R, 0, scene_epsilon, RT_DEFAULT_MAX);
      rtTrace(top_object, refl_ray, new_prd);
      result += p_reflectivity * new_prd.result * beer_attenuation;
      //result = make_float3(1, 1, 1);//new_prd.result;
      // result = p_reflectivity * new_prd.result;
    }
  }

  // Refraction
  if( fmaxf( p_Kt ) > 0.0f && prd.depth < max_depth ) {

    float3 t;// = ray.direction;
    if (my_refract(t, ray.direction, p_normal, 1.3f)) {

      // check for external or internal reflection
      float cos_theta = dot(ray.direction, p_normal);
      if (cos_theta < 0.0f)
        cos_theta = -cos_theta;
      else
        cos_theta = dot(t, p_normal);

      // float fresnel_exponent = 3.0f, fresnel_minimum = 1.0f, fresnel_maximum = 1.0f

      float reflection = 0.f;//fresnel_schlick(cos_theta, 4.0f, 0.1f, 1.f);
      PerRayData_radiance new_prd;
      new_prd.importance = prd.importance * (1.0f-reflection) * optix::luminance(p_Kt * beer_attenuation);//optix::luminance( p_Kt );
      new_prd.depth = prd.depth + 1;

      if (true/*new_prd.importance > 0.01f*/) {
        optix::Ray refl_ray = optix::make_Ray(bhp, t, 0, scene_epsilon, RT_DEFAULT_MAX);
        rtTrace(top_object, refl_ray, new_prd);
        if (new_prd.distance >= 0 && dot(p_normal, ray.direction) < 0) {
            result += (1.0f-reflection) * p_Kt * beer_attenuation * new_prd.result * lerp(make_float3(0.3f, 0.7f, 0.9f), make_float3(1.0f, 1.0f, 1.0f), exp(-5 * new_prd.distance));
        } else {
            result += (1.0f-reflection) * p_Kt * beer_attenuation * new_prd.result;
        }
        // result = new_prd.result;
      }
    } else { // full inner reflection
        // ray tree attenuation
        PerRayData_radiance new_prd;
        new_prd.importance = prd.importance * optix::luminance(p_Kt);
        new_prd.depth = prd.depth + 1;
        new_prd.result = make_float3(0.f, 0.f, 0.f);

        // reflection ray
        if( new_prd.importance > 0.01f) {
          float3 R = optix::reflect(ray.direction, p_normal);
          optix::Ray refl_ray = optix::make_Ray(hit_point, R, 0, scene_epsilon, RT_DEFAULT_MAX);
          rtTrace(top_object, refl_ray, new_prd);
          if (new_prd.distance >= 0 && dot(p_normal, ray.direction) < 0) {
            result += p_Kt * new_prd.result * beer_attenuation * lerp(make_float3(0.2f, 0.4f, 0.5f), make_float3(1.0f, 1.0f, 1.0f), exp(-5 * new_prd.distance));
          } else {
            result += p_Kt * new_prd.result * beer_attenuation;
          }

          //result = make_float3(1, 1, 1);//new_prd.result;
          // result = p_reflectivity * new_prd.result;
        }
    }
  }
  if (prd.depth >= max_depth && fmaxf(result) == 0) {
      result = make_float3(0.2f, 0.4f, 0.5f);
  }
/*
    // ray tree attenuation
    PerRayData_radiance new_prd;
    new_prd.importance = prd.importance * optix::luminance( p_Kt );
    new_prd.depth = prd.depth + 1;
    new_prd.result = make_float3(0.f, 0.f, 0.f);

    // refrection ray
    if( new_prd.importance > 0.01f ) {

      float3 R;
      if ( refract( R, ray.direction, p_normal, 1.2f) ) {
        optix::Ray refl_ray = optix::make_Ray(bhp, R, 0, scene_epsilon, RT_DEFAULT_MAX);
        rtTrace(top_object, refl_ray, new_prd);
        result += p_Kt * new_prd.result;
        // if (fminf(p_Kt * new_prd.result) < 0) {
        //    result = make_float3(1.f, 0.f, 0.f);
        //}
      }
    }
  }*/

/*  if (fminf(result) < 0) {
      result = make_float3(1.f, 0.f, 0.f);
  }
  // pass the color back up the tree
  prd.result = make_float3(1.f, 0.f, 0.f);//result;
  */
  prd.result = result;
}
#endif /* __CUDACC__ */
