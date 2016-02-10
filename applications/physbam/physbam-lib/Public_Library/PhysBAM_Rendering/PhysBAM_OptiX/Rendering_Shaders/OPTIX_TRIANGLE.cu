//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <optix_world.h>
#include "OPTIX_HELPERS.h"

using namespace optix;

rtDeclareVariable(float3, v1, , );
rtDeclareVariable(float3, v2, , );
rtDeclareVariable(float3, v3, , );

rtDeclareVariable(float3, geometric_normal, attribute geometric_normal, );
rtDeclareVariable(float3, shading_normal, attribute shading_normal, );
// rtDeclareVariable(float3, back_hit_point, attribute back_hit_point, );
// rtDeclareVariable(float3, front_hit_point, attribute front_hit_point, );
rtDeclareVariable(float, scene_epsilon, , );
rtDeclareVariable(optix::Ray, ray, rtCurrentRay, );

RT_PROGRAM void intersect(int) {
  // Intersect ray with triangle
  float3 n;
  float  t, beta, gamma;

  if(intersect_triangle(ray, v1, v2, v3, n, t, beta, gamma)) {
    if(rtPotentialIntersection(t)) {
      shading_normal = geometric_normal = normalize( n );

      // back_hit_point = front_hit_point = ray.origin + t * ray.direction + shading_normal * scene_epsilon;
      rtReportIntersection(0);
    }
  }
}

RT_PROGRAM void box_bounds (int, float result[6]) {
  optix::Aabb* aabb = (optix::Aabb*)result;

  aabb->m_min = fminf( fminf( v1, v2), v3 );
  aabb->m_max = fmaxf( fmaxf( v1, v2), v3 );
}
