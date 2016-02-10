//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <optix_world.h>
#include "OPTIX_HELPERS.h"
using namespace optix;

rtDeclareVariable(float3,center,,);
rtDeclareVariable(float,radius,,);
rtDeclareVariable(float3,geometric_normal,attribute geometric_normal,);
rtDeclareVariable(float3,shading_normal,attribute shading_normal,);
rtDeclareVariable(float,scene_epsilon,,);
rtDeclareVariable(optix::Ray,ray,rtCurrentRay,);

RT_PROGRAM void intersect(int) 
{
  float3 O=ray.origin-center;
  float3 D=ray.direction;

  float b=dot(O, D);
  float c=dot(O, O)-radius*radius;
  float disc=b*b-c;
  if(disc>0.0f){
    float sdisc=sqrtf(disc);
    float root1=(-b-sdisc);
    bool check_second=true;
    if(rtPotentialIntersection(root1)){
      shading_normal=geometric_normal=(O + root1*D)/radius;
      float3 hit_p=ray.origin+root1*ray.direction;
      float3 offset=shading_normal*scene_epsilon;
      if(rtReportIntersection(0))check_second = false;
    }
    if(check_second){
      float root2=(-b+sdisc);
      if(rtPotentialIntersection(root2)){
        shading_normal=geometric_normal=(O+root2*D)/radius;
        float3 hit_p=ray.origin+root2*ray.direction;
        float3 offset=shading_normal*scene_epsilon;
        float t=dot(shading_normal,ray.direction)<0?1:-1;
        rtReportIntersection(0);
      }
    }
  }
}

RT_PROGRAM void box_bounds(int,float result[6]) 
{
  const float3 rad=make_float3(radius);
  optix::Aabb* aabb=(optix::Aabb*)result;

  if(radius>0.0f&&!isinf(radius)){
    aabb->m_min=center-rad;
    aabb->m_max=center+rad;
  } 
  else{
    aabb->invalidate();
  }
}
