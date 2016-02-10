//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <optix_world.h>
#include "OPTIX_HELPERS.h"
using namespace optix;

rtDeclareVariable(float3,low_corner,,);
rtDeclareVariable(float3,up_corner,,);
rtDeclareVariable(optix::Ray,ray,rtCurrentRay,);
rtDeclareVariable(float3,geometric_normal,attribute geometric_normal,);
rtDeclareVariable(float3,shading_normal,attribute shading_normal,);

__device__ float3 boxnormal(float t) 
{
  float3 t0=(low_corner-ray.origin)/ray.direction;
  float3 t1=(up_corner-ray.origin)/ray.direction;
  float3 neg=make_float3(t==t0.x?1:0,t==t0.y?1:0,t==t0.z?1:0);
  float3 pos=make_float3(t==t1.x?1:0,t==t1.y?1:0,t==t1.z?1:0);
  return pos-neg;
}

RT_PROGRAM void intersect(int) 
{
    float3 t1=(up_corner-ray.origin)/ray.direction;
    float3 t2=(low_corner-ray.origin)/ray.direction;
    float3 T1=fminf(t1,t2);
    float3 T2=fmaxf(t1,t2);
    float Tnear=fmaxf(T1);
    float Tfar=fminf(T2);

    if(Tnear<=Tfar){
        bool checkAnother=true;
        if(rtPotentialIntersection(Tnear)){
            shading_normal=geometric_normal=boxnormal(Tnear);
            if(rtReportIntersection(0))checkAnother=false;
        }
        if(checkAnother&&rtPotentialIntersection(Tfar)){
            shading_normal=geometric_normal=boxnormal(Tfar);
            rtReportIntersection(0);
        }
    }
}

RT_PROGRAM void box_bounds(int,float result[6]) 
{
  optix::Aabb* aabb=(optix::Aabb*)result;
  aabb->m_min=low_corner;
  aabb->m_max=up_corner;
}

