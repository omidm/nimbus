//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <optix_world.h>
#include "OPTIX_HELPERS.h"
#include "OPTIX_RAY_STRUCTS.h"
using namespace optix;

rtDeclareVariable(float3,low_corner,,);
rtDeclareVariable(float3,up_corner,,);
rtDeclareVariable(float,scene_epsilon,,);
rtDeclareVariable(optix::Ray,ray,rtCurrentRay,);
rtDeclareVariable(float,front_box_hit_t,attribute front_box_hit_t,);
rtDeclareVariable(float,back_box_hit_t,attribute back_box_hit_t,);
rtDeclareVariable(float3,back_hit_point,attribute back_hit_point,);
rtDeclareVariable(float3,front_hit_point,attribute front_hit_point,);

__device__ bool operator<(float3 a,float3 b){return a.x<b.x&&a.y<b.y&&a.z<b.z;}
__device__ bool operator>(float3 a,float3 b){return a.x > b.x&&a.y>b.y&&a.z>b.z;}

RT_PROGRAM void intersect(int) 
{
    float3 t1=(up_corner-ray.origin)/ray.direction;
    float3 t2=(low_corner-ray.origin)/ray.direction;
    float3 T1=fminf(t1,t2);
    float3 T2=fmaxf(t1,t2);
    float Tnear=fmaxf(T1);
    float Tfar=fminf(T2);

    if (Tnear<Tfar&&Tfar>0){
        if (Tnear<0){Tnear=RAY_T_MIN*2;}
        if(rtPotentialIntersection(Tnear)){
            back_hit_point=front_hit_point=ray.origin+Tnear*ray.direction;
            front_box_hit_t=Tnear;
            back_box_hit_t=Tfar;
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

