//#####################################################################
// Copyright 2012, Bo Zhu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <optix_world.h>
#include "OPTIX_HELPERS.h"
using namespace optix;

rtDeclareVariable(float,scene_epsilon,,);
rtDeclareVariable(optix::Ray,ray,rtCurrentRay,);
rtDeclareVariable(float3,geometric_normal,attribute geometric_normal,);
rtDeclareVariable(float3,shading_normal,attribute shading_normal,);

RT_PROGRAM void intersect(int)
{
    ////unlimited scale groud
    float tmin=(-scene_epsilon-ray.origin.y)/ray.direction.y;
    if(rtPotentialIntersection(tmin)){
        shading_normal=geometric_normal=make_float3(0,1.0f,0);
        rtReportIntersection(0);
    }
}

RT_PROGRAM void box_bounds(int,float result[6]) 
{
    optix::Aabb* aabb=(optix::Aabb*)result;
    aabb->invalidate(); ////infinite area
}

