//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <optix_world.h>
#include <optix.h>
#include <optixu/optixu_math_namespace.h>
#include "OPTIX_RAY_STRUCTS.h"
using namespace optix;

rtDeclareVariable(float3,bg_color,,);

rtDeclareVariable(PerRayData_radiance,prd_radiance,rtPayload,);
rtDeclareVariable(optix::Ray,ray,rtCurrentRay,);

RT_PROGRAM void miss() 
{
    prd_radiance.result=bg_color;
}
