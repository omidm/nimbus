//#####################################################################
// Copyright 2012, Bo Zhu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_RAY_STRUCTS__
#define __OPTIX_RAY_STRUCTS__

#include <optix_world.h>
#include <optixu/optixu_math_namespace.h>
#include "OPTIX_HELPERS.h"
using namespace optix;

struct PerRayData_radiance 
{
    float3 result;
    float importance;
    float distance;
    int depth;
};

struct PerRayData_shadow 
{
    float3 attenuation;
};

struct PerRayData_photon 
{
    float result;
    float stop_distance;
};

rtDeclareVariable(uint,RAY_TYPE_RADIANCE,,)=0u;
rtDeclareVariable(uint,RAY_TYPE_SHADOW,,)=1u;
rtDeclareVariable(float,RAY_T_MIN,,)=1.0e-4f;
rtDeclareVariable(float,RAY_T_MAX,,)=RT_DEFAULT_MAX;

#endif
#endif