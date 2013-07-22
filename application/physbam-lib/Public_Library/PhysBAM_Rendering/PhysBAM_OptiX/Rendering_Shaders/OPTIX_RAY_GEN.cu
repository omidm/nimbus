//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <optix_world.h>
#include "OPTIX_HELPERS.h"
#include "OPTIX_RAY_STRUCTS.h"
using namespace optix;

rtDeclareVariable(float2,halfed_proj_metrics,,);
rtDeclareVariable(uint2,screen_metrics,,);
rtDeclareVariable(float3,up,,);
rtDeclareVariable(float3,right,,);
rtDeclareVariable(float3,dir,,);
rtDeclareVariable(float3,loc,,);
rtDeclareVariable(float3,from_loc_to_focal,,);
rtDeclareVariable(uint2,launch_index,rtLaunchIndex,);
rtDeclareVariable(uint2,launch_dim,rtLaunchDim,);
rtDeclareVariable(rtObject,top_object,,);
rtDeclareVariable(float,scene_epsilon,,);
rtBuffer<uchar4, 2> output_buffer;

RT_PROGRAM void ray_gen() 
{
    float2 d=make_float2(launch_index)/make_float2(launch_dim)*2.f-1.f;
    float3 ray_origin=loc;
    float3 ray_direction=normalize(d.x*right*halfed_proj_metrics.x+d.y*up*halfed_proj_metrics.y+from_loc_to_focal);
    optix::Ray ray=optix::make_Ray(ray_origin,ray_direction,0,scene_epsilon,RT_DEFAULT_MAX);

    PerRayData_radiance prd;
    prd.depth=0;
    prd.importance=1;

    rtTrace(top_object,ray,prd);
    output_buffer[launch_index]=make_color(prd.result);
}

RT_PROGRAM void exception() 
{
  const unsigned int code=rtGetExceptionCode();
  rtPrintExceptionDetails();
  rtPrintf("Caught exception 0x%X at launch index (%d,%d)\n",code,launch_index.x,launch_index.y);
  output_buffer[launch_index]=make_uchar4(1.0f,1.0f,1.0f,0);
}

