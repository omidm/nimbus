//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth, Bo Zhu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_RENDERING_SPOT_LIGHT__
#define __OPTIX_RENDERING_SPOT_LIGHT__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Shaders/OPTIX_COMMONSTRUCTS.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Lights/OPTIX_RENDERING_LIGHT.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD.h>
#include <optixu/optixpp_namespace.h>
#include <optixu/optixu_math_namespace.h>

namespace PhysBAM{
using optix::float3;
using optix::make_float3;
template<class T> class OPTIX_RENDER_WORLD;

template<class T>
class OPTIX_RENDERING_SPOT_LIGHT : public OPTIX_RENDERING_LIGHT<T>
{
    BasicLight optix_light;
public:
    int light_index;
    VECTOR<T,3> position;
    VECTOR<T,3> color;
    T brightness;
    int start_index;      ////index in the render_world light queue

    OPTIX_RENDERING_SPOT_LIGHT(VECTOR<T,3> position_input,VECTOR<T,3> color_input,T brightness_input) 
        :position(position_input),color(color_input),brightness(brightness_input),start_index(0)
    {
        optix_light.color=make_float3(color.x,color.y,color.z);
        optix_light.pos=make_float3(position.x,position.y,position.z);
        optix_light.casts_shadow=0;
        optix_light.light_on=0;
    }

    OPTIX_RENDERING_SPOT_LIGHT& operator=(const OPTIX_RENDERING_SPOT_LIGHT& l) 
    {
        light_index=l.light_index;
        position=l.position;
        color=l.color;
        brightness=l.brightness;
        optix_light.color=make_float3(l.color.x,l.color.y,l.color.z);
        optix_light.pos=make_float3(l.position.x,l.position.y,l.position.z);
        return *this;
    }

    BasicLight* Get_Basic_Light(int index=0){return &optix_light;}
    virtual void Add_To_Optix_Render_World(OPTIX_RENDER_WORLD<T>* optix_render_world){optix_render_world->Add_Spot_Light(this);}
};
}
#endif
#endif
