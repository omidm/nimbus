//#####################################################################
// Copyright 2012, Bo Zhu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_RENDERING_FIRE_LIGHT__
#define __OPTIX_RENDERING_FIRE_LIGHT__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Lights/OPTIX_RENDERING_LIGHT.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Shaders/OPTIX_COMMONSTRUCTS.h>
#include <string>
#include <optixu/optixpp_namespace.h>
#include <optixu/optixu_math_namespace.h>

namespace PhysBAM{
using namespace optix;
template<class T> class OPTIX_RENDER_WORLD;

template<class T>
class OPTIX_RENDERING_FIRE_LIGHT : public OPTIX_RENDERING_LIGHT<T>
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;typedef GRID<TV> T_GRID;
public:
    int spot_light_number;
    T_GRID spot_light_grid;
    T_GRID temperature_grid;
    int start_index;      ////index in the render_world light queue

    OPTIX_RENDERING_FIRE_LIGHT(const TV_INT& light_counts_input,const RANGE<TV>& light_domain_input,
        const TV_INT& temperature_counts_input,const RANGE<TV>& temperature_domain_input,const std::string& fire_color_table_file_input) 
        :spot_light_number(light_counts_input.Product()),spot_light_grid(light_counts_input,light_domain_input,true),
        temperature_grid(temperature_counts_input,temperature_domain_input,true),start_index(0)
    {Initialize_Fire_Color_Table(fire_color_table_file_input);Initialize_Spot_Lights();}

    virtual void Add_To_Optix_Render_World(OPTIX_RENDER_WORLD<T>* optix_render_world){optix_render_world->Add_Fire_Light(this);}
    BasicLight* Get_Basic_Light(int index){return &optix_lights.array(index);}
    void Update_Spot_Lights(ARRAY_VIEW<T,TV_INT>& temperature);  ////all the updates are done on CPU
private:
    ARRAY<BasicLight,TV_INT> optix_lights;
    ARRAY<RANGE<TV>,TV_INT> light_ranges;
    ARRAY<float4> fire_color_table;
    T temperature_max,temperature_min,temperature_step;

    void Initialize_Fire_Color_Table(const std::string& fire_color_table_file);
    void Initialize_Spot_Lights();
    void Max_Temperature_In_Range(T& max_t,TV& max_t_pos,const RANGE<TV>& range,const ARRAY_VIEW<T,TV_INT>& temperature);
    float3 Fire_Color(T t);
};
}
#endif
#endif
