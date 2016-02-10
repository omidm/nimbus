//#####################################################################
// Copyright 2012, Bo Zhu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifdef USE_OPTIX
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Lights/OPTIX_RENDERING_FIRE_LIGHT.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace PhysBAM;
using namespace optix;

template<class T>
void PhysBAM::OPTIX_RENDERING_FIRE_LIGHT<T>::Initialize_Fire_Color_Table(const std::string& fire_color_table_file)
{
    fire_color_table.Exact_Resize(4096);
    std::ifstream in(fire_color_table_file.c_str(),std::ios::binary);
    if(!in){PHYSBAM_FATAL_ERROR("cannot open file in OPTIX_RENDERING_FIRE_LIGHT<T>::Initialize_Fire_Color_Table()");}
    in.read((char*)fire_color_table.Get_Array_Pointer(),fire_color_table.m*sizeof(float4));
    in.close();
}

template<class T>
void PhysBAM::OPTIX_RENDERING_FIRE_LIGHT<T>::Initialize_Spot_Lights()
{
    optix_lights.Resize(RANGE<TV_INT>(TV_INT::Constant_Vector(1),spot_light_grid.counts));
    light_ranges.Resize(RANGE<TV_INT>(TV_INT::Constant_Vector(1),spot_light_grid.counts));

    for(int i=1;i<=spot_light_grid.counts.x;i++){
        for(int j=1;j<=spot_light_grid.counts.y;j++){
            for(int k=1;k<=spot_light_grid.counts.z;k++){
                TV center=spot_light_grid.Center(i,j,k);

                BasicLight light;
                light.pos=OPTIX_UTILITIES::Get_Float3(center);
                light.color=OPTIX_UTILITIES::Get_Float3(TV(0,0,0));
                light.casts_shadow=0;
                light.light_on=0;

                TV range_min=center-(T)0.5*spot_light_grid.dX;
                TV range_max=center+(T)0.5*spot_light_grid.dX;
                if(i==1)range_min.x=temperature_grid.domain.min_corner.x+(T)0.1*temperature_grid.dX.x;
                else if(i==spot_light_grid.counts.x)range_max.x=temperature_grid.domain.max_corner.x-(T)0.1*temperature_grid.dX.x;
                if(j==1)range_min.y=temperature_grid.domain.min_corner.y+(T)0.1*temperature_grid.dX.y;
                else if(j==spot_light_grid.counts.y)range_max.y=temperature_grid.domain.max_corner.y-(T)0.1*temperature_grid.dX.y;
                if(k==1)range_min.z=temperature_grid.domain.min_corner.z+(T)0.1*temperature_grid.dX.z;
                else if(k==spot_light_grid.counts.z)range_max.z=temperature_grid.domain.max_corner.z-(T)0.1*temperature_grid.dX.z;

                optix_lights(i,j,k)=light;
                light_ranges(i,j,k)=RANGE<TV>(range_min,range_max);
            }
        }
    }
}

template<class T>
void PhysBAM::OPTIX_RENDERING_FIRE_LIGHT<T>::Update_Spot_Lights( ARRAY_VIEW<T,TV_INT>& temperature )
{
    T max_t=0;
    for(int i=1;i<=temperature.Size().x;i++){
        for(int j=1;j<=temperature.Size().y;j++){
            for(int k=1;k<=temperature.Size().z;k++){
                if(temperature(i,j,k)>max_t){
                    max_t=temperature(i,j,k);
                }
            }
        }
    }

    float light_t_threshold=2500.0f;
    for(int i=1;i<=spot_light_number;i++){
        T max_t=0;
        TV max_t_pos;
        Max_Temperature_In_Range(max_t,max_t_pos,light_ranges.array(i),temperature);
        if(max_t>light_t_threshold){
            optix_lights.array(i).light_on=true;
            optix_lights.array(i).color=Fire_Color(max_t)*0.1f;
            optix_lights.array(i).pos=OPTIX_UTILITIES::Get_Float3(max_t_pos);
        }
        else{
            optix_lights.array(i).light_on=false;
        }
    }
}

template<class T>
float3 PhysBAM::OPTIX_RENDERING_FIRE_LIGHT<T>::Fire_Color( T t )
{
    if(t>=temperature_max)return make_float3(fire_color_table(fire_color_table.m));
    else if(t<=temperature_min)return make_float3(fire_color_table(1));

    int index=(int)((t-temperature_min)/temperature_step);
    float frac=t-(float)index*temperature_step;
    float4 c1=fire_color_table(index+1),c2=fire_color_table(index+2);
    return make_float3((1.0f-frac)*c1+frac*c2);
}

template<class T>
void PhysBAM::OPTIX_RENDERING_FIRE_LIGHT<T>::Max_Temperature_In_Range( T& max_t,TV& max_t_pos,const RANGE<TV>& range,const ARRAY_VIEW<T,TV_INT>& temperature )
{
    TV_INT index_min=TV_INT((range.min_corner-temperature_grid.domain.min_corner)*temperature_grid.one_over_dX)+1;
    TV_INT index_max=TV_INT((range.max_corner-temperature_grid.domain.min_corner)*temperature_grid.one_over_dX)+1;
    PHYSBAM_ASSERT(temperature_grid.Inside_Domain(index_min)&&temperature_grid.Inside_Domain(index_max));

    T max_t0=0;TV max_t_pos0;
    for(int i=index_min.x;i<=index_max.x;i++){
        for(int j=index_min.y;j<=index_max.y;j++){
            for(int k=index_min.z;k<=index_max.z;k++){
                T t=temperature(i,j,k);
                if(t>max_t0){max_t0=t;max_t_pos0=temperature_grid.Center(i,j,k);}
            }
        }
    }
    max_t=max_t0;max_t_pos=max_t_pos0;
}

template class OPTIX_RENDERING_FIRE_LIGHT<float>;
#endif
