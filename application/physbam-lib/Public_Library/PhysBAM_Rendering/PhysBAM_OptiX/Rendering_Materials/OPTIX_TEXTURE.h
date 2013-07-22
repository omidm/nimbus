//#####################################################################
// Copyright 2012, Bo Zhu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_TEXTURE
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_TEXTURE__
#define __OPTIX_TEXTURE__
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND_VIEW.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <optix_world.h>
#include <optixu/optixpp_namespace.h>
#include <string>

namespace PhysBAM{
using namespace optix;

template<class T,int DATA_N,int GRID_N>
class OPTIX_TEXTURE
{
    typedef VECTOR<int,GRID_N> TV_INT;
    typedef typename IF<DATA_N==1,float,typename IF<DATA_N==2,float2,typename IF<DATA_N==3,float3,float4>::TYPE>::TYPE>::TYPE OPTIX_DATA_TYPE;
public:
    const std::string name;
    
    OPTIX_TEXTURE(std::string name_input,Context rt_context_input,TV_INT size_input);
    OPTIX_DATA_TYPE* Map(){return (OPTIX_DATA_TYPE*)buffer->map();}
    void Unmap(){buffer->unmap();}
    TV_INT Size(){return size;}
    void Read_From_File(std::string file_name);
    void Copy_From(const OPTIX_DATA_TYPE* src_data);
    ARRAY_VIEW<OPTIX_DATA_TYPE,TV_INT>* Get_Texture_Data(){return data;}
private:
    Context rt_context;
    Buffer buffer;
    TV_INT size;
    ARRAY_VIEW<OPTIX_DATA_TYPE,TV_INT>* data;

    Buffer Create_Input_Buffer(RTformat format,VECTOR<int,1> size){return rt_context->createBuffer(RT_BUFFER_INPUT,format,size.x);}
    Buffer Create_Input_Buffer(RTformat format,VECTOR<int,2> size){return rt_context->createBuffer(RT_BUFFER_INPUT,format,size.x,size.y);}
    Buffer Create_Input_Buffer(RTformat format,VECTOR<int,3> size){return rt_context->createBuffer(RT_BUFFER_INPUT,format,size.x,size.y,size.z);}
};
}
#endif
#endif
