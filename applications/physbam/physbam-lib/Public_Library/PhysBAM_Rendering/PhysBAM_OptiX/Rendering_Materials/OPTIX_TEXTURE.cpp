//#####################################################################
// Copyright 2012, Bo Zhu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_MATERIAL_FIRE
//#####################################################################
#ifdef USE_OPTIX
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Materials/OPTIX_TEXTURE.h>
#include <iostream>
#include <fstream>
using namespace PhysBAM;
using namespace optix;


template<class T,int DATA_N,int GRID_N>
PhysBAM::OPTIX_TEXTURE<T, DATA_N, GRID_N>::OPTIX_TEXTURE( std::string name_input,Context rt_context_input,TV_INT size_input )
:name(name_input),rt_context(rt_context_input),size(size_input),data(0)
{
    PHYSBAM_ASSERT(typeid(T)==typeid(float));   ////now only supports float texture
    RTformat format=(DATA_N==1?RT_FORMAT_FLOAT:(DATA_N==2?RT_FORMAT_FLOAT2:(DATA_N==3?RT_FORMAT_FLOAT3:RT_FORMAT_FLOAT4)));
    buffer=Create_Input_Buffer(format,size);

    TextureSampler texture_sampler=rt_context->createTextureSampler();
    texture_sampler->setWrapMode(0,RT_WRAP_CLAMP_TO_EDGE);
    texture_sampler->setWrapMode(1,RT_WRAP_CLAMP_TO_EDGE);
    texture_sampler->setWrapMode(2,RT_WRAP_CLAMP_TO_EDGE);
    texture_sampler->setFilteringModes(RT_FILTER_LINEAR,RT_FILTER_LINEAR,RT_FILTER_NONE);
    texture_sampler->setIndexingMode(RT_TEXTURE_INDEX_NORMALIZED_COORDINATES);
    texture_sampler->setReadMode(RT_TEXTURE_READ_ELEMENT_TYPE);
    texture_sampler->setMaxAnisotropy(1.0f);
    texture_sampler->setMipLevelCount(1);
    texture_sampler->setArraySize(1);
    texture_sampler->setBuffer(0,0,buffer);

    rt_context[name.c_str()]->setTextureSampler(texture_sampler);
    data=new ARRAY_VIEW<OPTIX_DATA_TYPE,TV_INT>(RANGE<TV_INT>(TV_INT::Constant_Vector(1),size),Map());Unmap();
}

template<class T,int DATA_N,int GRID_N>
void PhysBAM::OPTIX_TEXTURE<T, DATA_N, GRID_N>::Read_From_File( std::string file_name )
{
    OPTIX_DATA_TYPE* data=Map();
    std::ifstream in(file_name.c_str(),std::ios::binary);
    if(!in){PHYSBAM_FATAL_ERROR("cannot open texture file");}
    in.read((char*)(data),size.Product()*sizeof(OPTIX_DATA_TYPE));
    in.close();Unmap();
}

template<class T,int DATA_N,int GRID_N>
void PhysBAM::OPTIX_TEXTURE<T, DATA_N, GRID_N>::Copy_From( const OPTIX_DATA_TYPE* src_data )
{
    //std::cout<<"copy from "<<size<<std::endl;
    OPTIX_DATA_TYPE* data=Map();
    memcpy((void*)data,(void*)src_data,size.Product()*sizeof(OPTIX_DATA_TYPE));
    //for(int i=0;i<size.Product();i++){
    //    data[i]=src_data[i]; std::cout<<i<<", "<<data[i]<<std::endl;
    //}
    Unmap();
}

template class OPTIX_TEXTURE<float,1,1>;
template class OPTIX_TEXTURE<float,1,2>;
template class OPTIX_TEXTURE<float,1,3>;
template class OPTIX_TEXTURE<float,4,1>;
template class OPTIX_TEXTURE<float,4,2>;
template class OPTIX_TEXTURE<float,4,3>;
#endif
