//#####################################################################
// Copyright 2012, Bo Zhu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_MATERIAL_FIRE
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_MATERIAL_FIRE__
#define __OPTIX_MATERIAL_FIRE__
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Materials/OPTIX_MATERIAL.h>
#include <PhysBAM_Rendering/PhysBAM_Optix/Rendering_Materials/OPTIX_TEXTURE.h>
#include <optix_world.h>
#include <optixu/optixpp_namespace.h>

namespace PhysBAM{
using namespace optix;

template<class T>
class OPTIX_MATERIAL_FIRE:public OPTIX_MATERIAL<T>
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,3> TV_INT;
public:
    TV_INT size;

    OPTIX_TEXTURE<T,4,1>* fire_color_texture;
    OPTIX_TEXTURE<T,1,3>* soot_texture;
    OPTIX_TEXTURE<T,1,3>* temperature_texture;

    OPTIX_MATERIAL_FIRE(const TV_INT& size_input,const std::string& fire_color_table_file_name_input)
        :OPTIX_MATERIAL<T>(false)/*fire is not opaque*/,size(size_input),world(NULL),fire_color_table_file_name(fire_color_table_file_name_input){}
    virtual Material getRTMaterial(){return rt_material;}
    virtual void ensureInitialized(OPTIX_RENDER_WORLD<T> *world_input)  
    {
        if (world&&world_input==world)return;
        PHYSBAM_ASSERT(world==NULL);world=world_input;
        Context rt_context=world->RTContext();
        Program closest_hit_program;
        Program any_hit_program;
        std::string path_to_ptx;

        path_to_ptx=world_input->shader_prefix+"_generated_OPTIX_PARTICIPATING_MEDIA.cu.ptx";
        closest_hit_program=rt_context->createProgramFromPTXFile(path_to_ptx,"closest_hit_radiance_fire");

        rt_material=rt_context->createMaterial();
        rt_material->setClosestHitProgram(world->RAY_TYPE_RADIANCE,closest_hit_program);

        fire_color_texture=new OPTIX_TEXTURE<T,4,1>("fire_color_texture",rt_context,VECTOR<int,1>(4096));
        fire_color_texture->Read_From_File(fire_color_table_file_name.c_str());

        soot_texture=new OPTIX_TEXTURE<T,1,3>("soot_texture",rt_context,size);
        temperature_texture=new OPTIX_TEXTURE<T,1,3>("temperature_texture",rt_context,size);
    }

    void Update_Texture_From_File(const std::string& soot_file,const std::string& temperature_file)
    {
        soot_texture->Read_From_File(soot_file.c_str());
        temperature_texture->Read_From_File(temperature_file.c_str());
    }
    void Update_Texture_From_Data(const T* soot,const T* temperature)
    {
        soot_texture->Copy_From(soot);
        temperature_texture->Copy_From(temperature);
    }
private:
    OPTIX_RENDER_WORLD<T>* world;
    Material rt_material;
    std::string fire_color_table_file_name;
};
}
#endif
#endif
