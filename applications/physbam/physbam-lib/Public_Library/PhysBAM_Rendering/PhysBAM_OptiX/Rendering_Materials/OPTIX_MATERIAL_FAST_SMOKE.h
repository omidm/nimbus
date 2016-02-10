//#####################################################################
// Copyright 2012, Bo Zhu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_MATERIAL_FAST_SMOKE
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_MATERIAL_FAST_SMOKE__
#define __OPTIX_MATERIAL_FAST_SMOKE__
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Materials/OPTIX_MATERIAL.h>
#include <PhysBAM_Rendering/PhysBAM_Optix/Rendering_Materials/OPTIX_TEXTURE.h>
#include <optix_world.h>
#include <optixu/optixpp_namespace.h>

namespace PhysBAM{
using namespace optix;

template<class T>
class OPTIX_MATERIAL_FAST_SMOKE:public OPTIX_MATERIAL<T>
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    TV_INT size;
    OPTIX_TEXTURE<T,1,3>* soot_texture;

    OPTIX_MATERIAL_FAST_SMOKE(const TV_INT& size_input)
        :OPTIX_MATERIAL<T>(false)/*smoke is not opaque*/,size(size_input),soot_texture(0),world(NULL){}
    virtual Material getRTMaterial(){return rt_material;}
    virtual void ensureInitialized(OPTIX_RENDER_WORLD<T> *world_input)
    {
        if (world&&world_input==world)return;
        PHYSBAM_ASSERT(world==NULL);world=world_input;
        Context rt_context=world->RTContext();
        Program closest_hit_program;
        Program closest_hit_program_shadow;
        std::string path_to_ptx;

        path_to_ptx=world_input->shader_prefix+"_generated_OPTIX_PARTICIPATING_MEDIA.cu.ptx";
        closest_hit_program=rt_context->createProgramFromPTXFile(path_to_ptx,"closest_hit_radiance_smoke");
        closest_hit_program_shadow=rt_context->createProgramFromPTXFile(path_to_ptx,"closest_hit_shadow_smoke");
        rt_material=rt_context->createMaterial();
        rt_material->setClosestHitProgram(world->RAY_TYPE_RADIANCE,closest_hit_program);
        rt_material->setClosestHitProgram(world->RAY_TYPE_SHADOW,closest_hit_program_shadow);
        soot_texture=new OPTIX_TEXTURE<T,1,3>("soot_texture",rt_context,size);
    }
    void Update_Texture_From_File(const std::string& soot_file)
    {
        soot_texture->Read_From_File(soot_file.c_str());
    }
    void Update_Texture_From_Data(const T* data)
    {
        soot_texture->Copy_From(data);
    }
private:
    OPTIX_RENDER_WORLD<T>* world;
    Material rt_material;
};
}
#endif
#endif