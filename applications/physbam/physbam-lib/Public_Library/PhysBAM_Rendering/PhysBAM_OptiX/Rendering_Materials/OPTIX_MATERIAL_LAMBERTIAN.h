//#####################################################################
// Copyright 2012, Bo Zhu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_MATERIAL_LAMBERTIAN
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_MATERIAL_LAMBERTIAN__
#define __OPTIX_MATERIAL_LAMBERTIAN__
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Materials/OPTIX_MATERIAL.h>
#include <optix_world.h>
#include <optixu/optixpp_namespace.h>

namespace PhysBAM{
using namespace optix;
template<class T>
class OPTIX_MATERIAL_LAMBERTIAN:public OPTIX_MATERIAL<T>
{
    typedef VECTOR<T,3> TV;
public:
    OPTIX_MATERIAL_LAMBERTIAN(const TV& kd_input=TV((T)0.2,(T)0.2,(T)0.2),const TV& ka_input=TV((T)0.8,(T)0.8,(T)0.8),
        bool use_constant_normal_input=false,bool use_shadow_input=false)
        :world(NULL),use_constant_normal(use_constant_normal_input),use_shadow(use_shadow_input),kd(kd_input),ka(ka_input){}
    virtual Material getRTMaterial(){return rt_material;}
    virtual void ensureInitialized(OPTIX_RENDER_WORLD<T> *world_input)  
    {
        if (world&&world_input==world)return;
        PHYSBAM_ASSERT(world==NULL);world=world_input;
        Context rt_context=world->RTContext();
        Program closest_hit_program;
        Program any_hit_program;
        std::string path_to_ptx;

        path_to_ptx=world_input->shader_prefix+"_generated_OPTIX_RADIANCE_CALC.cu.ptx";
        if(use_constant_normal&&use_shadow)closest_hit_program=rt_context->createProgramFromPTXFile(path_to_ptx,"closest_hit_radiance_floor_lambertian_with_shadow");
        else if(use_constant_normal&&!use_shadow)closest_hit_program=rt_context->createProgramFromPTXFile(path_to_ptx,"closest_hit_radiance_floor_lambertian");
        else if(!use_constant_normal&&use_shadow)closest_hit_program=rt_context->createProgramFromPTXFile(path_to_ptx,"closest_hit_radiance_lambertian_with_shadow");
        else/*!use_constant_normal&&!use_shadow*/closest_hit_program=rt_context->createProgramFromPTXFile(path_to_ptx,"closest_hit_radiance_lambertian");

        any_hit_program=rt_context->createProgramFromPTXFile(path_to_ptx,"any_hit_shadow_lambertian");

        rt_material=rt_context->createMaterial();
        rt_material->setClosestHitProgram(world->RAY_TYPE_RADIANCE,closest_hit_program);
        rt_material->setAnyHitProgram(world->RAY_TYPE_SHADOW,any_hit_program);

        rt_material["Ka"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(ka));
        rt_material["Kd"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(kd));
    }
private:
    OPTIX_RENDER_WORLD<T>* world;
    Material rt_material;
    bool use_constant_normal;
    bool use_shadow;
    TV kd,ka;
};
}
#endif
#endif
