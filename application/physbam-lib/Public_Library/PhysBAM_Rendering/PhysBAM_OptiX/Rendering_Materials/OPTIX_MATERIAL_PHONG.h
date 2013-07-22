//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_MATERIAL_PHONG
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_MATERIAL_PHONG__
#define __OPTIX_MATERIAL_PHONG__
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Materials/OPTIX_MATERIAL.h>
#include <optix_world.h>
#include <optixu/optixpp_namespace.h>

namespace PhysBAM{
using namespace optix;
template<class T>
class OPTIX_MATERIAL_PHONG:public OPTIX_MATERIAL<T> 
{
    typedef VECTOR<T,3> TV;
public:
    OPTIX_MATERIAL_PHONG():world(NULL),ka(0,0,0),kd((T)0.5,(T)0.5,(T)0.5),ks(0,0,0),
        kt(0,0,0),reflectivity(0,0,0),cutoff_color((T)1.0,(T)1.0,(T)1.0),phong_exp((T)1.0),refract_coef((T)1.3){}
    OPTIX_MATERIAL_PHONG(TV Ka,TV Kd,TV Ks,TV Kt,TV Reflectivity,T Phong_exp,TV Cutoff_color=TV((T)1.0,(T)1.0,(T)1.0),T Refract_coef=(T)1.3):
        world(NULL),ka(Ka),kd(Kd),ks(Ks),kt(Kt),reflectivity(Reflectivity),cutoff_color(Cutoff_color),phong_exp(Phong_exp),refract_coef(Refract_coef){}
    ~OPTIX_MATERIAL_PHONG(){}

    Material getRTMaterial(){return rt_material;}
    void ensureInitialized(OPTIX_RENDER_WORLD<T> *world_input)  
    {
        if (world&&world_input==world)return;
        // can not use one material in different render worlds (can't share same RTcontext)
        PHYSBAM_ASSERT(world==NULL);world=world_input;
        Context rt_context=world->RTContext();

        // TODO: add another ray type together with any hit program for rendering shadows
        Program closest_hit_program;
        std::string path_to_ptx;

        path_to_ptx=world_input->shader_prefix+"_generated_OPTIX_RADIANCE_CALC.cu.ptx";
        closest_hit_program=rt_context->createProgramFromPTXFile(path_to_ptx,"closest_hit_radiance_phong");
        Program any_hit_program=rt_context->createProgramFromPTXFile(path_to_ptx,"any_hit_shadow_phong");

        rt_material=rt_context->createMaterial();
        rt_material->setClosestHitProgram(0,closest_hit_program);
        rt_material->setAnyHitProgram(1,any_hit_program);

        rt_material["Ka"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(ka));
        rt_material["Ks"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(ks));
        rt_material["Kd"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(kd));
        rt_material["Kt"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(kt));
        rt_material["reflectivity"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(reflectivity));
        rt_material["phong_exp"]->setFloat(phong_exp);
        rt_material["refract_coef"]->setFloat(refract_coef);
        rt_material["transpar_atten_coef"]->setFloat(5.f);
        rt_material["cutoff_color"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(cutoff_color));
        float3 extinction = make_float3(.80f,.80f,.80f);
        rt_material["extinction_constant"]->setFloat(log(extinction.x),log(extinction.y),log(extinction.z));
    }
private:
    OPTIX_RENDER_WORLD<T>* world;
    Material rt_material;
    TV ka,kd,ks,kt,reflectivity,cutoff_color;
    T phong_exp,refract_coef;
};
}
#endif
#endif
