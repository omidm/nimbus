//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_MATERIAL_PHOTON_MAP
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_MATERIAL_PHOTON_MAP__
#define __OPTIX_MATERIAL_PHOTON_MAP__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Materials/OPTIX_MATERIAL.h>

#include <optix_world.h>
#include <optixu/optixpp_namespace.h>

namespace PhysBAM{

using namespace optix;

template<class T,int d> class VECTOR;
template<class T> class OPTIX_RENDER_WORLD;

template<class T>
class OPTIX_MATERIAL_PHOTON_MAP:public OPTIX_MATERIAL<T> {
private:
    typedef VECTOR<T,3> TV;
    OPTIX_RENDER_WORLD<T>* world;
    Material rt_material;

public:
    OPTIX_MATERIAL_PHOTON_MAP():world(NULL) {}
    ~OPTIX_MATERIAL_PHOTON_MAP() {}

    Material getRTMaterial() { return rt_material; }
    void ensureInitialized(OPTIX_RENDER_WORLD<T> *world_input)  {
        if (world && world_input == world) {
            return;
        }

        // can not use one material in different render worlds (can't share same RTcontext)
        PHYSBAM_ASSERT(world == NULL);
        world = world_input;
        Context rt_context = world->RTContext();

        // TODO: add another ray type together with any hit program for rendering shadows
        Program closest_hit_program;
        std::string path_to_ptx;

        path_to_ptx = world_input->shader_prefix + "_generated_OPTIX_PHOTON_ATTENUATE_CALC.cu.ptx";
        closest_hit_program = rt_context->createProgramFromPTXFile(path_to_ptx, "closest_hit_radiance");

        rt_material = rt_context->createMaterial();
        rt_material->setClosestHitProgram(0, closest_hit_program);

        float step = 0.01f, color_multiplier = 4.f;
        rt_material["step"]->setFloat(step);
        rt_material["color_multiplier"]->setFloat(color_multiplier);
        rt_material["exp_tex_resolution"]->setFloat(color_multiplier);
    }

    void setColorMultiplier(float color_multiplier) {
        rt_material["color_multiplier"]->setFloat(color_multiplier);
    }
};
}
#endif
#endif


