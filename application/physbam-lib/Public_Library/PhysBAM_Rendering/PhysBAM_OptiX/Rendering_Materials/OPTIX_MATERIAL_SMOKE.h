//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_MATERIAL_SMOKE
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_MATERIAL_SMOKE__
#define __OPTIX_MATERIAL_SMOKE__

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
class OPTIX_MATERIAL_SMOKE:public OPTIX_MATERIAL<T> {
private:
    typedef VECTOR<T,3> TV;
    OPTIX_RENDER_WORLD<T>* world;
    Material rt_material;

    int pow2(int n) {
      int x = 1;

      while(x < n) {
        x = x << 1;
      }

      return x;
    }

public:
    OPTIX_MATERIAL_SMOKE():OPTIX_MATERIAL(false),world(NULL){}
    ~OPTIX_MATERIAL_SMOKE(){}

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

        path_to_ptx = world_input->shader_prefix + "_generated_OPTIX_FOD_ATTENUATE_CALC.cu.ptx";
        closest_hit_program = rt_context->createProgramFromPTXFile(path_to_ptx, "closest_hit_radiance");

        rt_material = rt_context->createMaterial();
        rt_material->setClosestHitProgram(0, closest_hit_program);

        float step = 0.01f, density_multiplier = 4.f;
        rt_material["step"]->setFloat(step);
        rt_material["density_multiplier"]->setFloat(density_multiplier);
        rt_material["exp_tex_resolution"]->setFloat(density_multiplier);

        // putting exponent into a texture
        float exp_multiplier = step * 0.1f, exp_multiplier2 = 1.f;
        rt_material["exp_multiplier"]->setFloat(exp_multiplier);
        rt_material["exp_multiplier_2"]->setFloat(exp_multiplier2);
        /*
        int max_steps_number = log(exp_multiplier * 1e+3) / (step * density_multiplier);
        int number_in_each_step = 5;
        int texture_size = max(pow2(max_steps_number * number_in_each_step), 2048);
        Buffer tex_buffer = rt_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_FLOAT, texture_size);
        TextureSampler tex_sampler = rt_context->createTextureSampler();
        tex_sampler->setWrapMode(0, RT_WRAP_CLAMP_TO_EDGE);
        tex_sampler->setFilteringModes(RT_FILTER_LINEAR, RT_FILTER_LINEAR, RT_FILTER_NONE);
        tex_sampler->setIndexingMode(RT_TEXTURE_INDEX_ARRAY_INDEX); // RT_TEXTURE_INDEX_NORMALIZED_COORDINATES
        tex_sampler->setReadMode(RT_TEXTURE_READ_ELEMENT_TYPE);
        tex_sampler->setMaxAnisotropy(1.0f);
        tex_sampler->setMipLevelCount(1);
        tex_sampler->setArraySize(1);
        tex_sampler->setBuffer(0, 0, tex_buffer);

        std::cout << "Exponent texture\n";
        float* tex_buffer_data = static_cast<float*>(tex_buffer->map());
        for (int n = 0; n < texture_size; n++) {
            tex_buffer_data[n] = exp(-step * density_multiplier * n / (float)number_in_each_step);
            std::cout << tex_buffer_data[n] << " ";
        }
        std::cout << "\n";
        tex_buffer_data[texture_size - 1] = 0.f;
        tex_buffer->unmap();

        rt_material["exp_tex"]->setTextureSampler(tex_sampler);
        */

        // rt_material["Ka"]->setFloat(Util::Get_Float3<T>(ka));
    }

    void setDensityMultiplier(float density_multiplier) {
        rt_material["density_multiplier"]->setFloat(density_multiplier);
    }

    void setStep(float step) {
        rt_material["step"]->setFloat(step);
    }
};
}
#endif
#endif

