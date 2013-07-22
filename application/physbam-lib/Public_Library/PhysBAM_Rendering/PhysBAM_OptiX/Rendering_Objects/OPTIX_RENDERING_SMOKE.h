//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_RENDERING_SMOKE
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_RENDERING_SMOKE__
#define __OPTIX_RENDERING_SMOKE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>

#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Dynamics/Geometry/GENERAL_GEOMETRY_FORWARD.h>

#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD_KEY_LISTENER.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Objects/OPTIX_RENDERING_AABB_BOX.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Materials/OPTIX_MATERIAL_SMOKE.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Materials/OPTIX_MATERIAL_PHOTON_MAP.h>

#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>

#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Tessellation/IMPLICIT_OBJECT_TESSELLATION.h>

#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/SCALAR_COMPONENT.h>

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <iostream>
#include <string>

#include <optix_world.h>
#include <optixu/optixpp_namespace.h>

#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Utilities/OPTIX_UTILITIES.h>
namespace PhysBAM{

using namespace optix;

template<class T> class OPTIX_RENDERING_OBJECT;
template<class T> class OPTIX_RENDER_WORLD;

template<class T> class OPTIX_RENDERING_SMOKE : public OPTIX_RENDERING_GEOMETRY<T>, public OPTIX_RENDER_WORLD_KEY_LISTENER {
    typedef VECTOR<T,3> TV;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;

    bool needReinitialization;
    Program  rt_object_intersection_program;
    Program  rt_object_bounding_box_program;
    Acceleration rt_acceleration;

    OPTIX_RENDERING_AABB_BOX<T> smoke_bounded_box;
    // OPTIX_RENDERING_AABB_BOX<T> smoke_bounded_box2;
    OPTIX_MATERIAL_SMOKE<T> smoke_material;
    float density_multiplier;

    ARRAY<Transform> rt_transform_array;
    // VECTOR<unsigned int, 3> dimensions;
    OPTIX_RENDER_WORLD<T>* my_world_input;
    // Material rt_material;
    Buffer photon_map_tex_buffer;

    ARRAY<T,TV_INT>* my_densities_simulated;
    GRID<TV>* my_grid;

    float* densities_data;
    int densities_data_size;

    TextureSampler tex_sampler, photon_map_tex_sampler;
    Buffer tex_buffer;
    // MATERIAL_PHOTON_MAP<T> photon_map_material;
    float color_multiplier;
    Buffer rt_result_photon_buffer;

    float my_step;
    Program rt_ray_gen_program, exception_program;

    void SetupForRenderingPhotonTexture() {
        Context rt_context = my_world_input->RTContext();

        // might neef RT_BUFFER_INPUT_OUTPUT type for several light sources - photons will be accumulated then,
        // but may also pass light sources as an array
        // std::cout << "Photon buffer metrix(" << my_grid->counts.x << ", " << my_grid->counts.y << ", " << my_grid->counts.z << ")\n";
        rt_result_photon_buffer = rt_context->createBuffer(RT_BUFFER_OUTPUT, RT_FORMAT_FLOAT,
                                                          my_grid->counts.x, my_grid->counts.y, my_grid->counts.z);
        rt_context["photon_map"]->set(rt_result_photon_buffer);

        // add ray generation program
        std::string path_to_ptx = my_world_input->shader_prefix + "_generated_OPTIX_PHOTON_RAY_GEN.cu.ptx";

        // hanging photong mapping onto second entry point
        rt_ray_gen_program = rt_context->createProgramFromPTXFile(path_to_ptx, "ray_gen");
        rt_context->setRayGenerationProgram(1, rt_ray_gen_program);
        exception_program = rt_context->createProgramFromPTXFile(path_to_ptx, "exception");
        rt_context->setExceptionProgram(1, exception_program);

        // lights sources are set by OPTIX_RENDER_WORLD and are available through context variable
        // TODO: should think at this point, how volumetric objects might be organized in a group to pass it as a top_object
        rt_ray_gen_program["low_corner"]->setFloat(OPTIX_UTILITIES::Get_Float3(smoke_bounded_box.getLowCorner()));
        rt_ray_gen_program["up_corner"]->setFloat(OPTIX_UTILITIES::Get_Float3(smoke_bounded_box.getUpCorner()));
        rt_context["smoke_dencities_tex2"]->setTextureSampler(tex_sampler);

        // std::cout << "my_step = " << my_step << "\n";
        rt_ray_gen_program["step"]->setFloat(my_step);
        rt_ray_gen_program["color_multiplier"]->setFloat(color_multiplier);
        /*
        // setting smoke object
        Geometry rt_object = smoke_bounded_box.getGeometryObject();
        photon_map_material.ensureInitialized(my_world_input);

        // assigning material to object - creating instance
        Material rt_material=photon_map_material.getRTMaterial();
        GeometryInstance rt_instance = rt_context->createGeometryInstance();
        rt_instance->setMaterialCount(1);
        rt_instance->setMaterial(0, rt_material);
        rt_instance->setGeometry(rt_object);

        // putting instance into a group
        GeometryGroup rt_geometrygroup = rt_context->createGeometryGroup();
        rt_geometrygroup->setChildCount(1);
        rt_geometrygroup->setChild(0, rt_instance);

        Acceleration rt_acceleration = rt_context->createAcceleration("NoAccel","NoAccel");
        rt_geometrygroup->setAcceleration(rt_acceleration);
        rt_acceleration->markDirty();

        rt_ray_gen_program["smoke_object"]->set(rt_geometrygroup);
        rt_material["low_corner"]->setFloat(OPTIX_UTILITIES::Get_Float3(smoke_bounded_box.getLowCorner()));
        rt_material["up_corner"]->setFloat(OPTIX_UTILITIES::Get_Float3(smoke_bounded_box.getUpCorner()));
        rt_material["smoke_dencities_tex"]->setTextureSampler(tex_sampler);
        */

        // setting smoke object
        // GeometryGroup rt_geometrygroup = smoke_bounded_box2.getGeometryGroup();

        /*
        Group group = rt_context->createGroup();
        group->setChildCount(1);
        group->setChild(0, smoke_bounded_box2.getTransform());// objects(i)->getRTTransform());

        Acceleration group_acceleration = rt_context->createAcceleration("MedianBvh","Bvh");
        group->setAcceleration(group_acceleration);
        group_acceleration->markDirty();

        photon_map_material.ensureInitialized(my_world_input);

        // assigning material to object - creating instance
        Material rt_material=photon_map_material.getRTMaterial();

        // set smoke_object
        std::cout << "Set smoke object\n";
        rt_ray_gen_program["smoke_object"]->set(group);

        // setting dimensions
        smoke_bounded_box2.setLowCorner(my_grid->domain.min_corner);
        smoke_bounded_box2.setUpCorner(my_grid->domain.max_corner);
        rt_material["low_corner"]->setFloat(OPTIX_UTILITIES::Get_Float3(smoke_bounded_box.getLowCorner()));
        rt_material["up_corner"]->setFloat(OPTIX_UTILITIES::Get_Float3(smoke_bounded_box.getUpCorner()));
        rt_ray_gen_program["low_corner"]->setFloat(OPTIX_UTILITIES::Get_Float3(smoke_bounded_box.getLowCorner()));
        rt_ray_gen_program["up_corner"]->setFloat(OPTIX_UTILITIES::Get_Float3(smoke_bounded_box.getUpCorner()));


        rt_material["smoke_dencities_tex"]->setTextureSampler(tex_sampler);
        */

        // rt_context->validate();
        // rt_context->compile();
    }

    void LoadPhotonOuputBufferIntoTexture() {
        // std::cout << "Copy photon mapping\n";
        Context rt_context = my_world_input->RTContext();
        Buffer generated_texture = rt_context["photon_map"]->getBuffer();
        float* tex_buffer_data = static_cast<float*>(photon_map_tex_buffer->map());
        float* generated_tex_buffer_data = static_cast<float*>(generated_texture->map());
        memcpy(tex_buffer_data, generated_tex_buffer_data, densities_data_size);
        // memset(tex_buffer_data, 0, densities_data_size);
        /*
        for (int i = 0; i < my_grid->counts.x * my_grid->counts.y * my_grid->counts.z; i++) {
            tex_buffer_data[i] = 1.f;
        }*/
        generated_texture->unmap();
        photon_map_tex_buffer->unmap();
    }

    void CreateEmptySamplers() {
        // std::cout << "Buffer metrix(" << my_grid->counts.x << ", " << my_grid->counts.y << ", " << my_grid->counts.z << ")\n";

        Context context = my_world_input->RTContext();
        photon_map_tex_buffer = context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_FLOAT, my_grid->counts.x, my_grid->counts.y, my_grid->counts.z);
        photon_map_tex_sampler = context->createTextureSampler();

        photon_map_tex_sampler->setWrapMode(0, RT_WRAP_CLAMP_TO_EDGE);
        photon_map_tex_sampler->setWrapMode(1, RT_WRAP_CLAMP_TO_EDGE);
        photon_map_tex_sampler->setWrapMode(2, RT_WRAP_CLAMP_TO_EDGE);
        photon_map_tex_sampler->setFilteringModes(RT_FILTER_LINEAR, RT_FILTER_LINEAR, RT_FILTER_NONE);
        photon_map_tex_sampler->setIndexingMode(RT_TEXTURE_INDEX_NORMALIZED_COORDINATES);
        photon_map_tex_sampler->setReadMode(RT_TEXTURE_READ_ELEMENT_TYPE);
        photon_map_tex_sampler->setMaxAnisotropy(1.0f);
        photon_map_tex_sampler->setMipLevelCount(1);
        photon_map_tex_sampler->setArraySize(1);
        photon_map_tex_sampler->setBuffer(0, 0, photon_map_tex_buffer);

        tex_buffer = context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_FLOAT, my_grid->counts.x, my_grid->counts.y, my_grid->counts.z);
        tex_sampler = context->createTextureSampler();
        tex_sampler->setWrapMode(0, RT_WRAP_CLAMP_TO_EDGE);
        tex_sampler->setWrapMode(1, RT_WRAP_CLAMP_TO_EDGE);
        tex_sampler->setWrapMode(2, RT_WRAP_CLAMP_TO_EDGE);
        tex_sampler->setFilteringModes(RT_FILTER_LINEAR, RT_FILTER_LINEAR, RT_FILTER_NONE);
        tex_sampler->setIndexingMode(RT_TEXTURE_INDEX_NORMALIZED_COORDINATES);
        tex_sampler->setReadMode(RT_TEXTURE_READ_ELEMENT_TYPE);
        tex_sampler->setMaxAnisotropy(1.0f);
        tex_sampler->setMipLevelCount(1);
        tex_sampler->setArraySize(1);
        tex_sampler->setBuffer(0, 0, tex_buffer);

        /*float* tex_buffer_data = static_cast<float*>(tex_buffer->map());
        tex_buffer_data[0] = 1.0f;
        tex_buffer->unmap();*/

        smoke_material.getRTMaterial()["smoke_dencities_tex"]->setTextureSampler(tex_sampler);
        smoke_material.getRTMaterial()["photon_map_tex"]->setTextureSampler(photon_map_tex_sampler);
    }

    void SetupForRenderingSmoke() {
        Context context = my_world_input->RTContext();

        // std::cout << "Rendering photon texture with metrix(" << my_grid->counts.x << ", " << my_grid->counts.y << ", " << my_grid->counts.z << ")\n";
        // std::cout << "Launch attenutation\n";
        context->launch(1, my_grid->counts.x, my_grid->counts.y, my_grid->counts.z);
        // std::cout << "Attenutation finished\n";
        LoadPhotonOuputBufferIntoTexture();

        // photon_map_tex_buffer = context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_FLOAT, my_grid->counts.x, my_grid->counts.y, my_grid->counts.z);
        // photon_map_tex_sampler->setBuffer(0, 0, photon_map_tex_buffer);

        // Buffer tex_buffer = context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_FLOAT, my_grid->counts.x, my_grid->counts.y, my_grid->counts.z);
        float* tex_buffer_data = static_cast<float*>(tex_buffer->map());
        memcpy(tex_buffer_data, densities_data, densities_data_size);
        tex_buffer->unmap();
        // tex_sampler->setBuffer(0, 0, tex_buffer);

        // float* tex_buffer_data = static_cast<float*>(tex_buffer->map());
        // memset(tex_buffer_data, 0, sizeof(float) * dimensions.x*dimensions.y*dimensions.z);
        // tex_buffer->unmap();
        // SetupForRenderingPhotonTexture();
        // filling texture with data

        // Buffer generated_texture = rt_context["photon_map"]->getBuffer();

        // float* generated_tex_buffer_data = static_cast<float*>(generated_texture->map());
        // memcpy(tex_buffer_data, generated_tex_buffer_data, sizeof(float) * dimensions.z * dimensions.y * dimensions.x);
        /*
        for (int z = 0; z < dimensions.z; z++) {
            float vz = ((dimensions.z - 1 - z) /(float)(dimensions.z - 1) - 0.5f) * 2;
            for (int y = 0; y < dimensions.y; y++) {
                float vy = ((dimensions.y - 1 - y)/(float)(dimensions.y - 1) - 0.5f) * 2;
                for (int x = 0; x < dimensions.x; x++) {
                    TV v(((dimensions.x - 1 - x)/(float)(dimensions.x - 1) - 0.5f) * 2, vy, vz);
                    int index = x + dimensions.x * (y + z * dimensions.y);
                    float length = (float)v.Magnitude();
                    if (length > 1) {
                        tex_buffer_data[index] = 0;
                    } else {
                        tex_buffer_data[index] = 1.f - length;
                    }
                }
            }
        }
        */
        /*
        for (int i = 0; i < dimensions.x*dimensions.y*dimensions.z;i++) {
            tex_buffer_data[i] = 1.f;
        }*/
        // rt_material["smoke_dencities_tex"]->setTextureSampler(tex_sampler);
        // rt_material["photon_map_tex"]->setTextureSampler(photon_map_tex_sampler);
    }
public:
    OPTIX_RENDERING_SMOKE():OPTIX_RENDERING_GEOMETRY<T>("Smoke"),smoke_bounded_box(TV(0, 0 , 0), TV(1, 1, 1), &smoke_material)/*,
    smoke_bounded_box2(TV(0, 0 , 0), TV(1, 1, 1), &photon_map_material)*//*,dimensions(100, 100, 100)*/,
    densities_data(NULL), needReinitialization(false) {
        opaque=false;   
        my_step = 0.01f, color_multiplier = 4.f;
    }

    void Initialize(OPTIX_RENDER_WORLD<T>* world_input) {
        // std::cout << "Initialize\n";
        my_world_input = world_input;
        smoke_bounded_box.Initialize(world_input);
        // photon_map_material.ensureInitialized(world_input);
        //smoke_bounded_box2.Initialize(world_input);
        world_input->Add_Key_Listener(this);
        rt_transform_array.Resize(1);
        rt_transform_array(1) = smoke_bounded_box.getTransform();

        Material rt_material = smoke_material.getRTMaterial();

        /*
        std::cout << "low_corner " << smoke_bounded_box.getLowCorner().x << ", " << smoke_bounded_box.getLowCorner().y << ", " << smoke_bounded_box.getLowCorner().z << ")\n";
        std::cout << "up_corner " << smoke_bounded_box.getUpCorner().x << ", " << smoke_bounded_box.getUpCorner().y << ", " << smoke_bounded_box.getUpCorner().z << ")\n";
        std::cout << "my_grid->x " << my_grid->counts.x << "\n";
        std::cout << "my_grid->y " << my_grid->counts.y << "\n";
        std::cout << "my_grid->z " << my_grid->counts.z << "\n";
        */

        rt_material["low_corner"]->setFloat(OPTIX_UTILITIES::Get_Float3(smoke_bounded_box.getLowCorner()));
        rt_material["up_corner"]->setFloat(OPTIX_UTILITIES::Get_Float3(smoke_bounded_box.getUpCorner()));
        // rt_material["dimension"]->setUint(dimensions.x, dimensions.y, dimensions.z);
        rt_material["scene_epsilon"]->setFloat(0.01f);

        CreateEmptySamplers();
        SetupForRenderingPhotonTexture();
        /*
        if (needReinitialization) {
            SetupForRenderingSmoke();
        }
        needReinitialization = false;
        */
        color_multiplier = 0.f;
        // photon_map_material.setColorMultiplier(color_multiplier);
        density_multiplier = 15.f;
        smoke_material.setDensityMultiplier(density_multiplier);
        smoke_material.setStep(my_step);
    }

    void BeforeLaunch() {
    }

    // Group getGroup() { return rt_group; }
    ARRAY<Transform>& getTransforms() { return rt_transform_array; }

    void KeyPressed(unsigned char key, int x, int y) {
        switch(key) {
        case 'Q':
            density_multiplier -= 0.2f;
            smoke_material.setDensityMultiplier(density_multiplier);
            break;
        case 'q':
            density_multiplier += 0.2f;
            smoke_material.setDensityMultiplier(density_multiplier);
            break;
        case 'W':
            color_multiplier -= 0.2f;
            rt_ray_gen_program["color_multiplier"]->setFloat(color_multiplier);
            // photon_map_material.setColorMultiplier(color_multiplier);
            break;
        case 'w':
            color_multiplier += 0.2f;
            rt_ray_gen_program["color_multiplier"]->setFloat(color_multiplier);
            // photon_map_material.setColorMultiplier(color_multiplier);
            break;
        }
        // std::cout << "color_multiplier = " << color_multiplier << "\n";
        // std::cout << "density_multiplier = " << density_multiplier << "\n";
    }

    void PrepareForDraw() {
        if (needReinitialization) {
            SetupForRenderingSmoke();
        }
        needReinitialization = false;
    }

    void Update(T time) {
        // std::cout << "Update\n";
        /*
        if (needReinitialization) {
            SetupForRenderingSmoke();
        }
        needReinitialization = false;
        */
    }

    void Set_Densities(ARRAY<T,TV_INT>* densities_simulated, GRID<TV>* grid) {
        // std::cout << "Set_Densities\n";
        my_step = ((grid->domain.max_corner - grid->domain.min_corner) / VECTOR<T, 3>((T)grid->counts.x,(T)grid->counts.y,(T)grid->counts.z)).Min() * (T)0.5;
        // std::cout << "my_step = " << my_step << "\n";
        // smoke_bounded_box2.setLowCorner(grid.domain.min_corner);
        // smoke_bounded_box2.setUpCorner(grid.domain.max_corner);

        smoke_bounded_box.setLowCorner(grid->domain.min_corner);
        smoke_bounded_box.setUpCorner(grid->domain.max_corner);

        my_densities_simulated = densities_simulated;
        my_grid = grid;
        if (densities_data != NULL) {
            delete [] densities_data;
        }
        densities_data = new float [grid->counts.x * grid->counts.y * grid->counts.z];
        densities_data_size = grid->counts.x * grid->counts.y * grid->counts.z * sizeof(float);

        ReinitializeFromSimulation();
    }

    void ReinitializeFromSimulation() {
        // std::cout << "ReinitializeFromSimulation\n";
            for(int i=1;i<=my_grid->counts.x;i++)
            for(int j=1;j<=my_grid->counts.y;j++)
            for(int k=1;k<=my_grid->counts.z;k++){
                densities_data[i - 1 + (j - 1) * my_grid->counts.x + (k - 1) * my_grid->counts.x * my_grid->counts.y] =
                    (float)(*my_densities_simulated)(i,j,k);
                    // std::cout << (float)(*my_densities_simulated)(i,j,k) << " ";
            }

        needReinitialization = true;
    }
};
}
#endif
#endif


