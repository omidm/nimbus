//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_RENDERING_IMPLICIT_SURFACE
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_RENDERING_IMPLICIT_SURFACE__
#define __OPTIX_RENDERING_IMPLICIT_SURFACE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>

#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Dynamics/Geometry/GENERAL_GEOMETRY_FORWARD.h>

#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Objects/OPTIX_RENDERING_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>

#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>

#include <iostream>
#include <string>

#include <optix_world.h>
#include <optixu/optixpp_namespace.h>

#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Utilities/OPTIX_UTILITIES.h>

namespace PhysBAM{

using namespace optix;

template<class T> class OPTIX_RENDERING_IMPLICIT_SURFACE : public OPTIX_RENDERING_OBJECT<T> {
    typedef VECTOR<T,3> TV;
    GENERIC_PARSER<T> *parser;
    Geometry rt_object;
    Program  rt_object_intersection_program;
    Program  rt_object_bounding_box_program;
    Acceleration rt_acceleration;
    LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* surface;

    typedef OPTIX_RENDERING_OBJECT<T> BASE;
    // using BASE::Set_Transform;using BASE::Update_Transform;
public:
    OPTIX_RENDERING_IMPLICIT_SURFACE(const std::string filename, OPTIX_MATERIAL<T>* material):OPTIX_RENDERING_OBJECT<T>("Implicit Surface", material) {
        Initialize_Geometry_Particle();
        Initialize_Read_Write_General_Structures();

         try{
            FILE_UTILITIES::Create_From_File<T>(filename,surface);
            /*
            phi=&surface->levelset.phi;
            std::cout<<filename<<" statistics:"<<std::endl;
            std::cout<<"  grid = "<<surface->levelset.grid<<std::endl;
            std::cout<<"  phi array bounds = "<<phi->domain.min_corner.x<<" "<<phi->domain.max_corner.x<<", "<<phi->domain.min_corner.y<<" "<<phi->domain.max_corner.y<<", "<<phi->domain.min_corner.z<<" "<<phi->domain.max_corner.z<<std::endl;
            // checking that all the surface lay inside the domain
            for(int i=phi->domain.min_corner.x;i<=phi->domain.max_corner.x;i++)for(int j=phi->domain.min_corner.y;j<=phi->domain.max_corner.y;j++)if((*phi)(i,j,phi->domain.min_corner.z)<=0){std::cout<<"  phi<=0 on domain.min_corner.z"<<std::endl;goto check1_end;}check1_end:
            for(int i=phi->domain.min_corner.x;i<=phi->domain.max_corner.x;i++)for(int j=phi->domain.min_corner.y;j<=phi->domain.max_corner.y;j++)if((*phi)(i,j,phi->domain.max_corner.z)<=0){std::cout<<"  phi<=0 on domain.max_corner.z"<<std::endl;goto check2_end;}check2_end:
            for(int i=phi->domain.min_corner.x;i<=phi->domain.max_corner.x;i++)for(int k=phi->domain.min_corner.z;k<=phi->domain.max_corner.z;k++)if((*phi)(i,phi->domain.min_corner.y,k)<=0){std::cout<<"  phi<=0 on domain.min_corner.y"<<std::endl;goto check3_end;}check3_end:
            for(int i=phi->domain.min_corner.x;i<=phi->domain.max_corner.x;i++)for(int k=phi->domain.min_corner.z;k<=phi->domain.max_corner.z;k++)if((*phi)(i,phi->domain.max_corner.y,k)<=0){std::cout<<"  phi<=0 on domain.max_corner.y"<<std::endl;goto check4_end;}check4_end:
            for(int j=phi->domain.min_corner.y;j<=phi->domain.max_corner.y;j++)for(int k=phi->domain.min_corner.z;k<=phi->domain.max_corner.z;k++)if((*phi)(phi->domain.min_corner.x,j,k)<=0){std::cout<<"  phi<=0 on domain.min_corner.x"<<std::endl;goto check5_end;}check5_end:
            for(int j=phi->domain.min_corner.y;j<=phi->domain.max_corner.y;j++)for(int k=phi->domain.min_corner.z;k<=phi->domain.max_corner.z;k++)if((*phi)(phi->domain.max_corner.x,j,k)<=0){std::cout<<"  phi<=0 on domain.max_corner.x"<<std::endl;goto check6_end;}check6_end:
            int i = 0;
            */
        }
        catch(FILESYSTEM_ERROR&){}
    }

    ~OPTIX_RENDERING_IMPLICIT_SURFACE() {}

    Geometry getGeometryObject() {
        return rt_object;
    }

    Acceleration getAcceleration() {
        return rt_acceleration;
    }

    void InitializeGeometry(OPTIX_RENDER_WORLD<T>& world_input) {
        Context rt_context = world_input.RTContext();

        ARRAY<T,VECTOR<int,3> >* phi=&surface->levelset.phi;

        // creating 3D phi texture
        int dx = (phi->domain.max_corner.x - phi->domain.min_corner.x),
            dy = (phi->domain.max_corner.y - phi->domain.min_corner.y),
            dz = (phi->domain.max_corner.z - phi->domain.min_corner.z);
        Buffer phi_tex_buffer = rt_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_FLOAT, dx + 1, dy + 1, dz + 1);
        TextureSampler phi_map_tex_sampler = rt_context->createTextureSampler();

        phi_map_tex_sampler->setWrapMode(0, RT_WRAP_CLAMP_TO_EDGE);
        phi_map_tex_sampler->setWrapMode(1, RT_WRAP_CLAMP_TO_EDGE);
        phi_map_tex_sampler->setWrapMode(2, RT_WRAP_CLAMP_TO_EDGE);
        phi_map_tex_sampler->setFilteringModes(RT_FILTER_LINEAR, RT_FILTER_LINEAR, RT_FILTER_NONE);
        phi_map_tex_sampler->setIndexingMode(RT_TEXTURE_INDEX_NORMALIZED_COORDINATES);
        phi_map_tex_sampler->setReadMode(RT_TEXTURE_READ_ELEMENT_TYPE);
        phi_map_tex_sampler->setMaxAnisotropy(1.0f);
        phi_map_tex_sampler->setMipLevelCount(1);
        phi_map_tex_sampler->setArraySize(1);
        phi_map_tex_sampler->setBuffer(0, 0, phi_tex_buffer);

        float* phi_tex_buffer_data = static_cast<float*>(phi_tex_buffer->map());
        for(int k=0;k<=dz;k++)
            for(int j=0;j<=dy;j++) {
                for(int i=0;i<=dx;i++) {
                    int index = i + (j + k * (dy + 1)) * (dx + 1);
                    /*
                    float x = i / (float)dx * 2 - 1;
                    float y = j / (float)dy * 2 - 1;
                    float z = k / (float)dz * 2 - 1;
                    phi_tex_buffer_data[index] = x*x * 0.5 + y*y*0.5 + z*z*0.5 - 1;
                    */
                    phi_tex_buffer_data[index] = (*phi)(i + phi->domain.min_corner.x,j + phi->domain.min_corner.y, k + phi->domain.min_corner.z);
                    /*
                    if (k == dz / 2)
                        std::cout << std::setprecision(2) << phi_tex_buffer_data[index] << " ";
                        */
                }
                /*
                if (k == dz / 2)
                    std::cout << "\n";
                    */
            }
        // std::cout << "\n";
        phi_tex_buffer->unmap();

        rt_object = rt_context->createGeometry();
        std::string path_to_ptx = world_input.shader_prefix + "_generated_OPTIX_LEVELSET_SURFACE.cu.ptx";

        rt_object_intersection_program = rt_context->createProgramFromPTXFile(path_to_ptx, "intersect");
        rt_object->setIntersectionProgram(rt_object_intersection_program);
        rt_object_bounding_box_program = rt_context->createProgramFromPTXFile(path_to_ptx, "box_bounds");
        rt_object->setBoundingBoxProgram(rt_object_bounding_box_program);

        rt_object->setPrimitiveCount(1);
        rt_acceleration = rt_context->createAcceleration("Sbvh","Bvh");
        rt_object["phi_tex"]->setTextureSampler(phi_map_tex_sampler);
        rt_object["phi_tex_delta"]->setFloat(1.f / dx, 1.f / dy, 1.f / dz);
    }
};
}
#endif
#endif

