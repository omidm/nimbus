//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_RENDERING_TRIANGULATED_SURFACE
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_RENDERING_TRIANGULATED_SURFACE__
#define __OPTIX_RENDERING_TRIANGULATED_SURFACE__

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

#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Objects/OPTIX_RENDERING_OBJECT.h>

#include <iostream>
#include <string>

#include <optix_world.h>
#include <optixu/optixpp_namespace.h>

#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Utilities/OPTIX_UTILITIES.h>

namespace PhysBAM{

using namespace optix;

template<class T> class OPTIX_RENDERING_TRIANGULATED_SURFACE : public OPTIX_RENDERING_OBJECT<T> {
    typedef VECTOR<T,3> TV;
    GENERIC_PARSER<T> *parser;
    Geometry rt_object;
    Program  rt_object_intersection_program;
    Program  rt_object_bounding_box_program;
    TRIANGULATED_SURFACE<T>* surface;
    Acceleration rt_acceleration;
    bool my_interpolate_normal;

    typedef OPTIX_RENDERING_OBJECT<T> BASE;
    // using BASE::Set_Transform;using BASE::Update_Transform;
    void Init() {
    }
public:
    OPTIX_RENDERING_TRIANGULATED_SURFACE(TRIANGULATED_SURFACE<T>* _surface, OPTIX_MATERIAL<T>* material, bool interpolate_normal = false):
        OPTIX_RENDERING_OBJECT<T>("Triangulated Surface", material), my_interpolate_normal(interpolate_normal) {
        surface = _surface;
    }

    OPTIX_RENDERING_TRIANGULATED_SURFACE(std::string file_name, OPTIX_MATERIAL<T>* material):OPTIX_RENDERING_OBJECT<T>("Triangulated Surface", material) {
        // std::cout << "\nInitializing triangulated surface\n";
        Initialize_Geometry_Particle();
        Initialize_Read_Write_General_Structures();

        // std::cout << "Reading surface from file\n";
        FILE_UTILITIES::Create_From_File<T>(file_name,surface);

        for(int t=1;t<=surface->mesh.elements.m;t++){
            int i,j,k;
            surface->mesh.elements(t).Get(i,j,k);
            TV v1 = surface->particles.X(i), v2 = surface->particles.X(j), v3 = surface->particles.X(k);
        }
    }

    ~OPTIX_RENDERING_TRIANGULATED_SURFACE() {}

    Geometry getGeometryObject(){return rt_object;}

    Acceleration getAcceleration(){return rt_acceleration;}

    void InitializeGeometry(OPTIX_RENDER_WORLD<T>& world_input) 
    {
        Context rt_context = world_input.RTContext();

        rt_object = rt_context->createGeometry();

        std::string path_to_ptx = world_input.shader_prefix + "_generated_OPTIX_TRIANGULATED_SURFACE.cu.ptx";

        rt_object_intersection_program = rt_context->createProgramFromPTXFile(path_to_ptx, "intersect");
        rt_object->setIntersectionProgram(rt_object_intersection_program);

        rt_object_bounding_box_program = rt_context->createProgramFromPTXFile(path_to_ptx, "box_bounds");
        rt_object->setBoundingBoxProgram(rt_object_bounding_box_program);

        ARRAY<VECTOR<T,3> > vertex_normals;
        vertex_normals.Resize(surface->particles.array_collection->Size());
        ARRAYS_COMPUTATIONS::Fill(vertex_normals,VECTOR<T,3>());
        for(int t=1;t<=surface->mesh.elements.m;t++){
            int i,j,k;surface->mesh.elements(t).Get(i,j,k);
            VECTOR<T,3> normal=TRIANGLE_3D<T>::Normal(surface->particles.X(i),surface->particles.X(j),surface->particles.X(k));
            vertex_normals(i)+=normal;vertex_normals(j)+=normal;vertex_normals(k)+=normal;}
        for(int p=1;p<=surface->particles.array_collection->Size();p++)vertex_normals(p).Normalize();


        unsigned int vertexNumber = surface->particles.X.m;
        // unsigned int vertexNumber = surface->mesh.elements.m * 3;
        Buffer rt_vertex_buffer = rt_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_FLOAT3, vertexNumber);
        float3* vbuffer_data = static_cast<float3*>(rt_vertex_buffer->map());

        Buffer rt_normal_buffer = rt_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_FLOAT3, vertexNumber);
        float3* nbuffer_data = static_cast<float3*>(rt_normal_buffer->map());
        for (unsigned int i = 0; i < vertexNumber; i++) {
            TV v = surface->particles.X(i + 1);
            vbuffer_data[i] = make_float3(v.x, v.y, v.z);
            TV n = vertex_normals(i + 1);
            nbuffer_data[i] = make_float3(n.x, n.y, n.z);
        }

        //rt_object["vertex_buffer"]->set(rt_vertex_buffer);

        unsigned int trianglesNumber = surface->mesh.elements.m;
        // std::cout << "trianglesNumber = "  << trianglesNumber << "\n";
        Buffer rt_index_buffer = rt_context->createBuffer(RT_BUFFER_INPUT, RT_FORMAT_UNSIGNED_INT3, trianglesNumber);
        uint3* ibuffer_data = static_cast<uint3*>(rt_index_buffer->map());
        /*for (unsigned int i = 0; i < trianglesNumber; i++) {
            int t,j,k;
            surface->mesh.elements(i + 1).Get(t, j, k);
            ibuffer_data[i] = make_uint3(i*3, i*3+1, i*3+2);

            TV v0 = surface->particles.X(t);
            TV v1 = surface->particles.X(j);
            TV v2 = surface->particles.X(k);
            // std::cout << "(" << v.x << " " << v.y << " " << v.z << ") ";
            vbuffer_data[i * 3] = make_float3(v0.x, v0.y, v0.z);
            vbuffer_data[i * 3 + 1] = make_float3(v1.x, v1.y, v1.z);
            vbuffer_data[i * 3 + 2] = make_float3(v2.x, v2.y, v2.z);

            TV n0 = vertex_normals(t);
            TV n1 = vertex_normals(j);
            TV n2 = vertex_normals(k);
            nbuffer_data[i * 3] = make_float3(n0.x, n0.y, n0.z);
            nbuffer_data[i * 3 + 1] = make_float3(n1.x, n1.y, n1.z);
            nbuffer_data[i * 3 + 2] = make_float3(n2.x, n2.y, n2.z);
        }*/
        for (unsigned int i = 0; i < trianglesNumber; i++) {
            int t,j,k;
            surface->mesh.elements(i + 1).Get(t, j, k);
            ibuffer_data[i] = make_uint3(t - 1, j - 1, k - 1);
/*
            TV v0 = surface->particles.X(t);
            TV v1 = surface->particles.X(j);
            TV v2 = surface->particles.X(k);
            // std::cout << "(" << v.x << " " << v.y << " " << v.z << ") ";
            vbuffer_data[i * 3] = make_float3(v0.x, v0.y, v0.z);
            vbuffer_data[i * 3 + 1] = make_float3(v1.x, v1.y, v1.z);
            vbuffer_data[i * 3 + 2] = make_float3(v2.x, v2.y, v2.z);

            TV n0 = vertex_normals(t);
            TV n1 = vertex_normals(j);
            TV n2 = vertex_normals(k);
            nbuffer_data[i * 3] = make_float3(n0.x, n0.y, n0.z);
            nbuffer_data[i * 3 + 1] = make_float3(n1.x, n1.y, n1.z);
            nbuffer_data[i * 3 + 2] = make_float3(n2.x, n2.y, n2.z);
            */
        }
        rt_index_buffer->unmap();
        rt_vertex_buffer->unmap();
        rt_normal_buffer->unmap();

        rt_object["vertex_buffer"]->set(rt_vertex_buffer);
        rt_object["vindex_buffer"]->set(rt_index_buffer);
        rt_object["normal_buffer"]->set(rt_normal_buffer);
        rt_object["interpolate_normal"]->setInt(my_interpolate_normal);
        rt_object->setPrimitiveCount(surface->mesh.elements.m);

        // rt_acceleration = rt_context->createAcceleration("Sbvh","Bvh");
        rt_acceleration = rt_context->createAcceleration("MedianBvh","Bvh");
        // rt_acceleration = rt_context->createAcceleration("NoAccel","NoAccel");

        rt_acceleration->setProperty( "vertex_buffer_name", "vertex_buffer" );
        // rt_acceleration->setProperty( "vertex_buffer_stride", sizeof(float3));
        rt_acceleration->setProperty( "index_buffer_name", "vindex_buffer" );
        // rt_acceleration->setProperty( "index_buffer_stride", sizeof(uint3));


        /**************************************/
        /*
        if(!surface->mesh.incident_elements)
            surface->mesh.Initialize_Incident_Elements();
        for(int t=1;t<=(*surface->mesh.incident_elements)(0).m;t++){
            int triangle=(*surface->mesh.incident_elements)(0)(t);
            int i,j,k;surface->mesh.elements(triangle).Get(i,j,k);
            }
            */

    }
};
}
#endif
#endif

