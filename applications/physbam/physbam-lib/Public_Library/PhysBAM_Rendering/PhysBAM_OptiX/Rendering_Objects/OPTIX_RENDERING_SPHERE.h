//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_RENDERING_SPHERE
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_RENDERING_SPHERE__
#define __OPTIX_RENDERING_SPHERE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>

#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Objects/OPTIX_RENDERING_OBJECT.h>
#include <iostream>
#include <string>

#include <optix_world.h>
#include <optixu/optixpp_namespace.h>

#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Utilities/OPTIX_UTILITIES.h>

namespace PhysBAM{

using namespace optix;

template<class T> class OPTIX_RENDERING_OBJECT;
template<class T> class OPTIX_RENDER_WORLD;

template<class T> class OPTIX_RENDERING_SPHERE : public OPTIX_RENDERING_OBJECT<T> {
    typedef VECTOR<T,3> TV;
    TV my_center;
    T my_radius;
    Geometry rt_object;
    Program  rt_object_intersection_program;
    Program  rt_object_bounding_box_program;
public:

    OPTIX_RENDERING_SPHERE(TV center, T radius, OPTIX_MATERIAL<T>* material):
        OPTIX_RENDERING_OBJECT<T>("Sphere", material),my_center(center),my_radius(radius) {
    }

    ~OPTIX_RENDERING_SPHERE() {}

    Geometry getGeometryObject() {
        return rt_object;
    }

    void InitializeGeometry(OPTIX_RENDER_WORLD<T>& world_input) {
        Context rt_context = world_input.RTContext();

        rt_object = rt_context->createGeometry();
        rt_object->setPrimitiveCount(1u);

        std::string path_to_ptx = world_input.shader_prefix + "_generated_OPTIX_SPHERE.cu.ptx";

        rt_object_intersection_program = rt_context->createProgramFromPTXFile(path_to_ptx, "intersect");
        rt_object->setIntersectionProgram(rt_object_intersection_program);

        rt_object_bounding_box_program = rt_context->createProgramFromPTXFile(path_to_ptx, "box_bounds");
        rt_object->setBoundingBoxProgram(rt_object_bounding_box_program);

        rt_object["center"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(my_center));
        rt_object["radius"]->setFloat(my_radius);
    }
};
}
#endif
#endif

