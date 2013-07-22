//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth, Bo Zhu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_RENDERING_AABB_BOX
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_RENDERING_AABB_BOX__
#define __OPTIX_RENDERING_AABB_BOX__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Objects/OPTIX_RENDERING_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Utilities/OPTIX_UTILITIES.h>

#include <optix_world.h>
#include <optixu/optixpp_namespace.h>

namespace PhysBAM{
using namespace optix;
template<class T> class OPTIX_RENDERING_OBJECT;
template<class T> class OPTIX_RENDER_WORLD;

template<class T> 
class OPTIX_RENDERING_AABB_BOX : public OPTIX_RENDERING_OBJECT<T> 
{
    typedef VECTOR<T,3> TV;
public:
    OPTIX_RENDERING_AABB_BOX(TV lower_corner,TV upper_corner,OPTIX_MATERIAL<T>* material):
        OPTIX_RENDERING_OBJECT<T>("AABB_BOX",material),low_corner(lower_corner),up_corner(upper_corner),wasInitialized(false){}
    OPTIX_RENDERING_AABB_BOX(const RANGE<TV>& range,OPTIX_MATERIAL<T>* material):
        OPTIX_RENDERING_OBJECT<T>("AABB_BOX",material),low_corner(range.Minimum_Corner()),up_corner(range.Maximum_Corner()),wasInitialized(false){}
    ~OPTIX_RENDERING_AABB_BOX(){}

    virtual Geometry getGeometryObject(){return rt_object;}

    virtual void InitializeGeometry(OPTIX_RENDER_WORLD<T>& world_input) 
    {
        Context rt_context=world_input.RTContext();
        wasInitialized=true;

        rt_object=rt_context->createGeometry();
        rt_object->setPrimitiveCount(1u);

        std::string path_to_ptx=world_input.shader_prefix+"_generated_OPTIX_AABB_BOX.cu.ptx";

        rt_object_intersection_program=rt_context->createProgramFromPTXFile(path_to_ptx,"intersect");
        rt_object->setIntersectionProgram(rt_object_intersection_program);

        rt_object_bounding_box_program=rt_context->createProgramFromPTXFile(path_to_ptx,"box_bounds");
        rt_object->setBoundingBoxProgram(rt_object_bounding_box_program);

        
        rt_object["low_corner"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(low_corner));
        rt_object["up_corner"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(up_corner));

    }
    virtual void InitializeInstanceVariables()
    {
        rt_instance["low_corner"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(low_corner));
        rt_instance["up_corner"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(up_corner));
        rt_instance["one_over_box_size"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(TV::Constant_Vector(1.0)/(up_corner-low_corner)));
    }

    TV getLowCorner(){return low_corner;}
    TV getUpCorner(){return up_corner;}
    void setLowCorner(TV _low_corner){low_corner=_low_corner;if (wasInitialized){rt_object["low_corner"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(low_corner));}}
    void setUpCorner(TV _up_corner){up_corner=_up_corner;if(wasInitialized)rt_object["up_corner"]->setFloat(OPTIX_UTILITIES::Get_Float3<T>(up_corner));}
private:
    TV low_corner,up_corner;
    bool wasInitialized;
    Geometry rt_object;
    Program  rt_object_intersection_program;
    Program  rt_object_bounding_box_program;
};
}
#endif
#endif

