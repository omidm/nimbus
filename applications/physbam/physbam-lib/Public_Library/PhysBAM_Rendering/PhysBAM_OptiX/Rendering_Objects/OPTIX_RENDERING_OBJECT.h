//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth, Bo Zhu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_RENDERING_OBJECT
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_RENDERING_OBJECT__
#define __OPTIX_RENDERING_OBJECT__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Objects/OPTIX_RENDERING_GEOMETRY.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Materials/OPTIX_MATERIAL.h>
#include <optix_world.h>
#include <optixu/optixpp_namespace.h>

namespace PhysBAM{
using namespace optix;
template<class T> class OPTIX_RENDER_WORLD;

template<class T>
class OPTIX_RENDERING_OBJECT : public OPTIX_RENDERING_GEOMETRY<T>
{
    typedef VECTOR<T,3> TV;
public:
    std::string name;
    OPTIX_RENDERING_OBJECT(std::string object_name,OPTIX_MATERIAL<T>* material);
    virtual ~OPTIX_RENDERING_OBJECT(){}
    
    void Set_Transform(const MATRIX<T,4>& A){transform=A;inverse_transform=transform.Inverse();}
    ARRAY<Transform>& getTransforms(){return rt_transform_array;}
    GeometryGroup getGeometryGroup(){return rt_geometrygroup;}
    OPTIX_MATERIAL<T>* getMaterial(){return material_shader;}
    void Update_Transform(const MATRIX<T,4>& A){transform=A*transform;inverse_transform=transform.Inverse();}

    virtual void InitializeGeometry(OPTIX_RENDER_WORLD<T>& world_input)=0; // is called strictly by OPTIX_RENDER_WORLD, passing reference to itself
    virtual void InitializeInstanceVariables(){}
    virtual Acceleration getAcceleration(){return NULL;}
    virtual Geometry getGeometryObject()=0;

    virtual void Initialize(OPTIX_RENDER_WORLD<T>* world_input);
    virtual void Update(T time){if(!frozen)rt_transform->setMatrix(true,transform.x,inverse_transform.x);}
    virtual Transform getTransform(){return rt_transform;}

protected:
    bool frozen;    ////only apply transform when frozen=false
    MATRIX<T,4> transform,inverse_transform; // transforms from standard position into actual position
    OPTIX_MATERIAL<T>* material_shader;
    OPTIX_RENDER_WORLD<T>* render_world;
    Transform rt_transform;
    ARRAY<Transform> rt_transform_array;
    GeometryInstance rt_instance;
    GeometryGroup rt_geometrygroup;
    Acceleration rt_acceleration;
};
}
#endif
#endif
