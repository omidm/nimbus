//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_RENDERING_OBJECT
//#####################################################################
#ifdef USE_OPTIX
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD.h>
#include "OPTIX_RENDERING_OBJECT.h"

using namespace PhysBAM;
using namespace optix;

template<class T> OPTIX_RENDERING_OBJECT<T>::OPTIX_RENDERING_OBJECT(std::string object_name, OPTIX_MATERIAL<T>* material)
    :OPTIX_RENDERING_GEOMETRY<T>(object_name),frozen(false),transform(MATRIX<T,4>::Identity_Matrix()),inverse_transform(MATRIX<T,4>::Identity_Matrix()),material_shader(material),render_world(NULL)
{
}

template<class T> void OPTIX_RENDERING_OBJECT<T>::Initialize(OPTIX_RENDER_WORLD<T>* world_input) 
{
    render_world=world_input;
    InitializeGeometry(*render_world);
    Context rt_context=world_input->RTContext();
    Geometry rt_object=getGeometryObject();
    material_shader->ensureInitialized(render_world);

    // assigning material to object - creating instance
    Material rt_material=material_shader->getRTMaterial();

    rt_instance=rt_context->createGeometryInstance();
    rt_instance->setMaterialCount(1);
    rt_instance->setMaterial(0,rt_material);
    rt_instance->setGeometry(rt_object);
    InitializeInstanceVariables();
    // putting instance into a group
    rt_geometrygroup=rt_context->createGeometryGroup();
    rt_geometrygroup->setChildCount(1);
    rt_geometrygroup->setChild(0,rt_instance);

    // setting empty acceleration structure
    rt_acceleration=getAcceleration();
    if (rt_acceleration==NULL)rt_acceleration=rt_context->createAcceleration("NoAccel","NoAccel");
    rt_geometrygroup->setAcceleration(rt_acceleration);
    rt_acceleration->markDirty();

    // adding transform to a geometry group
    rt_transform=rt_context->createTransform();
    rt_transform->setChild(rt_geometrygroup);
    rt_transform->setMatrix(true,transform.x,inverse_transform.x);

    rt_transform_array.Resize(1);
    rt_transform_array(1)=rt_transform;
}

template class OPTIX_RENDERING_OBJECT<float>;
// template class OPTIX_RENDERING_OBJECT<double>;
#endif
