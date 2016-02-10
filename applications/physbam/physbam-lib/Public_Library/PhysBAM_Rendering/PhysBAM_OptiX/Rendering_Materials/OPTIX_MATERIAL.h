//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth, Bo Zhu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_MATERIAL
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_MATERIAL__
#define __OPTIX_MATERIAL__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD.h>

#include <optix_world.h>
#include <optixu/optixpp_namespace.h>

namespace PhysBAM{

using namespace optix;

template<class T,int d> class VECTOR;
template<class T> class OPTIX_RENDER_WORLD;

template<class T>
class OPTIX_MATERIAL:public NONCOPYABLE
{
public:
    bool opaque;
    OPTIX_MATERIAL(bool opaque_input=true):opaque(opaque_input){}
    virtual ~OPTIX_MATERIAL(){}

    virtual Material getRTMaterial()=0;
    virtual void ensureInitialized(OPTIX_RENDER_WORLD<T> *world_input)=0;
};
}
#endif
#endif
