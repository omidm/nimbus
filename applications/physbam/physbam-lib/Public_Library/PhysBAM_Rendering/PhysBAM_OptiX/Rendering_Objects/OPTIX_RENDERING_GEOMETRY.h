//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_RENDERING_GEOMETRY
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_RENDERING_GEOMETRY__
#define __OPTIX_RENDERING_GEOMETRY__

#include <string>
#include <optix_world.h>
#include <optixu/optixpp_namespace.h>

namespace PhysBAM{
using namespace optix;
template<class T> class OPTIX_RENDER_WORLD;

template<class T>
class OPTIX_RENDERING_GEOMETRY:public NONCOPYABLE {
public:
    std::string name;

    OPTIX_RENDERING_GEOMETRY(std::string name_input):name(name_input){}
    virtual ~OPTIX_RENDERING_GEOMETRY(){}
    virtual void Initialize(OPTIX_RENDER_WORLD<T>* world_input)=0; // is called strictly by OPTIX_RENDER_WORLD, passing reference to itself
    virtual void Update(T time)=0;
    virtual ARRAY<Transform>& getTransforms()=0;
    virtual void BeforeLaunch(){}
    virtual void PrepareForDraw(){}
};
}
#endif
#endif
