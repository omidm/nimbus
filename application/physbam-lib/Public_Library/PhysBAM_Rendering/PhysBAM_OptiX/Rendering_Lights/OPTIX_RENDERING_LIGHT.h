//#####################################################################
// Copyright 2012, Bo Zhu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_RENDERING_LIGHT__
#define __OPTIX_RENDERING_LIGHT__

namespace PhysBAM{
template<class T> class OPTIX_RENDER_WORLD;

template<class T>
class OPTIX_RENDERING_LIGHT 
{
public:
    OPTIX_RENDERING_LIGHT(){}
    virtual void Add_To_Optix_Render_World(OPTIX_RENDER_WORLD<T>* optix_render_world)=0;
};
}
#endif
#endif
