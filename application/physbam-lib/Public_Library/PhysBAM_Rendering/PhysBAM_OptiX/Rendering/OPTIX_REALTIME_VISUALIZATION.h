//#####################################################################
// Copyright 2012, Bo Zhu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPTIX_REALTIME_VISUALIZATION
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_REALTIME_VISUALIZATION__
#define __OPTIX_REALTIME_VISUALIZATION__
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering/OPTIX_RENDER_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WINDOW_GLUT.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Objects/OPTIX_RENDERING_AABB.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Objects/OPTIX_RENDERING_SPHERE.h>
#include <PhysBAM_Rendering/PhysBAM_Optix/Rendering_Objects/OPTIX_RENDERING_PLANE.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Materials/OPTIX_MATERIAL_LAMBERTIAN.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Materials/OPTIX_MATERIAL.h>
#include <PhysBAM_Rendering/PhysBAM_OptiX/Rendering_Lights/OPTIX_RENDERING_LIGHT.h>

namespace PhysBAM
{
template<class T>
class OPTIX_REALTIME_VISUALIZATION
{
public:
    OPTIX_REALTIME_VISUALIZATION(){}

    virtual void Initialize_Window()
    {
        OPTIX_RENDER_WORLD<T>::Instance()->window=new OPENGL_WINDOW_GLUT(*OPTIX_RENDER_WORLD<T>::Instance(),"PhysBAM OptiX",
            OPTIX_RENDER_WORLD<T>::Instance()->camera->screen_width,OPTIX_RENDER_WORLD<T>::Instance()->camera->screen_height);
    }

    virtual void Initialize_Scene(){}
    
    virtual void Initialize_Renderer(){OPTIX_RENDER_WORLD<T>::Instance()->Initialize();}
    
    virtual void Initialize()
    {
        Initialize_Window();
        Initialize_Scene();
        Initialize_Renderer();
    }

    virtual void Run(){OPTIX_RENDER_WORLD<T>::Instance()->Render_Frame();}
    virtual void Update_Data(){}
};
}

#endif
#endif