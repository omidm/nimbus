//#####################################################################
// Copyright 2011, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_POLICY 
//#####################################################################
#ifndef __OPENGL_COMPONENT_POLICY__
#define __OPENGL_COMPONENT_POLICY__

#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Deformable_Components/OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Deformable_Components/OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Rigids_Components/OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Rigids_Components/OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Rigids_Components/OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_SCALAR_FIELD_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_SCALAR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_3D.h>

namespace PhysBAM{

template<class TV> struct OPENGL_COMPONENT_POLICY;
struct UNUSABLE{};

//#####################################################################
// 0D
//#####################################################################
template<class T>
struct OPENGL_COMPONENT_POLICY<VECTOR<T,0> >
{
};
//#####################################################################
// 1D
//#####################################################################
template<class T>
struct OPENGL_COMPONENT_POLICY<VECTOR<T,1> >
{
    typedef UNUSABLE DEFORMABLE_BODY_COLLECTION;
    typedef OPENGL_COMPONENT_RIGID_BODY_COLLECTION_1D<T> RIGID_BODY_COLLECTION;
    typedef OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T> RIGID_GEOMETRY_COLLECTION;
    typedef OPENGL_COMPONENT_SCALAR_FIELD_1D<T> SCALAR_FIELD;
    typedef UNUSABLE MAC_VELOCITY_FIELD;
    typedef OPENGL_COMPONENT_LEVELSET_1D<T> LEVELSET;
};
//#####################################################################
// 2D
//#####################################################################
template<class T>
struct OPENGL_COMPONENT_POLICY<VECTOR<T,2> >
{
    typedef OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_2D<T> DEFORMABLE_BODY_COLLECTION;
    typedef OPENGL_COMPONENT_RIGID_BODY_COLLECTION_2D<T> RIGID_BODY_COLLECTION;
    typedef OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T> RIGID_GEOMETRY_COLLECTION;
    typedef OPENGL_COMPONENT_SCALAR_FIELD_2D<T> SCALAR_FIELD;
    typedef OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T> MAC_VELOCITY_FIELD;
    typedef OPENGL_COMPONENT_LEVELSET_2D<T> LEVELSET;
};
//#####################################################################
// 3D
//#####################################################################
template<class T>
struct OPENGL_COMPONENT_POLICY<VECTOR<T,3> >
{
    typedef OPENGL_COMPONENT_DEFORMABLE_BODY_COLLECTION_3D<T> DEFORMABLE_BODY_COLLECTION;
    typedef OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T> RIGID_BODY_COLLECTION;
    typedef OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T> RIGID_GEOMETRY_COLLECTION;
    typedef OPENGL_COMPONENT_SCALAR_FIELD_3D<T> SCALAR_FIELD;
    typedef OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T> MAC_VELOCITY_FIELD;
    typedef OPENGL_COMPONENT_LEVELSET_3D<T> LEVELSET;
};
//#####################################################################
}
#endif
