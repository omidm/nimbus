//#####################################################################
// Copyright 2002-2008, Ronald Fedkiw, Eran Guendelman, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_VISITOR
//##################################################################### 
#ifndef __PARTICLE_LEVELSET_VISITOR__
#define __PARTICLE_LEVELSET_VISITOR__

#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SPATIAL_ACCELERATION_FORWARD.h>
namespace PhysBAM{

template<class TV> struct RIGID_BODY_PARTICLE_INTERSECTION;
template<class TV> class IMPLICIT_OBJECT;
template<class TV>
struct PARTICLE_LEVELSET_VISITOR
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};

    ARRAY_VIEW<const RANGE<TV> > box_hierarchy; 
    ARRAY_VIEW<const T> box_radius;
    ARRAY_VIEW<const TV> X;
    ARRAY_VIEW<const ARRAY<int> > particles_in_group;
    const IMPLICIT_OBJECT<TV>& object_space_implicit_object;
    const MATRIX<T,d> rotation;
    const VECTOR<T,d> translation;
    const T contour_value;
    ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections;
    const int particle_body_id,levelset_body_id;

    PARTICLE_LEVELSET_VISITOR(const PARTICLE_HIERARCHY<TV>& particle_hierarchy,const IMPLICIT_OBJECT<TV>& object_space_implicit_object,const MATRIX<T,d>& rotation,
        const VECTOR<T,d>& translation,const T contour_value,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const int particle_body_id,
        const int levelset_body_id);

    ~PARTICLE_LEVELSET_VISITOR();

    bool Cull(const int box) const
    {TV center=rotation*box_hierarchy(box).Center()+translation;
    return object_space_implicit_object.Lazy_Outside_Extended_Levelset(center,box_radius(box)+contour_value);}

//#####################################################################
    void Store(const int box) const;
//#####################################################################
};
}
#endif
