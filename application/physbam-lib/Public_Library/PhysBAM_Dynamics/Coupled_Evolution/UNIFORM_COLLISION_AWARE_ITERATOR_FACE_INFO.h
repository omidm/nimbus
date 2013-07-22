//#####################################################################
// Copyright 2009, Nipun Kwatra, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO
//#####################################################################
#ifndef __UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO__
#define __UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO__
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COUPLING_CONSTRAINT_ID.h>

namespace PhysBAM{

template<class T_GRID> class GRID_BASED_COLLISION_GEOMETRY_UNIFORM;
template<class TV> class GRID;
template<class TV> class COLLISION_GEOMETRY;
template<class TV> class BOUNDARY_CONDITIONS_CALLBACKS;
template<class TV> class IMPLICIT_BOUNDARY_CONDITION_COLLECTION;

template<class TV>
struct COLLISION_FACE_INFO
{
    int axis;
    VECTOR<int,TV::dimension> index;
    int side;
    ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> > simplices;

    bool operator<(const COLLISION_FACE_INFO& cfi) const
    {
        if(axis!=cfi.axis) return axis<cfi.axis;
        for(int i=1;i<=TV::m;i++) if(index(i)!=cfi.index(i)) return index(i)<cfi.index(i);
        return false;
    }
};

template<class TV>
class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO:public NONCOPYABLE
{
    enum WORKAROUND {d=TV::dimension};
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::dimension> TV_INT;

    const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >& collision_bodies_affecting_fluid;

public:
    ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID> coupling_bodies;
    ARRAY<COLLISION_FACE_INFO<TV> > collision_face_info;
    GRID<TV>& grid;
    const ARRAY<bool,TV_INT>* outside_fluid;
    const bool& use_collision_face_neighbors;
    T iterator_rasterization_thickness;

    explicit UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO(const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >& collection);

//#####################################################################
    void Initialize_Collision_Aware_Face_Iterator(const ARRAY<bool,TV_INT>& outside_fluid_input,const ARRAY<bool,FACE_INDEX<d> >& kinematic_faces,int ghost_cells,const bool disable_thinshell);
    void Register_Neighbors_As_Collision_Faces();
//#####################################################################
};
}
#endif
