//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_BOUNDARY_CONDITION_COLLISIONS
//#####################################################################
#ifndef __IMPLICIT_BOUNDARY_CONDITION_COLLISIONS__
#define __IMPLICIT_BOUNDARY_CONDITION_COLLISIONS__
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_CALLBACKS.h>
namespace PhysBAM{
template<class TV>
class IMPLICIT_BOUNDARY_CONDITION_COLLISIONS:public IMPLICIT_BOUNDARY_CONDITION<TV>
{
    typedef VECTOR<int,TV::dimension> TV_INT;typedef typename TV::SCALAR T;
    typedef GRID<TV> T_GRID;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX T_SIMPLEX;

    COLLISION_GEOMETRY_COLLECTION<TV>& collision_geometry_collection;
    bool use_implicit_geometry;
    FLUIDS_PARAMETERS_CALLBACKS<T_GRID>& callbacks;

public:
    IMPLICIT_BOUNDARY_CONDITION_COLLISIONS(COLLISION_GEOMETRY_COLLECTION<TV>& collision_geometry_collection_input,
        const bool use_implicit_geometry_input, FLUIDS_PARAMETERS_CALLBACKS<T_GRID>& callbacks);

    virtual ~IMPLICIT_BOUNDARY_CONDITION_COLLISIONS();

//#####################################################################
    void Update_Boundary_Conditions(const GRID<TV>& grid,ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,ARRAY<T,TV_INT>& p,
        ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,const T time) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
