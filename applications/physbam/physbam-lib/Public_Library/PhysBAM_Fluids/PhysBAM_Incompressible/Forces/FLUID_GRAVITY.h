//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_GRAVITY
//#####################################################################
#ifndef __FLUID_GRAVITY__
#define __FLUID_GRAVITY__

#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/INCOMPRESSIBLE_FLUIDS_FORCES.h>
namespace PhysBAM{

template<class T_GRID>
class FLUID_GRAVITY:public INCOMPRESSIBLE_FLUIDS_FORCES<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
public:
    T gravity;
    TV gravity_direction;

    FLUID_GRAVITY(T gravity=(T)9.8,TV gravity_direction=-TV::Axis_Vector(TV::m==1?1:2))
        :gravity(gravity),gravity_direction(gravity_direction)
    {}

    virtual ~FLUID_GRAVITY()
    {}

//#####################################################################
    void Add_Explicit_Forces(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time) PHYSBAM_OVERRIDE
    {for(int axis=1;axis<=T_GRID::dimension;axis++) if(gravity_direction[axis]) face_velocities.Component(axis)+=dt*gravity*gravity_direction[axis];}
    void Add_Implicit_Forces_Projection(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Initialize_Grids(const T_GRID& grid) PHYSBAM_OVERRIDE {}
    T CFL(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& face_velocities) PHYSBAM_OVERRIDE
    {return abs(gravity)*(gravity_direction/grid.DX()).L1_Norm();}
//#####################################################################
};
}
#endif
