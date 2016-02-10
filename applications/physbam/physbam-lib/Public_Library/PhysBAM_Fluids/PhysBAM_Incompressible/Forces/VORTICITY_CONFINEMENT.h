//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VORTICITY_CONFINEMENT
//#####################################################################
#ifndef __VORTICITY_CONFINEMENT__
#define __VORTICITY_CONFINEMENT__

#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/INCOMPRESSIBLE_FLUIDS_FORCES.h>
namespace PhysBAM{
template<class T_GRID> class GRID_BASED_COLLISION_GEOMETRY_UNIFORM;

template<class T_GRID>
class VORTICITY_CONFINEMENT:public INCOMPRESSIBLE_FLUIDS_FORCES<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
public:
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>* collision_body_list;
    T_FACE_ARRAYS_BOOL* valid_mask;
    bool use_variable_vorticity_confinement;
    T vorticity_confinement;
    T_ARRAYS_SCALAR variable_vorticity_confinement;

    VORTICITY_CONFINEMENT(GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>* collision_body_list=0,T_FACE_ARRAYS_BOOL* valid_mask=0,const bool use_variable_vorticity_confinement=false,const T vorticity_confinement=.3);
    virtual ~VORTICITY_CONFINEMENT();

    void Set_Vorticity_Confinement(const T vorticity_confinement_input=.3)
    {vorticity_confinement=vorticity_confinement_input;}

//#####################################################################
    void Apply_Vorticity_Confinement_Force(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_velocities,T_ARRAYS_VECTOR& F);
    virtual void Compute_Vorticity_Confinement_Force(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,const T_FACE_ARRAYS_BOOL* valid_mask,T_ARRAYS_VECTOR& F);
    void Add_Explicit_Forces(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time) PHYSBAM_OVERRIDE;
    void Add_Implicit_Forces_Projection(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time) PHYSBAM_OVERRIDE;
    void Initialize_Grids(const T_GRID& grid) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
