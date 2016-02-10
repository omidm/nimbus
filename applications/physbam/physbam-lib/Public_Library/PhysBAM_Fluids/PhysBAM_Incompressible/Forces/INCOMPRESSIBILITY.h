//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBILITY
//#####################################################################
#ifndef __INCOMPRESSIBILITY__
#define __INCOMPRESSIBILITY__

#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/INCOMPRESSIBLE_FLUIDS_FORCES.h>
namespace PhysBAM{
template<class T_GRID> class PROJECTION_UNIFORM;

template<class T_GRID>
class INCOMPRESSIBILITY:public INCOMPRESSIBLE_FLUIDS_FORCES<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    PROJECTION_UNIFORM<T_GRID>& projection;
public:

    INCOMPRESSIBILITY(PROJECTION_UNIFORM<T_GRID>& projection_input);
    virtual ~INCOMPRESSIBILITY();

//#####################################################################
    void Add_Explicit_Forces(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_Implicit_Forces_Projection(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time) PHYSBAM_OVERRIDE;
    void Initialize_Grids(const T_GRID& grid) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
