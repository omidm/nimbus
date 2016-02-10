//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_FLUIDS_FORCES
//#####################################################################
#ifndef __INCOMPRESSIBLE_FLUIDS_FORCES__
#define __INCOMPRESSIBLE_FLUIDS_FORCES__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class T_GRID>
class INCOMPRESSIBLE_FLUIDS_FORCES:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
public:

    INCOMPRESSIBLE_FLUIDS_FORCES()
    {}

    virtual ~INCOMPRESSIBLE_FLUIDS_FORCES()
    {}

//#####################################################################
    virtual void Add_Explicit_Forces(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
    {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}

    virtual void Add_Implicit_Forces_Before_Projection(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
    {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}

    virtual void Add_Implicit_Forces_Projection(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_velocities_ghost,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
    {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}

    virtual void Initialize_Grids(const T_GRID& grid)
    {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    
    virtual T CFL(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& face_velocities)
    {return 0;}
//#####################################################################
};
}
#endif
