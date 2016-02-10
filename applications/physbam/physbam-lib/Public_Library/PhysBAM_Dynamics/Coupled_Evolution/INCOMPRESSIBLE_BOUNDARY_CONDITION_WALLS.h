//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS
//#####################################################################
#ifndef __INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS__
#define __INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS__
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS_BINARY_UNIFORM.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION.h>
namespace PhysBAM{
template<class TV>
class INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS:public IMPLICIT_BOUNDARY_CONDITION<TV>
{
    typedef VECTOR<int,TV::dimension> TV_INT;typedef typename TV::SCALAR T;
public:
    const VECTOR<VECTOR<bool,2>,TV::dimension>& walls;
    const VECTOR<VECTOR<bool,2>,TV::dimension> mpi_boundary;

    INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS(const VECTOR<VECTOR<bool,2>,TV::dimension>& walls_input,const VECTOR<VECTOR<bool,2>,TV::dimension>& mpi_boundary_input);
    virtual ~INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS();

    void Update_Boundary_Conditions(const GRID<TV>& grid,ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,ARRAY<T,TV_INT>& p,
        ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,const T time) PHYSBAM_OVERRIDE;
};
}
#endif
