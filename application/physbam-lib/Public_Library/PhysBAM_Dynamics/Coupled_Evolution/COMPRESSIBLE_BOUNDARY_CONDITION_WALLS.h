//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMPRESSIBLE_BOUNDARY_CONDITION_WALLS
//#####################################################################
#ifndef __COMPRESSIBLE_BOUNDARY_CONDITION_WALLS__
#define __COMPRESSIBLE_BOUNDARY_CONDITION_WALLS__
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS_BINARY_UNIFORM.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION.h>
namespace PhysBAM{

template<class T_GRID> class FLUIDS_PARAMETERS_UNIFORM;

template<class TV>
class COMPRESSIBLE_BOUNDARY_CONDITION_WALLS:public IMPLICIT_BOUNDARY_CONDITION<TV>
{
    typedef VECTOR<int,TV::dimension> TV_INT;typedef typename TV::SCALAR T;
    typedef GRID<TV> T_GRID;
public:
    const VECTOR<VECTOR<bool,2>,TV::dimension>& walls;
    const VECTOR<VECTOR<bool,2>,TV::dimension> mpi_boundary;
    const FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters;

    COMPRESSIBLE_BOUNDARY_CONDITION_WALLS(const VECTOR<VECTOR<bool,2>,TV::dimension>& walls_input,const VECTOR<VECTOR<bool,2>,TV::dimension>& mpi_boundary_input,const FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters_input);
    virtual ~COMPRESSIBLE_BOUNDARY_CONDITION_WALLS();

    void Update_Boundary_Conditions(const GRID<TV>& grid,ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,ARRAY<T,TV_INT>& p,
        ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,const T time) PHYSBAM_OVERRIDE;
};
}
#endif
