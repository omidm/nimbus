//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE<TV>::
INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE(const ARRAY<T,TV_INT>& phi_input)
    :phi(phi_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE<TV>::
~INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE()
{
}
//#####################################################################
// Function Update_Boundary_Conditions
//#####################################################################
template<class TV> void INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE<TV>::
Update_Boundary_Conditions(const GRID<TV>& grid,ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,ARRAY<T,TV_INT>& p,
    ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,const T time)
{
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next())
        if(phi(iterator.Cell_Index())>0){
            psi_D(iterator.Cell_Index())=true;
            p(iterator.Cell_Index())=0;}
}
template class INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE<VECTOR<float,1> >;
template class INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE<VECTOR<float,2> >;
template class INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE<VECTOR<double,1> >;
template class INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE<VECTOR<double,2> >;
template class INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE<VECTOR<double,3> >;
#endif
