//#####################################################################
// Copyright 2002-2011, Mridul Aanjaneya, Doug Enright, Ronald Fedkiw, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMPRESSIBLE_BOUNDARY_UNIFORM
//#####################################################################
#ifndef __COMPRESSIBLE_BOUNDARY_UNIFORM__
#define __COMPRESSIBLE_BOUNDARY_UNIFORM__

#include <PhysBAM_Tools/Parallel_Computation/MPI_GRID_POLICY.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class T_GRID>
class COMPRESSIBLE_BOUNDARY_UNIFORM : public BOUNDARY_UNIFORM<T_GRID,VECTOR<typename T_GRID::SCALAR,T_GRID::dimension+2> >
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT; typedef typename T_GRID::VECTOR_T TV;typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;
    typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,T_GRID::dimension> TV_SIDES;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;typedef typename T_ARRAYS_BASE::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_BASE;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef VECTOR<T,2*T_GRID::dimension> T_FACE_VECTOR;typedef VECTOR<TV,2*T_GRID::dimension> TV_FACE_VECTOR;
    typedef VECTOR<bool,2*T_GRID::dimension> T_FACE_VECTOR_BOOL;
    typedef VECTOR<TV_DIMENSION,2*T_GRID::dimension> TV_DIMENSION_FACE_VECTOR;
    typedef typename MPI_GRID_POLICY<T_GRID>::MPI_GRID T_MPI_GRID;
    enum {d=TV_DIMENSION::m};
  public:
    typedef BOUNDARY_UNIFORM<T_GRID,TV_DIMENSION> BASE;
    using BASE::Set_Constant_Extrapolation;using BASE::Constant_Extrapolation;using BASE::Fill_Single_Ghost_Region;

    EOS<T>* eos;
    VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>& eigensystem;
    T_MPI_GRID* mpi_grid;
    bool attenuate_using_riemann_invariants;
    T inflow_attenuation;
    bool always_attenuate;
    T_FACE_VECTOR linear_attenuations;
    T_FACE_VECTOR_BOOL linear_attenuation_faces;
    T_FACE_VECTOR S_far_field,iL_far_field,iR_far_field;
    TV_DIMENSION_FACE_VECTOR U_far_field;

    COMPRESSIBLE_BOUNDARY_UNIFORM(EOS<T>* eos_input,VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>& eigensystem_input,T_MPI_GRID* mpi_grid_input,
        const T_FACE_VECTOR rho_far_field,const T_FACE_VECTOR p_far_field,const TV_FACE_VECTOR velocity_far_field,const T inflow_attenuation_input=(T)1,
        const TV_SIDES& constant_extrapolation=TV_SIDES::Constant_Vector(TV_BOOL2::Constant_Vector(true)),
        const bool always_attenuate_input=false,const T_FACE_VECTOR linear_attenuations_input=T_FACE_VECTOR(),
        const T_FACE_VECTOR_BOOL linear_attenuation_faces_input=T_FACE_VECTOR_BOOL());
//#####################################################################
    void Attenuate_To_Far_Field_Values_Using_Riemann_Invariants(const T_ARRAYS_DIMENSION_BASE& u_ghost,const TV_INT& node_index,const int side,TV_DIMENSION &U,const T dt) const;
    void Attenuate_To_Far_Field_Values_Using_Characteristics(const T_ARRAYS_DIMENSION_BASE& u_ghost,const TV_INT& node_index,const int side,TV_DIMENSION &U,const T dt) const;
    void Attenuate_To_Far_Field_Values(const T_ARRAYS_DIMENSION_BASE& u_ghost,const TV_INT& node_index,const int side,TV_DIMENSION &U,const T dt) const;
    void Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_DIMENSION_BASE& u,T_ARRAYS_DIMENSION_BASE& u_ghost,const T dt,const T time,const int number_of_ghost_cells=3) PHYSBAM_OVERRIDE;
    void Fill_Single_Ghost_Region(const T_GRID& grid,T_ARRAYS_DIMENSION_BASE& u_ghost,const RANGE<TV_INT>& region,const int side,const T dt,const T time,const int number_of_ghost_cells) const PHYSBAM_OVERRIDE;
    void Apply_Boundary_Condition_Single_Side(const T_GRID& grid,T_ARRAYS_DIMENSION_BASE& u,const int side,const T time) const;
    void Apply_Boundary_Condition(const T_GRID& grid,T_ARRAYS_DIMENSION_BASE& u,const T time) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
