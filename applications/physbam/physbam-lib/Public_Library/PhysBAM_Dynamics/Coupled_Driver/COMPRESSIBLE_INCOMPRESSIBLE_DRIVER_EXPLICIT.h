//#####################################################################
// Copyright 2012, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT
//#####################################################################
#ifndef __COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT__
#define __COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT__

#include <PhysBAM_Dynamics/Coupled_Driver/COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/DRIVER.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM {

template<class TV>
class COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT:public DRIVER<TV>
{
    typedef DRIVER<TV> BASE;
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::dimension> TV_INT;
    typedef VECTOR<T,TV::dimension+2> TV_DIMENSION;typedef GRID<TV> T_GRID;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;
   
  public:
    using BASE::output_number;using BASE::current_frame;using BASE::time;using BASE::Write_Output_Files;using BASE::Write_Substep;using BASE::Read_Time;

    COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>& example;

    COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT(COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>& example_input);

//#####################################################################
    T Compute_Dt(const T time,const T target_time);
    void Initialize() PHYSBAM_OVERRIDE;
    void Run(RANGE<TV_INT>& domain,const T dt,const T time);
    void Advance_Levelset(const T dt,const T time);
    void Advance_One_Time_Step_Explicit_Part(const T dt,const T time);
    void Advance_One_Time_Step_Implicit_Part(const T dt,const T time);
    void Advance_To_Target_Time(const T target_time) PHYSBAM_OVERRIDE;
    void Simulate_To_Frame(const int frame) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif  // __COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT__
