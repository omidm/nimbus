//#####################################################################
// Copyright 2011, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMPRESSIBLE_DRIVER
//#####################################################################
#ifndef __COMPRESSIBLE_DRIVER__
#define __COMPRESSIBLE_DRIVER__

#include <PhysBAM_Fluids/PhysBAM_Compressible/COMPRESSIBLE_EXAMPLE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/DRIVER.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM {

template<class TV>
class COMPRESSIBLE_DRIVER:public DRIVER<TV>
{
    typedef DRIVER<TV> BASE;
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::dimension> TV_INT;
    typedef VECTOR<T,TV::dimension+2> TV_DIMENSION;typedef GRID<TV> T_GRID;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
   
  public:
    using BASE::output_number;using BASE::current_frame;using BASE::time;using BASE::Write_Output_Files;using BASE::Write_Substep;using BASE::Read_Time;

    COMPRESSIBLE_EXAMPLE<TV>& example;
    T restart_dt;
    bool reset_with_restart;

    COMPRESSIBLE_DRIVER(COMPRESSIBLE_EXAMPLE<TV>& example_input);

    void Preprocess_Frame(const int frame)
    {example.Preprocess_Frame(frame);}

    void Postprocess_Frame(const int frame)
    {example.Postprocess_Frame(frame);}

    void Preprocess_Substep(const int frame,const int substep)
    {example.Preprocess_Substep(frame,substep);}

    void Postprocess_Substep(const T dt,const T time)
    {example.Postprocess_Substep(dt,time);}

//#####################################################################
    T Compute_Dt(const T time,const T target_time,bool& done);
//#####################################################################
    void Initialize() PHYSBAM_OVERRIDE;
    void Update_Bodies(const T dt,const T time);
    void Calculate_Maximum_Allowable_dt(const T dt,T& min_dt,const int substep,RUNGEKUTTA<T_ARRAYS_DIMENSION_SCALAR>& rungekutta_u);
    void Advance_To_Target_Time(const T target_time) PHYSBAM_OVERRIDE;
    void Simulate_To_Frame(const int frame) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif  // __COMPRESSIBLE_DRIVER__
