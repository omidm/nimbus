//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLS_FSI_DRIVER
//#####################################################################
#ifndef __PLS_FSI_DRIVER__
#define __PLS_FSI_DRIVER__    

#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/DRIVER.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Driver/PLS_FSI_DRIVER.h>
#include <PhysBAM_Dynamics/Coupled_Driver/PLS_FSI_EXAMPLE.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_PARAMETERS.h>
namespace PhysBAM{

template<class TV>
class PLS_FSI_DRIVER:public DRIVER<TV>,public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef DRIVER<TV> BASE;
    typedef GRID<TV> T_GRID;typedef typename T_GRID::VECTOR_INT TV_INT;typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename LEVELSET_POLICY<T_GRID>::FAST_LEVELSET_T T_FAST_LEVELSET;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    using BASE::time;
public:
    using BASE::output_number;using BASE::Write_Output_Files;using BASE::Read_Time;using BASE::Write_Time;
    using BASE::Write_First_Frame;using BASE::Write_Last_Frame;using BASE::Write_Substep;

    PLS_FSI_EXAMPLE<TV>& example;
    int current_frame;

    PLS_FSI_DRIVER(PLS_FSI_EXAMPLE<TV>& example_input);
    virtual ~PLS_FSI_DRIVER();

//#####################################################################
    void Initialize() PHYSBAM_OVERRIDE;
    void Initialize_Fluids_Grids();
    void Advance_To_Target_Time(const T target_time) PHYSBAM_OVERRIDE;
    void First_Order_Time_Step(int substep,T dt);
    virtual void Postprocess_Frame(const int frame);
    virtual void Preprocess_Frame(const int frame);
    T Compute_Dt(const T time,const T target_time,bool& done);
    void Write_Output_Files(const int frame) PHYSBAM_OVERRIDE;
    void Advect_Fluid(const T dt,const int substep);
    void Execute_Main_Program() PHYSBAM_OVERRIDE;
    void Simulate_To_Frame(const int frame_input) PHYSBAM_OVERRIDE;
    void Delete_Particles_Inside_Objects(const T time);
    template<class T_PARTICLES> void Delete_Particles_Inside_Objects(ARRAY<T_PARTICLES*,TV_INT>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time);
    void Extrapolate_Velocity_Across_Interface(ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T_FAST_LEVELSET& phi,const T band_width);
    void Advance_Particles_With_PLS(T dt);
    void Extrapolate_Velocity_Across_Interface(T time,T dt);
//#####################################################################
};
}
#endif
