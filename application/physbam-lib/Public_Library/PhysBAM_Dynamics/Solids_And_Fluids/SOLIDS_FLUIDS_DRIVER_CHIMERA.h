//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Andrew Selle, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_DRIVER_CHIMERA
//#####################################################################
#ifndef __SOLIDS_FLUIDS_DRIVER_CHIMERA__
#define __SOLIDS_FLUIDS_DRIVER_CHIMERA__    

#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_CHIMERA.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_PARAMETERS.h>
namespace PhysBAM{

template<class T_GRID> class SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA;

template<class T_GRID>
class SOLIDS_FLUIDS_DRIVER_CHIMERA:public SOLIDS_FLUIDS_DRIVER<typename T_GRID::VECTOR_T>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename TV::SCALAR T;typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename LEVELSET_POLICY<T_GRID>::FAST_LEVELSET_T T_FAST_LEVELSET;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;

    typedef SOLIDS_FLUIDS_DRIVER<TV> BASE;
    using BASE::time;
public:
    using BASE::project_at_frame_boundaries;using BASE::current_frame;using BASE::next_dt;using BASE::next_done;using BASE::Write_Time;using BASE::Write_First_Frame;
    using BASE::Write_Last_Frame;using BASE::Write_Substep;

    SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>& example;
    T last_dt;
    T restart_dt;
    bool reset_with_restart;

    SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>* evolution;

    SOLIDS_FLUIDS_DRIVER_CHIMERA(SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>& example_input);
    virtual ~SOLIDS_FLUIDS_DRIVER_CHIMERA();

//#####################################################################
    void Initialize() PHYSBAM_OVERRIDE;
    void Advance_To_Target_Time(const T target_time) PHYSBAM_OVERRIDE;
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE;
    T Compute_Dt(const T time,const T target_time,bool& done);
    void Write_Output_Files(const int frame) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
