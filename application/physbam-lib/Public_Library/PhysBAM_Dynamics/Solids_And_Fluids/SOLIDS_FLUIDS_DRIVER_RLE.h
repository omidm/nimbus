#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_DRIVER_RLE
//#####################################################################
#ifndef __SOLIDS_FLUIDS_DRIVER_RLE__
#define __SOLIDS_FLUIDS_DRIVER_RLE__

#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_RLE.h>
namespace PhysBAM{

template<class T_GRID>
class SOLIDS_FLUIDS_DRIVER_RLE:public SOLIDS_FLUIDS_DRIVER<typename T_GRID::VECTOR_T>
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::BLOCK_ITERATOR BLOCK_ITERATOR;typedef typename T_GRID::VECTOR_HORIZONTAL TV_HORIZONTAL;
public:
    typedef SOLIDS_FLUIDS_DRIVER<TV> BASE;
    using BASE::time;using BASE::current_frame;using BASE::output_number;using BASE::Write_Time;using BASE::Write_Last_Frame;using BASE::Write_First_Frame;using BASE::Write_Substep;
    using BASE::next_dt;using BASE::next_done;

    SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>& example;

    SOLIDS_FLUIDS_DRIVER_RLE(SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>& example_input);
    virtual ~SOLIDS_FLUIDS_DRIVER_RLE();

//#####################################################################
    void Initialize() PHYSBAM_OVERRIDE;
    void Initialize_Grids();
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE;
    void Advance_To_Target_Time(const T target_time) PHYSBAM_OVERRIDE;
    void Compute_Next_Dt(const T next_time,const bool done_this_frame,T& next_dt,bool& next_done);
    void Rebuild_Grid(const T time,ARRAY<T>* V_ghost);
    void Write_Output_Files(const int frame) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
#endif 
