#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2006, Frank Losasso, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_DRIVER_DYADIC
//#####################################################################
#ifndef __SOLIDS_FLUIDS_DRIVER_DYADIC__
#define __SOLIDS_FLUIDS_DRIVER_DYADIC__    

#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_DYADIC.h>
namespace PhysBAM{

template<class T_GRID>
class SOLIDS_FLUIDS_DRIVER_DYADIC:public SOLIDS_FLUIDS_DRIVER<typename T_GRID::VECTOR_T>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::MAP_MESH MAP_MESH;typedef typename T_GRID::CELL T_CELL;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FLOOD_FILL FLOOD_FILL;
    typedef typename T_GRID::UNIFORM_GRID UNIFORM_GRID;
    typedef typename UNIFORM_GRID::CELL_ITERATOR UNIFORM_CELL_ITERATOR;
    typedef typename UNIFORM_GRID::NODE_ITERATOR UNIFORM_NODE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<UNIFORM_GRID>::ARRAYS_SCALAR T_UNIFORM_ARRAYS;
    typedef typename LEVELSET_POLICY<T_GRID>::PARTICLE_LEVELSET T_PARTICLE_LEVELSET;typedef typename TURBULENCE_POLICY<TV>::TURBULENCE T_TURBULENCE;
    typedef typename INCOMPRESSIBLE_POLICY<T_GRID>::INCOMPRESSIBLE T_INCOMPRESSIBLE;typedef typename LEVELSET_POLICY<T_GRID>::EXTRAPOLATION_SCALAR EXTRAPOLATION_SCALAR;
    typedef typename T_UNIFORM_ARRAYS::template REBIND<TV>::TYPE UNIFORM_ARRAYS_VECTOR;
    typedef typename LEVELSET_POLICY<UNIFORM_GRID>::LEVELSET UNIFORM_LEVELSET;typedef typename LEVELSET_POLICY<T_GRID>::EXTRAPOLATION_VECTOR EXTRAPOLATION_VECTOR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::INTERPOLATION_SCALAR T_INTERPOLATION_SCALAR;
    typedef typename T_GRID::CELL_HELPER CELL_HELPER;
public:
    typedef SOLIDS_FLUIDS_DRIVER<TV> BASE;
    using BASE::time;using BASE::current_frame;using BASE::next_dt;using BASE::next_done;using BASE::Write_Time;using BASE::Write_Last_Frame;using BASE::Write_Substep;
    using BASE::Write_First_Frame;

    SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>& example;
    ARRAY<bool> occupied_cell_for_coarsening;

    SOLIDS_FLUIDS_DRIVER_DYADIC(SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>& example_input);
    virtual ~SOLIDS_FLUIDS_DRIVER_DYADIC();

//#####################################################################
    void Initialize() PHYSBAM_OVERRIDE;
    void Initialize_Fluids_Grids();
    void Set_Ghost_Density_And_Temperature_Inside_Flame_Core();
    void Advance_To_Target_Time(const T target_time) PHYSBAM_OVERRIDE;
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE;
    void Compute_Next_Dt(const T next_time,const bool done_this_frame,T& next_dt,bool& next_done);
private:
    static void Find_Ghost_Cells_To_Refine(void* data,const T_CELL* cell1,const T_CELL* cell2,int axis);
public:
    bool Coarsen_Tree(T_CELL* cell,const ARRAY<typename T_CELL::REFINE_ACTION>& refine_action);
    void Coarsen_Ghost_Cells(T_CELL* cell);
    void Delete_All_Particles_In_Cell(T_CELL* cell);
    void Interpolate_To_Direct_Children(T_CELL* cell);
    void Update_Tree_Topology(const T dt,const T time);
    void Write_Output_Files(const int frame) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
#endif
