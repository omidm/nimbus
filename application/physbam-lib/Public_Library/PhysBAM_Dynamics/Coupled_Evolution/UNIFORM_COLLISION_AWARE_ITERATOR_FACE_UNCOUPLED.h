//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED
//#####################################################################
#ifndef __UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED__
#define __UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM{
template<class TV> class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO;

template<class TV>
class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED:public UNIFORM_GRID_ITERATOR_FACE<TV>
{
    typedef VECTOR<int,TV::m> TV_INT;
    const ARRAY<bool,TV_INT>& outside_fluid;
public:
    typedef typename TV::SCALAR T;
    typedef UNIFORM_GRID_ITERATOR_FACE<TV> BASE;typedef typename GRID<TV>::REGION T_REGION;
    using BASE::grid;using BASE::index;using BASE::region;using BASE::valid;using BASE::Reset;using BASE::current_region;using BASE::number_of_ghost_cells;using BASE::axis;
    using BASE::Add_Region;using BASE::Reset_Regions;using BASE::First_Cell_Index;using BASE::Second_Cell_Index;using BASE::Valid;using BASE::Full_Index;using BASE::region_type;

    int collision_index;
    const ARRAY<COLLISION_FACE_INFO<TV> >& collision_face_info;
    int scan_end;
    bool use_outside;

    // axis_input==0 means iterate through faces in all dimensions
    explicit UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info,const int number_of_ghost_cells_input=0,bool use_outside_input=true,
        const T_REGION& region_type_input=GRID<TV>::WHOLE_REGION,const int side_input=0,int axis_input=0);

    ~UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED();

    void Next()
    {index(TV::dimension)++;if(index(TV::dimension)>=scan_end) Next_Helper();}

    // TODO: Careful about ghost cells and the validity of outside_fluid.
    void Next_Fluid()
    {Next();if(use_outside) while(Valid() && outside_fluid(First_Cell_Index()) && outside_fluid(Second_Cell_Index())) Next();}

    void Next_Helper();
private:
    int Compare_Collision_Index() const;
};
}
#endif
