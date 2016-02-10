//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED
//#####################################################################
#ifndef __UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED__
#define __UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED__

#include <PhysBAM_Tools/Grids_Uniform/SIDED_FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>

namespace PhysBAM{

template<class TV> class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO;
template<class TV> struct COLLISION_FACE_INFO;

template<class TV>
class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED:public UNIFORM_GRID_ITERATOR_FACE<TV>
{
public:
    enum WORKAROUND {d=TV::dimension};
    typedef typename TV::SCALAR T;typedef VECTOR<int,d> TV_INT;typedef typename GRID<TV>::REGION T_REGION;typedef UNIFORM_GRID_ITERATOR_FACE<TV> BASE;
    using BASE::grid;using BASE::index;using BASE::region;using BASE::valid;using BASE::Reset;using BASE::current_region;using BASE::Add_Region;using BASE::Reset_Regions;using BASE::axis;
    using BASE::First_Cell_Index;using BASE::Second_Cell_Index;

    int collision_index;
    int side;
    const ARRAY<COLLISION_FACE_INFO<TV> >& collision_face_info;

    // TODO: Handle ghost cells, regions, side, and axis properly construction values properly.
    // axis_input==0 means iterate through faces in all dimensions
    explicit UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info,const int number_of_ghost_cells_input=0,
        const T_REGION& region_type_input=GRID<TV>::WHOLE_REGION,const int side_input=0,int axis_input=0);

    void Next()
    {if(++collision_index<=collision_face_info.m){index=collision_face_info(collision_index).index;axis=collision_face_info(collision_index).axis;side=collision_face_info(collision_index).side;}}

    bool Valid()
    {return collision_index<=collision_face_info.m;}

    const ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> >& Get_Simplices() const
    {return collision_face_info(collision_index).simplices;}

    TV_INT Real_Cell_Index() const
    {return side==1?First_Cell_Index():Second_Cell_Index();}

    TV_INT Ghost_Cell_Index() const
    {return side==2?First_Cell_Index():Second_Cell_Index();}

    //void Neighbor_Face_Stencil(VECTOR<SIDED_FACE_INDEX<d>,2*d-2>& neighbor_faces) const
    //{TV_INT real_cell=Real_Cell_Index();
    //for(int i=1;i<=d-1;i++){int a=i+(i>=axis);
    //    neighbor_faces(2*i-1)=SIDED_FACE_INDEX<d>(2,a,real_cell);
    //    neighbor_faces(2*i)=SIDED_FACE_INDEX<d>(1,a,real_cell+TV_INT::Axis_Vector(a));}}
};
}
#endif
