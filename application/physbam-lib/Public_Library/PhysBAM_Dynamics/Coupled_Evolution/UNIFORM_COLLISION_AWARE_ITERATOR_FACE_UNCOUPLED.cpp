//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED
//##################################################################### 
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<TV>::
UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info,const int number_of_ghost_cells_input,bool use_outside_input,
    const T_REGION& region_type_input,const int side_input,int axis_input)
    :BASE(info.grid,number_of_ghost_cells_input,region_type_input,side_input,axis_input),outside_fluid(*info.outside_fluid),
    collision_index(1),collision_face_info(info.collision_face_info),scan_end(INT_MIN),use_outside(use_outside_input)
{
    index(TV::dimension)-=2;
    Next_Fluid();
}
template<class TV> UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<TV>::
~UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED()
{
}
//#####################################################################
// Function Next_Helper
//#####################################################################
template<class TV> void UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<TV>::
Next_Helper()
{
    BASE::Next();
    for(int c;(c=Compare_Collision_Index())<=0;collision_index++) if(!c) BASE::Next();

    scan_end=grid.counts(TV::dimension)+(TV::dimension==axis)+number_of_ghost_cells+(region_type==GRID<TV>::WHOLE_REGION?1:0);
    if(collision_index<=collision_face_info.Size()){
        const COLLISION_FACE_INFO<TV>& cfi=collision_face_info(collision_index);
        if(axis==cfi.axis && index.Remove_Index(TV::dimension)==cfi.index.Remove_Index(TV::dimension))
            scan_end=cfi.index(TV::dimension);}
}
//#####################################################################
// Function Compare_Collision_Index
//#####################################################################
template<class TV> int UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<TV>::
Compare_Collision_Index() const
{
    if(collision_index>collision_face_info.Size()) return 1;
    const COLLISION_FACE_INFO<TV>& cfi=collision_face_info(collision_index);
    if(cfi.axis<axis) return -1;
    if(cfi.axis>axis) return 1;
    for(int i=1;i<=TV::dimension;i++){
        if(cfi.index(i)<index(i)) return -1;
        if(cfi.index(i)>index(i)) return 1;}
    return 0;
}
//#####################################################################
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<VECTOR<float,1> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<VECTOR<float,2> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<VECTOR<double,1> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<VECTOR<double,2> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<VECTOR<double,3> >;
#endif
