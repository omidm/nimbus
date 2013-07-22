//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED
//##################################################################### 
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<TV>::
UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info,const int number_of_ghost_cells_input,const T_REGION& region_type_input,
    const int side_input,int axis_input)
    :BASE(info.grid,number_of_ghost_cells_input,region_type_input,side_input,axis_input),collision_index(0),collision_face_info(info.collision_face_info)
{
    Next();
}
//#####################################################################
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<VECTOR<float,1> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<VECTOR<float,2> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<VECTOR<double,1> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<VECTOR<double,2> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<VECTOR<double,3> >;
#endif
