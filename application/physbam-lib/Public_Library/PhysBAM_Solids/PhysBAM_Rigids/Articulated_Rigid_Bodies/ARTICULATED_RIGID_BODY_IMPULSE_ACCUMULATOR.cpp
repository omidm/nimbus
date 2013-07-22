//#####################################################################
// Copyright 2008, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR.h>
using namespace PhysBAM;
//#####################################################################
// Function Reset
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>::
Reset()
{
    accumulated_impulse=TWIST<TV>();
}
//#####################################################################
// Function Add_Impulse
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>::
Add_Impulse(const TV& location,const TWIST<TV>& impulse)
{
    accumulated_impulse+=impulse;
    //{std::stringstream ss;ss<<"M: Adding impulse to "<<joint.id_number<<" with size "<<impulse<<" making total "<<accumulated_impulse<<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
//template class ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<VECTOR<float,1> >;
template class ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<VECTOR<float,2> >;
template class ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
//template class ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<VECTOR<double,1> >;
template class ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<VECTOR<double,2> >;
template class ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<VECTOR<double,3> >;
#endif
