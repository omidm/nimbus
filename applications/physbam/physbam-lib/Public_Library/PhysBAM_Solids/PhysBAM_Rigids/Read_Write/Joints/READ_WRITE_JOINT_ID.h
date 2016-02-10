//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_JOINT_ID
//#####################################################################
#ifndef __READ_WRITE_JOINT_ID__
#define __READ_WRITE_JOINT_ID__

#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_ELEMENT_ID.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_ID.h>
namespace PhysBAM{

#ifndef COMPILE_ID_TYPES_AS_INT
template<class RW>
class Read_Write<JOINT_ID,RW>:public Read_Write<ELEMENT_ID<JOINT_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T>,RW>
{
};
#endif
}
#endif
