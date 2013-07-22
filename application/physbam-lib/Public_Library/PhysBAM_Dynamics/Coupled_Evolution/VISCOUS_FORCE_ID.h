//#####################################################################
// Copyright 2010, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VISCOUS_FORCE_ID
//#####################################################################
#ifndef __VISCOUS_FORCE_ID__
#define __VISCOUS_FORCE_ID__

#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
namespace PhysBAM{

#ifndef COMPILE_ID_TYPES_AS_INT
PHYSBAM_DECLARE_ELEMENT_ID(VISCOUS_FORCE_ID,int,ELEMENT_ID_HELPER::for_loop|ELEMENT_ID_HELPER::logical|ELEMENT_ID_HELPER::add_T);
#else
typedef int VISCOUS_FORCE_ID;
#endif
}
#endif
