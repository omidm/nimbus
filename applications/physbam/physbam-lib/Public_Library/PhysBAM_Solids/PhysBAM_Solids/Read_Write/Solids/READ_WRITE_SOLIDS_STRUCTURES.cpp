//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_SOLIDS_STRUCTURES
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_FORWARD.h>
namespace PhysBAM{

void Initialize_Read_Write_Structures();
void Initialize_Bindings();
void Initialize_Rigid_Body_Binding();

void Initialize_Read_Write_Solids_Structures()
{
    static bool done=false;if(done) return;done=true;
    Initialize_Read_Write_Structures();
    Initialize_Bindings();
    Initialize_Rigid_Body_Binding();
}
}
