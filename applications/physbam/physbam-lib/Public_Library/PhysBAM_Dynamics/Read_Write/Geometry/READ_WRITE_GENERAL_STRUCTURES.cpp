//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_GENERAL_STRUCTURES
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_FORWARD.h>
#include <PhysBAM_Dynamics/Geometry/GENERAL_GEOMETRY_FORWARD.h>
namespace PhysBAM{

void Register_Read_Write_Embedded_Tetrahedralized_Volume_Boundary_Surface();
void Register_Read_Write_Embedded_Material_Surface();
void Register_Read_Write_Embedding();
void Initialize_Read_Write_Solids_Structures();

void Initialize_Read_Write_General_Structures()
{
    static bool done=false;if(done) return;done=true;
    Initialize_Read_Write_Solids_Structures();
    Register_Read_Write_Embedded_Tetrahedralized_Volume_Boundary_Surface();
    Register_Read_Write_Embedded_Material_Surface();
    Register_Read_Write_Embedding();
}
}
