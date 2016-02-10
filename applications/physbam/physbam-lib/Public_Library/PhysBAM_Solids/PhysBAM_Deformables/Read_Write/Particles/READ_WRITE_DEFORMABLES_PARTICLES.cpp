//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_DEFORMABLES_PARTICLES
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY_COLLECTION.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLES_PARTICLES_FORWARD.h>
namespace PhysBAM{

void Initialize_Deformables_Particles()
{
    static bool done=false;if(done) return;done=true;
    Register_Attribute_Name(ATTRIBUTE_ID_MASS,"mass");
    Register_Attribute_Name(ATTRIBUTE_ID_ONE_OVER_MASS,"one_over_mass");
    Register_Attribute_Name(ATTRIBUTE_ID_EFFECTIVE_MASS,"effective_mass");
    Register_Attribute_Name(ATTRIBUTE_ID_ONE_OVER_EFFECTIVE_MASS,"one_over_effective_mass");

    Initialize_Geometry_Particle();
}
}
