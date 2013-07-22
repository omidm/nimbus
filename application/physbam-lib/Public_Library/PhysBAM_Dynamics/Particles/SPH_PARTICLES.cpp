//#####################################################################
// Copyright 2004-2008, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Dynamics/Particles/PARTICLES_FORWARD.h>
#include <PhysBAM_Dynamics/Particles/SPH_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SPH_PARTICLES<TV>::
SPH_PARTICLES()
    :V(0,0)
{
    array_collection->Add_Array(ATTRIBUTE_ID_V,&V);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SPH_PARTICLES<TV>::
~SPH_PARTICLES()
{}
//#####################################################################
template class SPH_PARTICLES<VECTOR<float,1> >;
template class SPH_PARTICLES<VECTOR<float,2> >;
template class SPH_PARTICLES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SPH_PARTICLES<VECTOR<double,1> >;
template class SPH_PARTICLES<VECTOR<double,2> >;
template class SPH_PARTICLES<VECTOR<double,3> >;
#endif
}
