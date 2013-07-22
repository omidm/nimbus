//#####################################################################
// Copyright 2007-2009, Jon Gretarsson, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Dynamics/Particles/COMPRESSIBLE_FLUID_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLES_FORWARD.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
/*template<class TV> COMPRESSIBLE_FLUID_PARTICLES<TV>::
COMPRESSIBLE_FLUID_PARTICLES(ARRAY_COLLECTION* array_collection_input)
    :BASE(array_collection_input),rho(0,0),E(0,0),phi(0,0),grad_phi(0,0),V(0,0)
{
    array_collection->Add_Array(ATTRIBUTE_ID_RHO,&rho);
    array_collection->Add_Array(ATTRIBUTE_ID_E,&E);
    array_collection->Add_Array(ATTRIBUTE_ID_PHI,&phi);
    array_collection->Add_Array(ATTRIBUTE_ID_GRAD_PHI,&grad_phi);
    array_collection->Add_Array(ATTRIBUTE_ID_V,&V);
    }*/
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMPRESSIBLE_FLUID_PARTICLES<TV>::
COMPRESSIBLE_FLUID_PARTICLES()
    :rho(0,0),E(0,0),phi(0,0),grad_phi(0,0),V(0,0)
{
    array_collection->Add_Array(ATTRIBUTE_ID_RHO,&rho);
    array_collection->Add_Array(ATTRIBUTE_ID_E,&E);
    array_collection->Add_Array(ATTRIBUTE_ID_PHI,&phi);
    array_collection->Add_Array(ATTRIBUTE_ID_GRAD_PHI,&grad_phi);
    array_collection->Add_Array(ATTRIBUTE_ID_V,&V);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMPRESSIBLE_FLUID_PARTICLES<TV>::
~COMPRESSIBLE_FLUID_PARTICLES()
{}
//#####################################################################
template class COMPRESSIBLE_FLUID_PARTICLES<VECTOR<float,1> >;
template class COMPRESSIBLE_FLUID_PARTICLES<VECTOR<float,2> >;
template class COMPRESSIBLE_FLUID_PARTICLES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMPRESSIBLE_FLUID_PARTICLES<VECTOR<double,1> >;
template class COMPRESSIBLE_FLUID_PARTICLES<VECTOR<double,2> >;
template class COMPRESSIBLE_FLUID_PARTICLES<VECTOR<double,3> >;
#endif
}
