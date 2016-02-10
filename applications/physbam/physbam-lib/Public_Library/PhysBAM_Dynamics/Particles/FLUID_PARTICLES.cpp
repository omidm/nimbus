//#####################################################################
// Copyright 2011, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Dynamics/Particles/PARTICLES_FORWARD.h>
#include <PhysBAM_Dynamics/Particles/FLUID_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_PARTICLES<TV>::
FLUID_PARTICLES()
    :F(0,0),momentum_residual(0,0),color(0,0),volume(0,0),volume_residual(0,0),on_boundary(0,0),type(0,0),normal(0,0),node(0,0),age(0,0),phi(0,0)
{
    Store_Id();
    Store_Velocity();
    array_collection->Add_Array(ATTRIBUTE_ID_FORCE,&F);
    //array_collection->Add_Array(ATTRIBUTE_ID_MOMENTUM_RESIDUAL,&momentum_residual);
    array_collection->Add_Array(ATTRIBUTE_ID_COLOR,&color);
    array_collection->Add_Array(ATTRIBUTE_ID_VOLUME,&volume);
    array_collection->Add_Array(ATTRIBUTE_ID_VOLUME_RESIDUAL,&volume_residual);
    array_collection->Add_Array(ATTRIBUTE_ID_ON_BOUNDARY,&on_boundary);
    array_collection->Add_Array(ATTRIBUTE_ID_TYPE,&type);
    //array_collection->Add_Array(ATTRIBUTE_ID_NORMAL,&normal);
    array_collection->Add_Array(ATTRIBUTE_ID_NODE,&node);
    array_collection->Add_Array(ATTRIBUTE_ID_AGE,&age);
    //array_collection->Add_Array(ATTRIBUTE_ID_PHI,&phi);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> FLUID_PARTICLES<TV>::
~FLUID_PARTICLES()
{}
//#####################################################################
template class FLUID_PARTICLES<VECTOR<float,1> >;
template class FLUID_PARTICLES<VECTOR<float,2> >;
template class FLUID_PARTICLES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FLUID_PARTICLES<VECTOR<double,1> >;
template class FLUID_PARTICLES<VECTOR<double,2> >;
template class FLUID_PARTICLES<VECTOR<double,3> >;
#endif
}
