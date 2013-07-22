//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Nipun Kwatra, Frank Losasso, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLES_FORWARD.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
/*template<class TV> PARTICLE_LEVELSET_PARTICLES<TV>::
PARTICLE_LEVELSET_PARTICLES(ARRAY_COLLECTION& array_collection_input)
    :POINT_CLOUD<TV>(array_collection_input),quantized_collision_distance(0,0),age(0,0),radius(0,0),next(0)
{
    array_collection->Add_Array(ATTRIBUTE_ID_QUANTIZED_COLLISION_DISTANCE,&quantized_collision_distance);
    array_collection->Add_Array(ATTRIBUTE_ID_AGE,&age);
    array_collection->Add_Array(ATTRIBUTE_ID_RADIUS,&radius);
    }*/
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_LEVELSET_PARTICLES<TV>::
PARTICLE_LEVELSET_PARTICLES()
    :quantized_collision_distance(0,0),age(0,0),radius(0,0),next(0)
{
    array_collection->Add_Array(ATTRIBUTE_ID_QUANTIZED_COLLISION_DISTANCE,&quantized_collision_distance);
    array_collection->Add_Array(ATTRIBUTE_ID_AGE,&age);
    array_collection->Add_Array(ATTRIBUTE_ID_RADIUS,&radius);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_LEVELSET_PARTICLES<TV>::
~PARTICLE_LEVELSET_PARTICLES()
{}
//#####################################################################
template class PARTICLE_LEVELSET_PARTICLES<VECTOR<float,1> >;
template class PARTICLE_LEVELSET_PARTICLES<VECTOR<float,2> >;
template class PARTICLE_LEVELSET_PARTICLES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PARTICLE_LEVELSET_PARTICLES<VECTOR<double,1> >;
template class PARTICLE_LEVELSET_PARTICLES<VECTOR<double,2> >;
template class PARTICLE_LEVELSET_PARTICLES<VECTOR<double,3> >;
#endif
}
