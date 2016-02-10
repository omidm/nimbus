//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDING
//#####################################################################
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDING.h>
namespace PhysBAM{
//#####################################################################
// Register this class as read-write
//#####################################################################
namespace {
bool Register_Embedding(){
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<EMBEDDING<VECTOR<float,2> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<EMBEDDING<VECTOR<float,3> > >();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<EMBEDDING<VECTOR<double,2> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<EMBEDDING<VECTOR<double,3> > >();
#endif
    return true;
}
static bool registered=Register_Embedding();
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> EMBEDDING<TV>::
EMBEDDING(GEOMETRY_PARTICLES<TV>& particles_input)
    :particles(particles_input),material_surface(material_surface_mesh,particles),need_destroy_particles(false)
{
    PHYSBAM_ASSERT(registered);
}
//#####################################################################
template class EMBEDDING<VECTOR<float,2> >;
template class EMBEDDING<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EMBEDDING<VECTOR<double,2> >;
template class EMBEDDING<VECTOR<double,3> >;
#endif
}
