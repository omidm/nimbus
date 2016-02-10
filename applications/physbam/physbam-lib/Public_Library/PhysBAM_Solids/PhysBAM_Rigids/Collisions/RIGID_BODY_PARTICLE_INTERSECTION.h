//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Craig Schroeder, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_PARTICLE_INTERSECTION
//##################################################################### 
#ifndef __RIGID_BODY_PARTICLE_INTERSECTION__
#define __RIGID_BODY_PARTICLE_INTERSECTION__

namespace PhysBAM{

template<class TV>
struct RIGID_BODY_PARTICLE_INTERSECTION
{
    TV particle_location; // object space location relative to particle body
    int particle_index;
    int particle_body,levelset_body;

    RIGID_BODY_PARTICLE_INTERSECTION()
        :particle_index(0),particle_body(0),levelset_body(0)
    {}

    RIGID_BODY_PARTICLE_INTERSECTION(const TV& initial_particle_location,const int initial_particle_index,const int initial_particle_body,const int initial_levelset_body)
        :particle_location(initial_particle_location),particle_index(initial_particle_index),particle_body(initial_particle_body),levelset_body(initial_levelset_body)
    {}
};

}
#endif
