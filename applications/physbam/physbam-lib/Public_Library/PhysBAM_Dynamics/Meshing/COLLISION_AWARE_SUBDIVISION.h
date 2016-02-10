//#####################################################################
// Copyright 2003-2006, Robert Bridson, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COLLISION_AWARE_SUBDIVISION 
//##################################################################### 
#ifndef __COLLISION_AWARE_SUBDIVISION__
#define __COLLISION_AWARE_SUBDIVISION__

#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
namespace PhysBAM{

template<class T> class TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY;

template<class T>
class COLLISION_AWARE_SUBDIVISION:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;
public:
    TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry;
    COLLISION_GEOMETRY_COLLECTION<TV> collision_body_list;
    T collision_tolerance;
    
    COLLISION_AWARE_SUBDIVISION(TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry_input)
        :geometry(geometry_input),collision_tolerance((T)1e-6)
    {}
    
//#####################################################################
    bool Push_Surface_Outside_Of_Collision_Bodies(const T push_distance=0,const int max_number_of_attempts=3);
    void Subdivide(const int number_of_subdivisions,const bool push_outside_collision_bodies=true,const T push_distance=1e-6,const int push_attempts=3,const bool verbose=true);
    void Subdivide();
//#####################################################################
};
}
#endif
