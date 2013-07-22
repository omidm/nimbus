//#####################################################################
// Copyright 2008, Craig Schroeder, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRACTURE_PATTERN
//##################################################################### 
#ifndef __FRACTURE_PATTERN__
#define __FRACTURE_PATTERN__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>

namespace PhysBAM{

template<class TV> class RIGID_BODY;
template<class T> class FRACTURE_REGION;

template<class T>
class FRACTURE_PATTERN
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,3> TV_INT;
public:
    typedef int HAS_TYPED_READ_WRITE;

    ARRAY<FRACTURE_REGION<T>*> regions;
    bool use_particle_partitions;
    
    FRACTURE_PATTERN();
    virtual ~FRACTURE_PATTERN();

//#####################################################################
    void Read(TYPED_ISTREAM& input);
    void Write(TYPED_OSTREAM& output) const;
    void Intersect_With_Rigid_Body(const RIGID_BODY<VECTOR<T,1> >& body,const VECTOR<T,1>& point_of_impact,ARRAY<int>& added_bodies,const bool allow_refracture,const bool use_particle_optimization,const bool generate_object_tessellation=false);
    void Intersect_With_Rigid_Body(const RIGID_BODY<VECTOR<T,2> >& body,const VECTOR<T,2>& point_of_impact,ARRAY<int>& added_bodies,const bool allow_refracture,const bool use_particle_optimization,const bool generate_object_tessellation=false);
    void Intersect_With_Rigid_Body(const RIGID_BODY<TV>& body,const TV& point_of_impact,ARRAY<int>& added_bodies,const bool allow_refracture,const bool use_particle_optimization,const bool generate_object_tessellation=false);
//#####################################################################
};
}
#endif
