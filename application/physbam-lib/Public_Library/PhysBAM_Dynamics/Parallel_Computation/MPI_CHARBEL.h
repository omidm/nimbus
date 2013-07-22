//#####################################################################
// Copyright 2008, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_CHARBEL
//#####################################################################
#ifndef __MPI_CHARBEL__
#define __MPI_CHARBEL__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Dynamics/Particles/COMPRESSIBLE_FLUID_PARTICLES.h>
namespace MPI{class Group;class Intracomm;class Request;class Status;class Op;}
namespace PhysBAM{

class MPI_PACKAGE;

template<class T>
class MPI_CHARBEL:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;
public:
    int rank; // global rank
    int number_of_processes;
    MPI::Intracomm* comm;
    MPI::Group *group;

    TETRAHEDRALIZED_VOLUME<T>* local_tet_volume;
    COMPRESSIBLE_FLUID_PARTICLES<TV> physbam_particles;
private:
    mutable int current_tag;
    ARRAY<ARRAY<int> > interior_particles_to_recv;
    ARRAY<ARRAY<int> > interior_particles_to_send;
    ARRAY<ARRAY<int> > ghost_particles_to_recv;
    ARRAY<ARRAY<int> > ghost_particles_to_send;
    ARRAY<int> global_to_local_physbam_map;
    ARRAY<int> global_to_local_aerof_map;
public:

    MPI_CHARBEL();
    ~MPI_CHARBEL();

    int Number_Of_Processors() const
    {return number_of_processes;}

    int Get_Unique_Tag() const
    {return current_tag=max(32,(current_tag+1)&((1<<15)-1));}

public:
    void Reduce_Max(int& local_value) const;
    void Setup_AEROF_PhysBAM_Mapping(TETRAHEDRALIZED_VOLUME<T>& tet_volume,ARRAY<ARRAY<int> >& tets_to_send,ARRAY<int>& local_to_global_map,
                     const int& global_particle_count,const RANGE<TV>& domain);
    void Exchange_Compressible_Data(COMPRESSIBLE_FLUID_PARTICLES<TV>& particles_aerof);
    void Exchange_Back_Compressible_Data(COMPRESSIBLE_FLUID_PARTICLES<TV>& particles_aerof);
//#####################################################################
};
}
#endif
