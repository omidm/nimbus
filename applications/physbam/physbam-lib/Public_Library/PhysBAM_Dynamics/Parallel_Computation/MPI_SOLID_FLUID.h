//#####################################################################
// Copyright 2008, Avi Robinson-Mosher, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_SOLID_FLUID
//#####################################################################
#ifndef __MPI_SOLID_FLUID__
#define __MPI_SOLID_FLUID__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
namespace MPI{class Group;class Intracomm;class Request;class Status;class Op;}
namespace PhysBAM{

class MPI_PACKAGE;
template<class TV> class FLUID_SYSTEM_MPI;
template<class TV> class SOLID_SYSTEM_MPI;
template<class TV> class GENERALIZED_VELOCITY;
template<class TV> class SOLID_BODY_COLLECTION;
class SPARSE_MATRIX_PARTITION;

template<class TV>
class MPI_SOLID_FLUID:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    int rank; // global rank
    int number_of_processes;
    int number_of_solid_processes;
    int solid_node;
    MPI::Intracomm* comm;
    MPI::Group *group,*solid_group,*fluid_group;
    VECTOR_ND<int> solid_ranks,fluid_ranks;
private:
    mutable int current_tag;
public:

    MPI_SOLID_FLUID();
    ~MPI_SOLID_FLUID();

    int Number_Of_Processors() const
    {return number_of_processes;}

    int Get_Unique_Tag() const
    {return current_tag=max(32,(current_tag+1)&((1<<15)-1));}

public:
    void Exchange_Solid_Positions_And_Velocities(SOLID_BODY_COLLECTION<TV>& solid_body_collection) const;
    bool Fluid_Node() const;
    bool Solid_Node() const;
    void Create_Fluid_Comm_For_Solid_Nodes() const;
    template<class T2> void Reduce_Add(const T2& input,T2& output) const;
    T Reduce_Min(const T local_value) const;
    T Reduce_Max(const T local_value) const;
    void Parallel_Solve_Fluid_Part(FLUID_SYSTEM_MPI<TV>& fluid_system,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& x_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& b_array,
        KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& p_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& ap_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& ar_array,
        KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& r_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& z_array,KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > >& zaq_array,
        const int min_iterations,const int max_iterations,const T tolerance,const bool recompute_preconditioner,ARRAY<MPI::Intracomm>* fluid_comm,
        ARRAY<SPARSE_MATRIX_PARTITION>* partitions);
    void Parallel_Solve_Solid_Part(SOLID_SYSTEM_MPI<TV>& solid_system,GENERALIZED_VELOCITY<TV>& x_array,GENERALIZED_VELOCITY<TV>& b_array,GENERALIZED_VELOCITY<TV>& p_array,
        GENERALIZED_VELOCITY<TV>& ap_array,GENERALIZED_VELOCITY<TV>& ar_array,GENERALIZED_VELOCITY<TV>& r_array,GENERALIZED_VELOCITY<TV>& z_array,GENERALIZED_VELOCITY<TV>& zaq_array,
        const int min_iterations,const int max_iterations,const T tolerance);
    void Distribute_Lists_From_Solid_Node(GENERALIZED_VELOCITY<TV>& F) const;
    void Exchange_Coupled_Deformable_Particle_List(ARRAY<int>* fluid_list,ARRAY<ARRAY<int> >* results);
    void Aggregate_Lists_To_Solid_Node(GENERALIZED_VELOCITY<TV>& F);
//#####################################################################
};
}
#endif
