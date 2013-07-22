//#####################################################################
// Copyright 2005-2009, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace MPI_UTILITIES
//#####################################################################
#ifdef USE_MPI
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_REPULSION_UTILITIES.h>
namespace PhysBAM{
namespace MPI_UTILITIES{
//#####################################################################
// Datatype conversion for POINT_FACE_REPULSION_PAIR
//#####################################################################
template<class T,int d> MPI::Datatype
DATATYPE_HELPER<POINT_FACE_REPULSION_PAIR<VECTOR<T,d> > >::Datatype()
{
    static MPI::Datatype datatype=MPI::DATATYPE_NULL;
    if(datatype==MPI::DATATYPE_NULL){ // NOTE: We are sending only VECTOR<int,d+1> nodes and T distance
        static POINT_FACE_REPULSION_PAIR<VECTOR<T,d> > zero_entry;static int lengths[3]={d+1,1,1};
        static MPI::Aint displacements[3]={0,(char*)&zero_entry.distance-(char*)&zero_entry,sizeof(zero_entry)};
        static MPI::Datatype old_types[3]={MPI::INT,MPI_UTILITIES::Datatype<T>(),MPI::UB};
        datatype=MPI::Datatype::Create_struct(3,lengths,displacements,old_types);
        datatype.Commit();}
    return datatype;
}
//#####################################################################
// Datatype conversion for EDGE_EDGE_REPULSION_PAIR
//#####################################################################
template<class T,int d> MPI::Datatype
DATATYPE_HELPER<EDGE_EDGE_REPULSION_PAIR<VECTOR<T,d> > >::Datatype()
{
    static MPI::Datatype datatype=MPI::DATATYPE_NULL;
    if(datatype==MPI::DATATYPE_NULL){ // NOTE: We are sending only VECTOR<int,2*d-2> nodes, VECTOR<T,d> collision_free_normal, and T distance
        static EDGE_EDGE_REPULSION_PAIR<VECTOR<T,d> > zero_entry;static int lengths[3]={2*d-2,d+1,1};
        static MPI::Aint displacements[3]={0,(char*)&zero_entry.collision_free_normal-(char*)&zero_entry,sizeof(zero_entry)};
        static MPI::Datatype old_types[3]={MPI::INT,MPI_UTILITIES::Datatype<T>(),MPI::UB};
        datatype=MPI::Datatype::Create_struct(3,lengths,displacements,old_types);
        datatype.Commit();}
    return datatype;
}
//#####################################################################
#define INSTANTIATION_HELPER_REPULSIONS(T,dim) \
    template struct DATATYPE_HELPER<POINT_FACE_REPULSION_PAIR<VECTOR<T,dim> > >; \
    template struct DATATYPE_HELPER<EDGE_EDGE_REPULSION_PAIR<VECTOR<T,dim> > >;
#define INSTANTIATION_HELPER(T) \
    INSTANTIATION_HELPER_REPULSIONS(T,2); \
    INSTANTIATION_HELPER_REPULSIONS(T,3);
INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
}
}
#endif
