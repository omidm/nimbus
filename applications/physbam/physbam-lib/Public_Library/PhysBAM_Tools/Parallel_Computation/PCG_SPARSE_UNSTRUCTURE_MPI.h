//#####################################################################
// Copyright 2011, Linhai Qiu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PCG_SPARSE_UNSTRUCTURE_MPI
//#####################################################################
#ifndef __PCG_SPARSE_UNSTRUCTURE_MPI__
#define __PCG_SPARSE_UNSTRUCTURE_MPI__

#ifdef USE_MPI

#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_MPI.h>
namespace PhysBAM{

class SPARSE_MATRIX_PARTITION;

template<class T_GRID>
class PCG_SPARSE_UNSTRUCTURE_MPI:public PCG_SPARSE_MPI<T_GRID>
{
    typedef PCG_SPARSE_MPI<T_GRID> BASE;
public:
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    using BASE::Parallel_Solve;using BASE::comm;using BASE::partition;

    ARRAY<ARRAY<bool> >& is_boundary_grid;
    ARRAY<int>& local_to_global_grid_index_map;
    ARRAY<ARRAY<MPI::Datatype> > boundary_datatypes_array,ghost_datatypes_array;
    int n_local_grids;

    PCG_SPARSE_UNSTRUCTURE_MPI(PCG_SPARSE<T>& pcg_input,MPI::Intracomm& comm_input,SPARSE_MATRIX_PARTITION& partition_input,ARRAY<ARRAY<bool> >& is_boundary_grid_input,ARRAY<int>& local_to_global_grid_index_map_input)
        :PCG_SPARSE_MPI<T_GRID>(pcg_input,comm_input,partition_input),is_boundary_grid(is_boundary_grid_input),local_to_global_grid_index_map(local_to_global_grid_index_map_input)
    {n_local_grids=is_boundary_grid.m;boundary_datatypes_array.Resize(n_local_grids);ghost_datatypes_array.Resize(n_local_grids);}

    void Fill_Ghost_Cells(VECTOR_ND<T>& v)
    {
        //LOG::SCOPE scope1("Fill_Ghost_Cells construct");
        ARRAY<MPI::Request> requests;requests.Preallocate(2*n_local_grids*partition.number_of_sides);
        for(int local_grid_index=1;local_grid_index<=n_local_grids;local_grid_index++){ 
            int this_grid_index=local_to_global_grid_index_map(local_grid_index);
            for(int s=1;s<=partition.number_of_sides;s++) if(boundary_datatypes_array(local_grid_index)(s)!=MPI::DATATYPE_NULL){
                requests.Append(comm.Isend(v.x-1,1,boundary_datatypes_array(local_grid_index)(s),partition.neighbor_ranks(s),Get_Send_Tag(this_grid_index,s)));}
            for(int s=1;s<=partition.number_of_sides;s++) if(ghost_datatypes_array(local_grid_index)(s)!=MPI::DATATYPE_NULL){
                requests.Append(comm.Irecv(v.x-1,1,ghost_datatypes_array(local_grid_index)(s),partition.neighbor_ranks(s),Get_Recv_Tag(s,this_grid_index)));}}
        MPI_UTILITIES::Wait_All(requests);
    }

    int Get_Send_Tag(int send_grid,int recv_grid,int increment=0) const
    {
        STATIC_ASSERT(TV_INT::m<=3);int tag=0;
        tag=(send_grid-1)*partition.number_of_sides+recv_grid-1;
        return tag+increment;
    }

    int Get_Recv_Tag(int send_grid,int recv_grid,int increment=0) const
    {
        return Get_Send_Tag(send_grid,recv_grid,increment);
    }
    
    void Initialize_Datatypes()
    {
        for(int local_grid_index=1;local_grid_index<=n_local_grids;local_grid_index++){
            MPI_UTILITIES::Free_Elements_And_Clean_Memory(boundary_datatypes_array(local_grid_index));
            MPI_UTILITIES::Free_Elements_And_Clean_Memory(ghost_datatypes_array(local_grid_index));
            boundary_datatypes_array(local_grid_index).Resize(partition.number_of_sides);
            ghost_datatypes_array(local_grid_index).Resize(partition.number_of_sides);}

        for(int local_grid_index=1;local_grid_index<=n_local_grids;local_grid_index++) for(int s=1;s<=partition.number_of_sides;s++) if(is_boundary_grid(local_grid_index)(s)){
            if(partition.boundary_indices(s).m){
                const ARRAY<int>& displacements=partition.boundary_indices(s);
                ARRAY<int> block_lengths(displacements.m,false);ARRAYS_COMPUTATIONS::Fill(block_lengths,1);
                boundary_datatypes_array(local_grid_index)(s)=MPI_UTILITIES::Datatype<T>().Create_indexed(displacements.m,&block_lengths(1),&displacements(1)); // TODO: collapse consecutive elements into blocks
                boundary_datatypes_array(local_grid_index)(s).Commit();}
            int ghost_indices_length=partition.ghost_indices(s).Size()+1;
            if(ghost_indices_length){
                ghost_datatypes_array(local_grid_index)(s)=MPI_UTILITIES::Datatype<T>().Create_indexed(1,&ghost_indices_length,&partition.ghost_indices(s).min_corner);
                ghost_datatypes_array(local_grid_index)(s).Commit();}}
    }
};
}
#endif
#endif
