//#####################################################################
// Copyright 2011, Linhai Qiu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CHIMERA_GRID
//#####################################################################
#ifndef __CHIMERA_GRID__
#define __CHIMERA_GRID__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GRID.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Parallel_Computation/SPARSE_MATRIX_PARTITION.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Grid_Based_Fields/INCOMPRESSIBLE_FLUID_CONTAINER.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>

namespace MPI{class Group;class Intracomm;class Request;class Status;class Op;}
namespace PhysBAM{

class MPI_PACKAGE;
template<class T_GRID,class T2> class BOUNDARY_CHIMERA;

template<class T_GRID,class T2=typename T_GRID::SCALAR>
class CHIMERA_GRID:public NONCOPYABLE
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::VECTOR_T TV;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;typedef typename T_ARRAYS_BASE::template REBIND<T2>::TYPE T_ARRAYS_T2;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename REBIND<T_ARRAYS_SCALAR,T2>::TYPE T_ARRAYS_SCALAR_T2;
    typedef typename REBIND<T_ARRAYS_SCALAR,bool>::TYPE T_ARRAYS_BOOL;typedef typename REBIND<T_ARRAYS_SCALAR,int>::TYPE T_ARRAYS_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename REBIND<T_FACE_ARRAYS_SCALAR,T2>::TYPE T_FACE_ARRAYS_T2;
    typedef typename REBIND<T_FACE_ARRAYS_SCALAR,int>::TYPE T_FACE_ARRAYS_INT;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_BOOL::template REBIND<VECTOR<bool,T_GRID::dimension> >::TYPE T_FACE_ARRAYS_BOOL_DIMENSION;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef PAIR<ARRAY<TV_INT>,ARRAY<T> > CELL_PACKAGE_TYPE;typedef PAIR<ARRAY<FACE_INDEX<TV::dimension> >,ARRAY<T> > FACE_PACKAGE_TYPE;
public:
    typedef FACE_INDEX<TV::dimension> D_FACE_INDEX;
    typedef PAIR<int,D_FACE_INDEX> GRID_FACE_INDEX;
    typedef PAIR<int,TV_INT> GRID_CELL_INDEX;
    typedef VECTOR<GRID_CELL_INDEX,2> VORONOI_FACE_INDICES;
    typedef PAIR<VORONOI_FACE_INDICES,T> VORONOI_FACE;

    RIGID_GRID_COLLECTION<T_GRID>& local_grid;
    RIGID_GRID_COLLECTION<T_GRID> global_grid;
    ARRAY<T_GRID> initial_grids;
    ARRAY<RANGE<TV> > interpolation_domains;
    ARRAY<ORIENTED_BOX<TV> > previous_oriented_boxes;
    ARRAY<ORIENTED_BOX<TV> > oriented_boxes;
    int number_of_global_grids;
    int number_of_local_grids;
    int number_of_ghost_cells;
    bool use_mpi;
    int rank;
    int number_of_processors;
    ARRAY<int> mpi_grid_indices;
    ARRAY<int> grid_ranks;
    ARRAY<int> local_to_global_grid_index_map;
    ARRAY<int> global_to_local_grid_index_map;
    ARRAY<MPI_UNIFORM_GRID<T_GRID>*> mpi_split_grids;
    ARRAY<ARRAY<int> > side_neighbor_ranks_chimera;
    ARRAY<ARRAY<int> > all_neighbor_ranks_chimera;
    ARRAY<TV_INT> all_coordinates;
    ARRAY<T_ARRAYS_INT> mpi_grid_coordinate_to_global_grid_index_map;
    ARRAY<int> grids_sorted;//from finest to coarsest
    MPI::Intracomm* comm;
    ARRAY<MPI_PACKAGE>* mpi_packages_buffers;
    ARRAY<MPI::Request>* mpi_requests_buffers;
    bool use_previous_grid_frame_cell;
    bool use_previous_grid_frame_face;
    PAIR<bool,int> boundary_regions_found_cell;
    PAIR<bool,int> boundary_regions_found_face;
    PAIR<bool,int> overlap_regions_found_cell;
    PAIR<bool,int> overlap_regions_found_face;

    ARRAY<PAIR<CELL_PACKAGE_TYPE,int> > boundary_packages_cell_send;
    ARRAY<PAIR<CELL_PACKAGE_TYPE,int> > boundary_packages_cell_recv;
    ARRAY<PAIR<FACE_PACKAGE_TYPE,int> > boundary_packages_face_send;
    ARRAY<PAIR<FACE_PACKAGE_TYPE,int> > boundary_packages_face_recv;
    ARRAY<PAIR<CELL_PACKAGE_TYPE,int> > overlap_packages_cell_send;
    ARRAY<PAIR<CELL_PACKAGE_TYPE,int> > overlap_packages_cell_recv;
    ARRAY<PAIR<FACE_PACKAGE_TYPE,int> > overlap_packages_face_send;
    ARRAY<PAIR<FACE_PACKAGE_TYPE,int> > overlap_packages_face_recv;
    ARRAY<T_ARRAYS_SCALAR> scalar_field_buffer;
    ARRAY<T_FACE_ARRAYS_SCALAR> face_scalar_field_buffer;

    ARRAY<INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>*>* incompressible_fluid_containers;

    CHIMERA_GRID(RIGID_GRID_COLLECTION<T_GRID>& global_grid_input,ARRAY<T_GRID>& grids_array,ARRAY<FRAME<TV> >& grid_frames,int number_of_ghost_cells_input,bool use_mpi_input=false,ARRAY<INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>*>* incompressible_fluid_containers_input=0,ARRAY<int>* number_divisions_per_grid_input=0);
    ~CHIMERA_GRID();

    static int Number_Of_Processors();

    int Get_Send_Tag(int send_grid,int recv_grid,int increment=0) const
    {
        STATIC_ASSERT(TV_INT::m<=3);int tag=0;
        tag=(send_grid-1)*number_of_global_grids+recv_grid-1;
        return tag+increment;
    }

    int Get_Recv_Tag(int send_grid,int recv_grid,int increment=0) const
    {
        return Get_Send_Tag(send_grid,recv_grid,increment);
    }

    void Get_Send_Recv_Grid(int tag,int& send_grid,int& recv_grid,int increment=0) const
    {
        if(tag-increment<0){send_grid=0;recv_grid=0;return;}
        send_grid=int((tag-increment)/number_of_global_grids)+1;
        recv_grid=(tag-increment)%number_of_global_grids+1;
        if(send_grid>number_of_global_grids || send_grid<1 || recv_grid>number_of_global_grids || recv_grid<1){
            send_grid=0;recv_grid=0;}
    }

    void Initialize_Grid_Index_Maps()
    {
        //difine local_to_global_grid_index_map
        for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++){
            if(grid_ranks(grid_index)==rank) local_to_global_grid_index_map.Append(grid_index);}
        //difine global_to_local_grid_index_map
        global_to_local_grid_index_map.Resize(number_of_global_grids);
        ARRAYS_COMPUTATIONS::Fill(global_to_local_grid_index_map,0);
        for(int local_grid_index=1;local_grid_index<=local_to_global_grid_index_map.Size();local_grid_index++)
            global_to_local_grid_index_map(local_to_global_grid_index_map(local_grid_index))=local_grid_index;
    }

    void Sort_Grid_Sizes()
    {
        ARRAY<PAIR<int,T> > grid_sizes_to_sort(number_of_local_grids);
        for(int i=1;i<=number_of_local_grids;i++){
            grid_sizes_to_sort(i)=PAIR<int,T>(i,local_grid.Rigid_Grid(i).grid.dX.Product());}
        for(int i=number_of_local_grids;i>=2;i--){
            for(int j=1;j<=i-1;j++){
                if(grid_sizes_to_sort(j).y>grid_sizes_to_sort(j+1).y){
                    PAIR<int,T> tmp=grid_sizes_to_sort(j+1);
                    grid_sizes_to_sort(j+1)=grid_sizes_to_sort(j);
                    grid_sizes_to_sort(j)=tmp;}}}
        grids_sorted.Resize(number_of_local_grids);
        for(int i=1;i<=number_of_local_grids;i++){
            grids_sorted(i)=grid_sizes_to_sort(i).x;}
    }

    void Synchronize_Grid_Position()
    {
        for(int grid_index=1;grid_index<=local_grid.particles.array_collection->Size();grid_index++){
            local_grid.Rigid_Grid(grid_index).X()=global_grid.Rigid_Grid(local_to_global_grid_index_map(grid_index)).X();
            local_grid.Rigid_Grid(grid_index).Rotation()=global_grid.Rigid_Grid(local_to_global_grid_index_map(grid_index)).Rotation();
            local_grid.Rigid_Grid(grid_index).previous_state=global_grid.Rigid_Grid(local_to_global_grid_index_map(grid_index)).previous_state;}
    }

    bool Boundary_Regions_Found_Cell(int number_of_ghost_cells_new)
    {
        return boundary_regions_found_cell.x && (boundary_regions_found_cell.y==number_of_ghost_cells_new);
    }

    bool Boundary_Regions_Found_Face(int number_of_ghost_cells_new)
    {
        return boundary_regions_found_face.x && (boundary_regions_found_face.y==number_of_ghost_cells_new);
    }

    bool Overlap_Regions_Found_Cell(int number_of_ghost_cells_new)
    {
        return overlap_regions_found_cell.x && (overlap_regions_found_cell.y==number_of_ghost_cells_new);
    }

    bool Overlap_Regions_Found_Face(int number_of_ghost_cells_new)
    {
        return overlap_regions_found_face.x && (overlap_regions_found_face.y==number_of_ghost_cells_new);
    }

    void Set_To_Use_Previous_Grid_Frame_Cell()
    {
        use_previous_grid_frame_cell=true;
        Clean_Buffers_Cell();
    }
    void Set_To_Use_Previous_Grid_Frame_Face()
    {
        use_previous_grid_frame_face=true;
        Clean_Buffers_Face();
    }

    void Set_To_Use_Current_Grid_Frame_Cell()
    {
        use_previous_grid_frame_cell=false;
        Clean_Buffers_Cell();
    }
    void Set_To_Use_Current_Grid_Frame_Face()
    {
        use_previous_grid_frame_face=false;
        Clean_Buffers_Face();
    }

    void Adjustment_After_Updating_Grids()
    {
        Set_To_Use_Previous_Grid_Frame_Cell();
        Set_To_Use_Previous_Grid_Frame_Face();
    }

    void Adjustment_After_Advecting_Scalar_Fields()
    {
        Set_To_Use_Current_Grid_Frame_Cell();
    }

    void Adjustment_After_Advecting_Velocities()
    {
        Set_To_Use_Current_Grid_Frame_Face();
    }

    void Copy_Data_To_Buffers_Cell(ARRAY<T_ARRAYS_SCALAR*> scalar_array)
    {
        scalar_field_buffer.Resize(scalar_array.m);
        for(int grid_index=1;grid_index<=scalar_array.m;grid_index++){
            scalar_field_buffer(grid_index)=*(scalar_array(grid_index));}
    }

    void Copy_Data_To_Buffers_Face(ARRAY<T_FACE_ARRAYS_SCALAR*> face_scalar_array)
    {
        face_scalar_field_buffer.Resize(face_scalar_array.m);
        for(int grid_index=1;grid_index<=face_scalar_array.m;grid_index++){
            face_scalar_field_buffer(grid_index)=*(face_scalar_array(grid_index));}
    }

    void Get_Pointer_Of_Incompressible_Fluid_Containers(ARRAY<INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>*>* incompressible_fluid_containers_input)
    {
        incompressible_fluid_containers=incompressible_fluid_containers_input;
    }

    int Find_Current_Grid_Local_Index(const T_GRID& grid)
    {
        int current_grid_index=0;
        for(int grid_index=1;grid_index<=number_of_local_grids;grid_index++) if(&grid==&(local_grid.Rigid_Grid(grid_index).grid)){
            current_grid_index=grid_index;
            break;}
        if(current_grid_index) return current_grid_index;
        for(int grid_index=1;grid_index<=number_of_local_grids;grid_index++) if(&grid==&((*incompressible_fluid_containers)(grid_index)->particle_levelset_evolution.grid)){
            current_grid_index=grid_index;
            break;}
        return current_grid_index;
    }

    bool Local_Grid(const int grid_index) const
    {
        return global_to_local_grid_index_map(grid_index)==0?false:true;
    }

//######################################################################
    static void Determine_Divisions_Of_Grids(ARRAY<T_GRID>& grids,ARRAY<int>& divisions);
    void Initialize_Grids(ARRAY<T_GRID>& grids_array,ARRAY<FRAME<TV> >& grid_frames,ARRAY<int>& number_divisions_per_grid);
    void Initialize_Global_Grid(ARRAY<T_GRID>& new_grids,ARRAY<FRAME<TV> >& new_frames);
    void Restrict_Local_Grid();
    void Update_Oriented_Boxes();
    bool Grid_Intersect_Local_Domain(int one_grid_index);
    int Find_Finest_Grid(const TV& location,bool use_previous_state=false);
    void Clean_Buffers_Cell();
    void Clean_Buffers_Face();
    void Clean_Buffers_Boundary_Cell();
    void Clean_Buffers_Overlap_Cell();
    void Clean_Buffers_Boundary_Face();
    void Clean_Buffers_Overlap_Face();
    void Find_Boundary_Regions(ARRAY<RANGE<VECTOR<int,1> > >& regions,const RANGE<VECTOR<int,1> >& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,const bool include_ghost_regions,const GRID<VECTOR<T,1> >& grid,int global_grid_index) const;
    void Find_Boundary_Regions(ARRAY<RANGE<VECTOR<int,2> > >& regions,const RANGE<VECTOR<int,2> >& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,const bool include_ghost_regions,const GRID<VECTOR<T,2> >& grid,int global_grid_index) const;
    void Find_Boundary_Regions(ARRAY<RANGE<VECTOR<int,3> > >& regions,const RANGE<VECTOR<int,3> >& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,const bool include_ghost_regions,const GRID<VECTOR<T,3> >& grid,int global_grid_index) const;
    void Find_Boundary_Regions_Cell(const int number_of_ghost_cells);
    void Find_Boundary_Regions_Face(const int number_of_ghost_cells);
    void Find_Overlap_Regions_Cell(const int& number_of_ghost_cells,ARRAY<T_ARRAYS_BOOL>* fixed_mask=0);
    void Find_Overlap_Regions_Face(const int& number_of_ghost_cells,ARRAY<T_FACE_ARRAYS_BOOL>* fixed_mask=0);
    void Exchange_Overlap_Cell_Data_Recursively(ARRAY<T_ARRAYS_SCALAR*>& data,const int& number_of_ghost_cells,bool only_overwrite_nan=false,ARRAY<T_ARRAYS_BOOL>* fixed_mask=0);
    void Exchange_Overlap_Face_Data_Recursively(ARRAY<T_FACE_ARRAYS_T2*>& data,const int& number_of_ghost_cells,ARRAY<T_FACE_ARRAYS_BOOL>* fixed_mask=0);
    void Exchange_Overlap_Cell_Data(ARRAY<T_ARRAYS_SCALAR*>& data,const int& number_of_ghost_cells,bool only_overwrite_nan=false,ARRAY<T_ARRAYS_BOOL>* fixed_mask=0);
    void Exchange_Overlap_Face_Data(ARRAY<T_FACE_ARRAYS_T2*>& data,const int& number_of_ghost_cells,ARRAY<T_FACE_ARRAYS_BOOL>* fixed_mask=0);
    void Exchange_Boundary_Cell_Data(ARRAY<T_ARRAYS_SCALAR*>& data,const int& number_of_ghost_cells,bool only_overwrite_nan=false);
    void Exchange_Boundary_Face_Data(ARRAY<T_FACE_ARRAYS_T2*>& data,const int& number_of_ghost_cells);
    void Pack_Send_Data_Cell(CELL_PACKAGE_TYPE& cell_package,T_ARRAYS_T2& data,int send_grid,int recv_grid);
    void Pack_Send_Data_Face(FACE_PACKAGE_TYPE& face_package,T_FACE_ARRAYS_T2& data,int send_grid,int recv_grid);
    void Send_Data_To_Recv_Grid_Cell(CELL_PACKAGE_TYPE& cell_package,int send_grid,int recv_grid,int tag);
    void Send_Data_To_Recv_Grid_Face(FACE_PACKAGE_TYPE& face_package,int send_grid,int recv_grid,int tag);
    void Recv_Data_From_Send_Grid_Cell(ARRAY<PAIR<CELL_PACKAGE_TYPE,int> >& packages_cell_send,CELL_PACKAGE_TYPE& cell_package,int send_grid,int recv_grid,int tag);
    void Recv_Data_From_Send_Grid_Face(ARRAY<PAIR<FACE_PACKAGE_TYPE,int> >& packages_face_send,FACE_PACKAGE_TYPE& face_package,int send_grid,int recv_grid,int tag);
    void Unpack_Recv_Data_Cell(CELL_PACKAGE_TYPE& cell_package,T_ARRAYS_T2& data,bool only_overwrite_nan=false);
    void Unpack_Recv_Data_Face(FACE_PACKAGE_TYPE& face_package,T_FACE_ARRAYS_T2& data);
    void Initialize_Communicator();
    void Synchronize_Dt(T& dt) const;
    void Synchronize_Done(bool& done) const;
    void Synchronize_Fix_Masks(ARRAY<T_ARRAYS_BOOL>& fixed_mask);
    void Exchange_Boundary_Cell_Arrays(ARRAY<ARRAY<TV_INT> >& boundary_cell_indices,ARRAY<ARRAY<bool> >& is_boundary_grid);
    void Exchange_Boundary_Voronoi_Faces(ARRAY<ARRAY<VORONOI_FACE> >& boundary_cell_voronoi_faces_to_send,ARRAY<ARRAY<VORONOI_FACE> >& boundary_cell_voronoi_faces_to_receive);
    void Exchange_Voronoi_Face_Positions_For_Interpolation(ARRAY<ARRAY<PAIR<TV,TV> > >& voronoi_face_positions_interpolated_from_others,ARRAY<ARRAY<PAIR<TV,TV> > >& voronoi_face_positions_interpolated_here);
    void Exchange_Voronoi_Face_Interpolated_Values(ARRAY<ARRAY<T> >& voronoi_face_values_interpolated_here,ARRAY<ARRAY<T> >& voronoi_face_values_interpolated_from_others);
    template<class T3> void Exchange_Boundary_Scalar_Values(ARRAY<ARRAY<T3> >& boundary_scalar_values,const ARRAY<ARRAY<bool> >& is_boundary_grid);
    void Exchange_Overlap_Packages_Cell_Send(ARRAY<CELL_PACKAGE_TYPE>& boundary_packages);
    void Exchange_Overlap_Packages_Face_Send(ARRAY<FACE_PACKAGE_TYPE>& boundary_packages);
    void Parallel_Solve(SPARSE_MATRIX_FLAT_NXN<T>& A_pcg,VECTOR_ND<T>& x_pcg,VECTOR_ND<T>& rhs_pcg,T tolerance,SPARSE_MATRIX_PARTITION& partition,ARRAY<ARRAY<bool> >& is_boundary_grid,int maximum_iterations,bool use_incomplete_cholesky);
    GRID_CELL_INDEX Find_Real_Grid_Cell_Index(const GRID_CELL_INDEX& index);
    TV_INT Unsplit_Cell_Index(const GRID_CELL_INDEX& grid_cell_index);
    TV_INT Split_Cell_Index(const int grid_index,const TV_INT& unsplit_cell_index);
    int Parent_Grid(const int grid_index);
//######################################################################
};
}
#endif
