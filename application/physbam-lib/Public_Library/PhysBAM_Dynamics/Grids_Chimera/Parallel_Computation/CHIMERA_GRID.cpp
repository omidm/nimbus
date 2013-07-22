//#####################################################################
// Copyright 2011, Linhai Qiu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CHIMERA_GRID
//#####################################################################
#include <PhysBAM_Dynamics/Grids_Chimera/Parallel_Computation/CHIMERA_GRID.h>
#include <PhysBAM_Dynamics/Grids_Chimera/Grids_Chimera_Boundaries/BOUNDARY_CHIMERA.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_UNSTRUCTURE_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Math_Tools/Is_NaN.h>
using namespace PhysBAM;

//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T2> CHIMERA_GRID<T_GRID,T2>::
CHIMERA_GRID(RIGID_GRID_COLLECTION<T_GRID>& global_grid_input,ARRAY<T_GRID>& grids_array,ARRAY<FRAME<TV> >& grid_frames,int number_of_ghost_cells_input,bool use_mpi_input,ARRAY<INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>*>* incompressible_fluid_containers_input,ARRAY<int>* number_divisions_per_grid_input)
    :local_grid(global_grid_input),initial_grids(grids_array),number_of_ghost_cells(number_of_ghost_cells_input),use_mpi(use_mpi_input),rank(0),use_previous_grid_frame_cell(false),use_previous_grid_frame_face(false),incompressible_fluid_containers(incompressible_fluid_containers_input)
{
#ifndef USE_MPI
    LOG::cout << "USE_MPI and USE_PTHREADS should be set" << std::endl;
    exit(0);
#endif
#ifndef USE_PTHREADS
    LOG::cout << "USE_MPI and USE_PTHREADS should be set" << std::endl;
    exit(0);
#endif

    if(!use_mpi){
        number_of_processors=1;
        Initialize_Grids(initial_grids,grid_frames,*number_divisions_per_grid_input);
        Clean_Buffers_Cell();
        Clean_Buffers_Face();
        Sort_Grid_Sizes();
        return;}
#ifdef USE_MPI
    else{
        number_of_processors=MPI::COMM_WORLD.Get_size();
        Initialize_Communicator();
        assert(number_divisions_per_grid_input);
        Initialize_Grids(initial_grids,grid_frames,*number_divisions_per_grid_input);
        mpi_packages_buffers=new ARRAY<MPI_PACKAGE>;
        mpi_requests_buffers=new ARRAY<MPI::Request>;
        Clean_Buffers_Cell();
        Clean_Buffers_Face();
        Sort_Grid_Sizes();}
#endif
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID,class T2> CHIMERA_GRID<T_GRID,T2>::
~CHIMERA_GRID()
{
#ifdef USE_MPI
    if(comm){comm->Free();delete comm;}
#endif
}
//#####################################################################
// Determine_Divisions_Of_Grids
//#####################################################################
bool Is_Prime_Number(int number)
{
    if(number<=1) return false;
    for(int i=2;i<=sqrt((float)number);i++){
        if(!(number%i)) return false;}
    return true;
}
void Minimize_2D_Surface_Area(const int n_processes,const int m,const int n,int& count_x)
{
    int min_area=INT_MAX;
    for(int try_count_x=1;try_count_x<=n_processes;try_count_x++) if(n_processes%try_count_x==0){
        int try_area=n*(try_count_x-1)+m*(n_processes/try_count_x-1);
        if(try_area < min_area){count_x=try_count_x;min_area=try_area;}}
}
template<class T>
void Split_Single_Grid(const GRID<VECTOR<T,1> >& grid,const int n_processes,VECTOR<int,1>& processes_per_dimension)
{}
template<class T>
void Split_Single_Grid(const GRID<VECTOR<T,2> >& grid,const int n_processes,VECTOR<int,2>& processes_per_dimension)
{
    Minimize_2D_Surface_Area(n_processes,grid.counts.x,grid.counts.y,processes_per_dimension.x);
    processes_per_dimension.y=n_processes/processes_per_dimension.x;
}
template<class T>
void Split_Single_Grid(const GRID<VECTOR<T,3> >& grid,const int n_processes,VECTOR<int,3>& processes_per_dimension)
{
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;
    int minimum_surface_area=INT_MAX;VECTOR<int,3> try_count;
    for(try_count.z=1;try_count.z<=n_processes;try_count.z++) if(n_processes%try_count.z==0){
        Minimize_2D_Surface_Area(n_processes/try_count.z,m,n,try_count.x);
        try_count.y=n_processes/(try_count.x*try_count.z);
        int surface_area=(try_count.x-1)*(n*mn)+(try_count.y-1)*(m*mn)+(try_count.z-1)*(m*n);
        if(surface_area<minimum_surface_area){processes_per_dimension=try_count;minimum_surface_area=surface_area;}}
}
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Determine_Divisions_Of_Grids(ARRAY<T_GRID>& grids,ARRAY<int>& divisions)
{
    int free_ranks=Number_Of_Processors()-grids.Size();
    while(free_ranks>0){
        int grid_index_maximum=1;
        for(int grid_index=1;grid_index<=grids.Size();grid_index++){
            if((grids(grid_index).Counts().Product()/divisions(grid_index))>(grids(grid_index_maximum).Counts().Product()/divisions(grid_index_maximum)))
                grid_index_maximum=grid_index;}
        divisions(grid_index_maximum)++;
        free_ranks--;}
    for(int grid_index=1;grid_index<=grids.Size();grid_index++){
        if(divisions(grid_index)%2!=0 && Is_Prime_Number(divisions(grid_index))){
            if(free_ranks>0){
                divisions(grid_index)++;
                free_ranks--;}
            else{
                divisions(grid_index)--;
                free_ranks++;}}}
}
//#####################################################################
// Initialize_Grids
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Initialize_Grids(ARRAY<T_GRID>& grids_array,ARRAY<FRAME<TV> >& grid_frames,ARRAY<int>& number_divisions_per_grid)
{
    ARRAY<T_GRID> new_grids;
    ARRAY<FRAME<TV> > new_frames;
    ARRAY<bool> contains_splitted_grid(number_of_processors);contains_splitted_grid.Fill(false);
    //split grids
    if(use_mpi){
#ifdef USE_MPI
        int current_rank=0;
        mpi_split_grids.Resize(grids_array.m);
        mpi_grid_coordinate_to_global_grid_index_map.Resize(grids_array.m);
        for(int grid_index=1;grid_index<=grids_array.m;grid_index++){
            MPI_UNIFORM_GRID<T_GRID>*& mpi_grid=mpi_split_grids(grid_index);
            mpi_grid=new MPI_UNIFORM_GRID<T_GRID>(grids_array(grid_index),number_of_ghost_cells,true);
            mpi_grid->number_of_processes=number_divisions_per_grid(grid_index);
            mpi_grid->global_grid=mpi_grid->local_grid.Get_Regular_Grid();
            TV_INT processes_per_dimension=TV_INT();
            Split_Single_Grid<T>(grids_array(grid_index),number_divisions_per_grid(grid_index),processes_per_dimension);
            mpi_grid->Split_Grid(processes_per_dimension);
            mpi_grid->process_ranks.Resize(mpi_grid->process_grid.Domain_Indices(1));ARRAYS_COMPUTATIONS::Fill(mpi_grid->process_ranks.array,MPI::PROC_NULL);
            bool this_mpi_grid_belongs_to_this_process=false;
            for(NODE_ITERATOR iterator(mpi_grid->process_grid);iterator.Valid();iterator.Next()){
                if(contains_splitted_grid(current_rank+1)) PHYSBAM_FATAL_ERROR("Two subgrids from splitting cannot be in the same process!");
                TV_INT node_index=iterator.Node_Index();
                T_GRID grid=mpi_grid->Restrict_Grid(node_index);
                mpi_grid->process_ranks(node_index)=current_rank;
                grid_ranks.Append(current_rank);
                if(mpi_grid->process_grid.counts.Product()>1) contains_splitted_grid(current_rank+1)=true;
                new_grids.Append(grid);
                new_frames.Append(grid_frames(grid_index));
                mpi_grid_indices.Append(grid_index);
                if(current_rank==rank){this_mpi_grid_belongs_to_this_process=true;mpi_grid->coordinates=node_index;mpi_grid->local_grid=grid;mpi_grid->comm=comm;}
                current_rank=(current_rank+1)%number_of_processors;}
            if(this_mpi_grid_belongs_to_this_process) mpi_grid->Initialize_Topology();
            mpi_grid_coordinate_to_global_grid_index_map(grid_index).Resize(mpi_grid->process_grid.Domain_Indices());
            for(NODE_ITERATOR iterator(mpi_grid->process_grid);iterator.Valid();iterator.Next()){
                TV_INT node_index=iterator.Node_Index();
                ARRAY<int> side_neighbor_ranks(T_GRID::number_of_neighbors_per_node);
                for(int n=1;n<=T_GRID::number_of_neighbors_per_node;n++)
                    side_neighbor_ranks(n)=mpi_grid->process_ranks(T_GRID::Node_Neighbor(node_index,n));
                side_neighbor_ranks_chimera.Append(side_neighbor_ranks);
                ARRAY<int> all_neighbor_ranks(T_GRID::number_of_one_ring_neighbors_per_cell);
                for(int n=1;n<=T_GRID::number_of_one_ring_neighbors_per_cell;n++)
                    all_neighbor_ranks(n)=mpi_grid->process_ranks(T_GRID::One_Ring_Neighbor(node_index,n));
                all_neighbor_ranks_chimera.Append(all_neighbor_ranks);
                all_coordinates.Append(node_index);
                mpi_grid_coordinate_to_global_grid_index_map(grid_index)(node_index)=all_coordinates.Size();}}
#endif
    }
    else{
        new_grids=grids_array;new_frames=grid_frames;
        for(int grid_index=1;grid_index<=grids_array.m;grid_index++) grid_ranks.Append(0);}
    //initialize rigid grid collection
    Initialize_Global_Grid(new_grids,new_frames);
    //initialize local and global grids index maps
    Initialize_Grid_Index_Maps();
    //Restrict local grid
    Restrict_Local_Grid();
}
//#####################################################################
// Initialize_Global_Grid
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Initialize_Global_Grid(ARRAY<T_GRID>& new_grids,ARRAY<FRAME<TV> >& new_frames)
{
    for(int grid_index=1;grid_index<=new_grids.m;grid_index++){
        global_grid.Add_Rigid_Grid(new_frames(grid_index).t,new_frames(grid_index).r);
        global_grid.Rigid_Grid(grid_index).grid=new_grids(grid_index);
        global_grid.Rigid_Grid(grid_index).previous_state=global_grid.Rigid_Grid(grid_index).Frame();}
    number_of_global_grids=global_grid.particles.array_collection->Size();

    interpolation_domains.Resize(number_of_global_grids);
    previous_oriented_boxes.Resize(number_of_global_grids);
    oriented_boxes.Resize(number_of_global_grids);

    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++){
        interpolation_domains(grid_index)=global_grid.Rigid_Grid(grid_index).grid.Domain();
        TV dx=global_grid.Rigid_Grid(grid_index).grid.DX();
        if(!use_mpi) interpolation_domains(grid_index).Change_Size((T)-.5*dx);
#ifdef USE_MPI
        else{
            for(int n=1;n<=T_GRID::number_of_neighbors_per_node;n++)
                if(side_neighbor_ranks_chimera(grid_index)(n)==MPI::PROC_NULL){
                    TV_INT direction=T_GRID::Node_Neighbor(TV_INT(),n);
                    int axis=direction.Dominant_Axis();
                    if(direction(axis)==1) interpolation_domains(grid_index).max_corner(axis)-=(T).5*dx(axis);
                    else interpolation_domains(grid_index).min_corner(axis)+=(T).5*dx(axis);}}
#endif
        const FRAME<TV>& previous_frame=global_grid.Rigid_Grid(grid_index).previous_state.frame;
        const FRAME<TV>& frame=global_grid.Rigid_Grid(grid_index).Frame();
        previous_oriented_boxes(grid_index)=ORIENTED_BOX<TV>(interpolation_domains(grid_index),previous_frame);
        oriented_boxes(grid_index)=ORIENTED_BOX<TV>(interpolation_domains(grid_index),frame);
    }
}
//#####################################################################
// Update_Oriented_Boxes
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Update_Oriented_Boxes()
{
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++){
        const FRAME<TV>& previous_frame=global_grid.Rigid_Grid(grid_index).previous_state.frame;
        const FRAME<TV>& frame=global_grid.Rigid_Grid(grid_index).Frame();
        previous_oriented_boxes(grid_index)=ORIENTED_BOX<TV>(interpolation_domains(grid_index),previous_frame);
        oriented_boxes(grid_index)=ORIENTED_BOX<TV>(interpolation_domains(grid_index),frame);}
}
//#####################################################################
// Restrict_Local_Grid
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Restrict_Local_Grid()
{
    number_of_local_grids=0;
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++) if(global_to_local_grid_index_map(grid_index)){
        if(incompressible_fluid_containers){
            INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>* new_container=new INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>(local_grid,number_of_ghost_cells);
            incompressible_fluid_containers->Append(new_container);number_of_local_grids++;
            RIGID_GRID<T_GRID>& local_rigid_grid=new_container->rigid_grid;
            RIGID_GRID<T_GRID>& global_rigid_grid=global_grid.Rigid_Grid(grid_index);
            local_rigid_grid.X()=global_rigid_grid.X();local_rigid_grid.Rotation()=global_rigid_grid.Rotation();
            local_rigid_grid.grid=global_rigid_grid.grid;
            local_rigid_grid.previous_state=local_rigid_grid.Frame();}
        else{
            RIGID_GRID<T_GRID>& global_rigid_grid=global_grid.Rigid_Grid(grid_index);
            local_grid.Add_Rigid_Grid(global_rigid_grid.X(),global_rigid_grid.Rotation());
            number_of_local_grids++;
            local_grid.Rigid_Grid(number_of_local_grids).grid=global_rigid_grid.grid;
            local_grid.Rigid_Grid(number_of_local_grids).previous_state=local_grid.Rigid_Grid(grid_index).Frame();}}
}
//#####################################################################
// Grid_Intersect_Local_Domain
//#####################################################################
template<class T_GRID,class T2> bool CHIMERA_GRID<T_GRID,T2>::
Grid_Intersect_Local_Domain(int one_grid_index)
{
    ORIENTED_BOX<TV> one_grid_box(global_grid.Rigid_Grid(one_grid_index).grid.Domain(),global_grid.Rigid_Grid(one_grid_index).Frame());
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++)
        if(Local_Grid(grid_index)){
            if(grid_index==one_grid_index) return true;
            ORIENTED_BOX<TV> grid_box(global_grid.Rigid_Grid(grid_index).grid.Domain(),global_grid.Rigid_Grid(grid_index).Frame());
            if(grid_index<one_grid_index?grid_box.Intersection(one_grid_box):one_grid_box.Intersection(grid_box)) return true;}
    return false;
}
//#####################################################################
// Find_Finest_Grid
//#####################################################################
template<class T_GRID,class T2> int CHIMERA_GRID<T_GRID,T2>::
Find_Finest_Grid(const TV& location,bool use_previous_state)
{
    int finest_grid=0;
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++){
        const ORIENTED_BOX<TV>& box=use_previous_state?previous_oriented_boxes(grid_index):oriented_boxes(grid_index);
        bool inside=box.Lazy_Inside(location);
        if(inside && (!finest_grid || (global_grid.Rigid_Grid(grid_index).grid.Cell_Size()<global_grid.Rigid_Grid(finest_grid).grid.Cell_Size()) || (global_grid.Rigid_Grid(grid_index).grid.Cell_Size()==global_grid.Rigid_Grid(finest_grid).grid.Cell_Size() && grid_index<finest_grid) ) )
            finest_grid=grid_index;}
    return finest_grid;
}
//#####################################################################
// Clean_Buffers_Cell
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Clean_Buffers_Cell()
{
    Clean_Buffers_Boundary_Cell();
    Clean_Buffers_Overlap_Cell();
}
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Clean_Buffers_Boundary_Cell()
{
    boundary_regions_found_cell.x=false;
    boundary_packages_cell_send.Clean_Memory();
    boundary_packages_cell_recv.Clean_Memory();
}
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Clean_Buffers_Overlap_Cell()
{
    overlap_regions_found_cell.x=false;
    overlap_packages_cell_send.Clean_Memory();
    overlap_packages_cell_recv.Clean_Memory();
}

//#####################################################################
// Clean_Buffers_Face
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Clean_Buffers_Face()
{
    Clean_Buffers_Boundary_Face();
    Clean_Buffers_Overlap_Face();
}
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Clean_Buffers_Boundary_Face()
{
    boundary_regions_found_face.x=false;
    boundary_packages_face_send.Clean_Memory();
    boundary_packages_face_recv.Clean_Memory();
}
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Clean_Buffers_Overlap_Face()
{
    overlap_regions_found_face.x=false;
    overlap_packages_face_send.Clean_Memory();
    overlap_packages_face_recv.Clean_Memory();
}
//#####################################################################
// Function Find_Boundary_Regions
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Find_Boundary_Regions(ARRAY<RANGE<VECTOR<int,1> > >& regions,const RANGE<VECTOR<int,1> >& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,
    const bool include_ghost_regions,const GRID<VECTOR<T,1> >& grid,int global_grid_index) const
{
    RANGE<VECTOR<int,1> > box=RANGE<VECTOR<int,1> >(1,grid.numbers_of_cells.x)+sentinels;
    RANGE<VECTOR<int,1> > band_x=band;
    if(skip_common_boundary && band.max_corner.x>0){
        if(band.min_corner.x<0) PHYSBAM_NOT_IMPLEMENTED();
        if(sentinels.max_corner.x) band_x+=VECTOR<int,1>(1);}
    regions.Clean_Memory();
    regions.Append(RANGE<VECTOR<int,1> >(box.min_corner.x,box.min_corner.x)+band_x);
    regions.Append(RANGE<VECTOR<int,1> >(box.max_corner.x,box.max_corner.x)-band_x);
}
//#####################################################################
// Function Find_Boundary_Regions
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Find_Boundary_Regions(ARRAY<RANGE<VECTOR<int,2> > >& regions,const RANGE<VECTOR<int,2> >& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,
    const bool include_ghost_regions,const GRID<VECTOR<T,2> >& grid,int global_grid_index) const
{
    RANGE<VECTOR<int,2> > box=RANGE<VECTOR<int,2> >(1,grid.numbers_of_cells.x,1,grid.numbers_of_cells.y)+sentinels;
    RANGE<VECTOR<int,2> > band_x(band.min_corner.x,band.max_corner.x,0,0),band_y(0,0,band.min_corner.x,band.max_corner.x);
    if(skip_common_boundary && band.max_corner.x>0){
        if(band.min_corner.x<0) PHYSBAM_NOT_IMPLEMENTED();
        if(sentinels.max_corner.x) band_x+=VECTOR<int,2>(1,0);
        if(sentinels.max_corner.y) band_y+=VECTOR<int,2>(0,1);}
    regions.Clean_Memory();
    if(!include_corners){ // add in side order
        regions.Append(RANGE<VECTOR<int,2> >(box.min_corner.x,box.min_corner.x,box.min_corner.y,box.max_corner.y)+band_x);
        regions.Append(RANGE<VECTOR<int,2> >(box.max_corner.x,box.max_corner.x,box.min_corner.y,box.max_corner.y)-band_x);
        regions.Append(RANGE<VECTOR<int,2> >(box.min_corner.x,box.max_corner.x,box.min_corner.y,box.min_corner.y)+band_y);
        regions.Append(RANGE<VECTOR<int,2> >(box.min_corner.x,box.max_corner.x,box.max_corner.y,box.max_corner.y)-band_y);}
    else{ // add in one ring neighbors order
        regions.Append(RANGE<VECTOR<int,2> >(box.min_corner.x,box.min_corner.x,box.min_corner.y,box.min_corner.y)+band_x+band_y);
        regions.Append(RANGE<VECTOR<int,2> >(box.min_corner.x,box.max_corner.x,box.min_corner.y,box.min_corner.y)+band_y);
        regions.Append(RANGE<VECTOR<int,2> >(box.max_corner.x,box.max_corner.x,box.min_corner.y,box.min_corner.y)-band_x+band_y);
        regions.Append(RANGE<VECTOR<int,2> >(box.min_corner.x,box.min_corner.x,box.min_corner.y,box.max_corner.y)+band_x);
        regions.Append(RANGE<VECTOR<int,2> >(box.max_corner.x,box.max_corner.x,box.min_corner.y,box.max_corner.y)-band_x);
        regions.Append(RANGE<VECTOR<int,2> >(box.min_corner.x,box.min_corner.x,box.max_corner.y,box.max_corner.y)+band_x-band_y);
        regions.Append(RANGE<VECTOR<int,2> >(box.min_corner.x,box.max_corner.x,box.max_corner.y,box.max_corner.y)-band_y);
        regions.Append(RANGE<VECTOR<int,2> >(box.max_corner.x,box.max_corner.x,box.max_corner.y,box.max_corner.y)-band_x-band_y);}
}
//#####################################################################
// Function Find_Boundary_Regions
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Find_Boundary_Regions(ARRAY<RANGE<VECTOR<int,3> > >& regions,const RANGE<VECTOR<int,3> >& sentinels,const bool skip_common_boundary,const RANGE<VECTOR<int,1> >& band,const bool include_corners,
    const bool include_ghost_regions,const GRID<VECTOR<T,3> >& grid,int global_grid_index) const
{
    RANGE<VECTOR<int,3> > box=RANGE<VECTOR<int,3> >(1,grid.numbers_of_cells.x,1,grid.numbers_of_cells.y,1,grid.numbers_of_cells.z)+sentinels;
    RANGE<VECTOR<int,3> > band_x(band.min_corner.x,band.max_corner.x,0,0,0,0),band_y(0,0,band.min_corner.x,band.max_corner.x,0,0),band_z(0,0,0,0,band.min_corner.x,band.max_corner.x);
    if(skip_common_boundary && band.max_corner.x>0){
        if(band.min_corner.x<0) PHYSBAM_NOT_IMPLEMENTED();
        if(sentinels.max_corner.x) band_x+=VECTOR<int,3>(1,0,0);
        if(sentinels.max_corner.y) band_y+=VECTOR<int,3>(0,1,0);
        if(sentinels.max_corner.z) band_z+=VECTOR<int,3>(0,0,1);}
    regions.Clean_Memory();
    if(!include_corners){ // add in side order
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.min_corner.x,box.min_corner.y,box.max_corner.y,box.min_corner.z,box.max_corner.z)+band_x);
        regions.Append(RANGE<VECTOR<int,3> >(box.max_corner.x,box.max_corner.x,box.min_corner.y,box.max_corner.y,box.min_corner.z,box.max_corner.z)-band_x);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.max_corner.x,box.min_corner.y,box.min_corner.y,box.min_corner.z,box.max_corner.z)+band_y);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.max_corner.x,box.max_corner.y,box.max_corner.y,box.min_corner.z,box.max_corner.z)-band_y);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.max_corner.x,box.min_corner.y,box.max_corner.y,box.min_corner.z,box.min_corner.z)+band_z);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.max_corner.x,box.min_corner.y,box.max_corner.y,box.max_corner.z,box.max_corner.z)-band_z);}
    else{ // add in one ring neighbors order
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.min_corner.x,box.min_corner.y,box.min_corner.y,box.min_corner.z,box.min_corner.z)+band_x+band_y+band_z);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.max_corner.x,box.min_corner.y,box.min_corner.y,box.min_corner.z,box.min_corner.z)+band_y+band_z);
        regions.Append(RANGE<VECTOR<int,3> >(box.max_corner.x,box.max_corner.x,box.min_corner.y,box.min_corner.y,box.min_corner.z,box.min_corner.z)-band_x+band_y+band_z);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.min_corner.x,box.min_corner.y,box.max_corner.y,box.min_corner.z,box.min_corner.z)+band_x+band_z);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.max_corner.x,box.min_corner.y,box.max_corner.y,box.min_corner.z,box.min_corner.z)+band_z); 
        regions.Append(RANGE<VECTOR<int,3> >(box.max_corner.x,box.max_corner.x,box.min_corner.y,box.max_corner.y,box.min_corner.z,box.min_corner.z)-band_x+band_z);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.min_corner.x,box.max_corner.y,box.max_corner.y,box.min_corner.z,box.min_corner.z)+band_x-band_y+band_z);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.max_corner.x,box.max_corner.y,box.max_corner.y,box.min_corner.z,box.min_corner.z)-band_y+band_z);
        regions.Append(RANGE<VECTOR<int,3> >(box.max_corner.x,box.max_corner.x,box.max_corner.y,box.max_corner.y,box.min_corner.z,box.min_corner.z)-band_x-band_y+band_z);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.min_corner.x,box.min_corner.y,box.min_corner.y,box.min_corner.z,box.max_corner.z)+band_x+band_y);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.max_corner.x,box.min_corner.y,box.min_corner.y,box.min_corner.z,box.max_corner.z)+band_y);
        regions.Append(RANGE<VECTOR<int,3> >(box.max_corner.x,box.max_corner.x,box.min_corner.y,box.min_corner.y,box.min_corner.z,box.max_corner.z)-band_x+band_y);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.min_corner.x,box.min_corner.y,box.max_corner.y,box.min_corner.z,box.max_corner.z)+band_x);
        regions.Append(RANGE<VECTOR<int,3> >(box.max_corner.x,box.max_corner.x,box.min_corner.y,box.max_corner.y,box.min_corner.z,box.max_corner.z)-band_x);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.min_corner.x,box.max_corner.y,box.max_corner.y,box.min_corner.z,box.max_corner.z)+band_x-band_y);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.max_corner.x,box.max_corner.y,box.max_corner.y,box.min_corner.z,box.max_corner.z)-band_y);
        regions.Append(RANGE<VECTOR<int,3> >(box.max_corner.x,box.max_corner.x,box.max_corner.y,box.max_corner.y,box.min_corner.z,box.max_corner.z)-band_x-band_y);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.min_corner.x,box.min_corner.y,box.min_corner.y,box.max_corner.z,box.max_corner.z)+band_x+band_y-band_z);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.max_corner.x,box.min_corner.y,box.min_corner.y,box.max_corner.z,box.max_corner.z)+band_y-band_z);
        regions.Append(RANGE<VECTOR<int,3> >(box.max_corner.x,box.max_corner.x,box.min_corner.y,box.min_corner.y,box.max_corner.z,box.max_corner.z)-band_x+band_y-band_z);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.min_corner.x,box.min_corner.y,box.max_corner.y,box.max_corner.z,box.max_corner.z)+band_x-band_z);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.max_corner.x,box.min_corner.y,box.max_corner.y,box.max_corner.z,box.max_corner.z)-band_z); 
        regions.Append(RANGE<VECTOR<int,3> >(box.max_corner.x,box.max_corner.x,box.min_corner.y,box.max_corner.y,box.max_corner.z,box.max_corner.z)-band_x-band_z);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.min_corner.x,box.max_corner.y,box.max_corner.y,box.max_corner.z,box.max_corner.z)+band_x-band_y-band_z);
        regions.Append(RANGE<VECTOR<int,3> >(box.min_corner.x,box.max_corner.x,box.max_corner.y,box.max_corner.y,box.max_corner.z,box.max_corner.z)-band_y-band_z);
        regions.Append(RANGE<VECTOR<int,3> >(box.max_corner.x,box.max_corner.x,box.max_corner.y,box.max_corner.y,box.max_corner.z,box.max_corner.z)-band_x-band_y-band_z);}
}
//#####################################################################
// Find_Boundary_Regions_Cell
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Find_Boundary_Regions_Cell(const int number_of_ghost_cells)
{
#ifdef USE_MPI
    Clean_Buffers_Boundary_Cell();
    MPI_UNIFORM_GRID<T_GRID>* mpi_grid=0;
    for(int boundary_grid_index=1;boundary_grid_index<=number_of_global_grids;boundary_grid_index++){
        if(mpi_split_grids.m) mpi_grid=mpi_split_grids(mpi_grid_indices(boundary_grid_index));
        ARRAY<CELL_PACKAGE_TYPE> boundary_packages(number_of_global_grids);
        RIGID_GRID<T_GRID>& boundary_grid=global_grid.Rigid_Grid(boundary_grid_index);
        ARRAY<RANGE<TV_INT> > boundary_regions;Find_Boundary_Regions(boundary_regions,RANGE<TV_INT>::Zero_Box(),false,RANGE<VECTOR<int,1> >(-number_of_ghost_cells,-1),true,true,boundary_grid.grid,boundary_grid_index);
        for(int side=1;side<=T_GRID::number_of_one_ring_neighbors_per_cell;side++) if(!mpi_grid || (mpi_grid && all_neighbor_ranks_chimera(boundary_grid_index)(side)==MPI::PROC_NULL)) for(CELL_ITERATOR iterator(boundary_grid.grid,boundary_regions(side));iterator.Valid();iterator.Next()){
            TV cell_location;
            if(use_previous_grid_frame_cell)
                cell_location=boundary_grid.previous_state.World_Space_Point(iterator.Location());
            else 
                cell_location=boundary_grid.World_Space_Point(iterator.Location());
            int background_grid_index=Find_Finest_Grid(cell_location,use_previous_grid_frame_cell);
            if(background_grid_index==0 || background_grid_index==boundary_grid_index) continue;
            boundary_packages(background_grid_index).x.Append(iterator.Cell_Index());
            boundary_packages(background_grid_index).y.Append((T)0);}

        for(int backgrid=1;backgrid<=number_of_global_grids;backgrid++) 
            if(boundary_packages(backgrid).x.Size() && (Local_Grid(backgrid) || Local_Grid(boundary_grid_index))){
                boundary_packages_cell_send.Append(PAIR<CELL_PACKAGE_TYPE,int>(boundary_packages(backgrid),Get_Send_Tag(backgrid,boundary_grid_index)));
                boundary_packages_cell_recv.Append(PAIR<CELL_PACKAGE_TYPE,int>(boundary_packages(backgrid),Get_Recv_Tag(backgrid,boundary_grid_index)));}}
    boundary_regions_found_cell.x=true;
    boundary_regions_found_cell.y=number_of_ghost_cells;
#endif
}
//#####################################################################
// Find_Boundary_Regions_Face
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Find_Boundary_Regions_Face(const int number_of_ghost_cells)
{
#ifdef USE_MPI
    Clean_Buffers_Boundary_Face();
    MPI_UNIFORM_GRID<T_GRID>* mpi_grid=0;
    for(int boundary_grid_index=1;boundary_grid_index<=number_of_global_grids;boundary_grid_index++){
        if(mpi_split_grids.m) mpi_grid=mpi_split_grids(mpi_grid_indices(boundary_grid_index));
        ARRAY<FACE_PACKAGE_TYPE> boundary_packages(number_of_global_grids);
        RIGID_GRID<T_GRID>& boundary_grid=global_grid.Rigid_Grid(boundary_grid_index);
        ARRAY<RANGE<TV_INT> > boundary_regions;
        for(int face_axis=1;face_axis<=T_GRID::dimension;face_axis++){
            T_GRID face_grid=boundary_grid.grid.Get_Face_Grid(face_axis);
            Find_Boundary_Regions(boundary_regions,RANGE<TV_INT>(TV_INT(),TV_INT::Axis_Vector(face_axis)),true,RANGE<VECTOR<int,1> >(-number_of_ghost_cells,-1),true,true,boundary_grid.grid,boundary_grid_index);
            for(int side=1;side<=T_GRID::number_of_one_ring_neighbors_per_cell;side++) if(!mpi_grid || (mpi_grid && all_neighbor_ranks_chimera(boundary_grid_index)(side)==MPI::PROC_NULL)) for(NODE_ITERATOR iterator(face_grid,boundary_regions(side));iterator.Valid();iterator.Next()){
                TV face_location;
                if(use_previous_grid_frame_face)
                    face_location=boundary_grid.previous_state.World_Space_Point(iterator.Location());
                else 
                    face_location=boundary_grid.World_Space_Point(iterator.Location());
                int background_grid_index=Find_Finest_Grid(face_location,use_previous_grid_frame_face);
                if(background_grid_index==0 || background_grid_index==boundary_grid_index) continue;
                boundary_packages(background_grid_index).x.Append(FACE_INDEX<TV::dimension>(face_axis,iterator.Node_Index()));
                boundary_packages(background_grid_index).y.Append((T)0);}}
     
        for(int backgrid=1;backgrid<=number_of_global_grids;backgrid++) 
            if(boundary_packages(backgrid).x.Size() && (Local_Grid(backgrid) || Local_Grid(boundary_grid_index))){
                boundary_packages_face_send.Append(PAIR<FACE_PACKAGE_TYPE,int>(boundary_packages(backgrid),Get_Send_Tag(backgrid,boundary_grid_index)));
                boundary_packages_face_recv.Append(PAIR<FACE_PACKAGE_TYPE,int>(boundary_packages(backgrid),Get_Recv_Tag(backgrid,boundary_grid_index)));}}
    boundary_regions_found_face.x=true;
    boundary_regions_found_face.y=number_of_ghost_cells;
#endif
}
//#####################################################################
// Find_Overlap_Regions_Cell
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Find_Overlap_Regions_Cell(const int& number_of_ghost_cells,ARRAY<T_ARRAYS_BOOL>* fixed_mask)
{
    if(!number_of_local_grids) return;
    Clean_Buffers_Overlap_Cell();
    ARRAY<CELL_PACKAGE_TYPE> boundary_packages;
    ARRAY<bool> recv_from(number_of_global_grids);recv_from.Fill(false);
    for(int boundary_grid_index=number_of_global_grids;boundary_grid_index>=1;boundary_grid_index--) if(Local_Grid(boundary_grid_index)){
        boundary_packages.Clean_Memory();boundary_packages.Resize(number_of_global_grids);
        RIGID_GRID<T_GRID>& boundary_grid=global_grid.Rigid_Grid(boundary_grid_index);
        for(CELL_ITERATOR iterator(boundary_grid.grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){
            if(fixed_mask && (*fixed_mask)(boundary_grid_index)(iterator.Cell_Index())) continue;
            TV cell_location;
            if(use_previous_grid_frame_cell)
                cell_location=boundary_grid.previous_state.World_Space_Point(iterator.Location());
            else 
                cell_location=boundary_grid.World_Space_Point(iterator.Location());
            int background_grid_index=Find_Finest_Grid(cell_location,use_previous_grid_frame_cell);
            if(background_grid_index==0 ||  background_grid_index==boundary_grid_index || (use_mpi && mpi_grid_indices(background_grid_index)==mpi_grid_indices(boundary_grid_index))) continue;
            RIGID_GRID<T_GRID>& background_grid=global_grid.Rigid_Grid(background_grid_index);
            if(boundary_grid.grid.Cell_Size()<background_grid.grid.Cell_Size() || (boundary_grid.grid.Cell_Size()==background_grid.grid.Cell_Size() && boundary_grid_index<background_grid_index) ) continue;
            boundary_packages(background_grid_index).x.Append(iterator.Cell_Index());
            boundary_packages(background_grid_index).y.Append((T)0);}
        for(int backgrid=1;backgrid<=number_of_global_grids;backgrid++) 
            if(boundary_packages(backgrid).x.Size() && (Local_Grid(backgrid) || Local_Grid(boundary_grid_index))){
                if(!use_mpi) overlap_packages_cell_send.Append(PAIR<CELL_PACKAGE_TYPE,int>(boundary_packages(backgrid),Get_Send_Tag(backgrid,boundary_grid_index,number_of_global_grids*number_of_global_grids)));
                recv_from(backgrid)=true;
                overlap_packages_cell_recv.Append(PAIR<CELL_PACKAGE_TYPE,int>(boundary_packages(backgrid),Get_Recv_Tag(backgrid,boundary_grid_index,number_of_global_grids*number_of_global_grids)));}}
    
    if(use_mpi){
        int local_grid_global_index=local_to_global_grid_index_map(1);
        Exchange_Overlap_Packages_Cell_Send(boundary_packages);
        for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++) if(boundary_packages(grid_index).x.Size() && !recv_from(grid_index))
            overlap_packages_cell_send.Append(PAIR<CELL_PACKAGE_TYPE,int>(boundary_packages(grid_index),Get_Send_Tag(local_grid_global_index,grid_index,number_of_global_grids*number_of_global_grids)));}

    overlap_regions_found_cell.x=true;
    overlap_regions_found_cell.y=number_of_ghost_cells;
}
//#####################################################################
// Find_Overlap_Regions_Face
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Find_Overlap_Regions_Face(const int& number_of_ghost_cells,ARRAY<T_FACE_ARRAYS_BOOL>* fixed_mask)
{
    if(!number_of_local_grids) return;
    Clean_Buffers_Overlap_Face();
    ARRAY<FACE_PACKAGE_TYPE> boundary_packages;
    ARRAY<bool> recv_from(number_of_global_grids);recv_from.Fill(false);
    for(int boundary_grid_index=number_of_global_grids;boundary_grid_index>=1;boundary_grid_index--) if(Local_Grid(boundary_grid_index)){
        boundary_packages.Clean_Memory();boundary_packages.Resize(number_of_global_grids);
        RIGID_GRID<T_GRID>& boundary_grid=global_grid.Rigid_Grid(boundary_grid_index);
        for(FACE_ITERATOR iterator(boundary_grid.grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){
            if(fixed_mask && (*fixed_mask)(boundary_grid_index)(iterator.Full_Index())) continue;
            TV face_location;
            if(use_previous_grid_frame_face)
                face_location=boundary_grid.previous_state.World_Space_Point(iterator.Location());
            else 
                face_location=boundary_grid.World_Space_Point(iterator.Location());
            int background_grid_index=Find_Finest_Grid(face_location,use_previous_grid_frame_face);
            // Note: The split grids belonging to the same mpi_grid should not try to interpolate shared faces from neighboring grids
            // Otherwise, due to slight numerical errors, the split grid with lower index can try to receive from the higher index grid and causing deadlock when recursively filling
            // coupling regions
            if(background_grid_index==0 ||  background_grid_index==boundary_grid_index || (use_mpi && mpi_grid_indices(background_grid_index)==mpi_grid_indices(boundary_grid_index))) continue;
            RIGID_GRID<T_GRID>& background_grid=global_grid.Rigid_Grid(background_grid_index);
            if(boundary_grid.grid.Cell_Size()<background_grid.grid.Cell_Size() || (boundary_grid.grid.Cell_Size()==background_grid.grid.Cell_Size() && boundary_grid_index<background_grid_index) ) continue;
            boundary_packages(background_grid_index).x.Append(iterator.Full_Index());
            boundary_packages(background_grid_index).y.Append((T)0);}
        for(int backgrid=1;backgrid<=number_of_global_grids;backgrid++) 
            if(boundary_packages(backgrid).x.Size() && (Local_Grid(backgrid) || Local_Grid(boundary_grid_index))){
                if(!use_mpi) overlap_packages_face_send.Append(PAIR<FACE_PACKAGE_TYPE,int>(boundary_packages(backgrid),Get_Send_Tag(backgrid,boundary_grid_index,number_of_global_grids*number_of_global_grids)));
                recv_from(backgrid)=true;
                overlap_packages_face_recv.Append(PAIR<FACE_PACKAGE_TYPE,int>(boundary_packages(backgrid),Get_Recv_Tag(backgrid,boundary_grid_index,number_of_global_grids*number_of_global_grids)));}}

    if(use_mpi){
        int local_grid_global_index=local_to_global_grid_index_map(1);
        Exchange_Overlap_Packages_Face_Send(boundary_packages);
        for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++) if(boundary_packages(grid_index).x.Size() && !recv_from(grid_index)){
            overlap_packages_face_send.Append(PAIR<FACE_PACKAGE_TYPE,int>(boundary_packages(grid_index),Get_Send_Tag(local_grid_global_index,grid_index,number_of_global_grids*number_of_global_grids)));}}
    
    overlap_regions_found_face.x=true;
    overlap_regions_found_face.y=number_of_ghost_cells;
}
//#####################################################################
// Exchange_Overlap_Cell_Data_Recursively
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Exchange_Overlap_Cell_Data_Recursively(ARRAY<T_ARRAYS_SCALAR*>& data,const int& number_of_ghost_cells,bool only_overwrite_nan,ARRAY<T_ARRAYS_BOOL>* fixed_mask)
{
    if(!Overlap_Regions_Found_Cell(number_of_ghost_cells)) Find_Overlap_Regions_Cell(number_of_ghost_cells,fixed_mask);
    for(int grid_index=1;grid_index<=number_of_local_grids;grid_index++){
        int current_grid=grids_sorted(grid_index);
        //recv information from other grids
        for(int pack=1;pack<=overlap_packages_cell_recv.Size();pack++){
            int send_grid=0,recv_grid=0,tag=overlap_packages_cell_recv(pack).y;Get_Send_Recv_Grid(tag,send_grid,recv_grid,number_of_global_grids*number_of_global_grids);
            if(local_to_global_grid_index_map(current_grid)==recv_grid){
                Recv_Data_From_Send_Grid_Cell(overlap_packages_cell_send,overlap_packages_cell_recv(pack).x,send_grid,recv_grid,tag);}}
#ifdef USE_MPI
        if(use_mpi){MPI_UTILITIES::Wait_All(*mpi_requests_buffers);MPI_PACKAGE::Free_All(*mpi_packages_buffers);}
#endif
        for(int pack=1;pack<=overlap_packages_cell_recv.Size();pack++){
            int send_grid=0,recv_grid=0,tag=overlap_packages_cell_recv(pack).y;Get_Send_Recv_Grid(tag,send_grid,recv_grid,number_of_global_grids*number_of_global_grids);
            if(local_to_global_grid_index_map(current_grid)==recv_grid){
                Unpack_Recv_Data_Cell(overlap_packages_cell_recv(pack).x,*data(current_grid),only_overwrite_nan);}}
        //send information to other grids
        for(int pack=1;pack<=overlap_packages_cell_send.Size();pack++){
            int send_grid=0,recv_grid=0,tag=overlap_packages_cell_send(pack).y;Get_Send_Recv_Grid(tag,send_grid,recv_grid,number_of_global_grids*number_of_global_grids);
            if(local_to_global_grid_index_map(current_grid)==send_grid){
                Pack_Send_Data_Cell(overlap_packages_cell_send(pack).x,*data(current_grid),send_grid,recv_grid);
                if(use_mpi) Send_Data_To_Recv_Grid_Cell(overlap_packages_cell_send(pack).x,send_grid,recv_grid,tag);}}}
#ifdef USE_MPI
    if(use_mpi){MPI_UTILITIES::Wait_All(*mpi_requests_buffers);MPI_PACKAGE::Free_All(*mpi_packages_buffers);}
#endif
}
//#####################################################################
// Exchange_Overlap_Face_Data_Recursively
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Exchange_Overlap_Face_Data_Recursively(ARRAY<T_FACE_ARRAYS_T2*>& data,const int& number_of_ghost_cells,ARRAY<T_FACE_ARRAYS_BOOL>* fixed_mask)
{
    if(!Overlap_Regions_Found_Face(number_of_ghost_cells)) Find_Overlap_Regions_Face(number_of_ghost_cells,fixed_mask);
    for(int grid_index=1;grid_index<=number_of_local_grids;grid_index++){
        int current_grid=grids_sorted(grid_index);
        //recv information from other grids
        for(int pack=1;pack<=overlap_packages_face_recv.Size();pack++){
            int send_grid=0,recv_grid=0,tag=overlap_packages_face_recv(pack).y;Get_Send_Recv_Grid(tag,send_grid,recv_grid,number_of_global_grids*number_of_global_grids);
            if(local_to_global_grid_index_map(current_grid)==recv_grid){
                Recv_Data_From_Send_Grid_Face(overlap_packages_face_send,overlap_packages_face_recv(pack).x,send_grid,recv_grid,tag);}}
#ifdef USE_MPI
        if(use_mpi){MPI_UTILITIES::Wait_All(*mpi_requests_buffers);MPI_PACKAGE::Free_All(*mpi_packages_buffers);}
#endif
        for(int pack=1;pack<=overlap_packages_face_recv.Size();pack++){
            int send_grid=0,recv_grid=0,tag=overlap_packages_face_recv(pack).y;Get_Send_Recv_Grid(tag,send_grid,recv_grid,number_of_global_grids*number_of_global_grids);
            if(local_to_global_grid_index_map(current_grid)==recv_grid){
                Unpack_Recv_Data_Face(overlap_packages_face_recv(pack).x,*data(current_grid));}}
        //send information to other grids
        for(int pack=1;pack<=overlap_packages_face_send.Size();pack++){
            int send_grid=0,recv_grid=0,tag=overlap_packages_face_send(pack).y;Get_Send_Recv_Grid(tag,send_grid,recv_grid,number_of_global_grids*number_of_global_grids);
            if(local_to_global_grid_index_map(current_grid)==send_grid){ 
                Pack_Send_Data_Face(overlap_packages_face_send(pack).x,*data(current_grid),send_grid,recv_grid);
                if(use_mpi) Send_Data_To_Recv_Grid_Face(overlap_packages_face_send(pack).x,send_grid,recv_grid,tag);}}}
#ifdef USE_MPI
    if(use_mpi){MPI_UTILITIES::Wait_All(*mpi_requests_buffers);MPI_PACKAGE::Free_All(*mpi_packages_buffers);}
#endif
}
//#####################################################################
// Exchange_Overlap_Cell_Data
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Exchange_Overlap_Cell_Data(ARRAY<T_ARRAYS_SCALAR*>& data,const int& number_of_ghost_cells,bool only_overwrite_nan,ARRAY<T_ARRAYS_BOOL>* fixed_mask)
{
    if(!Overlap_Regions_Found_Cell(number_of_ghost_cells)) Find_Overlap_Regions_Cell(number_of_ghost_cells,fixed_mask);
    //send information to other grids
    for(int grid_index=1;grid_index<=number_of_local_grids;grid_index++){
        for(int pack=1;pack<=overlap_packages_cell_send.Size();pack++){
            int send_grid=0,recv_grid=0,tag=overlap_packages_cell_send(pack).y;Get_Send_Recv_Grid(tag,send_grid,recv_grid,number_of_global_grids*number_of_global_grids);
            if(local_to_global_grid_index_map(grid_index)==send_grid){
                Pack_Send_Data_Cell(overlap_packages_cell_send(pack).x,*data(grid_index),send_grid,recv_grid);
                if(use_mpi) Send_Data_To_Recv_Grid_Cell(overlap_packages_cell_send(pack).x,send_grid,recv_grid,tag);}}}
    //receive information from other grids
    for(int grid_index=1;grid_index<=number_of_local_grids;grid_index++){
        for(int pack=1;pack<=overlap_packages_cell_recv.Size();pack++){
            int send_grid=0,recv_grid=0,tag=overlap_packages_cell_recv(pack).y;Get_Send_Recv_Grid(tag,send_grid,recv_grid,number_of_global_grids*number_of_global_grids);
            if(local_to_global_grid_index_map(grid_index)==recv_grid){
                Recv_Data_From_Send_Grid_Cell(overlap_packages_cell_send,overlap_packages_cell_recv(pack).x,send_grid,recv_grid,tag);}}}
#ifdef USE_MPI
    if(use_mpi){MPI_UTILITIES::Wait_All(*mpi_requests_buffers);MPI_PACKAGE::Free_All(*mpi_packages_buffers);}
#endif
    for(int grid_index=1;grid_index<=number_of_local_grids;grid_index++){
        for(int pack=1;pack<=overlap_packages_cell_recv.Size();pack++){
            int send_grid=0,recv_grid=0,tag=overlap_packages_cell_recv(pack).y;Get_Send_Recv_Grid(tag,send_grid,recv_grid,number_of_global_grids*number_of_global_grids);
            if(local_to_global_grid_index_map(grid_index)==recv_grid){
                Unpack_Recv_Data_Cell(overlap_packages_cell_recv(pack).x,*data(grid_index),only_overwrite_nan);}}}
}
//#####################################################################
// Exchange_Overlap_Face_Data
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Exchange_Overlap_Face_Data(ARRAY<T_FACE_ARRAYS_T2*>& data,const int& number_of_ghost_cells,ARRAY<T_FACE_ARRAYS_BOOL>* fixed_mask)
{
    if(!Overlap_Regions_Found_Face(number_of_ghost_cells)) Find_Overlap_Regions_Face(number_of_ghost_cells,fixed_mask);
    //send information to other grids
    for(int grid_index=1;grid_index<=number_of_local_grids;grid_index++){
        for(int pack=1;pack<=overlap_packages_face_send.Size();pack++){
            int send_grid=0,recv_grid=0,tag=overlap_packages_face_send(pack).y;Get_Send_Recv_Grid(tag,send_grid,recv_grid,number_of_global_grids*number_of_global_grids);
            if(local_to_global_grid_index_map(grid_index)==send_grid){ 
                Pack_Send_Data_Face(overlap_packages_face_send(pack).x,*data(grid_index),send_grid,recv_grid);
                if(use_mpi) Send_Data_To_Recv_Grid_Face(overlap_packages_face_send(pack).x,send_grid,recv_grid,tag);}}}
    //recv information from other grids
    for(int grid_index=1;grid_index<=number_of_local_grids;grid_index++){
        for(int pack=1;pack<=overlap_packages_face_recv.Size();pack++){
            int send_grid=0,recv_grid=0,tag=overlap_packages_face_recv(pack).y;Get_Send_Recv_Grid(tag,send_grid,recv_grid,number_of_global_grids*number_of_global_grids);
            if(local_to_global_grid_index_map(grid_index)==recv_grid){
                Recv_Data_From_Send_Grid_Face(overlap_packages_face_send,overlap_packages_face_recv(pack).x,send_grid,recv_grid,tag);}}}
#ifdef USE_MPI
    if(use_mpi){MPI_UTILITIES::Wait_All(*mpi_requests_buffers);MPI_PACKAGE::Free_All(*mpi_packages_buffers);}
#endif
    for(int grid_index=1;grid_index<=number_of_local_grids;grid_index++){
        for(int pack=1;pack<=overlap_packages_face_recv.Size();pack++){
            int send_grid=0,recv_grid=0,tag=overlap_packages_face_recv(pack).y;Get_Send_Recv_Grid(tag,send_grid,recv_grid,number_of_global_grids*number_of_global_grids);
            if(local_to_global_grid_index_map(grid_index)==recv_grid){
                Unpack_Recv_Data_Face(overlap_packages_face_recv(pack).x,*data(grid_index));}}}
}
//#####################################################################
// Exchange_Boundary_Cell_Data
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Exchange_Boundary_Cell_Data(ARRAY<T_ARRAYS_SCALAR*>& data,const int& number_of_ghost_cells,bool only_overwrite_nan)
{
    if(!Boundary_Regions_Found_Cell(number_of_ghost_cells)){Find_Boundary_Regions_Cell(number_of_ghost_cells);}
    //send information to other grids
    for(int grid_index=1;grid_index<=number_of_local_grids;grid_index++){
        for(int pack=1;pack<=boundary_packages_cell_send.Size();pack++){
            int send_grid=0,recv_grid=0,tag=boundary_packages_cell_send(pack).y;Get_Send_Recv_Grid(tag,send_grid,recv_grid);
            if(local_to_global_grid_index_map(grid_index)==send_grid){
                Pack_Send_Data_Cell(boundary_packages_cell_send(pack).x,*data(grid_index),send_grid,recv_grid);
                if(use_mpi) Send_Data_To_Recv_Grid_Cell(boundary_packages_cell_send(pack).x,send_grid,recv_grid,tag);}}}
    //recv information from other grids
    for(int grid_index=1;grid_index<=number_of_local_grids;grid_index++){
        for(int pack=1;pack<=boundary_packages_cell_recv.Size();pack++){
            int send_grid=0,recv_grid=0,tag=boundary_packages_cell_recv(pack).y;Get_Send_Recv_Grid(tag,send_grid,recv_grid);
            if(local_to_global_grid_index_map(grid_index)==recv_grid){
                Recv_Data_From_Send_Grid_Cell(boundary_packages_cell_send,boundary_packages_cell_recv(pack).x,send_grid,recv_grid,tag);}}}
#ifdef USE_MPI
    if(use_mpi){MPI_UTILITIES::Wait_All(*mpi_requests_buffers);MPI_PACKAGE::Free_All(*mpi_packages_buffers);}
#endif
    for(int grid_index=1;grid_index<=number_of_local_grids;grid_index++){
        for(int pack=1;pack<=boundary_packages_cell_recv.Size();pack++){
            int send_grid=0,recv_grid=0,tag=boundary_packages_cell_recv(pack).y;Get_Send_Recv_Grid(tag,send_grid,recv_grid);
            if(local_to_global_grid_index_map(grid_index)==recv_grid){
                Unpack_Recv_Data_Cell(boundary_packages_cell_recv(pack).x,*data(grid_index),only_overwrite_nan);}}}
}
//#####################################################################
// Exchange_Boundary_Face_Data
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Exchange_Boundary_Face_Data(ARRAY<T_FACE_ARRAYS_T2*>& data,const int& number_of_ghost_cells)
{
    if(!Boundary_Regions_Found_Face(number_of_ghost_cells)){Find_Boundary_Regions_Face(number_of_ghost_cells);}
    //send information to other grids
    for(int grid_index=1;grid_index<=number_of_local_grids;grid_index++){
        for(int pack=1;pack<=boundary_packages_face_send.Size();pack++){
            int send_grid=0,recv_grid=0,tag=boundary_packages_face_send(pack).y;Get_Send_Recv_Grid(tag,send_grid,recv_grid);
            if(local_to_global_grid_index_map(grid_index)==send_grid){ 
                Pack_Send_Data_Face(boundary_packages_face_send(pack).x,*data(grid_index),send_grid,recv_grid);
                if(use_mpi) Send_Data_To_Recv_Grid_Face(boundary_packages_face_send(pack).x,send_grid,recv_grid,tag);}}}
    //recv information from other grids
    for(int grid_index=1;grid_index<=number_of_local_grids;grid_index++){
        for(int pack=1;pack<=boundary_packages_face_recv.Size();pack++){
            int send_grid=0,recv_grid=0,tag=boundary_packages_face_recv(pack).y;Get_Send_Recv_Grid(tag,send_grid,recv_grid);
            if(local_to_global_grid_index_map(grid_index)==recv_grid){
                Recv_Data_From_Send_Grid_Face(boundary_packages_face_send,boundary_packages_face_recv(pack).x,send_grid,recv_grid,tag);}}}
#ifdef USE_MPI
    if(use_mpi){MPI_UTILITIES::Wait_All(*mpi_requests_buffers);MPI_PACKAGE::Free_All(*mpi_packages_buffers);}
#endif
    for(int grid_index=1;grid_index<=number_of_local_grids;grid_index++){
        for(int pack=1;pack<=boundary_packages_face_recv.Size();pack++){
            int send_grid=0,recv_grid=0,tag=boundary_packages_face_recv(pack).y;Get_Send_Recv_Grid(tag,send_grid,recv_grid);
            if(local_to_global_grid_index_map(grid_index)==recv_grid){
                Unpack_Recv_Data_Face(boundary_packages_face_recv(pack).x,*data(grid_index));}}}
}
//#####################################################################
// Pack_Send_Data_Cell
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Pack_Send_Data_Cell(CELL_PACKAGE_TYPE& cell_package,T_ARRAYS_T2& data,int send_grid,int recv_grid)
{
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2> interpolation;
    for(int cell=1;cell<=cell_package.x.Size();cell++){
        TV cell_location_in_recv_grid=global_grid.Rigid_Grid(recv_grid).grid.X(cell_package.x(cell));
        TV cell_location_in_current_grid;
        if(use_previous_grid_frame_cell)
            cell_location_in_current_grid=global_grid.Rigid_Grid(send_grid).previous_state.Object_Space_Point(global_grid.Rigid_Grid(recv_grid).previous_state.World_Space_Point(cell_location_in_recv_grid));
        else
            cell_location_in_current_grid=global_grid.Rigid_Grid(send_grid).Object_Space_Point(global_grid.Rigid_Grid(recv_grid).World_Space_Point(cell_location_in_recv_grid));
        T value=interpolation.Clamped_To_Array(global_grid.Rigid_Grid(send_grid).grid,data,cell_location_in_current_grid);
        cell_package.y(cell)=value;}
}
//#####################################################################
// Pack_Send_Data_Face
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Pack_Send_Data_Face(FACE_PACKAGE_TYPE& face_package,T_FACE_ARRAYS_T2& data,int send_grid,int recv_grid)
{
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2> interpolation;
    for(int face=1;face<=face_package.x.Size();face++){
        TV face_location_in_recv_grid=global_grid.Rigid_Grid(recv_grid).grid.Face(face_package.x(face).axis,face_package.x(face).index);
        TV face_location_in_current_grid;
        if(use_previous_grid_frame_face)
            face_location_in_current_grid=global_grid.Rigid_Grid(send_grid).previous_state.Object_Space_Point(global_grid.Rigid_Grid(recv_grid).previous_state.World_Space_Point(face_location_in_recv_grid));
        else
            face_location_in_current_grid=global_grid.Rigid_Grid(send_grid).Object_Space_Point(global_grid.Rigid_Grid(recv_grid).World_Space_Point(face_location_in_recv_grid));
        TV interpolated_vector=interpolation.Clamped_To_Array_Face(global_grid.Rigid_Grid(send_grid).grid,data,face_location_in_current_grid);
        if(use_previous_grid_frame_face)
            interpolated_vector=global_grid.Rigid_Grid(recv_grid).previous_state.Object_Space_Vector(global_grid.Rigid_Grid(send_grid).previous_state.World_Space_Vector(interpolated_vector));
        else
            interpolated_vector=global_grid.Rigid_Grid(recv_grid).Object_Space_Vector(global_grid.Rigid_Grid(send_grid).World_Space_Vector(interpolated_vector));
        face_package.y(face)=interpolated_vector(face_package.x(face).axis);}
}
//#####################################################################
// Send_Data_To_Recv_Grid_Cell
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Send_Data_To_Recv_Grid_Cell(CELL_PACKAGE_TYPE& cell_package,int send_grid,int recv_grid,int tag)
{
#ifdef USE_MPI
    if(use_mpi && !global_to_local_grid_index_map(recv_grid)){
        MPI_PACKAGE package=MPI_PACKAGE(cell_package.y);
        mpi_packages_buffers->Append(package);
        mpi_requests_buffers->Append(package.Isend(*comm,grid_ranks(recv_grid),tag));}
#endif
}
//#####################################################################
// Send_Data_To_Recv_Grid_Face
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Send_Data_To_Recv_Grid_Face(FACE_PACKAGE_TYPE& face_package,int send_grid,int recv_grid,int tag)
{
#ifdef USE_MPI
    if(use_mpi && !global_to_local_grid_index_map(recv_grid)){
        MPI_PACKAGE package=MPI_PACKAGE(face_package.y);
        mpi_packages_buffers->Append(package);
        mpi_requests_buffers->Append(package.Isend(*comm,grid_ranks(recv_grid),tag));}
#endif
}
//#####################################################################
// Recv_Data_From_Send_Grid_Cell
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Recv_Data_From_Send_Grid_Cell(ARRAY<PAIR<CELL_PACKAGE_TYPE,int> >& packages_cell_send,CELL_PACKAGE_TYPE& cell_package,int send_grid,int recv_grid,int tag)
{
    if(!use_mpi || global_to_local_grid_index_map(send_grid)) for(int pack=1;pack<=packages_cell_send.Size();pack++){
        if(tag==packages_cell_send(pack).y){
            for(int cell=1;cell<=cell_package.x.Size();cell++){
                cell_package.y(cell)=(packages_cell_send(pack).x).y(cell);}}}
#ifdef USE_MPI
    else{
        MPI_PACKAGE package=MPI_PACKAGE(cell_package.y);
        mpi_packages_buffers->Append(package);
        mpi_requests_buffers->Append(package.Irecv(*comm,grid_ranks(send_grid),tag));}
#endif
}
//#####################################################################
// Recv_Data_From_Send_Grid
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Recv_Data_From_Send_Grid_Face(ARRAY<PAIR<FACE_PACKAGE_TYPE,int> >& packages_face_send,FACE_PACKAGE_TYPE& face_package,int send_grid,int recv_grid,int tag)
{
    if(!use_mpi || global_to_local_grid_index_map(send_grid)) for(int pack=1;pack<=packages_face_send.Size();pack++){
        if(tag==packages_face_send(pack).y){
            for(int face=1;face<=face_package.x.Size();face++){
                face_package.y(face)=(packages_face_send(pack).x).y(face);}}}
#ifdef USE_MPI
    else{
        MPI_PACKAGE package=MPI_PACKAGE(face_package.y);
        mpi_packages_buffers->Append(package);
        mpi_requests_buffers->Append(package.Irecv(*comm,grid_ranks(send_grid),tag));}
#endif
}
//#####################################################################
// Unpack_Recv_Data_Cell
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Unpack_Recv_Data_Cell(CELL_PACKAGE_TYPE& cell_package,T_ARRAYS_T2& data,bool only_overwrite_nan)
{
    for(int cell=1;cell<=cell_package.x.Size();cell++){
        bool is_nan=Is_NaN(data(cell_package.x(cell)));
        if(!only_overwrite_nan || is_nan)
            data(cell_package.x(cell))=cell_package.y(cell);}
}
//#####################################################################
// Unpack_Recv_Data_Face
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Unpack_Recv_Data_Face(FACE_PACKAGE_TYPE& face_package,T_FACE_ARRAYS_T2& data)
{
    for(int face=1;face<=face_package.x.Size();face++){
        data(face_package.x(face))=face_package.y(face);}
}
//#####################################################################
// Function Initialize_Communicator
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Initialize_Communicator()
{
#ifdef USE_MPI
    comm=new MPI::Intracomm;*comm=MPI::COMM_WORLD.Dup();
    LOG::filecout("After create comm\n");
    rank=comm->Get_rank();
#endif
}
//#####################################################################
// Function Synchronize_Dt
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Synchronize_Dt(T& dt) const
{
#ifdef USE_MPI
    T dt_local=dt;
    MPI_UTILITIES::Reduce(dt_local,dt,MPI::MIN,*comm);
#endif
}
//#####################################################################
// Function Synchronize_Done
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Synchronize_Done(bool& done) const
{
#ifdef USE_MPI
    bool done_local=done;
    comm->Allreduce(&done_local,&done,1,MPI_UTILITIES::Datatype<bool>(),MPI::LAND);
#endif
}
//#####################################################################
// Function Synchronize_Fix_Masks
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Synchronize_Fix_Masks(ARRAY<T_ARRAYS_BOOL>& fixed_mask)
{
    if(!use_mpi) return;
    #ifdef USE_MPI
    //send to other grids
    for(int local_grid_index=1;local_grid_index<=number_of_local_grids;local_grid_index++){
        int this_grid_index=local_to_global_grid_index_map(local_grid_index);
        for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++){
            if(global_to_local_grid_index_map(grid_index)) continue;
            MPI_PACKAGE package(fixed_mask(this_grid_index),fixed_mask(this_grid_index).domain);
            mpi_packages_buffers->Append(package);
            mpi_requests_buffers->Append(package.Isend(*comm,grid_ranks(grid_index),Get_Send_Tag(this_grid_index,grid_index)));}}
    //receive from other grids
    for(int local_grid_index=1;local_grid_index<=number_of_local_grids;local_grid_index++){
        int this_grid_index=local_to_global_grid_index_map(local_grid_index);
        for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++){
            if(global_to_local_grid_index_map(grid_index)) continue;
            MPI_PACKAGE package(fixed_mask(grid_index),fixed_mask(grid_index).domain);
            mpi_packages_buffers->Append(package);
            mpi_requests_buffers->Append(package.Irecv(*comm,grid_ranks(grid_index),Get_Recv_Tag(grid_index,this_grid_index)));}}
    MPI_UTILITIES::Wait_All(*mpi_requests_buffers);MPI_PACKAGE::Free_All(*mpi_packages_buffers);
    #endif
}
//#####################################################################
// Function Exchange_Boundary_Cell_Arrays
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Exchange_Boundary_Cell_Arrays(ARRAY<ARRAY<TV_INT> >& boundary_cell_indices,ARRAY<ARRAY<bool> >& is_boundary_grid)
{
#ifdef USE_MPI
    if(!use_mpi) return;
    //Exchange sizes of boundary_arrays
    ARRAY<int> boundary_cell_indices_sizes(number_of_global_grids);boundary_cell_indices_sizes.Fill(0);
    for(int local_grid_index=1;local_grid_index<=number_of_local_grids;local_grid_index++){
        int this_grid_index=local_to_global_grid_index_map(local_grid_index);
        boundary_cell_indices_sizes(this_grid_index)=boundary_cell_indices(this_grid_index).Size();
        for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++) if(is_boundary_grid(local_grid_index)(grid_index))
            mpi_requests_buffers->Append(comm->Isend(&boundary_cell_indices_sizes(this_grid_index),1,MPI_UTILITIES::Datatype<int>(),grid_ranks(grid_index),Get_Send_Tag(this_grid_index,grid_index)));
        for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++) if(is_boundary_grid(local_grid_index)(grid_index))
            mpi_requests_buffers->Append(comm->Irecv(&boundary_cell_indices_sizes(grid_index),1,MPI_UTILITIES::Datatype<int>(),grid_ranks(grid_index),Get_Recv_Tag(grid_index,this_grid_index)));}
    MPI_UTILITIES::Wait_All(*mpi_requests_buffers);
    //Resize boundary_cell_indices arrays for receiving
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++){
        boundary_cell_indices(grid_index).Resize(boundary_cell_indices_sizes(grid_index));}
    //Exchange boundary_cell_indices
    for(int local_grid_index=1;local_grid_index<=number_of_local_grids;local_grid_index++){
        int this_grid_index=local_to_global_grid_index_map(local_grid_index);
        for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++) if(is_boundary_grid(local_grid_index)(grid_index)){
            MPI_PACKAGE package(boundary_cell_indices(this_grid_index));
            mpi_packages_buffers->Append(package);
            mpi_requests_buffers->Append(package.Isend(*comm,grid_ranks(grid_index),Get_Send_Tag(this_grid_index,grid_index)));}
        for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++) if(is_boundary_grid(local_grid_index)(grid_index)){
            MPI_PACKAGE package(boundary_cell_indices(grid_index));
            mpi_packages_buffers->Append(package);
            mpi_requests_buffers->Append(package.Irecv(*comm,grid_ranks(grid_index),Get_Recv_Tag(grid_index,this_grid_index)));}}
    MPI_UTILITIES::Wait_All(*mpi_requests_buffers);MPI_PACKAGE::Free_All(*mpi_packages_buffers);
#endif
}
//#####################################################################
// Function Exchange_Boundary_Voronoi_Faces
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Exchange_Boundary_Voronoi_Faces(ARRAY<ARRAY<VORONOI_FACE> >& boundary_cell_voronoi_faces_to_send,ARRAY<ARRAY<VORONOI_FACE> >& boundary_cell_voronoi_faces_to_receive)
{
#ifdef USE_MPI
    if(!use_mpi) return;
    //Exchange sizes of boundary_arrays
    ARRAY<int> boundary_sizes_to_send(number_of_global_grids),boundary_sizes_to_receive(number_of_global_grids);
    boundary_sizes_to_receive.Fill(0);
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++) boundary_sizes_to_send(grid_index)=boundary_cell_voronoi_faces_to_send(grid_index).Size();
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++){
        if(!Local_Grid(grid_index)){
            mpi_requests_buffers->Append(comm->Isend(&boundary_sizes_to_send(grid_index),1,MPI_UTILITIES::Datatype<int>(),grid_ranks(grid_index),0));
            mpi_requests_buffers->Append(comm->Irecv(&boundary_sizes_to_receive(grid_index),1,MPI_UTILITIES::Datatype<int>(),grid_ranks(grid_index),0));}}
    MPI_UTILITIES::Wait_All(*mpi_requests_buffers);
    //Resize boundary_cell_voronoi_faces arrays for receiving
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++)
        boundary_cell_voronoi_faces_to_receive(grid_index).Resize(boundary_sizes_to_receive(grid_index));
    //Exchange boundary_cell_voronoi_faces
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++){
        if(boundary_sizes_to_send(grid_index)){
            MPI_PACKAGE package(boundary_cell_voronoi_faces_to_send(grid_index));
            mpi_packages_buffers->Append(package);
            mpi_requests_buffers->Append(package.Isend(*comm,grid_ranks(grid_index),0));}
        if(boundary_sizes_to_receive(grid_index)){
            MPI_PACKAGE package(boundary_cell_voronoi_faces_to_receive(grid_index));
            mpi_packages_buffers->Append(package);
            mpi_requests_buffers->Append(package.Irecv(*comm,grid_ranks(grid_index),0));}}
    MPI_UTILITIES::Wait_All(*mpi_requests_buffers);MPI_PACKAGE::Free_All(*mpi_packages_buffers);
#endif
}
//#####################################################################
// Function Exchange_Voronoi_Face_Positions_For_Interpolation
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Exchange_Voronoi_Face_Positions_For_Interpolation(ARRAY<ARRAY<PAIR<TV,TV> > >& voronoi_face_positions_interpolated_from_others,ARRAY<ARRAY<PAIR<TV,TV> > >& voronoi_face_positions_interpolated_here)
{
#ifdef USE_MPI
    if(!use_mpi) return;
    //Exchange sizes of arrays
    ARRAY<int> sizes_to_send(number_of_global_grids),sizes_to_receive(number_of_global_grids);sizes_to_receive.Fill(0);
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++) sizes_to_send(grid_index)=voronoi_face_positions_interpolated_from_others(grid_index).Size();
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++){
        if(!Local_Grid(grid_index)){
            mpi_requests_buffers->Append(comm->Isend(&sizes_to_send(grid_index),1,MPI_UTILITIES::Datatype<int>(),grid_ranks(grid_index),0));
            mpi_requests_buffers->Append(comm->Irecv(&sizes_to_receive(grid_index),1,MPI_UTILITIES::Datatype<int>(),grid_ranks(grid_index),0));}}
    MPI_UTILITIES::Wait_All(*mpi_requests_buffers);
    //Resize arrays for receiving
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++)
        voronoi_face_positions_interpolated_here(grid_index).Resize(sizes_to_receive(grid_index));
    //Exchange 
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++){
        if(sizes_to_send(grid_index)){
            MPI_PACKAGE package(voronoi_face_positions_interpolated_from_others(grid_index));
            mpi_packages_buffers->Append(package);
            mpi_requests_buffers->Append(package.Isend(*comm,grid_ranks(grid_index),0));}
        if(sizes_to_receive(grid_index)){
            MPI_PACKAGE package(voronoi_face_positions_interpolated_here(grid_index));
            mpi_packages_buffers->Append(package);
            mpi_requests_buffers->Append(package.Irecv(*comm,grid_ranks(grid_index),0));}}
    MPI_UTILITIES::Wait_All(*mpi_requests_buffers);MPI_PACKAGE::Free_All(*mpi_packages_buffers);
#endif
}
//#####################################################################
// Function Exchange_Voronoi_Face_Interpolated_Values
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Exchange_Voronoi_Face_Interpolated_Values(ARRAY<ARRAY<T> >& voronoi_face_values_interpolated_here,ARRAY<ARRAY<T> >& voronoi_face_values_interpolated_from_others)
{
#ifdef USE_MPI
    if(!use_mpi) return;
    //Exchange sizes of arrays
    ARRAY<int> sizes_to_send(number_of_global_grids),sizes_to_receive(number_of_global_grids);sizes_to_receive.Fill(0);
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++) sizes_to_send(grid_index)=voronoi_face_values_interpolated_here(grid_index).Size();
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++){
        if(!Local_Grid(grid_index)){
            mpi_requests_buffers->Append(comm->Isend(&sizes_to_send(grid_index),1,MPI_UTILITIES::Datatype<int>(),grid_ranks(grid_index),0));
            mpi_requests_buffers->Append(comm->Irecv(&sizes_to_receive(grid_index),1,MPI_UTILITIES::Datatype<int>(),grid_ranks(grid_index),0));}}
    MPI_UTILITIES::Wait_All(*mpi_requests_buffers);
    //Resize arrays for receiving
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++)
        voronoi_face_values_interpolated_from_others(grid_index).Resize(sizes_to_receive(grid_index));
    //Exchange 
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++){
        if(sizes_to_send(grid_index)){
            MPI_PACKAGE package(voronoi_face_values_interpolated_here(grid_index));
            mpi_packages_buffers->Append(package);
            mpi_requests_buffers->Append(package.Isend(*comm,grid_ranks(grid_index),0));}
        if(sizes_to_receive(grid_index)){
            MPI_PACKAGE package(voronoi_face_values_interpolated_from_others(grid_index));
            mpi_packages_buffers->Append(package);
            mpi_requests_buffers->Append(package.Irecv(*comm,grid_ranks(grid_index),0));}}
    MPI_UTILITIES::Wait_All(*mpi_requests_buffers);MPI_PACKAGE::Free_All(*mpi_packages_buffers);
#endif
}
//#####################################################################
// Function Exchange_Boundary_Scalar_Values
//#####################################################################
template<class T_GRID,class T2> template<class T3> void CHIMERA_GRID<T_GRID,T2>::
Exchange_Boundary_Scalar_Values(ARRAY<ARRAY<T3> >& boundary_scalar_values,const ARRAY<ARRAY<bool> >& is_boundary_grid)
{
#ifdef USE_MPI
    if(!use_mpi) return;
    for(int local_grid_index=1;local_grid_index<=number_of_local_grids;local_grid_index++){
        int this_grid_index=local_to_global_grid_index_map(local_grid_index);
        for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++) if(is_boundary_grid(local_grid_index)(grid_index)){
            MPI_PACKAGE package(boundary_scalar_values(this_grid_index));
            mpi_packages_buffers->Append(package);
            mpi_requests_buffers->Append(package.Isend(*comm,grid_ranks(grid_index),Get_Send_Tag(this_grid_index,grid_index)));}}
    for(int local_grid_index=1;local_grid_index<=number_of_local_grids;local_grid_index++){
        int this_grid_index=local_to_global_grid_index_map(local_grid_index);
        for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++) if(is_boundary_grid(local_grid_index)(grid_index)){
            MPI_PACKAGE package(boundary_scalar_values(grid_index));
            mpi_packages_buffers->Append(package);
            mpi_requests_buffers->Append(package.Irecv(*comm,grid_ranks(grid_index),Get_Recv_Tag(grid_index,this_grid_index)));}}
    MPI_UTILITIES::Wait_All(*mpi_requests_buffers);MPI_PACKAGE::Free_All(*mpi_packages_buffers);
#endif
}
//#####################################################################
// Function Exchange_Overlap_Packages_Cell_Send
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Exchange_Overlap_Packages_Cell_Send(ARRAY<CELL_PACKAGE_TYPE>& boundary_packages)
{
#ifdef USE_MPI
    if(!use_mpi) return;
    int local_grid_global_index=local_to_global_grid_index_map(1);
    //Exchange sizes of boundary_arrays
    ARRAY<int> cell_indices_sizes_send(number_of_global_grids);
    ARRAY<int> cell_indices_sizes_recv(number_of_global_grids);cell_indices_sizes_recv.Fill(0);
    for(int recv_grid_index=1;recv_grid_index<=number_of_global_grids;recv_grid_index++) 
        cell_indices_sizes_send(recv_grid_index)=boundary_packages(recv_grid_index).x.Size();
    for(int recv_grid_index=1;recv_grid_index<=number_of_global_grids;recv_grid_index++) if(Grid_Intersect_Local_Domain(recv_grid_index) && !Local_Grid(recv_grid_index))
        mpi_requests_buffers->Append(comm->Isend(&cell_indices_sizes_send(recv_grid_index),1,MPI_UTILITIES::Datatype<int>(),grid_ranks(recv_grid_index),Get_Send_Tag(local_grid_global_index,recv_grid_index)));
    for(int send_grid_index=1;send_grid_index<=number_of_global_grids;send_grid_index++) if(Grid_Intersect_Local_Domain(send_grid_index) && !Local_Grid(send_grid_index))
        mpi_requests_buffers->Append(comm->Irecv(&cell_indices_sizes_recv(send_grid_index),1,MPI_UTILITIES::Datatype<int>(),grid_ranks(send_grid_index),Get_Recv_Tag(send_grid_index,local_grid_global_index)));
    MPI_UTILITIES::Wait_All(*mpi_requests_buffers);
    //Resize boundary_cell_indices arrays for receiving
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++) if(!Local_Grid(grid_index) && cell_indices_sizes_send(grid_index)==0){
        boundary_packages(grid_index).x.Resize(cell_indices_sizes_recv(grid_index));
        boundary_packages(grid_index).y.Resize(cell_indices_sizes_recv(grid_index));}
    //Exchange boundary_cell_indices
    for(int recv_grid_index=1;recv_grid_index<=number_of_global_grids;recv_grid_index++) if(cell_indices_sizes_send(recv_grid_index)){
        MPI_PACKAGE package(boundary_packages(recv_grid_index).x);
        mpi_packages_buffers->Append(package);
        mpi_requests_buffers->Append(package.Isend(*comm,grid_ranks(recv_grid_index),Get_Send_Tag(local_grid_global_index,recv_grid_index)));}
    for(int send_grid_index=1;send_grid_index<=number_of_global_grids;send_grid_index++) if(cell_indices_sizes_recv(send_grid_index)){
        MPI_PACKAGE package(boundary_packages(send_grid_index).x);
        mpi_packages_buffers->Append(package);
        mpi_requests_buffers->Append(package.Irecv(*comm,grid_ranks(send_grid_index),Get_Recv_Tag(send_grid_index,local_grid_global_index)));}
    MPI_UTILITIES::Wait_All(*mpi_requests_buffers);MPI_PACKAGE::Free_All(*mpi_packages_buffers);
#endif
}
//#####################################################################
// Function Exchange_Overlap_Packages_Face_Send
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Exchange_Overlap_Packages_Face_Send(ARRAY<FACE_PACKAGE_TYPE>& boundary_packages)
{
#ifdef USE_MPI
    if(!use_mpi) return;
    int local_grid_global_index=local_to_global_grid_index_map(1);
    //Exchange sizes of boundary_arrays
    ARRAY<int> face_indices_sizes_send(number_of_global_grids);
    ARRAY<int> face_indices_sizes_recv(number_of_global_grids);face_indices_sizes_recv.Fill(0);
    for(int recv_grid_index=1;recv_grid_index<=number_of_global_grids;recv_grid_index++) 
        face_indices_sizes_send(recv_grid_index)=boundary_packages(recv_grid_index).x.Size();
    for(int recv_grid_index=1;recv_grid_index<=number_of_global_grids;recv_grid_index++) if(Grid_Intersect_Local_Domain(recv_grid_index) && !Local_Grid(recv_grid_index))
        mpi_requests_buffers->Append(comm->Isend(&face_indices_sizes_send(recv_grid_index),1,MPI_UTILITIES::Datatype<int>(),grid_ranks(recv_grid_index),Get_Send_Tag(local_grid_global_index,recv_grid_index)));
    for(int send_grid_index=1;send_grid_index<=number_of_global_grids;send_grid_index++) if(Grid_Intersect_Local_Domain(send_grid_index) && !Local_Grid(send_grid_index))
        mpi_requests_buffers->Append(comm->Irecv(&face_indices_sizes_recv(send_grid_index),1,MPI_UTILITIES::Datatype<int>(),grid_ranks(send_grid_index),Get_Recv_Tag(send_grid_index,local_grid_global_index)));
    MPI_UTILITIES::Wait_All(*mpi_requests_buffers);
    //Resize boundary_face_indices arrays for receiving
    for(int grid_index=1;grid_index<=number_of_global_grids;grid_index++) if(!Local_Grid(grid_index) && face_indices_sizes_send(grid_index)==0){
        boundary_packages(grid_index).x.Resize(face_indices_sizes_recv(grid_index));
        boundary_packages(grid_index).y.Resize(face_indices_sizes_recv(grid_index));}
    //Exchange boundary_face_indices
    for(int recv_grid_index=1;recv_grid_index<=number_of_global_grids;recv_grid_index++) if(face_indices_sizes_send(recv_grid_index)){
        MPI_PACKAGE package(boundary_packages(recv_grid_index).x);
        mpi_packages_buffers->Append(package);
        mpi_requests_buffers->Append(package.Isend(*comm,grid_ranks(recv_grid_index),Get_Send_Tag(local_grid_global_index,recv_grid_index)));}
    for(int send_grid_index=1;send_grid_index<=number_of_global_grids;send_grid_index++) if(face_indices_sizes_recv(send_grid_index)){
        MPI_PACKAGE package(boundary_packages(send_grid_index).x);
        mpi_packages_buffers->Append(package);
        mpi_requests_buffers->Append(package.Irecv(*comm,grid_ranks(send_grid_index),Get_Recv_Tag(send_grid_index,local_grid_global_index)));}
    MPI_UTILITIES::Wait_All(*mpi_requests_buffers);MPI_PACKAGE::Free_All(*mpi_packages_buffers);
#endif
}
//#####################################################################
// Function Solve_Parallel
//#####################################################################
template<class T_GRID,class T2> void CHIMERA_GRID<T_GRID,T2>::
Parallel_Solve(SPARSE_MATRIX_FLAT_NXN<T>& A_pcg,VECTOR_ND<T>& x_pcg,VECTOR_ND<T>& rhs_pcg,T tolerance,SPARSE_MATRIX_PARTITION& partition,ARRAY<ARRAY<bool> >& is_boundary_grid,int maximum_iterations,bool use_incomplete_cholesky)
{
#ifdef USE_MPI
    PCG_SPARSE<T> pcg_sparse;pcg_sparse.Use_Conjugate_Gradient();
    PCG_SPARSE_UNSTRUCTURE_MPI<T_GRID> pcg_sparse_mpi(pcg_sparse,*comm,partition,is_boundary_grid,local_to_global_grid_index_map);
    pcg_sparse.show_results=true;
    pcg_sparse.maximum_iterations=maximum_iterations;
    if(use_incomplete_cholesky) pcg_sparse.Use_Incomplete_Cholesky();
    pcg_sparse_mpi.Parallel_Solve(A_pcg,x_pcg,rhs_pcg,tolerance);
#endif
}
//#####################################################################
// Function Unsplit_Cell_Index
//#####################################################################
template<class T_GRID,class T2> typename T_GRID::VECTOR_INT CHIMERA_GRID<T_GRID,T2>::Unsplit_Cell_Index(const GRID_CELL_INDEX& grid_cell_index)
{
    if(!use_mpi) return grid_cell_index.y;
    const TV_INT& coordinate=all_coordinates(grid_cell_index.x);
    TV_INT unsplit_cell_index;
    int mpi_grid_index=mpi_grid_indices(grid_cell_index.x);
    for(int axis=1;axis<=TV::dimension;axis++)
        unsplit_cell_index(axis)=grid_cell_index.y(axis)+mpi_split_grids(mpi_grid_index)->boundaries(axis)(coordinate(axis))-1;
    return unsplit_cell_index;
}
//#####################################################################
// Function Split_Cell_Index
//#####################################################################
template<class T_GRID,class T2> typename T_GRID::VECTOR_INT CHIMERA_GRID<T_GRID,T2>::Split_Cell_Index(const int grid_index,const TV_INT& unsplit_cell_index)
{
    if(!use_mpi) return unsplit_cell_index;
    const TV_INT& coordinate=all_coordinates(grid_index);
    TV_INT split_cell_index;
    int mpi_grid_index=mpi_grid_indices(grid_index);
    for(int axis=1;axis<=TV::dimension;axis++)
        split_cell_index(axis)=unsplit_cell_index(axis)-mpi_split_grids(mpi_grid_index)->boundaries(axis)(coordinate(axis))+1;
    return split_cell_index;
}
//#####################################################################
// Function Split_Cell_Index
//#####################################################################
template<class T_GRID,class T2> int CHIMERA_GRID<T_GRID,T2>::Parent_Grid(const int grid_index)
{
    if(!use_mpi) return grid_index;
    return mpi_grid_indices(grid_index);
}
//#####################################################################
// Find_Real_Grid_Cell_Index
//#####################################################################
template<class T_GRID,class T2> PAIR<int,typename T_GRID::VECTOR_INT> CHIMERA_GRID<T_GRID,T2>::
Find_Real_Grid_Cell_Index(const GRID_CELL_INDEX& index)
{
    if(!use_mpi) return index;
    int mpi_grid_index=mpi_grid_indices(index.x);
    if(mpi_split_grids(mpi_grid_index)->number_of_processes<=1) return index;

    TV_INT unsplit_cell_index=Unsplit_Cell_Index(index);
    TV_INT real_coordinate;
    TV_INT split_cell_index;
    for(int axis=1;axis<=TV::dimension;axis++){
        real_coordinate(axis)=1;
        int count=mpi_split_grids(mpi_grid_index)->boundaries(axis).Size()-1;
        while(real_coordinate(axis)<count && unsplit_cell_index(axis)>=mpi_split_grids(mpi_grid_index)->boundaries(axis)(real_coordinate(axis)+1))
            real_coordinate(axis)++;
        split_cell_index(axis)=unsplit_cell_index(axis)-mpi_split_grids(mpi_grid_index)->boundaries(axis)(real_coordinate(axis))+1;}

    return GRID_CELL_INDEX(mpi_grid_coordinate_to_global_grid_index_map(mpi_grid_index)(real_coordinate),split_cell_index);
}
//#####################################################################
// Function Solve_Parallel
//#####################################################################
template<class T_GRID,class T2> int CHIMERA_GRID<T_GRID,T2>::
Number_Of_Processors()
{
#ifdef USE_MPI
    return MPI::COMM_WORLD.Get_size();
#else
    return (int) 1;
#endif
}
//#####################################################################
#define INSTANTIATION_HELPER(T,D)                      \
    template class CHIMERA_GRID<GRID<VECTOR<T,D> >,T>; \
    template void CHIMERA_GRID<GRID<VECTOR<T,D> >,T>::Exchange_Boundary_Scalar_Values(ARRAY<ARRAY<T> >& scalar_values,const ARRAY<ARRAY<bool> >& is_boundary_grid); \
    template void CHIMERA_GRID<GRID<VECTOR<T,D> >,T>::Exchange_Boundary_Scalar_Values(ARRAY<ARRAY<bool> >& scalar_values,const ARRAY<ARRAY<bool> >& is_boundary_grid);

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif
