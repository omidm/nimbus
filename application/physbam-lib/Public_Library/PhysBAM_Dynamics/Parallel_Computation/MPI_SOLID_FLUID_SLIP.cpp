// Copyright 2008-2009, Avi Robinson-Mosher, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_SOLID_FLUID_SLIP
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID_SLIP.h>
#include <PhysBAM_Dynamics/Parallel_Computation/SYMMQMR_SPARSE_MPI.h>
#ifdef USE_MPI
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#endif
namespace PhysBAM{

#ifdef USE_MPI

//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPI_SOLID_FLUID_SLIP<TV>::
MPI_SOLID_FLUID_SLIP()
    :number_of_solid_processes(1),solid_node(0),fluid_comm(0),mpi_grid(0),poisson(0),current_tag(0)
{
    LOG::SCOPE scope("MPI INITIALIZE","Initializing MPI");
    number_of_processes=MPI::COMM_WORLD.Get_size();
    std::stringstream ss;
    ss<<"number of processes = "<<number_of_processes<<std::endl;

    comm=new MPI::Intracomm;
    *comm=MPI::COMM_WORLD.Dup();
    rank=comm->Get_rank();
    ss<<"global rank = "<<rank<<std::endl;
    LOG::filecout(ss.str());

    group=new MPI::Group(comm->Get_group());
    solid_ranks.Resize(number_of_solid_processes);
    for(int i=1;i<=number_of_solid_processes;i++)solid_ranks(i)=i-1;
    solid_group=new MPI::Group(group->Incl(solid_ranks.n,&solid_ranks(1)));
    fluid_ranks.Resize(number_of_processes-number_of_solid_processes);
    for(int i=1;i<=fluid_ranks.n;i++)fluid_ranks(i)=i+number_of_solid_processes-1;
    fluid_group=new MPI::Group(group->Incl(fluid_ranks.n,&fluid_ranks(1)));

}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPI_SOLID_FLUID_SLIP<TV>::
~MPI_SOLID_FLUID_SLIP()
{
    if(comm){comm->Free();delete comm;}
    if(group){group->Free();delete group;}
    if(fluid_comm){fluid_comm->Free();delete fluid_comm;}
    if(solid_group){solid_group->Free();delete solid_group;}
    if(fluid_group){fluid_group->Free();delete fluid_group;}
}
//#####################################################################
// Function Exchange_Solid_Positions_And_Velocities
//#####################################################################
template<> void MPI_SOLID_FLUID_SLIP<VECTOR<float,1> >::Exchange_Solid_Positions_And_Velocities(SOLID_BODY_COLLECTION<VECTOR<float,1> >& solid_body_collection) const {PHYSBAM_NOT_IMPLEMENTED();}
template<> void MPI_SOLID_FLUID_SLIP<VECTOR<double,1> >::Exchange_Solid_Positions_And_Velocities(SOLID_BODY_COLLECTION<VECTOR<double,1> >& solid_body_collection) const {PHYSBAM_NOT_IMPLEMENTED();}
template<class TV> void MPI_SOLID_FLUID_SLIP<TV>::
Exchange_Solid_Positions_And_Velocities(SOLID_BODY_COLLECTION<TV>& solid_body_collection) const
{
    int tag=Get_Unique_Tag();
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particle;
    if(Solid_Node()){
        ARRAY<ARRAY<char> > send_buffers(fluid_ranks.n);ARRAY<MPI::Request> requests;
        for(int i=1;i<=fluid_ranks.n;i++){
            int buffer_size=MPI_UTILITIES::Pack_Size(particles.X,rigid_body_particles.X,rigid_body_particles.rotation,particles.V,rigid_body_particles.angular_velocity,
                rigid_body_particles.angular_momentum,*comm)+1;
            send_buffers(i).Resize(buffer_size);int position=0;
            MPI_UTILITIES::Pack(particles.X,rigid_body_particles.X,rigid_body_particles.rotation,particles.V,rigid_body_particles.angular_velocity,rigid_body_particles.angular_momentum,
                send_buffers(i),position,*comm);
            requests.Append(comm->Isend(&(send_buffers(i)(1)),position,MPI::PACKED,fluid_ranks(i),tag));}
        MPI_UTILITIES::Wait_All(requests);}
    else{
        int buffer_size=MPI_UTILITIES::Pack_Size(particles.X,rigid_body_particles.X,rigid_body_particles.rotation,particles.V,rigid_body_particles.angular_velocity,
            rigid_body_particles.angular_momentum,*comm)+1;
        ARRAY<char> buffer(buffer_size);int position=0;
        comm->Recv(&buffer(1),buffer_size,MPI::PACKED,solid_node,tag);
        MPI_UTILITIES::Unpack(particles.X,rigid_body_particles.X,rigid_body_particles.rotation,particles.V,rigid_body_particles.angular_velocity,rigid_body_particles.angular_momentum,
            buffer,position,*comm);}
}
//#####################################################################
// Function Fluid_Node
//#####################################################################
template<class TV> bool MPI_SOLID_FLUID_SLIP<TV>::
Fluid_Node() const
{
    return fluid_group->Get_rank()!=MPI::UNDEFINED;
}
//#####################################################################
// Function Solid_Node
//#####################################################################
template<class TV> bool MPI_SOLID_FLUID_SLIP<TV>::
Solid_Node() const
{
    return solid_group->Get_rank()!=MPI::UNDEFINED;
}
//#####################################################################
// Function Create_Fluid_Comm_For_Solid_Nodes
//#####################################################################
template<class TV> void MPI_SOLID_FLUID_SLIP<TV>::
Create_Fluid_Comm_For_Solid_Nodes() const
{
    MPI::COMM_WORLD.Create(*fluid_group);
}
//#####################################################################
// Function Reduce_Add
//#####################################################################
template<class TV> template<class T2> void MPI_SOLID_FLUID_SLIP<TV>::
Reduce_Add(const T2& input,T2& output) const
{
    MPI_UTILITIES::Reduce(input,output,MPI::SUM,*comm);
}
//#####################################################################
// Function Reduce_Min
//#####################################################################
template<class TV> typename TV::SCALAR MPI_SOLID_FLUID_SLIP<TV>::
Reduce_Min(const T local_value) const
{
    T global_value;
    MPI_UTILITIES::Reduce(local_value,global_value,MPI::MIN,*comm);
    return global_value;
}
//#####################################################################
// Function Parallel_Solve_Fluid_Part
//#####################################################################
template<class TV> void MPI_SOLID_FLUID_SLIP<TV>::
Set_MPI_Grid(T_MPI_GRID& mpi_grid_input)
{
    if(!mpi_grid){
        mpi_grid=&mpi_grid_input;
        // allocate partitions and compute neighbor ranks
        if(fluid_comm){fluid_comm->Free();delete fluid_comm;}
        fluid_comm=new MPI::Intracomm;
        *fluid_comm=mpi_grid->comm->Create(*fluid_group);
        partition.Set_Number_Of_Sides(T_PARALLEL_GRID::number_of_faces_per_cell);
        for(int s=1;s<=T_PARALLEL_GRID::number_of_faces_per_cell;s++){
            int global_rank=mpi_grid->side_neighbor_ranks(s);
            if(global_rank==MPI::PROC_NULL) partition.neighbor_ranks(s)=MPI::PROC_NULL;
            else MPI::Group::Translate_ranks(*mpi_grid->group,1,&global_rank,*fluid_group,&partition.neighbor_ranks(s));}
    }
}
//#####################################################################
// Function Parallel_Solve_Fluid_Part
//#####################################################################
template<class TV> void MPI_SOLID_FLUID_SLIP<TV>::
Parallel_Solve_Fluid_Part(FLUID_SYSTEM_MPI_SLIP<TV>& fluid_system,VECTOR_ND<T>& x_array,VECTOR_ND<T>& b_array,VECTOR_ND<T>& p_array,VECTOR_ND<T>& ap_array,
    VECTOR_ND<T>& ar_array,VECTOR_ND<T>& r_array,VECTOR_ND<T>& z_array,VECTOR_ND<T>& zaq_array,const int min_iterations,const int max_iterations,
    const T tolerance,const bool recompute_preconditioner,const bool leakproof_solve)
{
    SYMMQMR_SPARSE_MPI<TV> symmqmr_mpi(leakproof_solve?*fluid_comm:*comm,fluid_comm,&partition);
    symmqmr_mpi.nullspace_tolerance=(T)0;
    symmqmr_mpi.restart_iterations=100;
    //symmqmr_mpi.maximum_iterations=solids_parameters.cg_iterations;
    //LOG::cout<<"Before solve"<<std::endl;
    KRYLOV_VECTOR_WRAPPER<T,VECTOR_ND<T>&> xw(x_array),bw(b_array),pw(p_array),apw(ap_array),arw(ar_array),rw(r_array),zw(z_array),zaqw(zaq_array);
    symmqmr_mpi.Parallel_Solve_Fluid_Part(fluid_system,xw,bw,pw,apw,arw,rw,zw,zaqw,min_iterations,max_iterations,tolerance,recompute_preconditioner);
}
//#####################################################################
// Function Parallel_Solve_Solid_Part
//#####################################################################
template<class TV> void MPI_SOLID_FLUID_SLIP<TV>::
Parallel_Solve_Solid_Part(SOLID_SYSTEM_MPI_SLIP<TV>& solid_system,GENERALIZED_VELOCITY<TV>& x_array,GENERALIZED_VELOCITY<TV>& b_array,GENERALIZED_VELOCITY<TV>& p_array,GENERALIZED_VELOCITY<TV>& ap_array,
    GENERALIZED_VELOCITY<TV>& ar_array,GENERALIZED_VELOCITY<TV>& r_array,GENERALIZED_VELOCITY<TV>& z_array,GENERALIZED_VELOCITY<TV>& zaq_array,const int min_iterations,const int max_iterations,const T tolerance)
{
    SYMMQMR_SPARSE_MPI<TV> symmqmr_mpi(*comm,0,0);
    symmqmr_mpi.nullspace_tolerance=(T)0;
    symmqmr_mpi.restart_iterations=100;
    //symmqmr_mpi.maximum_iterations=solids_parameters.cg_iterations;
    symmqmr_mpi.Parallel_Solve_Solid_Part(solid_system,x_array,b_array,p_array,ap_array,ar_array,r_array,z_array,zaq_array,min_iterations,max_iterations,tolerance);
}
//#####################################################################
// Function Exchange_Fluid_Domain_To_Solid
//#####################################################################
template<class TV> void MPI_SOLID_FLUID_SLIP<TV>::
Exchange_Coupled_Deformable_Particle_List(ARRAY<int>* fluid_list,ARRAY<ARRAY<int> >* results)
{
    int tag=Get_Unique_Tag();
    if(Solid_Node()){
        for(int i=1;i<=fluid_ranks.n;i++){
            MPI::Status status;
            comm->Probe(MPI::ANY_SOURCE,tag,status);
            int source=status.Get_source();
            ARRAY<char> buffer(status.Get_count(MPI::PACKED));int position=0;
            comm->Recv(&buffer(1),buffer.m,MPI::PACKED,source,tag);
            MPI_UTILITIES::Unpack((*results)(source),buffer,position,*comm);}}
    else{
        int buffer_size=MPI_UTILITIES::Pack_Size(*fluid_list,*comm)+1;
        ARRAY<char> buffer(buffer_size);int position=0;
        MPI_UTILITIES::Pack(*fluid_list,buffer,position,*comm);
        comm->Send(&buffer(1),buffer_size,MPI::PACKED,solid_node,tag);}
}
//#####################################################################
// Function Find_Matrix_Indices
//#####################################################################
template<class TV> void MPI_SOLID_FLUID_SLIP<TV>::
Find_Matrix_Indices(const GRID<TV>& local_grid,const T_ARRAYS_BOOL& valid_divergence_cells,T_ARRAYS_INT& cell_index_to_matrix_index,T_FACE_ARRAYS_INT& face_ghost_cell_index,
    ARRAY<int>& face_ghost_cell_index_map,T_FACE_ARRAYS_INT& face_lambdas,INTERVAL<int>& divergence_indices,int& cell_count)
{
    assert(local_grid.Is_MAC_Grid());
    cell_count=0;
    RANGE<TV_INT> domain_indices=local_grid.Domain_Indices();
    Find_Matrix_Indices_In_Region(local_grid,valid_divergence_cells,0,domain_indices,cell_index_to_matrix_index,face_ghost_cell_index,face_ghost_cell_index_map,face_lambdas,divergence_indices,cell_count);
    INTERVAL<int> dummy_box;
    for(int axis=1;axis<=TV::dimension;axis++){
        RANGE<TV_INT> region(domain_indices);
        region.min_corner(axis)=region.max_corner(axis)=domain_indices.min_corner(axis)-1;
        Find_Matrix_Indices_In_Region(local_grid,valid_divergence_cells,2*(axis-1)+1,region,cell_index_to_matrix_index,face_ghost_cell_index,face_ghost_cell_index_map,face_lambdas,dummy_box,cell_count);
        region.min_corner(axis)=region.max_corner(axis)=domain_indices.max_corner(axis)+1;
        Find_Matrix_Indices_In_Region(local_grid,valid_divergence_cells,2*(axis-1)+2,region,cell_index_to_matrix_index,face_ghost_cell_index,face_ghost_cell_index_map,face_lambdas,dummy_box,cell_count);
    }
    for(int axis=1;axis<=TV::dimension;axis++){
        RANGE<TV_INT> region(domain_indices);
        region.min_corner(axis)=region.max_corner(axis)=domain_indices.min_corner(axis);
        Find_Boundary_Indices_In_Region(local_grid,valid_divergence_cells,2*(axis-1)+1,region,cell_index_to_matrix_index,face_ghost_cell_index,face_lambdas);
        region.min_corner(axis)=region.max_corner(axis)=domain_indices.max_corner(axis);
        Find_Boundary_Indices_In_Region(local_grid,valid_divergence_cells,2*(axis-1)+2,region,cell_index_to_matrix_index,face_ghost_cell_index,face_lambdas);
    }
}
//#####################################################################
// Function Find_Matrix_Indices_In_Region
//#####################################################################
template<class TV> void MPI_SOLID_FLUID_SLIP<TV>::
Find_Matrix_Indices_In_Region(const GRID<TV>& local_grid,const T_ARRAYS_BOOL& valid_divergence_cells,const int region_index,const RANGE<TV_INT>& region,
    T_ARRAYS_INT& cell_index_to_matrix_index,T_FACE_ARRAYS_INT& face_ghost_cell_index,ARRAY<int>& face_ghost_cell_index_map,T_FACE_ARRAYS_INT& face_lambdas,INTERVAL<int>& divergence_indices,int& cell_count)
{
    if(region_index)
        partition.ghost_indices(region_index).min_corner=cell_count+1;
    else
        partition.interior_indices.min_corner=cell_count+1;

    divergence_indices.min_corner=cell_count+1;
#ifdef BRICK
    int count_before=cell_count;
#endif
    for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(!poisson->psi_D(cell_index) && valid_divergence_cells(cell_index)){
            cell_index_to_matrix_index(cell_index)=++cell_count;
        }
    }
    divergence_indices.max_corner=cell_count;

    for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(!poisson->psi_D(cell_index) && !valid_divergence_cells(cell_index)){
            cell_index_to_matrix_index(cell_index)=++cell_count;
        }
    }

#ifdef BRICK
    std::stringstream ss;
    ss<<cell_count-count_before<<" real cells in region "<<region_index<<std::endl;
    LOG::filecout(ss.str());
    count_before=cell_count;
#endif

    // TODO: set the region/axis properly.  Probably, expand by one cell, see if ghost cell is inside region
    /*typename FACE_ITERATOR::T_REGION region_type;
    int side=0;
    if(region_index){
        region_type=GRID<TV>::BOUNDARY_REGION;
        side=region_index;
    }
    else
        region_type=GRID<TV>::WHOLE_REGION;*/

    RANGE<TV_INT> expanded_domain_indices=local_grid.Domain_Indices().Thickened(1);
    for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        for(int axis=1;axis<=TV::dimension;axis++){
            TV_INT axis_vector=TV_INT::Axis_Vector(axis);
            int first_ghost_cell_index=face_ghost_cell_index(1,axis,cell_index); // looking IN to this cell
            if(first_ghost_cell_index/* && expanded_domain_indices.Lazy_Inside(cell_index-axis_vector)*/){
                //if(region_index)
                    //    DEBUG_UTILITIES::Debug_Breakpoint();
#ifdef BRICK
                if(!region_index){
                    std::stringstream ss;
                    ss<<"Found boundary ghost cell across face (1, "<<axis<<", "<<cell_index<<") with temp index "<<first_ghost_cell_index<<std::endl;
                    LOG::filecout(ss.str());}
#endif
                if(!face_ghost_cell_index_map(first_ghost_cell_index))
                    face_ghost_cell_index_map(first_ghost_cell_index)=++cell_count;
            }
            int second_ghost_cell_index=face_ghost_cell_index(2,axis,cell_index+axis_vector);
            if(second_ghost_cell_index/* && expanded_domain_indices.Lazy_Inside(cell_index+axis_vector)*/){
#ifdef BRICK
                if(!region_index){
                    std::stringstream ss;
                    ss<<"Found boundary ghost cell across face (2, "<<axis<<", "<<cell_index+axis_vector<<") with temp index "<<second_ghost_cell_index<<std::endl;
                    LOG::filecout(ss.str());}
#endif
                if(!face_ghost_cell_index_map(second_ghost_cell_index))
                    face_ghost_cell_index_map(second_ghost_cell_index)=++cell_count;}}}

    for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        for(int axis=1;axis<=TV::dimension;axis++){
            TV_INT axis_vector=TV_INT::Axis_Vector(axis);
            int first_ghost_cell_index=face_ghost_cell_index(1,axis,cell_index);
            if(first_ghost_cell_index){
#ifdef BRICK
                if(!region_index){
                    std::stringstream ss;
                    ss<<"Found boundary ghost cell across face (1, "<<axis<<", "<<cell_index<<") with index "<<face_ghost_cell_index_map(first_ghost_cell_index)<<std::endl;
                    LOG::filecout(ss.str());}
#endif
                face_ghost_cell_index(1,axis,cell_index)=face_ghost_cell_index_map(first_ghost_cell_index);
            }
            int second_ghost_cell_index=face_ghost_cell_index(2,axis,cell_index+axis_vector);
            if(second_ghost_cell_index){
#ifdef BRICK
                if(!region_index){
                    std::stringstream ss;
                    ss<<"Found boundary ghost cell across face (2, "<<axis<<", "<<cell_index+axis_vector<<") with index"<<face_ghost_cell_index_map(second_ghost_cell_index)<<std::endl;
                    LOG::filecout(ss.str());}
#endif
                face_ghost_cell_index(2,axis,cell_index+axis_vector)=face_ghost_cell_index_map(second_ghost_cell_index);}}}

#ifdef BRICK
    std::stringstream ss;
    ss<<cell_count-count_before<<" ghost cells in region "<<region_index<<std::endl;
    LOG::filecout(ss.str());
    count_before=cell_count;
#endif
    
    /*for(FACE_ITERATOR iterator(local_grid,0,region_type,side);iterator.Valid();iterator.Next()){
        int axis=iterator.Axis();
        TV_INT face_index=iterator.Face_Index();
        for(int side=1;side<=2;side++){
            TV_INT cell_index=face_index+(1-side)*TV_INT::Axis_Vector(axis);
            int ghost_cell_index=face_ghost_cell_index(side,axis,face_index);
            if(ghost_cell_index && region.Lazy_Inside(cell_index)){
                if(!face_ghost_cell_index_map(ghost_cell_index))
                    face_ghost_cell_index_map(ghost_cell_index)=++cell_count;}}}

    for(FACE_ITERATOR iterator(local_grid,0,region_type,side);iterator.Valid();iterator.Next()){
        int axis=iterator.Axis();
        TV_INT face_index=iterator.Face_Index();
        for(int side=1;side<=2;side++){
            TV_INT cell_index=face_index+(1-side)*TV_INT::Axis_Vector(axis);
            int ghost_cell_index=face_ghost_cell_index(side,axis,face_index);
            if(ghost_cell_index && region.Lazy_Inside(cell_index))
            face_ghost_cell_index(side,axis,face_index)=face_ghost_cell_index_map(ghost_cell_index);}}*/

    //PHYSBAM_FATAL_ERROR("pretty sure don't want Lazy_Inside");
    for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        //if(cell_index(1)==0 && cell_index(2)==1)
        //DEBUG_UTILITIES::Debug_Breakpoint();
        for(int axis=1;axis<=TV::dimension;axis++){
            TV_INT axis_vector=TV_INT::Axis_Vector(axis);
            if(expanded_domain_indices.Lazy_Inside(cell_index-axis_vector) && face_lambdas(2,axis,cell_index)!=0){
#ifdef BRICK
                std::stringstream ss;
                ss<<"Lambda in ghost region at (2, "<<axis<<", "<<cell_index<<std::endl;
                LOG::filecout(ss.str());
#endif
                face_lambdas(2,axis,cell_index)=++cell_count;
            }
            if(expanded_domain_indices.Lazy_Inside(cell_index+axis_vector) && face_lambdas(1,axis,cell_index+TV_INT::Axis_Vector(axis))!=0){
#ifdef BRICK
                std::stringstream ss;ss<<"Lambda in ghost region at (1, "<<axis<<", "<<cell_index+axis_vector<<std::endl;LOG::filecout(ss.str());
#endif
                face_lambdas(1,axis,cell_index+TV_INT::Axis_Vector(axis))=++cell_count;}}}

#ifdef BRICK
    std::stringstream ss;ss<<cell_count-count_before<<" face lambdas in region "<<region_index<<std::endl;LOG::filecout(ss.str());
#endif

    if(region_index)
        partition.ghost_indices(region_index).max_corner=cell_count;
    else
        partition.interior_indices.max_corner=cell_count;
}
//#####################################################################
// Function Find_Matrix_Indices_In_Region
//#####################################################################
template<class TV> void MPI_SOLID_FLUID_SLIP<TV>::
Find_Boundary_Indices_In_Region(const GRID<TV>& local_grid,const T_ARRAYS_BOOL& valid_divergence_cells,const int domain_side,const RANGE<TV_INT>& region,const T_ARRAYS_INT& cell_index_to_matrix_index,const T_FACE_ARRAYS_INT& face_ghost_cell_index,const T_FACE_ARRAYS_INT& face_lambdas)
{
    int axis=(domain_side-1)/2+1,cell_side=domain_side&1;
    RANGE<TV_INT> face_region=region;if(!cell_side) face_region+=TV_INT::Axis_Vector(axis);
    //int side=domain_side%2?2:1;
    TV_INT face_offset=cell_side?TV_INT():TV_INT::Axis_Vector(axis);
    int boundary_cell_count=0;
    for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(/*!poisson->psi_N(axis,cell_index+face_offset) && */!poisson->psi_D(cell_index))
            //(face_ghost_cell_index(side,axis,cell_index+face_offset) || (!face_ghost_cell_index(side,axis,cell_index+face_offset) && !poisson->psi_D(cell_index))))
            ++boundary_cell_count;}

    RANGE<TV_INT> domain_indices_shrink_axis=local_grid.Domain_Indices();domain_indices_shrink_axis.Change_Size(-TV_INT::Axis_Vector(axis));
    ARRAY<int> boundary_ghost_cells;
    for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        for(int second_axis=1;second_axis<=TV::dimension;second_axis++){
            TV_INT axis_vector=TV_INT::Axis_Vector(second_axis);
            if(!poisson->psi_N(second_axis,cell_index) && face_ghost_cell_index(1,second_axis,cell_index) && domain_indices_shrink_axis.Lazy_Outside(cell_index-axis_vector)){
#ifdef BRICK
                std::stringstream ss;
                ss<<"Found boundary ghost cell across face (1, "<<second_axis<<", "<<cell_index<<") with index "<<face_ghost_cell_index(1,second_axis,cell_index)<<std::endl;
                LOG::filecout(ss.str());
#endif
                boundary_ghost_cells.Append_Unique(face_ghost_cell_index(1,second_axis,cell_index));
            }
            if(!poisson->psi_N(second_axis,cell_index+axis_vector) && face_ghost_cell_index(2,second_axis,cell_index+axis_vector) && domain_indices_shrink_axis.Lazy_Outside(cell_index+axis_vector)){
#ifdef BRICK
                std::stringstream ss;
                ss<<"Found boundary ghost cell across face (2, "<<second_axis<<", "<<cell_index+axis_vector<<") with index "<<face_ghost_cell_index(2,second_axis,cell_index+axis_vector)<<std::endl;
                LOG::filecout(ss.str());
#endif
                boundary_ghost_cells.Append_Unique(face_ghost_cell_index(2,second_axis,cell_index+axis_vector));}}}
    boundary_cell_count+=boundary_ghost_cells.m;
    
    RANGE<TV_INT> contracted_domain_indices=local_grid.Domain_Indices();contracted_domain_indices.Change_Size(-TV_INT::Axis_Vector(axis));
    for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        for(int second_axis=1;second_axis<=TV::dimension;second_axis++){
            TV_INT axis_vector=TV_INT::Axis_Vector(second_axis);
            if(contracted_domain_indices.Lazy_Outside(cell_index-axis_vector) && face_lambdas(2,second_axis,cell_index)!=0)
                ++boundary_cell_count;
            if(contracted_domain_indices.Lazy_Outside(cell_index+axis_vector) && face_lambdas(1,second_axis,cell_index+TV_INT::Axis_Vector(second_axis))!=0)
                ++boundary_cell_count;}}

    partition.boundary_indices(domain_side).Resize(boundary_cell_count);
    boundary_cell_count=0;

#ifdef BRICK
    int count_before=boundary_cell_count;
#endif
    for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        //if(!poisson->psi_N(axis,cell_index+face_offset)){
        if(/*!face_ghost_cell_index(side,axis,cell_index+face_offset) &&*/ !poisson->psi_D(cell_index) && valid_divergence_cells(cell_index)){
                partition.boundary_indices(domain_side)(++boundary_cell_count)=cell_index_to_matrix_index(cell_index);
        }
            /*else if(face_ghost_cell_index(side,axis,cell_index+face_offset)){
                DEBUG_UTILITIES::Debug_Breakpoint();
                partition.boundary_indices(domain_side)(++boundary_cell_count)=face_ghost_cell_index(side,axis,cell_index+face_offset);
                }*/
        //}
    }

    for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        //if(!poisson->psi_N(axis,cell_index+face_offset)){
        if(/*!face_ghost_cell_index(side,axis,cell_index+face_offset) &&*/ !poisson->psi_D(cell_index) && !valid_divergence_cells(cell_index)){
                partition.boundary_indices(domain_side)(++boundary_cell_count)=cell_index_to_matrix_index(cell_index);
        }
            /*else if(face_ghost_cell_index(side,axis,cell_index+face_offset)){
                DEBUG_UTILITIES::Debug_Breakpoint();
                partition.boundary_indices(domain_side)(++boundary_cell_count)=face_ghost_cell_index(side,axis,cell_index+face_offset);
                }*/
        //}
    }

#ifdef BRICK
    std::stringstream ss;
    ss<<boundary_cell_count-count_before<<" real cells in boundary region "<<domain_side<<std::endl;
    LOG::filecout(ss.str());
    count_before=boundary_cell_count;
#endif


    /*for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        for(int second_axis=1;second_axis<=TV::dimension;second_axis++){
            TV_INT second_axis_vector=TV_INT::Axis_Vector(second_axis);
            if(!poisson->psi_N(second_axis,cell_index) && face_ghost_cell_index(1,second_axis,cell_index) && domain_indices_shrink_axis.Lazy_Outside(cell_index-second_axis_vector))
                partition.boundary_indices(domain_side)(++boundary_cell_count)=face_ghost_cell_index(1,second_axis,cell_index);
            if(!poisson->psi_N(second_axis,cell_index+second_axis_vector) && face_ghost_cell_index(2,second_axis,cell_index+second_axis_vector) && domain_indices_shrink_axis.Lazy_Outside(cell_index+second_axis_vector))
            partition.boundary_indices(domain_side)(++boundary_cell_count)=face_ghost_cell_index(2,second_axis,cell_index+second_axis_vector);}}*/
    for(int i=1;i<=boundary_ghost_cells.m;i++){
        partition.boundary_indices(domain_side)(++boundary_cell_count)=boundary_ghost_cells(i);
    }

#ifdef BRICK
    std::stringstream ss;
    ss<<boundary_cell_count-count_before<<" ghost cells in boundary region "<<domain_side<<std::endl;
    LOG::filecout(ss.str());
    count_before=boundary_cell_count;
#endif

    for(CELL_ITERATOR iterator(local_grid,region);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        for(int second_axis=1;second_axis<=TV::dimension;second_axis++){
            TV_INT axis_vector=TV_INT::Axis_Vector(second_axis);
            if(contracted_domain_indices.Lazy_Outside(cell_index-axis_vector) && face_lambdas(2,second_axis,cell_index)!=0){
#ifdef BRICK
                std::stringstream ss;ss<<"Lambda in boundary region at (2, "<<second_axis<<", "<<cell_index<<std::endl;LOG::filecout(ss.str());
#endif
                partition.boundary_indices(domain_side)(++boundary_cell_count)=face_lambdas(2,second_axis,cell_index);
            }
            if(contracted_domain_indices.Lazy_Outside(cell_index+axis_vector) && face_lambdas(1,second_axis,cell_index+axis_vector)!=0){
#ifdef BRICK
                std::stringstream ss;ss<<"Lambda in boundary region at (1, "<<second_axis<<", "<<cell_index+axis_vector<<std::endl;LOG::filecout(ss.str());
#endif
                partition.boundary_indices(domain_side)(++boundary_cell_count)=face_lambdas(1,second_axis,cell_index+axis_vector);}}}

#ifdef BRICK
    std::stringstream ss;ss<<boundary_cell_count-count_before<<" face lambdas in boundary region "<<domain_side<<std::endl;LOG::filecout(ss.str());
#endif

}
//#####################################################################
#else
//#####################################################################
template<class TV> MPI_SOLID_FLUID_SLIP<TV>::MPI_SOLID_FLUID_SLIP(){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> MPI_SOLID_FLUID_SLIP<TV>::~MPI_SOLID_FLUID_SLIP(){}
template<class TV> void MPI_SOLID_FLUID_SLIP<TV>::Set_MPI_Grid(T_MPI_GRID& mpi_grid_input){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_SOLID_FLUID_SLIP<TV>::Exchange_Solid_Positions_And_Velocities(SOLID_BODY_COLLECTION<TV>& solid_body_collection) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool MPI_SOLID_FLUID_SLIP<TV>::Fluid_Node() const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> bool MPI_SOLID_FLUID_SLIP<TV>::Solid_Node() const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_SOLID_FLUID_SLIP<TV>::Create_Fluid_Comm_For_Solid_Nodes() const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> template<class T2> void MPI_SOLID_FLUID_SLIP<TV>::Reduce_Add(const T2& input,T2& output) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> typename TV::SCALAR MPI_SOLID_FLUID_SLIP<TV>::Reduce_Min(const T local_value) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_SOLID_FLUID_SLIP<TV>::Parallel_Solve_Fluid_Part(FLUID_SYSTEM_MPI_SLIP<TV>& fluid_system,VECTOR_ND<T>& x_array,VECTOR_ND<T>& b_array,VECTOR_ND<T>& p_array,VECTOR_ND<T>& ap_array,
    VECTOR_ND<T>& ar_array,VECTOR_ND<T>& r_array,VECTOR_ND<T>& z_array,VECTOR_ND<T>& zaq_array,const int min_iterations,const int max_iterations,
    const T tolerance,const bool recompute_preconditioner,const bool leakproof_solve){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_SOLID_FLUID_SLIP<TV>::Parallel_Solve_Solid_Part(SOLID_SYSTEM_MPI_SLIP<TV>& solid_system,GENERALIZED_VELOCITY<TV>& x_array,GENERALIZED_VELOCITY<TV>& b_array,GENERALIZED_VELOCITY<TV>& p_array,GENERALIZED_VELOCITY<TV>& ap_array,
    GENERALIZED_VELOCITY<TV>& ar_array,GENERALIZED_VELOCITY<TV>& r_array,GENERALIZED_VELOCITY<TV>& z_array,GENERALIZED_VELOCITY<TV>& zaq_array,const int min_iterations,const int max_iterations,const T tolerance){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_SOLID_FLUID_SLIP<TV>::Exchange_Coupled_Deformable_Particle_List(ARRAY<int>* fluid_list,ARRAY<ARRAY<int> >* results){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_SOLID_FLUID_SLIP<TV>::Find_Matrix_Indices(const GRID<TV>& local_grid,const T_ARRAYS_BOOL& valid_divergence_cells,T_ARRAYS_INT& cell_index_to_matrix_index,T_FACE_ARRAYS_INT& face_ghost_cell_index,ARRAY<int>& face_ghost_cell_index_map,T_FACE_ARRAYS_INT& face_lambdas,INTERVAL<int>& divergence_indices,int &cell_count){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_SOLID_FLUID_SLIP<TV>::Find_Boundary_Indices_In_Region(const GRID<TV>& local_grid,const T_ARRAYS_BOOL& valid_divergence_cells,const int domain_side,const RANGE<TV_INT>& region,const T_ARRAYS_INT& cell_index_to_matrix_index,const T_FACE_ARRAYS_INT& face_ghost_cell_index,const T_FACE_ARRAYS_INT& face_lambdas){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
template<class TV> void MPI_SOLID_FLUID_SLIP<TV>::Find_Matrix_Indices_In_Region(const GRID<TV>& local_grid,const T_ARRAYS_BOOL& valid_divergence_cells,const int region_index,const RANGE<TV_INT>& region,
    T_ARRAYS_INT& cell_index_to_matrix_index,T_FACE_ARRAYS_INT& face_ghost_cell_index,ARRAY<int>& face_ghost_cell_index_map,T_FACE_ARRAYS_INT& face_lambdas,INTERVAL<int>& divergence_indices,int& cell_count){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
#endif
#define INSTANTIATION_HELPER_UNIFORM(T,d) \
    template class MPI_SOLID_FLUID_SLIP<VECTOR<T,d> >; \
    template void MPI_SOLID_FLUID_SLIP<VECTOR<T,d> >::Reduce_Add(const ARRAY<T>&,ARRAY<T>&) const;
#define INSTANTIATION_HELPER(T) \
    INSTANTIATION_HELPER_UNIFORM(T,1);INSTANTIATION_HELPER_UNIFORM(T,2);INSTANTIATION_HELPER_UNIFORM(T,3);
INSTANTIATION_HELPER(float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double);
#endif
}
