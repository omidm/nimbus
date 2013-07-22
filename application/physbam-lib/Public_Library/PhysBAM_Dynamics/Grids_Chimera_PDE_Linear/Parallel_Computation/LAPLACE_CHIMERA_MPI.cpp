//#####################################################################
// Copyright 2012
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_CHIMERA_MPI
//#####################################################################
#include <PhysBAM_Dynamics/Grids_Chimera_PDE_Linear/Parallel_Computation/LAPLACE_CHIMERA_MPI.h>
#include <PhysBAM_Dynamics/Grids_Chimera_PDE_Linear/Parallel_Computation/LAPLACE_CHIMERA_GRID_MPI.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Utilities/INTERRUPTS.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Dynamics/Meshing/POLYGONAL_TRIANGULATION.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/POLYGON_HYPERPLANE_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/POLYGON.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Utilities/TIMER.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
using namespace PhysBAM;
//#####################################################################
// Macros for chimera grid and data access
//#####################################################################
#define PHI(I) (*chimera_grid.incompressible_fluid_containers)(I)->particle_levelset_evolution.phi
#define FACE_VELOCITIES(I) (*chimera_grid.incompressible_fluid_containers)(I)->face_velocities
#define PRESSURE(I) (*chimera_grid.incompressible_fluid_containers)(I)->pressure
#define PSI_N(I) (*chimera_grid.incompressible_fluid_containers)(I)->psi_N
#define PSI_D(I) (*chimera_grid.incompressible_fluid_containers)(I)->psi_D
//#####################################################################
// Psi_D
//#####################################################################
template<class T_GRID> bool LAPLACE_CHIMERA_MPI<T_GRID>::Psi_D(const int grid_index,const TV_INT& cell_index,T& pressure) const
{
    if(laplace_grid.Local_Grid(grid_index)){
        pressure=PRESSURE(laplace_grid.Local_Grid_Index(grid_index))(cell_index);
        return PSI_D(laplace_grid.Local_Grid_Index(grid_index))(cell_index);}
    else if(boundary_cell_psi_D(grid_index)(laplace_grid.boundary_cell_indices_to_linear_index(grid_index).Get(cell_index)).x){
        pressure=boundary_cell_psi_D(grid_index)(laplace_grid.boundary_cell_indices_to_linear_index(grid_index).Get(cell_index)).y;
        return true;}
    return false;
}
//#####################################################################
// Phi
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LAPLACE_CHIMERA_MPI<T_GRID>::Phi(const int grid_index,const TV_INT& cell_index) const
{
    //ASSUME THAT GHOST PHI IS VALID
    //ONE GHOST CELL IS NEEDED FOR COMPUTING BOUNDARY FACE VOLUME FRACTIONS
    //EETODO: RETURN PHI FOR BOUNDARY CELLS
    if(laplace_grid.Local_Grid(grid_index)) return PHI(laplace_grid.Local_Grid_Index(grid_index))(cell_index);
    else return boundary_cell_phi(grid_index)(laplace_grid.boundary_cell_indices_to_linear_index(grid_index).Get(cell_index));
}
//#####################################################################
// Neumann_Pocket
//#####################################################################
template<class T_GRID> bool LAPLACE_CHIMERA_MPI<T_GRID>::Neumann_Pocket(const int grid_index,const TV_INT& cell_index) const
{
    if(laplace_grid.Local_Grid(grid_index))
        for(int face_axis=1;face_axis<=TV::dimension;face_axis++)
            for(int face_side=1;face_side<=2;face_side++){
                D_FACE_INDEX face_index(face_axis,cell_index+(face_side-1)*TV_INT::Axis_Vector(face_axis));
                if(laplace_grid.Chimera_Face(grid_index,face_index))
                    if(!PSI_N(laplace_grid.Local_Grid_Index(grid_index))(face_index))
                        return false;}
    
    //LOG::cout << "neumann pocket " << grid_index << " " << cell_index << std::endl;
    if(laplace_grid.Boundary_Cell(grid_index,cell_index)){
        GRID_CELL_INDEX grid_cell_index(grid_index,cell_index);
        if(laplace_grid.cell_indices_to_incident_voronoi_face_indices.Contains(grid_cell_index)){
            const ARRAY<int>& incident_voronoi_face_indices=laplace_grid.cell_indices_to_incident_voronoi_face_indices.Get(grid_cell_index);
            for(int incident_index=1;incident_index<=incident_voronoi_face_indices.Size();incident_index++){
                const VORONOI_FACE_INDICES& indices=laplace_grid.voronoi_faces(incident_voronoi_face_indices(incident_index)).x;
                if((laplace_grid.Local_Grid(indices(1).x) || laplace_grid.Local_Grid(indices(2).x)) && !voronoi_face_psi_N_velocities(incident_voronoi_face_indices(incident_index)).x)
                    return false;}}}
    
    return true;
}
//#####################################################################
// Dual_Cell_Volume_Fraction
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LAPLACE_CHIMERA_MPI<T_GRID>::
Dual_Cell_Volume_Fraction(const int grid_index_1,const TV_INT& cell_index_1,const int grid_index_2,const TV_INT& cell_index_2) const
{
    return 1;
    /*T minimum_dual_cell_volume_fractions=1e-5;

    if(!example.fluids_parameters.water)
        return 1;
    
    T phi_1=Phi(grid_index_1,cell_index_1);
    T phi_2=Phi(grid_index_2,cell_index_2);
    
    if(phi_1>0 && phi_2>0)
        return 0;
    
    if(phi_1<=0 && phi_2<=0)
        return 1;
    
    T min_phi=min(phi_1,phi_2);
    T max_phi=max(phi_1,phi_2);
    return max(minimum_dual_cell_volume_fractions,min((T)1,min_phi/(min_phi-max_phi)));*/
}

//#####################################################################
// Matrix_Cell_Index
//#####################################################################
template<class T_GRID> int LAPLACE_CHIMERA_MPI<T_GRID>::Matrix_Cell_Index(const int grid_index,const TV_INT& cell_index) const
{
    if(laplace_grid.Local_Grid(grid_index)){
        if(laplace_grid.Grid(grid_index).Inside_Domain(cell_index)) return cell_indices_to_matrix_cell_indices(grid_index)(cell_index);}
    else return boundary_cell_indices_to_matrix_cell_indices(grid_index).Get_Default(cell_index,0);
    return 0;
}
//#####################################################################
// Matrix_Face_Index
//#####################################################################
template<class T_GRID> int LAPLACE_CHIMERA_MPI<T_GRID>::Matrix_Face_Index(const int grid_index,const D_FACE_INDEX& face_index) const
{
    if(laplace_grid.Local_Grid(grid_index)) return face_indices_to_matrix_face_indices(grid_index)(face_index);
    else return 0;
}
//#####################################################################
// Matrix_Face_Index
//#####################################################################
template<class T_GRID> int LAPLACE_CHIMERA_MPI<T_GRID>::
Matrix_Face_Index(const int voronoi_face_index) const
{
    return voronoi_face_indices_to_matrix_face_indices(voronoi_face_index);
}
//#####################################################################
// Compute_Local_And_Boundary_Neumann_Pockets
//#####################################################################
/*template<class T_GRID> void LAPLACE_CHIMERA_MPI<T_GRID>::
Compute_Local_And_Boundary_Neumann_Pockets(const T time,const T dt)
{
    cell_indices_to_valid_pressure_cells.Resize(n_global_grids);
    face_invalid_indices.Resize(n_local_grids);
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
        cell_indices_to_valid_pressure_cells(grid_index).Resize(Grid(grid_index).Domain_Indices(1),true,false,true);
        if(Local_Grid(grid_index)) face_invalid_indices(Local_Grid_Index(grid_index)).Resize(TV::dimension);
        for(CELL_ITERATOR iterator(Grid(grid_index));iterator.Valid();iterator.Next()){
            TV_INT cell_index=iterator.Cell_Index();
            if(Local_Grid(grid_index) || Boundary_Cell(grid_index,cell_index))
                cell_indices_to_valid_pressure_cells(grid_index)(cell_index)=!(Neumann_Pocket(grid_index,cell_index));
        }}
}*/

//#####################################################################
// Active_Cell
//#####################################################################
template<class T_GRID> bool LAPLACE_CHIMERA_MPI<T_GRID>::Active_Cell_Compute(const int grid_index,const TV_INT& cell_index) const
{
    if(!laplace_grid.Chimera_Cell(grid_index,cell_index)) return false;
    //if(example.fluids_parameters.water && Phi(grid_index,cell_index)>0) return false;
    T pressure;
    if(Psi_D(grid_index,cell_index,pressure)) return false; // Dirichlet cells
    //LOG::cout << "cell not dirichlet " << grid_index << " " << cell_index << std::endl;
    return !Neumann_Pocket(grid_index,cell_index); // Neumann pockets
}
//#####################################################################
// Function Compute_Matrix_Indices
//#####################################################################
template<class T_GRID> void LAPLACE_CHIMERA_MPI<T_GRID>::Construct_Matrix_Indices(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time)
{
    LOG::SCOPE scope("Compute_Matrix_Indices");

    //LOG::SCOPE scope("Compute_Boundary_Conditions");
    
    voronoi_face_psi_N_velocities.Resize(laplace_grid.voronoi_faces.Size());;

    //QUERY EXAMPLE FOR NEUMANN BOUNDARY CONDITIONS
    T normal_gradient;
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index))
            for(FACE_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                D_FACE_INDEX face_index=iterator.Full_Index();
                TV normal=laplace_grid.Rigid_Grid(grid_index).Frame().r.Rotate(TV::Axis_Vector(iterator.Axis()));
                TV location=laplace_grid.Rigid_Grid(grid_index).Frame()*iterator.Location();
                PSI_N(laplace_grid.Local_Grid_Index(grid_index))(face_index)=callbacks.Get_Neumann_Boundary_Condition(location,normal,normal_gradient,time);}
    voronoi_face_psi_N_velocities.Resize(laplace_grid.voronoi_faces.Size());
    for(int voronoi_face_index=1;voronoi_face_index<=laplace_grid.voronoi_faces.Size();voronoi_face_index++){
        GRID_CELL_INDEX& grid_cell_index_1=laplace_grid.voronoi_faces(voronoi_face_index).x(1);
        GRID_CELL_INDEX& grid_cell_index_2=laplace_grid.voronoi_faces(voronoi_face_index).x(2);
        TV location_1=laplace_grid.Rigid_Grid(grid_cell_index_1.x).Frame()*laplace_grid.Grid(grid_cell_index_1.x).X(grid_cell_index_1.y);
        TV location_2=laplace_grid.Rigid_Grid(grid_cell_index_2.x).Frame()*laplace_grid.Grid(grid_cell_index_2.x).X(grid_cell_index_2.y);
        PAIR<bool,T>& psi_N_velocity=voronoi_face_psi_N_velocities(voronoi_face_index);
        psi_N_velocity.x=callbacks.Get_Neumann_Boundary_Condition((location_1+location_2)/2,(location_1-location_2).Normalized(),normal_gradient,time);}
    
    //QUERY EXAMPLE FOR DIRICHLET BOUNDARY_CONDITIONS
    boundary_cell_psi_D.Resize(n_global_grids);
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
        if(laplace_grid.Local_Grid(grid_index))
            for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index),1);iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                TV location=laplace_grid.Frame(grid_index)*iterator.Location();
                T value;
                PSI_D(laplace_grid.Local_Grid_Index(grid_index))(cell_index)=callbacks.Get_Dirichlet_Boundary_Condition(location,value,time);}
        else if(laplace_grid.Boundary_Grid(grid_index))
            for(int boundary_cell_index=1;boundary_cell_index<=laplace_grid.boundary_cell_indices(grid_index).Size();boundary_cell_index++){
                TV_INT cell_index=laplace_grid.boundary_cell_indices(grid_index)(boundary_cell_index);
                TV location=laplace_grid.Frame(grid_index)*laplace_grid.Grid(grid_index).X(cell_index);
                T value;
                bool psi_d=callbacks.Get_Dirichlet_Boundary_Condition(location,value,time);
                boundary_cell_psi_D(grid_index).Append(PAIR<bool,T>(psi_d,0));}}
    
    //boundary_cell_indices_to_active_cells.Clean_Memory();
    /*boundary_cell_active.Resize(n_global_grids);
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
        if(Local_Grid(grid_index))
        {
            boundary_cell_active(grid_index).Resize(boundary_cell_indices(grid_index).Size());
            for(int boundary_cell_index=1;boundary_cell_index<=boundary_cell_indices(grid_index).Size();boundary_cell_index++)
                boundary_cell_active(grid_index)(boundary_cell_index)=Active_Cell_Compute(grid_index,boundary_cell_indices(grid_index)(boundary_cell_index));
        }
        else boundary_cell_active(grid_index).Resize(boundary_cell_indices(grid_index).Size());
    }

    chimera_grid.Exchange_Active_Boundary_Cells(boundary_cell_active,is_boundary_grid);*/
    
    /*boundary_cell_indices_to_active_cells.Clean_Memory();
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(Boundary_Grid(grid_index))
            for(int boundary_cell_index=1;boundary_cell_index<=boundary_cell_indices(grid_index).Size();boundary_cell_index++){
            boundary_cell_indices_to_active_cells.Insert(GRID_CELL_INDEX(grid_index,boundary_cell_indices(grid_index)(boundary_cell_index)),boundary_cell_active(grid_index)(boundary_cell_index));}*/

    //DEBUG OUTPUT
    //voronoi.laplace_grid.voronoi_faces=laplace_grid.voronoi_faces;
    //voronoi.voronoi_psi_N=voronoi_face_psi_N_velocities;
    
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("after compute boundary conditions"),2,2);

    n_matrix_cells=0;
    n_matrix_faces=0;
    
    //ALLOCATE MATRIX CELL INDICES FOR LOCAL GRIDS FIRST
    cell_indices_to_matrix_cell_indices.Resize(n_global_grids);
    boundary_cell_active.Resize(n_global_grids);
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
        if(laplace_grid.Local_Grid(grid_index)){
            cell_indices_to_matrix_cell_indices(grid_index).Resize(laplace_grid.Grid(grid_index).Domain_Indices());
            boundary_cell_active(grid_index).Remove_All();
            for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                if(Active_Cell_Compute(grid_index,cell_index)){
                    n_matrix_cells++;
                    cell_indices_to_matrix_cell_indices(grid_index)(cell_index)=n_matrix_cells;}
                else
                    cell_indices_to_matrix_cell_indices(grid_index)(cell_index)=0;
                if(laplace_grid.Boundary_Cell(grid_index,cell_index))
                    boundary_cell_active(grid_index).Append(cell_indices_to_matrix_cell_indices(grid_index)(cell_index)==0?false:true);}}
        else
            boundary_cell_active(grid_index).Resize(laplace_grid.boundary_cell_indices(grid_index).Size());}
    
    n_local_matrix_cells=n_matrix_cells;

    LOG::cout << "n_matrix_cells after local grids " << n_matrix_cells << std::endl;

    chimera_grid.Exchange_Boundary_Scalar_Values(boundary_cell_active,laplace_grid.is_boundary_grid);

    //ALLOCATE MATRIX CELL INDICES FOR BOUNDARY GRIDS SECOND
    boundary_cell_indices_to_matrix_cell_indices.Resize(n_global_grids);
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Boundary_Grid(grid_index)){
            boundary_cell_indices_to_matrix_cell_indices(grid_index).Remove_All();
            for(int boundary_cell_index=1;boundary_cell_index<=laplace_grid.boundary_cell_indices(grid_index).Size();boundary_cell_index++){
                TV_INT& cell_index=laplace_grid.boundary_cell_indices(grid_index)(boundary_cell_index);
                if(boundary_cell_active(grid_index)(boundary_cell_index) && laplace_grid.Boundary_Cell_Incident_To_Local_Grid(grid_index,cell_index)){
                    n_matrix_cells++;
                    //LOG::cout << "allocating matrix cell index " << grid_index << " " << cell_index << " " << n_matrix_cells << std::endl;
                    boundary_cell_indices_to_matrix_cell_indices(grid_index).Insert(cell_index,n_matrix_cells);}
                else
                    boundary_cell_indices_to_matrix_cell_indices(grid_index).Insert(cell_index,0);}}
    
    LOG::cout << "n_matrix_cells after boundary grids " << n_matrix_cells << std::endl;

    face_indices_to_matrix_face_indices.Resize(n_global_grids);
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index)){
            face_indices_to_matrix_face_indices(grid_index).Resize(laplace_grid.Grid(grid_index).Domain_Indices(1));
            for(FACE_ITERATOR iterator(laplace_grid.Grid(grid_index),1);iterator.Valid();iterator.Next())
                face_indices_to_matrix_face_indices(grid_index)(iterator.Full_Index())=0;
            for(FACE_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                TV_INT first_cell_index=iterator.First_Cell_Index();
                TV_INT second_cell_index=iterator.Second_Cell_Index();
                D_FACE_INDEX face_index=iterator.Full_Index();
                if(laplace_grid.Chimera_Face(grid_index,face_index) && (Matrix_Cell_Index(grid_index,first_cell_index) || Matrix_Cell_Index(grid_index,second_cell_index)) && !PSI_N(laplace_grid.Local_Grid_Index(grid_index))(face_index)){
                    n_matrix_faces++;
                    face_indices_to_matrix_face_indices(grid_index)(face_index)=n_matrix_faces;}
                else
                    face_indices_to_matrix_face_indices(grid_index)(face_index)=0;}}
    
    voronoi_face_indices_to_matrix_face_indices.Resize(laplace_grid.voronoi_faces.Size(),true,false,0);
    for(int voronoi_face_index=1;voronoi_face_index<=laplace_grid.voronoi_faces.Size();voronoi_face_index++){
        GRID_CELL_INDEX& index_1=laplace_grid.voronoi_faces(voronoi_face_index).x(1);
        GRID_CELL_INDEX& index_2=laplace_grid.voronoi_faces(voronoi_face_index).x(2);
        if(((laplace_grid.Local_Grid(index_1.x) && Matrix_Cell_Index(index_1.x,index_1.y)) || (laplace_grid.Local_Grid(index_2.x) && Matrix_Cell_Index(index_2.x,index_2.y))) && !voronoi_face_psi_N_velocities(voronoi_face_index).x){
            n_matrix_faces++;
            voronoi_face_indices_to_matrix_face_indices(voronoi_face_index)=n_matrix_faces;}
        else
            voronoi_face_indices_to_matrix_face_indices(voronoi_face_index)=0;}
  
    // LOG::cout << "after compute matrix indices sizes " << laplace_grid.voronoi_faces.Size() << " " << voronoi_face_indices_to_matrix_face_indices.Size() << " " <<
    // voronoi_face_psi_N_velocities.Size() << " " << n_matrix_faces << " " << laplace_grid.n_local_voronoi_faces << " " << chimera_grid.rank << std::endl;

    //DEBUG OUTPUT
    //voronoi.cell_indices_to_matrix_cell_indices=cell_indices_to_matrix_cell_indices;
    //voronoi.face_indices_to_matrix_face_indices=face_indices_to_matrix_face_indices;
    //voronoi.voronoi_face_indices_to_matrix_face_indices=voronoi_face_indices_to_matrix_face_indices;
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("after compute matrix indices"),2,2);
}
//#####################################################################
// Function Construct_Matrix_Partitions
//#####################################################################
template<class T_GRID> void LAPLACE_CHIMERA_MPI<T_GRID>::Construct_Matrix_Partitions()
{
    LOG::SCOPE scope("Construct_Matrix_Partitions");

    if(!chimera_grid.use_mpi) return;

    partition.boundary_indices.Clean_Memory();

    partition.interior_indices.min_corner=1;partition.interior_indices.max_corner=n_local_matrix_cells;
    partition.boundary_indices.Resize(n_global_grids);partition.ghost_indices.Resize(n_global_grids);
    partition.neighbor_ranks.Resize(n_global_grids);partition.neighbor_ranks.Fill(-2);

    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index))
            for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                if(laplace_grid.Boundary_Cell(grid_index,cell_index)){
                    int matrix_cell_index=Matrix_Cell_Index(grid_index,cell_index);
                    if(matrix_cell_index){
                        GRID_CELL_INDEX grid_cell_index(grid_index,cell_index);
                        ARRAY<int>& incident_voronoi_faces=laplace_grid.cell_indices_to_incident_voronoi_face_indices.Get(grid_cell_index);
                        HASHTABLE<int> added_destination;//added_destination.Clean_Memory();
                        for(int incident_voronoi_face_index=1;incident_voronoi_face_index<=incident_voronoi_faces.Size();incident_voronoi_face_index++){
                            VORONOI_FACE_INDICES& voronoi_face_indices=laplace_grid.voronoi_faces(incident_voronoi_faces(incident_voronoi_face_index)).x;
                            GRID_CELL_INDEX& destination_grid_cell_index=(voronoi_face_indices(1)==grid_cell_index)?voronoi_face_indices(2):voronoi_face_indices(1);
                            if(laplace_grid.Boundary_Grid(destination_grid_cell_index.x) && !added_destination.Contains(destination_grid_cell_index.x)){
                                added_destination.Set(destination_grid_cell_index.x);
                                partition.boundary_indices(destination_grid_cell_index.x).Append(matrix_cell_index);}}}}}
    
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
        if(laplace_grid.Boundary_Grid(grid_index)){
            int first_matrix_cell_index=0;
            int last_matrix_cell_index=-1;
            for(int boundary_cell_index=1;boundary_cell_index<=laplace_grid.boundary_cell_indices(grid_index).Size();boundary_cell_index++){
                int matrix_cell_index=Matrix_Cell_Index(grid_index,laplace_grid.boundary_cell_indices(grid_index)(boundary_cell_index));
                if(matrix_cell_index){
                    if(!first_matrix_cell_index)
                        first_matrix_cell_index=matrix_cell_index;
                    last_matrix_cell_index=matrix_cell_index;}}
            partition.ghost_indices(grid_index).min_corner=first_matrix_cell_index;
            partition.ghost_indices(grid_index).max_corner=last_matrix_cell_index;
            partition.neighbor_ranks(grid_index)=chimera_grid.grid_ranks(grid_index);}}
    partition.number_of_sides=n_global_grids;
    
    LOG::cout << "after construct matrix partitions sizes " << laplace_grid.voronoi_faces.Size() << " " << voronoi_face_indices_to_matrix_face_indices.Size() << " " << voronoi_face_psi_N_velocities.Size() << std::endl;

////////////////////////////////////////////////DEBUG
    /*LOG::cout<<"lqiu debug rank:"<<chimera_grid.rank<<",interior indices:"<<partition.interior_indices.min_corner<<","<<partition.interior_indices.max_corner<<std::endl;
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
        LOG::cout<<"lqiu debug grid_index"<<grid_index<<",rank:"<<chimera_grid.rank<<", boundary indices size:"<<partition.boundary_indices(grid_index).Size()<<",ghost region:"<<partition.ghost_indices(grid_index).min_corner<<","<<partition.ghost_indices(grid_index).max_corner<<",boundary cell indices size:"<<laplace_grid.boundary_cell_indices(grid_index).Size()<<std::endl;
    }*/
////////////////////////////////////////////////////
}
//#####################################################################
// Function Construct_Laplacian
//#####################################################################
template<class T_GRID> void LAPLACE_CHIMERA_MPI<T_GRID>::Construct_Laplacian(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks)
{
    Construct_Divergence();
    Construct_Inverse_Mass(callbacks);
    
    SPARSE_MATRIX_FLAT_MXN<T> divergence_matrix_transpose;divergence_matrix.Transpose(divergence_matrix_transpose);

    //LOG::cout << "dual_celL_inverse_mass " << dual_cell_inverse_mass << std::endl;
    
    //should store cell sizes and then compute dual cell inverse mass on the fly here

    system_matrix=divergence_matrix.Times_Diagonal_Times(dual_cell_inverse_mass,divergence_matrix_transpose).Create_NXN_Matrix();
}
//#####################################################################
// Function Build_Divergence_Matrix
//#####################################################################
template<class T_GRID> void LAPLACE_CHIMERA_MPI<T_GRID>::Construct_Divergence()
{
    LOG::SCOPE scope("Build_Divergence_Matrix");

    LOG::cout << "n_matrix_cells " << n_matrix_cells << std::endl;

    //CALCULATE DIVERGENCE MATRIX ROW LENGTHS
    ARRAY<int> divergence_matrix_row_lengths(n_matrix_cells);
    ARRAYS_COMPUTATIONS::Fill(divergence_matrix_row_lengths,0);
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index))
            for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                int matrix_cell_index=Matrix_Cell_Index(grid_index,cell_index);
                if(matrix_cell_index){
                    int n_cell_faces=0;
                    for(int axis=1;axis<=TV::dimension;axis++)
                        for(int side=1;side<=2;side++)
                            if(Matrix_Face_Index(grid_index,D_FACE_INDEX(axis,cell_index+(side-1)*TV_INT::Axis_Vector(axis))))
                                n_cell_faces++;
                    divergence_matrix_row_lengths(matrix_cell_index)=n_cell_faces;}}
    for(int voronoi_face_index=1;voronoi_face_index<=laplace_grid.voronoi_faces.Size();voronoi_face_index++){
        int matrix_face_index=Matrix_Face_Index(voronoi_face_index);
        if(matrix_face_index){
            for(int side=1;side<=2;side++){
                GRID_CELL_INDEX& grid_cell_index=laplace_grid.voronoi_faces(voronoi_face_index).x(side);
                int matrix_cell_index=Matrix_Cell_Index(grid_cell_index.x,grid_cell_index.y);
                if(matrix_cell_index)
                    divergence_matrix_row_lengths(matrix_cell_index)+=1;}}}
    
    //BUILD DIVERGENCE MATRIX
    divergence_matrix.Reset(0);
    divergence_matrix.Set_Row_Lengths(divergence_matrix_row_lengths);
    divergence_matrix.n=n_matrix_faces;
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index)){
            TV face_sizes=laplace_grid.Grid(grid_index).Face_Sizes();
            for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                int matrix_cell_index=Matrix_Cell_Index(grid_index,cell_index);
                if(matrix_cell_index){
                    for(int axis=1;axis<=TV::dimension;axis++)
                        for(int side=1;side<=2;side++){
                            int matrix_face_index=face_indices_to_matrix_face_indices(grid_index).Component(axis)(cell_index+(side-1)*TV_INT::Axis_Vector(axis));
                            if(matrix_face_index)
                                divergence_matrix.Set_Element(matrix_cell_index,matrix_face_index,(2*side-3)*face_sizes(axis));}}}}
    for(int voronoi_face_index=1;voronoi_face_index<=laplace_grid.voronoi_faces.Size();voronoi_face_index++){
        int matrix_face_index=Matrix_Face_Index(voronoi_face_index);
        if(matrix_face_index){
            T face_size=laplace_grid.voronoi_faces(voronoi_face_index).y;
            for(int side=1;side<=2;side++){
                GRID_CELL_INDEX& grid_cell_index=laplace_grid.voronoi_faces(voronoi_face_index).x(side);
                int matrix_cell_index=Matrix_Cell_Index(grid_cell_index.x,grid_cell_index.y);
                if(matrix_cell_index)
                    divergence_matrix.Set_Element(matrix_cell_index,matrix_face_index,(3-2*side)*face_size);}}}
    
    /////////////////////////////DEBUG
    /*t cell=29;
    for(int i=divergence_matrix.offsets(cell);i<divergence_matrix.offsets(cell+1);i++)
    LOG::cout << "cell-face coefficient " << cell << " " << divergence_matrix.A(i).j << " " << divergence_matrix.A(i).a << std::endl;*/
}
//#####################################################################
// Function Inverse_Mass
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LAPLACE_CHIMERA_MPI<T_GRID>::Inverse_Mass(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const int grid_index,const D_FACE_INDEX& face_index)
{
    TV location=laplace_grid.Frame(grid_index)*laplace_grid.Grid(grid_index).Face(face_index.axis,face_index.index);
    T density=callbacks.Get_Density(location);
    T size=laplace_grid.Face_Size(grid_index,face_index);
    T mass=size*density*Dual_Cell_Volume_Fraction(grid_index,face_index.First_Cell_Index(),grid_index,face_index.Second_Cell_Index());
    return (T)1.0/mass;
}
template<class T_GRID> typename T_GRID::SCALAR LAPLACE_CHIMERA_MPI<T_GRID>::Inverse_Mass(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const int voronoi_face_index)
{
    GRID_CELL_INDEX& grid_cell_index_1=laplace_grid.voronoi_faces(voronoi_face_index).x(1);
    GRID_CELL_INDEX& grid_cell_index_2=laplace_grid.voronoi_faces(voronoi_face_index).x(2);
    TV location_1=laplace_grid.Frame(grid_cell_index_1.x)*laplace_grid.Grid(grid_cell_index_1.x).X(grid_cell_index_1.y);
    TV location_2=laplace_grid.Frame(grid_cell_index_2.x)*laplace_grid.Grid(grid_cell_index_2.x).X(grid_cell_index_2.y);
    TV location=(location_1+location_2)/2;
    T density=callbacks.Get_Density(location);
    T size=laplace_grid.Face_Size(voronoi_face_index);
    T mass=size*density*Dual_Cell_Volume_Fraction(grid_cell_index_1.x,grid_cell_index_1.y,grid_cell_index_2.x,grid_cell_index_2.y);
    return (T)1.0/mass;
}
//#####################################################################
// Function Construct_Inverse_mass
//#####################################################################
template<class T_GRID> void LAPLACE_CHIMERA_MPI<T_GRID>::Construct_Inverse_Mass(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks)
{
    LOG::SCOPE scope("Build_Inverse_Mass_Matrix");
    
    dual_cell_inverse_mass.Resize(n_matrix_faces);
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
      if(laplace_grid.Local_Grid(grid_index))
        for(FACE_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
            int matrix_face_index=face_indices_to_matrix_face_indices(grid_index)(iterator.Full_Index());
            if(matrix_face_index)
                dual_cell_inverse_mass(matrix_face_index)=Inverse_Mass(callbacks,grid_index,iterator.Full_Index());}
    for(int voronoi_face_index=1;voronoi_face_index<=laplace_grid.voronoi_faces.Size();voronoi_face_index++){
        int matrix_face_index=voronoi_face_indices_to_matrix_face_indices(voronoi_face_index);
        if(matrix_face_index)
            dual_cell_inverse_mass(matrix_face_index)=Inverse_Mass(callbacks,voronoi_face_index);}
}
//#####################################################################
// Function Compute_Neumann_Divergence
//#####################################################################
template<class T_GRID> void LAPLACE_CHIMERA_MPI<T_GRID>::Compute_Neumann_Divergence(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time,VECTOR_ND<T>& divergence)
{
    divergence.Fill(0);
    
    T normal_gradient;
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index)){
            TV face_sizes=laplace_grid.Grid(grid_index).Face_Sizes();
            for(FACE_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next())
                if(laplace_grid.Chimera_Face(grid_index,iterator.Full_Index())){
                    if(callbacks.Get_Neumann_Boundary_Condition(laplace_grid.Frame(grid_index)*iterator.Location(),laplace_grid.Frame(grid_index).r.Rotate(TV::Axis_Vector(iterator.Axis())),normal_gradient,time)){
                        //LOG::cout << "neumann bc " << grid_index << " " << iterator.Full_Index() << " " << normal_gradient << std::endl;
                        //normal_gradient*=Inverse_Mass(callbacks,grid_index,iterator.Full_Index());
                        for(int side=1;side<=2;side++){
                            TV_INT cell_index=iterator.Face_Index()+(side-2)*TV_INT::Axis_Vector(iterator.Axis());
                            int matrix_cell_index=Matrix_Cell_Index(grid_index,cell_index);
                            if(matrix_cell_index)
                                divergence(matrix_cell_index)+=(3-2*side)*face_sizes(iterator.Axis())*normal_gradient;}}}}
    
    for(int voronoi_face_index=1;voronoi_face_index<=laplace_grid.voronoi_faces.Size();voronoi_face_index++){
        GRID_CELL_INDEX& index_1=laplace_grid.voronoi_faces(voronoi_face_index).x(1);
        GRID_CELL_INDEX& index_2=laplace_grid.voronoi_faces(voronoi_face_index).x(2);
        TV location_1=laplace_grid.Frame(index_1.x)*laplace_grid.Grid(index_1.x).X(index_1.y);
        TV location_2=laplace_grid.Frame(index_2.x)*laplace_grid.Grid(index_2.x).X(index_2.y);
        TV normal=(location_2-location_1).Normalized();
        TV location=(location_2+location_1)*(T)0.5;
        if(callbacks.Get_Neumann_Boundary_Condition(location,normal,normal_gradient,time)){
            //LOG::cout << "neumann bc " << voronoi_face_index << " " << normal_gradient << std::endl;
            //normal_gradient*=Inverse_Mass(callbacks,voronoi_face_index);
            T face_size=laplace_grid.voronoi_faces(voronoi_face_index).y;
            for(int side=1;side<=2;side++){
                GRID_CELL_INDEX& grid_cell_index=laplace_grid.voronoi_faces(voronoi_face_index).x(side);
                int matrix_cell_index=Matrix_Cell_Index(grid_cell_index.x,grid_cell_index.y);
                if(matrix_cell_index)
                    divergence(matrix_cell_index)+=(T)(3-2*side)*face_size*normal_gradient;}}}
}
//#####################################################################
// Function Compute_Dirichlet_Divergence
//#####################################################################
template<class T_GRID> void LAPLACE_CHIMERA_MPI<T_GRID>::Compute_Dirichlet_Divergence(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time,VECTOR_ND<T>& divergence)
{
    VECTOR_ND<T> face_gradients(n_matrix_faces);
    
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
        if(laplace_grid.Local_Grid(grid_index)){
            TV face_sizes=laplace_grid.Grid(grid_index).Face_Sizes();
            for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index),1);iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                if(laplace_grid.Chimera_Cell(grid_index,cell_index)){
                    TV location=laplace_grid.Frame(grid_index)*laplace_grid.Grid(grid_index).X(cell_index);
                    T p;
                    if(callbacks.Get_Dirichlet_Boundary_Condition(location,p,time)){
                        for(int axis=1;axis<=TV::dimension;axis++)
                            for(int side=1;side<=2;side++){
                                TV_INT face_index=cell_index+(side-1)*TV_INT::Axis_Vector(axis);
                                int matrix_face_index=Matrix_Face_Index(grid_index,D_FACE_INDEX(axis,face_index));
                                if(matrix_face_index) face_gradients(matrix_face_index)+=(3-2*side)*dual_cell_inverse_mass(matrix_face_index)*face_sizes(axis)*p;}
                        GRID_CELL_INDEX grid_cell_index(grid_index,cell_index);
                        if(laplace_grid.cell_indices_to_incident_voronoi_face_indices.Contains(grid_cell_index)){
                            const ARRAY<int>& incident_voronoi_face_indices=laplace_grid.cell_indices_to_incident_voronoi_face_indices.Get(grid_cell_index);
                            for(int incident_index=1;incident_index<=incident_voronoi_face_indices.Size();incident_index++){
                                const VORONOI_FACE_INDICES& indices=laplace_grid.voronoi_faces(incident_voronoi_face_indices(incident_index)).x;
                                int side=indices(1).x==grid_index;
                                int matrix_face_index=Matrix_Face_Index(incident_voronoi_face_indices(incident_index));
                                if(matrix_face_index) face_gradients(matrix_face_index)+=(3-2*side)*dual_cell_inverse_mass(matrix_face_index)*laplace_grid.voronoi_faces(incident_voronoi_face_indices(incident_index)).y*p;}}}}}}
        else if(laplace_grid.Boundary_Grid(grid_index))
            for(int boundary_cell_index=1;boundary_cell_index<=laplace_grid.boundary_cell_indices(grid_index).Size();boundary_cell_index++){
                TV_INT cell_index=laplace_grid.boundary_cell_indices(grid_index)(boundary_cell_index);
                TV location=laplace_grid.Frame(grid_index)*laplace_grid.Grid(grid_index).X(cell_index);
                T p;
                if(callbacks.Get_Dirichlet_Boundary_Condition(location,p,time)){
                    GRID_CELL_INDEX grid_cell_index(grid_index,cell_index);
                    if(laplace_grid.cell_indices_to_incident_voronoi_face_indices.Contains(grid_cell_index)){
                        const ARRAY<int>& incident_voronoi_face_indices=laplace_grid.cell_indices_to_incident_voronoi_face_indices.Get(grid_cell_index);
                        for(int incident_index=1;incident_index<=incident_voronoi_face_indices.Size();incident_index++){
                            const VORONOI_FACE_INDICES& indices=laplace_grid.voronoi_faces(incident_voronoi_face_indices(incident_index)).x;
                            int side=indices(1).x==grid_index;
                            int matrix_face_index=Matrix_Face_Index(incident_voronoi_face_indices(incident_index));
                            if(matrix_face_index) face_gradients(matrix_face_index)+=(3-2*side)*dual_cell_inverse_mass(matrix_face_index)*laplace_grid.voronoi_faces(incident_voronoi_face_indices(incident_index)).y*p;}}}}}

    //LOG::cout << "face gradients " << face_gradients << std::endl;
    
    divergence_matrix.Times(face_gradients,divergence);
}
//#####################################################################
// Function Construct_System
//#####################################################################
template<class T_GRID> void LAPLACE_CHIMERA_MPI<T_GRID>::Construct_System(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time,const T dt,const bool diffusion)
{
    n_global_grids=chimera_grid.number_of_global_grids;
    n_local_grids=chimera_grid.number_of_local_grids;
    
    Construct_Matrix_Indices(callbacks,time);
    Construct_Matrix_Partitions();
    Construct_Laplacian(callbacks);
    
    if(diffusion){
        system_matrix*=dt;
        for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
            if(laplace_grid.Local_Grid(grid_index))
                for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                    TV_INT cell_index=iterator.Cell_Index();
                    int matrix_cell_index=Matrix_Cell_Index(grid_index,cell_index);
                    if(matrix_cell_index)
                        system_matrix.Add_Element(matrix_cell_index,matrix_cell_index,laplace_grid.Cell_Size(grid_index,cell_index));}}

    //LOG::cout << "matrix " << diffusion << " " << dt << std::endl << system_matrix << std::endl;
}
//#####################################################################
// Function Multiply_By_Volume_Weighted_Divergence
//#####################################################################
template<class T_GRID> void LAPLACE_CHIMERA_MPI<T_GRID>::Multiply_By_Volume_Weighted_Divergence(const VECTOR_ND<T>& x,VECTOR_ND<T>& result)
{
    divergence_matrix.Times(x,result);
}
//#####################################################################
// Function Multiply_By_Gradient
//#####################################################################
template<class T_GRID> void LAPLACE_CHIMERA_MPI<T_GRID>::Multiply_By_Gradient(const VECTOR_ND<T>& x,VECTOR_ND<T>& result)
{
    divergence_matrix.Transpose_Times(x,result);
    result.Negate();
}
//#####################################################################
// Function Multiply_By_Mass_Inverse
//#####################################################################
template<class T_GRID> void LAPLACE_CHIMERA_MPI<T_GRID>::Multiply_By_Inverse_Mass(const VECTOR_ND<T>& x,VECTOR_ND<T>& result)
{
    for(int i=1;i<=n_matrix_faces;i++)
        result(i)=x(i)*dual_cell_inverse_mass(i);
}
//#####################################################################
// Function Multiply_By_Cell_Size
//#####################################################################
template<class T_GRID> void LAPLACE_CHIMERA_MPI<T_GRID>::Multiply_By_Cell_Size(const VECTOR_ND<T>& x,VECTOR_ND<T>& result)
{
    /*for(int i=1;i<=n_cell_faces;i++)
        result(i)=x(i)*laplace_grid.Cell_Size((i);*/
}
//#####################################################################
// Function Multiply_By_Cell_Size
//#####################################################################
template<class T_GRID> void LAPLACE_CHIMERA_MPI<T_GRID>::Multiply_By_Dual_Cell_Size(const VECTOR_ND<T>& x,VECTOR_ND<T>& result)
{
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index))
            for(FACE_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                int matrix_face_index=Matrix_Face_Index(grid_index,iterator.Full_Index());
                if(matrix_face_index)
                    result(matrix_face_index)=x(matrix_face_index)*laplace_grid.Face_Size(grid_index,iterator.Full_Index());}
    
    for(int voronoi_face_index=1;voronoi_face_index<=laplace_grid.voronoi_faces.Size();voronoi_face_index++){
        int matrix_face_index=Matrix_Face_Index(voronoi_face_index);
        if(matrix_face_index)
            result(matrix_face_index)=x(matrix_face_index)*laplace_grid.Face_Size(voronoi_face_index);}
}
//#####################################################################

#define INSTANTIATION_HELPER(T,D)                               \
    template class LAPLACE_CHIMERA_MPI<GRID<VECTOR<T,D> > >;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif
