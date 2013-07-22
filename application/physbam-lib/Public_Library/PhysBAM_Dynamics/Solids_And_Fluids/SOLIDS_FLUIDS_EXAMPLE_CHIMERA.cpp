//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_BINARY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Fluids/Coupled_Evolution/COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/FLUID_COLLISION_BODY_INACCURATE_UNION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
#include <PhysBAM_Dynamics/Forces_And_Torques/EULER_FLUID_FORCES.h>
#include <PhysBAM_Dynamics/Incompressible_Flows/INCOMPRESSIBLE_MULTIPHASE_UNIFORM.h>
#include <PhysBAM_Dynamics/Interpolation/FACE_LOOKUP_FIRE_MULTIPHASE_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_MATRIX_1X1.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_1D.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_2D.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_3D.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_ND.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_HASHTABLE.h>
#include <PhysBAM_Geometry/Read_Write/Topology/READ_WRITE_SIMPLEX_MESH.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/Find_Type.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_1D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_CHIMERA.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_SEGMENT_2D.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_POLYGON.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
SOLIDS_FLUIDS_EXAMPLE_CHIMERA(const STREAM_TYPE stream_type,const int number_of_regions,const typename FLUIDS_PARAMETERS<T_GRID>::TYPE type,const int array_collection_type)
    :SOLIDS_FLUIDS_EXAMPLE<TV>(stream_type,array_collection_type),rigid_grid_collection(*new RIGID_GRID_COLLECTION<T_GRID>()),fluids_parameters(number_of_regions,type,*new INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>()),cfl_alpha(2.),cfl_beta(1.),boundary_chimera(0),boundary_phi_water(0),resolution(0),simulate_solids(false),kinematic_grid(false),run_advection_test(false),solve_poisson_equation(false),solve_heat_equation(false),solve_heat_equation_on_faces(false),solve_heat_equation_on_faces_coupled(false),use_flip_face_update(true),use_flip_face_update_consistent(true),use_flip_cell_update(false),use_flip_update_overlap(false),enforce_second_order_maccormack_scalar(false),enforce_second_order_maccormack_vector(false),clamp_extrema_maccormack_scalar(true),clamp_extrema_maccormack_vector(true),fill_overlap_region_recursively(true),use_maccormack_advection_scalar(false),use_maccormack_advection_vector(false),use_collidable_advection(false),revalidate_velocities_in_projection(true),define_slip_boundary(false)
{
    fluids_parameters.implicit_viscosity=true;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
~SOLIDS_FLUIDS_EXAMPLE_CHIMERA()
{
    if(dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(solids_evolution)) fluids_parameters.projection=0;
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Register_Options()
{
    BASE::Register_Options();
    this->parse_args->Add_Option_Argument("-skip_debug_data","turn off file io for debug data");
    this->parse_args->Add_Option_Argument("-use_fmm_extrapolation","use fast marching to extrapolate compressible flow data into solid state.");
    this->parse_args->Add_Integer_Argument("-resolution",1);
    this->parse_args->Add_Double_Argument("-vc",0.3,"vorticity confinement");
    this->parse_args->Add_Option_Argument("-write_voronoi");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Parse_Options()
{
    BASE::Parse_Options();
    if(this->parse_args->Is_Value_Set("-skip_debug_data")) fluids_parameters.write_debug_data=false;
    if(this->parse_args->Is_Value_Set("-use_fmm_extrapolation")) fluids_parameters.euler_solid_fluid_coupling_utilities->use_fast_marching=true;
    if(this->parse_args->Is_Value_Set("-vc")){fluids_parameters.use_vorticity_confinement=true;fluids_parameters.confinement_parameter=(T)(this->parse_args->Get_Double_Value("-vc"));}
    else fluids_parameters.use_vorticity_confinement=false;
    resolution=this->parse_args->Get_Integer_Value("-resolution");
    write_voronoi=this->parse_args->Get_Option_Value("-write_voronoi");
}
//#####################################################################
// Function Advect_Scalar_Field
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Advect_Scalar_Field(ARRAY<T_ARRAYS_SCALAR*>& u_array,ARRAY<T_ARRAYS_SCALAR*>& u_ghost_array,ARRAY<T_FACE_ARRAYS_SCALAR*>& face_velocities_ghost,const T dt,const T time)
{
    if(!use_maccormack_advection_scalar) for(int grid_index=1;grid_index<=rigid_grid_collection.particles.array_collection->Size();grid_index++){
        advection(grid_index)->Update_Advection_Equation_Cell(rigid_grid_collection.Rigid_Grid(grid_index).grid,*u_array(grid_index),*u_ghost_array(grid_index),*face_velocities_ghost(grid_index),*boundary_chimera,dt,time);}
    else{
        advection_maccormack_scalar->Update_Advection_Equation_Cell_Chimera(u_array,u_ghost_array,face_velocities_ghost,*boundary_chimera,dt,time);}
}
//#####################################################################
// Function Advect_Velocities
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Advect_Vector_Field(ARRAY<T_FACE_ARRAYS_SCALAR*>& u_array,ARRAY<T_FACE_ARRAYS_SCALAR*>& u_ghost_array,ARRAY<T_FACE_ARRAYS_SCALAR*>& face_velocities_ghost,const T dt,const T time)
{
    if(!use_maccormack_advection_vector) for(int grid_index=1;grid_index<=rigid_grid_collection.particles.array_collection->Size();grid_index++){
        advection(grid_index)->Update_Advection_Equation_Face(rigid_grid_collection.Rigid_Grid(grid_index).grid,*u_array(grid_index),*u_ghost_array(grid_index),*face_velocities_ghost(grid_index),*boundary_chimera,dt,time);}
    else{
        advection_maccormack_vector->Update_Advection_Equation_Face_Chimera(u_array,u_ghost_array,face_velocities_ghost,*boundary_chimera,dt,time);}
}
//#####################################################################
// Function Extrapolate_Scalar_Fields_Into_Objects
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Extrapolate_Scalar_Field_Into_Objects(ARRAY<T_ARRAYS_SCALAR*>& u_array,const T dt,const T time)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    for(int i=1;i<=rigid_grid_collection.particles.array_collection->Size();i++){
        //int grid_index=chimera_grid->local_to_global_grid_index_map(i);
        T_GRID& grid=rigid_grid_collection.Rigid_Grid(i).grid;
        T_ARRAYS_SCALAR phi_object(grid.Domain_Indices(fluids_parameters.number_of_ghost_cells));
        for(CELL_ITERATOR iterator(grid,fluids_parameters.number_of_ghost_cells);iterator.Valid();iterator.Next()){
            T phi=std::numeric_limits<T>::infinity();
            TV location=incompressible_fluid_containers(i)->rigid_grid.Frame()*iterator.Location();
            for(int id=1;id<=rigid_body_collection.rigid_body_particle.array_collection->Size();id++)
                phi=min(phi,rigid_body_collection.Rigid_Body(id).Implicit_Geometry_Extended_Value(location));
            phi_object(iterator.Cell_Index())=-phi;}
        
        EXTRAPOLATION_UNIFORM<T_GRID,T> extrapolate(grid,phi_object,*u_array(i),fluids_parameters.number_of_ghost_cells);
        extrapolate.Set_Band_Width((T)fluids_parameters.number_of_ghost_cells);
        extrapolate.Extrapolate(0,false);}
}
//#####################################################################
// Function Fill_Ghost_Cells_Chimera
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Fill_Ghost_Cells_Chimera(int number_of_ghost_cells,ARRAY<T_ARRAYS_SCALAR*>& u_array,ARRAY<T_ARRAYS_SCALAR*>& u_ghost_array,bool only_overwrite_nan)
{
    for(int grid_index=1;grid_index<=rigid_grid_collection.particles.array_collection->Size();grid_index++){
        int unsplit_grid_index=chimera_grid->use_mpi?chimera_grid->mpi_grid_indices(chimera_grid->local_to_global_grid_index_map(grid_index)):grid_index;
        if(only_overwrite_nan && unsplit_grid_index!=1) continue;
        T_GRID& grid=rigid_grid_collection.Rigid_Grid(grid_index).grid;
        boundary_chimera->boundary.Fill_Ghost_Cells(grid,*u_array(grid_index),*u_ghost_array(grid_index),0,0,number_of_ghost_cells);}
    if(chimera_grid->mpi_split_grids.m) for(int grid_index=1;grid_index<=rigid_grid_collection.particles.array_collection->Size();grid_index++){
        int mpi_grid_index=chimera_grid->mpi_grid_indices(chimera_grid->local_to_global_grid_index_map(grid_index));
        MPI_UNIFORM_GRID<T_GRID>* mpi_grid=chimera_grid->mpi_split_grids(mpi_grid_index);
        mpi_grid->Exchange_Boundary_Cell_Data(*u_ghost_array(grid_index),number_of_ghost_cells,true);}
    chimera_grid->Exchange_Boundary_Cell_Data(u_ghost_array,number_of_ghost_cells,only_overwrite_nan);
}
//#####################################################################
// Function Fill_Ghost_Cells_Face_Chimera
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Fill_Ghost_Cells_Face_Chimera(int number_of_ghost_cells,ARRAY<T_FACE_ARRAYS_SCALAR*>& u_array,ARRAY<T_FACE_ARRAYS_SCALAR*>& u_ghost_array)
{
    for(int grid_index=1;grid_index<=rigid_grid_collection.particles.array_collection->Size();grid_index++){
        T_GRID& grid=rigid_grid_collection.Rigid_Grid(grid_index).grid;
        boundary_chimera->boundary.Fill_Ghost_Cells_Face(grid,*u_array(grid_index),*u_ghost_array(grid_index),0,number_of_ghost_cells);}
    if(chimera_grid->mpi_split_grids.m) for(int grid_index=1;grid_index<=rigid_grid_collection.particles.array_collection->Size();grid_index++){
        int mpi_grid_index=chimera_grid->mpi_grid_indices(chimera_grid->local_to_global_grid_index_map(grid_index));
        MPI_UNIFORM_GRID<T_GRID>* mpi_grid=chimera_grid->mpi_split_grids(mpi_grid_index);
        mpi_grid->Exchange_Boundary_Face_Data(*u_ghost_array(grid_index),number_of_ghost_cells);}
    chimera_grid->Exchange_Boundary_Face_Data(u_ghost_array,number_of_ghost_cells);
}
//#####################################################################
// Function Coupling_Overlap_Regions_Cell
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::Exchange_Split_Grid_Ghost_Cells(ARRAY<T_ARRAYS_SCALAR*>& u_array)
{
    if(chimera_grid->mpi_split_grids.m) for(int grid_index=1;grid_index<=rigid_grid_collection.particles.array_collection->Size();grid_index++){
        int mpi_grid_index=chimera_grid->mpi_grid_indices(chimera_grid->local_to_global_grid_index_map(grid_index));
        MPI_UNIFORM_GRID<T_GRID>* mpi_grid=chimera_grid->mpi_split_grids(mpi_grid_index);
        mpi_grid->Exchange_Boundary_Cell_Data(*u_array(grid_index),1,true);}
}
//#####################################################################
// Function Coupling_Overlap_Regions_Cell
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Coupling_Overlap_Regions_Cell(int number_of_ghost_cells,ARRAY<T_ARRAYS_SCALAR*>& u_array,bool only_overwrite_nan)
{
    Exchange_Split_Grid_Ghost_Cells(u_array);

    if(fill_overlap_region_recursively){
        chimera_grid->Exchange_Overlap_Cell_Data_Recursively(u_array,0,only_overwrite_nan);}
    else{
        chimera_grid->Exchange_Overlap_Cell_Data(u_array,0,only_overwrite_nan);}

    Fill_Ghost_Cells_Chimera(number_of_ghost_cells,u_array,u_array,only_overwrite_nan);
}
//#####################################################################
// Function Coupling_Overlap_Regions_Face
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Coupling_Overlap_Regions_Face(int number_of_ghost_cells,ARRAY<T_FACE_ARRAYS_SCALAR*>& u_array)
{
    if(chimera_grid->mpi_split_grids.m) for(int grid_index=1;grid_index<=rigid_grid_collection.particles.array_collection->Size();grid_index++){
        int mpi_grid_index=chimera_grid->mpi_grid_indices(chimera_grid->local_to_global_grid_index_map(grid_index));
        MPI_UNIFORM_GRID<T_GRID>* mpi_grid=chimera_grid->mpi_split_grids(mpi_grid_index);
        mpi_grid->Exchange_Boundary_Face_Data(*u_array(grid_index),1);}

    if(fill_overlap_region_recursively){
        chimera_grid->Exchange_Overlap_Face_Data_Recursively(u_array,0);}
    else{
        chimera_grid->Exchange_Overlap_Face_Data(u_array,0);}

    Fill_Ghost_Cells_Face_Chimera(number_of_ghost_cells,u_array,u_array);
}
//#####################################################################
// Function Add_Volumetric_Body_To_Fluid_Simulation
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Add_Volumetric_Body_To_Fluid_Simulation(int grid_index,RIGID_BODY<TV>& rigid_body,bool add_collision,bool add_coupling)
{
    RIGID_COLLISION_GEOMETRY<TV>* collision_geometry=new RIGID_COLLISION_GEOMETRY<TV>(rigid_body);
    if(add_collision){
        incompressible_fluid_containers(grid_index)->collision_bodies_affecting_fluid.collision_geometry_collection.Add_Body(collision_geometry,rigid_body.particle_index,true);}
    //if(add_coupling)
    //  if(SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(solids_evolution))
    //      coupled_evolution->iterator_info.coupling_bodies.Append(collision_geometry);
    if(rigid_body.simplicial_object) rigid_body.simplicial_object->Initialize_Hierarchy();
}
//#####################################################################
// Function Add_Thin_Shell_To_Fluid_Simulation
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Add_Thin_Shell_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,bool add_collision,bool add_coupling)
{
    RIGID_COLLISION_GEOMETRY<TV>* collision_geometry=new RIGID_COLLISION_GEOMETRY<TV>(rigid_body);
    if(add_collision) fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(collision_geometry,rigid_body.particle_index,true);
    if(add_coupling)
        if(SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(solids_evolution))
            coupled_evolution->iterator_info.coupling_bodies.Append(collision_geometry);
    rigid_body.thin_shell=true;
    if(collision_geometry->collision_thickness<minimum_collision_thickness) collision_geometry->Set_Collision_Thickness(minimum_collision_thickness);
    if(rigid_body.simplicial_object) rigid_body.simplicial_object->Initialize_Hierarchy();
}
//#####################################################################
// Function Add_To_Fluid_Simulation
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Add_To_Fluid_Simulation(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions,bool add_collision,bool add_coupling)
{
    if(add_collision) fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(&deformable_collisions,0,false);
    if(add_coupling)
        if(SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(solids_evolution))
            coupled_evolution->iterator_info.coupling_bodies.Append(&deformable_collisions);
    if(deformable_collisions.collision_thickness<minimum_collision_thickness) deformable_collisions.Set_Collision_Thickness(minimum_collision_thickness);
    deformable_collisions.Initialize_For_Thin_Shells_Fluid_Coupling();
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Initialize_Phi()
{
}
//#####################################################################
// Function Initialize_MPI
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Initialize_MPI()
{
    fluids_parameters.mpi_grid->Initialize(fluids_parameters.domain_walls);
    if(fluids_parameters.number_of_regions==1)fluids_parameters.phi_boundary=new BOUNDARY_MPI<T_GRID>(fluids_parameters.mpi_grid,*fluids_parameters.phi_boundary);
    else if(fluids_parameters.number_of_regions>=2)for(int i=1;i<=fluids_parameters.number_of_regions;i++)
        fluids_parameters.phi_boundary_multiphase(i)=new BOUNDARY_MPI<T_GRID>(fluids_parameters.mpi_grid,*fluids_parameters.phi_boundary_multiphase(i));
    fluids_parameters.fluid_boundary=new BOUNDARY_MPI<T_GRID>(fluids_parameters.mpi_grid,*fluids_parameters.fluid_boundary);
    if(fluids_parameters.use_soot){
        if(fluids_parameters.soot_boundary) fluids_parameters.soot_boundary=new BOUNDARY_MPI<T_GRID>(fluids_parameters.mpi_grid,*fluids_parameters.soot_boundary);
        fluids_parameters.soot_container.Set_Custom_Boundary(*new BOUNDARY_MPI<T_GRID>(fluids_parameters.mpi_grid,*fluids_parameters.soot_container.boundary));
        fluids_parameters.soot_fuel_container.Set_Custom_Boundary(*new BOUNDARY_MPI<T_GRID>(fluids_parameters.mpi_grid,*fluids_parameters.soot_fuel_container.boundary));}
    if(fluids_parameters.use_density){
        if(fluids_parameters.density_boundary) fluids_parameters.density_boundary=new BOUNDARY_MPI<T_GRID>(fluids_parameters.mpi_grid,*fluids_parameters.density_boundary);
        fluids_parameters.density_container.Set_Custom_Boundary(*new BOUNDARY_MPI<T_GRID>(fluids_parameters.mpi_grid,*fluids_parameters.density_container.boundary));}
    if(fluids_parameters.use_temperature){
        if(fluids_parameters.temperature_boundary) fluids_parameters.temperature_boundary=new BOUNDARY_MPI<T_GRID>(fluids_parameters.mpi_grid,*fluids_parameters.temperature_boundary);
        fluids_parameters.temperature_container.Set_Custom_Boundary(*new BOUNDARY_MPI<T_GRID>(fluids_parameters.mpi_grid,*fluids_parameters.temperature_container.boundary));}
    for(int i=1;i<=fluids_parameters.number_of_regions;i++) fluids_parameters.particle_levelset_evolution->Particle_Levelset(i).mpi_grid=fluids_parameters.mpi_grid;
    if(fluids_parameters.use_strain||fluids_parameters.use_multiphase_strain.Count_Matches(0)<fluids_parameters.number_of_regions)
        fluids_parameters.strain_boundary=new BOUNDARY_MPI<T_GRID,SYMMETRIC_MATRIX>(fluids_parameters.mpi_grid,*fluids_parameters.strain_boundary);
    if(fluids_parameters.incompressible){
        fluids_parameters.incompressible->mpi_grid=fluids_parameters.mpi_grid;
        fluids_parameters.incompressible->projection.elliptic_solver->mpi_grid=fluids_parameters.mpi_grid;}
    if(fluids_parameters.compressible){
        fluids_parameters.compressible_boundary=new BOUNDARY_MPI<T_GRID,VECTOR<typename T_GRID::SCALAR,T_GRID::dimension+2> >(fluids_parameters.mpi_grid,*fluids_parameters.compressible_boundary);
        fluids_parameters.euler->euler_projection.pressure_boundary=new BOUNDARY_MPI<T_GRID>(fluids_parameters.mpi_grid,*fluids_parameters.euler->euler_projection.pressure_boundary);
        fluids_parameters.euler->mpi_grid=fluids_parameters.mpi_grid;
        fluids_parameters.euler->euler_projection.elliptic_solver->mpi_grid=fluids_parameters.mpi_grid;}
    if(!restart && fluids_parameters.store_particle_ids){ // TODO: this should be fixed so that id's are scalable
        for(int i=1;i<=fluids_parameters.number_of_regions;i++)
            fluids_parameters.particle_levelset_evolution->Particle_Levelset(i).last_unique_particle_id=fluids_parameters.mpi_grid->rank*30000000;}       
}
//#####################################################################
// Function Initialize_Solid_Fluid_Coupling
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization()
{
/*    fluids_parameters.use_poisson=true;
    if(fluids_parameters.use_slip){
        SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=new SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>(solids_parameters,solid_body_collection,fluids_parameters,solids_fluids_parameters,
            *new INCOMPRESSIBLE_FLUID_CONTAINER());
        delete solids_evolution;
        solids_evolution=coupled_evolution;
        fluids_parameters.Set_Projection(coupled_evolution);}
    else{
        delete solids_evolution;
        solids_evolution=new SOLID_FLUID_COUPLED_EVOLUTION<TV>(solids_parameters,solid_body_collection,fluids_parameters,*new INCOMPRESSIBLE_FLUID_CONTAINER(),solids_fluids_parameters);}
        // TODO: set up anything that needs to be set up for solid_affects_fluid only case*/
}
//#####################################################################
// Function Initialize_Solid_Fluid_Coupling
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Initialize_Solid_Fluid_Coupling_After_Grid_Initialization()
{
    fluids_parameters.collision_bodies_affecting_fluid->Initialize_Grids();

#if 0
    bool initialize_valid_value_mask=!restart; // If we're restarting we'll be reading in the values anyway
    // TODO: update for new advections!
    thin_shells_semi_lagrangian_density.Initialize(fluids_parameters.grid,initialize_valid_value_mask);
    thin_shells_semi_lagrangian_density.interpolation->Set_Default_Replacement_Value(fluids_parameters.density_container.ambient_density);
    thin_shells_semi_lagrangian_temperature.Initialize(fluids_parameters.grid,initialize_valid_value_mask);
    thin_shells_semi_lagrangian_temperature.interpolation->Set_Default_Replacement_Value(fluids_parameters.temperature_container.ambient_temperature);
    thin_shells_semi_lagrangian_velocity.Initialize(fluids_parameters.grid,initialize_valid_value_mask);
    if(fluids_parameters.water||fluids_parameters.fire) thin_shells_semi_lagrangian_phi.Initialize(fluids_parameters.grid,initialize_valid_value_mask);
    if(fluids_parameters.water) thin_shells_semi_lagrangian_phi.interpolation->Set_Default_Replacement_Value((T)1e-5);
    if(fluids_parameters.fire){ // special treatment for fire because we advect different quantities using different velocity fields
        thin_shells_semi_lagrangian_phi.interpolation->Set_Default_Replacement_Value((T)1e-5);} // bias to negative phi for fire
#endif

    if(fluids_parameters.compressible){
        EULER_UNIFORM<T_GRID>& euler=*fluids_parameters.euler;
        SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES<TV>& euler_solid_fluid_coupling_utilities=*fluids_parameters.euler_solid_fluid_coupling_utilities;
        fluids_parameters.euler_solid_fluid_coupling_utilities->Initialize_Solid_Fluid_Coupling(fluids_parameters.collision_bodies_affecting_fluid);
        if(fluids_parameters.fluid_affects_solid && !euler.timesplit){
            EULER_FLUID_FORCES<T_GRID>* euler_fluid_forces=new EULER_FLUID_FORCES<T_GRID>(euler.grid,euler_solid_fluid_coupling_utilities.pressure_at_faces,euler_solid_fluid_coupling_utilities.solid_fluid_face_time_n,
                euler_solid_fluid_coupling_utilities.cells_inside_fluid_time_n,euler_solid_fluid_coupling_utilities.collision_bodies_affecting_fluid,
                solid_body_collection.deformable_body_collection.particles,solid_body_collection.rigid_body_collection);
            solid_body_collection.Add_Force(euler_fluid_forces);}}

    if(SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(solids_evolution))
        coupled_evolution->Setup_Boundary_Condition_Collection();
}
//#####################################################################
// Function Initialize_Compressible_Incompressible_Coupling
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Initialize_Compressible_Incompressible_Coupling()
{
    if(fluids_parameters.compressible && fluids_parameters.number_of_regions){
        fluids_parameters.euler->euler_projection.Set_Incompressible_Coupling_Callbacks(fluids_parameters.compressible_incompressible_coupling_utilities);}
}
//#####################################################################
// Function Set_Ghost_Density_And_Temperature_Inside_Flame_Core
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Set_Ghost_Density_And_Temperature_Inside_Flame_Core()
{
    for(CELL_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
        int region=fluids_parameters.particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Inside_Region(iterator.Cell_Index());
        if(fluids_parameters.fuel_region(region)){
            fluids_parameters.density_container.density(iterator.Cell_Index())=fluids_parameters.density;
            fluids_parameters.temperature_container.temperature(iterator.Cell_Index())=fluids_parameters.temperature_products;}}
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Set_Dirichlet_Boundary_Conditions(const T time)
{
    if(fluids_parameters.smoke||fluids_parameters.fire||fluids_parameters.water){
        if(fluids_parameters.number_of_regions>=2)
            fluids_parameters.incompressible_multiphase->Set_Dirichlet_Boundary_Conditions(fluids_parameters.particle_levelset_evolution_multiple->phis,fluids_parameters.dirichlet_regions);
        else if(fluids_parameters.number_of_regions==1) fluids_parameters.incompressible->Set_Dirichlet_Boundary_Conditions(&fluids_parameters.particle_levelset_evolution->phi,0);
        else fluids_parameters.incompressible->Set_Dirichlet_Boundary_Conditions(0,0);} // 1 phase
    else if(fluids_parameters.sph) for(CELL_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
        fluids_parameters.incompressible->projection.elliptic_solver->psi_D(iterator.Cell_Index())=true;fluids_parameters.incompressible->projection.p(iterator.Cell_Index())=0;}
    else if(fluids_parameters.compressible){
        fluids_parameters.euler->euler_projection.Set_Dirichlet_Boundary_Conditions(time);
        if(fluids_parameters.number_of_regions) fluids_parameters.incompressible->Set_Dirichlet_Boundary_Conditions(&fluids_parameters.particle_levelset_evolution->phi,
            fluids_parameters.compressible_incompressible_coupling_utilities->p_dirichlet_incompressible);}

    if(fluids_parameters.solid_affects_fluid && (fluids_parameters.use_slip || fluids_parameters.fluid_affects_solid)){
        if(fluids_parameters.use_slip){
            // slip sets up its own boundary conditions
            //SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(solids_parameters.solids_evolution);
            //coupled_evolution->Set_Dirichlet_Boundary_Conditions(incompressible_fluid_containers(chimera_grid->current_grid)->face_velocities,time);
        }
        else{
            SOLID_FLUID_COUPLED_EVOLUTION<TV>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION<TV>*>(solids_evolution);
            coupled_evolution->Set_Dirichlet_Boundary_Conditions(time);}}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Get_Source_Velocities_Chimera(const T time,const T dt)
{
    for(int i=1;i<=rigid_grid_collection.particles.array_collection->Size();i++)
        for(FACE_ITERATOR iterator(incompressible_fluid_containers(i)->rigid_grid.grid);iterator.Valid();iterator.Next())
        {
            TV_INT face_index=iterator.Face_Index();
            int face_axis=iterator.Axis();
            TV location=incompressible_fluid_containers(i)->rigid_grid.Frame()*iterator.Location();
            TV normal=incompressible_fluid_containers(i)->rigid_grid.World_Space_Vector(TV::Axis_Vector(face_axis));
            T& velocity=incompressible_fluid_containers(i)->face_velocities.Component(face_axis)(face_index);
            incompressible_fluid_containers(i)->psi_N.Component(face_axis)(face_index)=Get_Source_Velocity(location,normal,velocity,time,dt);
        }
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Get_Source_Velocities(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity)
{
    //for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) if(source.Lazy_Inside(world_to_source.Homogeneous_Times(iterator.Location()))){
    //    int axis=iterator.Axis();fluids_parameters.incompressible->projection.elliptic_solver->psi_N.Component(axis)(iterator.Face_Index())=true;
    //    incompressible_fluid_containers(chimera_grid->current_grid)->face_velocities.Component(axis)(iterator.Face_Index())=constant_source_velocity[axis];}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Get_Source_Velocities(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity,const T_FACE_ARRAYS_BOOL& invalid_mask)
{
    //for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){const int axis=iterator.Axis();const TV_INT face_index=iterator.Face_Index();
    //    if(!invalid_mask(axis,face_index) && source.Lazy_Inside(world_to_source.Homogeneous_Times(iterator.Location()))){
    //        fluids_parameters.incompressible->projection.elliptic_solver->psi_N.Component(axis)(face_index)=true;
    //        incompressible_fluid_containers(chimera_grid->current_grid)->face_velocities.Component(axis)(face_index)=constant_source_velocity[axis];}}
}
//#####################################################################
// Adjust_Phi_With_Sources (yuey)
//#####################################################################
template<class T_GRID> bool SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Adjust_Phi_With_Sources(const T time)
{
/*    RANGE<TV> range=RANGE<TV>::Zero_Box();
    range.min_corner(2)=0.45;range.max_corner(2)=0.55;range.min_corner(1)=-0.01;
    if(TV::dimension==3)range.min_corner(3)=-0.05;
    BOX<TV> source_box(range);
    T_ARRAYS_SCALAR& phi1=incompressible_fluid_containers(1)->particle_levelset_evolution.phi,phi2=incompressible_fluid_containers(2)->particle_levelset_evolution.phi;
    FRAME<TV> frame1=incompressible_fluid_containers(1)->rigid_grid.Frame(),frame2=incompressible_fluid_containers(2)->rigid_grid.Frame();
    for(typename GRID<TV>::CELL_ITERATOR iterator(incompressible_fluid_containers(1)->rigid_grid.grid);iterator.Valid();iterator.Next()){
        const TV &X=iterator.Location();
        phi1(iterator.Cell_Index())=min(phi1(iterator.Cell_Index()),source_box.Signed_Distance(X));}
    for(typename GRID<TV>::CELL_ITERATOR iterator(incompressible_fluid_containers(2)->rigid_grid.grid);iterator.Valid();iterator.Next()){
        const TV &X=iterator.Location();
        phi2(iterator.Cell_Index())=min(phi2(iterator.Cell_Index()),source_box.Signed_Distance(frame1.Inverse()*(frame2*X)));}*/
    return true;
}
//#####################################################################
// Function Adjust_Phi_With_Source
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Adjust_Phi_With_Source(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source)
{
    for(CELL_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){TV source_X=world_to_source.Homogeneous_Times(iterator.Location());
        if(source.Lazy_Inside(source_X)) 
            fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())=min(fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index()),
                                                                                          source.Signed_Distance(source_X));}
}
//#####################################################################
// Function Adjust_Phi_With_Source
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Adjust_Phi_With_Source(const GEOMETRY& source,const int region,const T_TRANSFORMATION_MATRIX& world_to_source)
{
    T bandwidth=3*fluids_parameters.grid->Minimum_Edge_Length();
    ARRAY<T_ARRAYS_SCALAR>& phis=fluids_parameters.particle_levelset_evolution_multiple->phis;
    for(CELL_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
        TV source_X=world_to_source.Homogeneous_Times(iterator.Location());
        if(source.Inside(source_X,-bandwidth)){
            T source_signed_distance=source.Signed_Distance(source_X);
            for(int i=1;i<=fluids_parameters.number_of_regions;i++){
                if(i==region) phis(i)(iterator.Cell_Index())=min(phis(i)(iterator.Cell_Index()),source_signed_distance);
                else phis(i)(iterator.Cell_Index())=max(phis(i)(iterator.Cell_Index()),-source_signed_distance);}}}
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Get_Source_Reseed_Mask(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,T_ARRAYS_BOOL*& cell_centered_mask,const bool reset_mask)
{
    if(reset_mask){if(cell_centered_mask) delete cell_centered_mask;cell_centered_mask=new T_ARRAYS_BOOL(fluids_parameters.grid->Domain_Indices(1));}
    T padding=3*fluids_parameters.grid->dX.Max();
    for(CELL_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) 
        if(!source.Outside(world_to_source.Homogeneous_Times(iterator.Location()),padding)) (*cell_centered_mask)(iterator.Cell_Index())=true;
}
//#####################################################################
// Function Adjust_Density_And_Temperature_With_Sources
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Adjust_Density_And_Temperature_With_Sources(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const T source_density,const T source_temperature)
{
    for(CELL_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) if(source.Lazy_Inside(world_to_source.Homogeneous_Times(iterator.Location()))){
        if(fluids_parameters.use_density) fluids_parameters.density_container.density(iterator.Cell_Index())=source_density;
        if(fluids_parameters.use_temperature) fluids_parameters.temperature_container.temperature(iterator.Cell_Index())=source_temperature;}
}
//#####################################################################
// Function Revalidate_Fluid_Scalars
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Revalidate_Fluid_Scalars()
{
    //for(int i=1;i<=fluids_parameters.number_of_regions;i++){
    //  T_FAST_LEVELSET& levelset=fluids_parameters.particle_levelset_evolution->Levelset(i);
    //  T_FAST_LEVELSET_ADVECTION& levelset_advection=fluids_parameters.particle_levelset_evolution->Levelset_Advection(i);
    //  int sign=1;if(fluids_parameters.number_of_regions>=2&&fluids_parameters.dirichlet_regions(i))sign=-1;
    //  if(levelset_advection.nested_semi_lagrangian_collidable)
    //      levelset_advection.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(*fluids_parameters.grid,sign*fluids_parameters.collidable_phi_replacement_value,levelset.phi);}
    //if(fluids_parameters.use_density){
    //  if(fluids_parameters.density_container.nested_semi_lagrangian_collidable)
    //      fluids_parameters.density_container.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(fluids_parameters.density_container.grid,fluids_parameters.density_container.ambient_density,
    //      fluids_parameters.density_container.density);}
    //if(fluids_parameters.use_temperature){
    //  if(fluids_parameters.temperature_container.nested_semi_lagrangian_collidable)
    //      fluids_parameters.temperature_container.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(fluids_parameters.temperature_container.grid,
    //      fluids_parameters.temperature_container.ambient_temperature,fluids_parameters.temperature_container.temperature);}
}
//#####################################################################
// Function Revalidate_Phi_After_Modify_Levelset
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Revalidate_Phi_After_Modify_Levelset()
{
    /* for(int i=1;i<=fluids_parameters.number_of_regions;i++){
        T_FAST_LEVELSET& levelset=fluids_parameters.particle_levelset_evolution->Levelset(i);
        T_FAST_LEVELSET_ADVECTION& levelset_advection=fluids_parameters.particle_levelset_evolution->Levelset_Advection(i);
        int sign=1;if(fluids_parameters.number_of_regions>=2&&fluids_parameters.dirichlet_regions(i))sign=-1;
        if(levelset_advection.nested_semi_lagrangian_collidable){
            levelset_advection.nested_semi_lagrangian_collidable->cell_valid_points_current=levelset_advection.nested_semi_lagrangian_collidable->cell_valid_points_next;
            levelset_advection.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(*fluids_parameters.grid,sign*fluids_parameters.collidable_phi_replacement_value,levelset.phi);}}*/
}
//#####################################################################
// Function Revalidate_Fluid_Velocity
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Revalidate_Fluid_Velocity(T_FACE_ARRAYS_SCALAR& face_velocities)
{

    //if(fluids_parameters.incompressible->nested_nested_semi_lagrangian_fire_multiphase_collidable) 
    //  fluids_parameters.incompressible->nested_nested_semi_lagrangian_fire_multiphase_collidable->Average_To_Invalidated_Face(*fluids_parameters.grid,face_velocities);
    //if(fluids_parameters.incompressible->nested_semi_lagrangian_collidable) 
    //  fluids_parameters.incompressible->nested_semi_lagrangian_collidable->Average_To_Invalidated_Face(*fluids_parameters.grid,face_velocities);
    //if(fluids_parameters.incompressible->nested_semi_lagrangian_collidable_slip)
    //    fluids_parameters.incompressible->nested_semi_lagrangian_collidable_slip->Average_To_Invalidated_Face(*fluids_parameters.grid,face_velocities);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Get_Object_Velocities(LAPLACE_UNIFORM<T_GRID>* elliptic_solver,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    if(fluids_parameters.solid_affects_fluid && (fluids_parameters.fluid_affects_solid || fluids_parameters.use_slip)){
        if(!fluids_parameters.use_slip){
           SOLID_FLUID_COUPLED_EVOLUTION<TV>& coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION<TV>&>(*solids_evolution);
           coupled_evolution.Apply_Solid_Boundary_Conditions(time,false,face_velocities);}
        else{
            SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>& coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>&>(*solids_evolution);
            coupled_evolution.Get_Coupled_Faces_And_Interpolated_Solid_Velocities(time,elliptic_solver->psi_N,face_velocities);}}
    else fluids_parameters.collision_bodies_affecting_fluid->Compute_Psi_N(elliptic_solver->psi_N,&face_velocities);
}
//#####################################################################
// Function Get_Levelset_Velocity
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Get_Levelset_Velocity(const T_GRID& grid,T_LEVELSET& levelset,T_FACE_ARRAYS_SCALAR& V_levelset,const T time) const
{
    //if(!fluids_parameters.use_reacting_flow)
    int current_grid=chimera_grid->Find_Current_Grid_Local_Index(grid);
    T_FACE_ARRAYS_SCALAR::Put(incompressible_fluid_containers(current_grid)->face_velocities,V_levelset);
    /*if(fluids_parameters.pseudo_dirichlet_regions.Number_True()>0){
        T_ARRAYS_SCALAR phi_for_pseudo_dirichlet_regions;T_GRID grid_temp(grid);T_FAST_LEVELSET levelset_for_pseudo_dirichlet_regions(grid_temp,phi_for_pseudo_dirichlet_regions);
        fluids_parameters.particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Get_Single_Levelset(fluids_parameters.pseudo_dirichlet_regions,levelset_for_pseudo_dirichlet_regions,false);
        fluids_parameters.incompressible->Extrapolate_Velocity_Across_Interface(V_levelset,phi_for_pseudo_dirichlet_regions);}*/
}
//#####################################################################
// Function Initialize_Swept_Occupied_Blocks_For_Advection
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Initialize_Swept_Occupied_Blocks_For_Advection(const T dt,const T time,const bool include_removed_particle_velocities)
{
    for(int grid_index=1;grid_index<=rigid_grid_collection.particles.array_collection->Size();grid_index++){
        T_GRID& grid=rigid_grid_collection.Rigid_Grid(grid_index).grid;
        T maximum_particle_speed=0;
        T maximum_fluid_speed=incompressible_fluid_containers(grid_index)->face_velocities.Maxabs().Max();
        //if(fluids_parameters.fire){
        //  if(fluids_parameters.number_of_regions==1){
        //      Get_Levelset_Velocity(*fluids_parameters.grid,fluids_parameters.particle_levelset_evolution->particle_levelset.levelset,
        //          fluids_parameters.particle_levelset_evolution->V,time);
        //      maximum_fluid_speed=max(maximum_fluid_speed,fluids_parameters.particle_levelset_evolution->V.Maxabs().Max());}
        //  else{
        //      Get_Levelset_Velocity(*fluids_parameters.grid,fluids_parameters.particle_levelset_evolution_multiple->Levelset_Multiple(),
        //          fluids_parameters.particle_levelset_evolution_multiple->V,time);
        //      maximum_fluid_speed=max(maximum_fluid_speed,fluids_parameters.particle_levelset_evolution_multiple->V.Maxabs().Max());}}
        if(include_removed_particle_velocities){
            PARTICLE_LEVELSET_UNIFORM<T_GRID>& particle_levelset=incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset;
            if(particle_levelset.use_removed_negative_particles) for(CELL_ITERATOR iterator(particle_levelset.levelset.grid);iterator.Valid();iterator.Next()){
                PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles=particle_levelset.removed_negative_particles(iterator.Cell_Index());
                if(particles) maximum_particle_speed=max(maximum_particle_speed,ARRAYS_COMPUTATIONS::Maximum_Magnitude(particles->V));}
            if(particle_levelset.use_removed_positive_particles) for(CELL_ITERATOR iterator(particle_levelset.levelset.grid);iterator.Valid();iterator.Next()){
                PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles=particle_levelset.removed_positive_particles(iterator.Cell_Index());
                if(particles) maximum_particle_speed=max(maximum_particle_speed,ARRAYS_COMPUTATIONS::Maximum_Magnitude(particles->V));}}
        T max_particle_collision_distance=0;
        
        max_particle_collision_distance=max(max_particle_collision_distance,incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.max_collision_distance_factor*grid.dX.Max());
        incompressible_fluid_containers(grid_index)->collision_bodies_affecting_fluid.Compute_Occupied_Blocks(true,
            dt*max(maximum_fluid_speed,maximum_particle_speed)+2*max_particle_collision_distance+(T).5*rigid_grid_collection.Rigid_Grid(grid_index).grid.dX.Max(),10);
        
        //if(fluids_parameters.use_maccormack_semi_lagrangian_advection && fluids_parameters.use_maccormack_compute_mask){
        //  typedef typename T_GRID::VECTOR_INT TV_INT;
        //  fluids_parameters.maccormack_cell_mask.Resize(grid.Domain_Indices(fluids_parameters.number_of_ghost_cells),false,false);
        //  fluids_parameters.maccormack_face_mask.Resize(grid,fluids_parameters.number_of_ghost_cells);
            // don't use maccormack near domain boundary conditions
        //  for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
        //      fluids_parameters.maccormack_cell_mask(iterator.Cell_Index())=TV_INT::Componentwise_Min(iterator.Cell_Index()-grid.Domain_Indices().Minimum_Corner(),
        //          grid.Domain_Indices().Maximum_Corner()-iterator.Cell_Index()).Min()>fluids_parameters.cfl;
        //  VECTOR<T_GRID,T_GRID::dimension> grids;
        //  for(int i=1;i<=T_GRID::dimension;i++) grids[i]=grid.Get_Face_Grid(i);
        //  for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        //      int axis=iterator.Axis();TV_INT index=iterator.Face_Index();
        //      fluids_parameters.maccormack_face_mask(axis,index)
        //          =TV_INT::Componentwise_Min(index-grids[axis].Domain_Indices().Minimum_Corner(),grids[axis].Domain_Indices().Maximum_Corner()-index).Min()>fluids_parameters.cfl;}
            //if(fluids_parameters.mpi_grid){ // if mpi check turned off domain boundary regions to make sure they are global domain boundaries
                //  RANGE<TV> global_domain=fluids_parameters.mpi_grid->global_grid.domain;T double_cfl_dx=(T)fluids_parameters.cfl*fluids_parameters.p_grid.dX.Max();
                //for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
            //   if(!fluids_parameters.maccormack_cell_mask(iterator.Cell_Index()) && global_domain.Inside(iterator.Location(),(T)double_cfl_dx))
            //          fluids_parameters.maccormack_cell_mask(iterator.Cell_Index())=true;
            //  for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            //      int axis=iterator.Axis();TV_INT index=iterator.Face_Index();
            //      if(!fluids_parameters.maccormack_face_mask(axis,index) && global_domain.Inside(iterator.Location(),(T)double_cfl_dx))
            //          fluids_parameters.maccormack_face_mask(axis,index)=true;}}
            // don't use maccormack near the interface
            //if(fluids_parameters.bandwidth_without_maccormack_near_interface){
            //  T dx_times_band=grid.Maximum_Edge_Length()*fluids_parameters.bandwidth_without_maccormack_near_interface;
            //  for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next())if(fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())>-dx_times_band){
            //      fluids_parameters.maccormack_cell_mask(iterator.Cell_Index())=false;
            //      for(int axis=1;axis<=T_GRID::dimension;axis++){
            //          fluids_parameters.maccormack_face_mask(axis,iterator.First_Face_Index(axis))=false;fluids_parameters.maccormack_face_mask(axis,iterator.Second_Face_Index(axis))=false;}}}
            // turn off maccormack near objects
            //for(typename T_GRID::NODE_ITERATOR node_iterator(grid,2);node_iterator.Valid();node_iterator.Next()) 
        //  if(fluids_parameters.collision_bodies_affecting_fluid->occupied_blocks(node_iterator.Node_Index())){
        //          TV_INT block_index=node_iterator.Node_Index();BLOCK_UNIFORM<T_GRID> block(grid,block_index);
        //          for(int cell_index=0;cell_index<T_GRID::number_of_cells_per_block;cell_index++) fluids_parameters.maccormack_cell_mask(block.Cell(cell_index))=false;
        //          for(int axis=1;axis<=T_GRID::dimension;axis++) for(int face=0;face<T_GRID::number_of_incident_faces_per_block;face++)
        //              fluids_parameters.maccormack_face_mask(axis,block.Incident_Face(axis,face))=false;}}}
    }
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time)
{}
//#####################################################################
// Function Read_Output_Files_Fluids
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Read_Output_Files_Fluids(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    std::string filename;
    for(int i=1;i<=incompressible_fluid_containers.Size();i++){
        std::string g=chimera_grid->use_mpi?STRING_UTILITIES::string_sprintf("%d",chimera_grid->local_to_global_grid_index_map(i)):STRING_UTILITIES::string_sprintf("%d",i);
        filename=output_directory+"/"+f+"/chimera_grid_"+g;
        if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading chimera_grid_ "<<i<<std::endl;LOG::filecout(ss.str());
            FILE_UTILITIES::Read_From_File(stream_type,filename,rigid_grid_collection.Rigid_Grid(i).grid);}
        filename=output_directory+"/"+f+"/rigid_grid_frame_"+g;
        if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading rigid_grid_frame_"<<i<<std::endl;LOG::filecout(ss.str());
            FRAME<TV> grid_frame;
            FILE_UTILITIES::Read_From_File(stream_type,filename,grid_frame);rigid_grid_collection.Rigid_Grid(i).Set_Frame(grid_frame);}
        if(fluids_parameters.smoke){
            filename=output_directory+"/"+f+"/density_"+g;
            if(FILE_UTILITIES::File_Exists(filename)){
                std::stringstream ss;ss<<"Reading density_"<<i<<std::endl;LOG::filecout(ss.str());
                FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible_fluid_containers(i)->density_container.density);}}
        filename=output_directory+"/"+f+"/mac_velocities_"+g;
        if(FILE_UTILITIES::File_Exists(filename)){
            std::stringstream ss;ss<<"Reading mac_velocities_"<<i<<std::endl;LOG::filecout(ss.str());
            if(solve_heat_equation_on_faces) FILE_UTILITIES::Read_From_File(stream_type,filename,face_temperatures(i));
            else FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible_fluid_containers(i)->face_velocities);}
        filename=output_directory+"/"+f+"/psi_N_"+g;
        if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading psi_N_"<<i<<std::endl;LOG::filecout(ss.str());
            FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible_fluid_containers(i)->psi_N);}
        filename=output_directory+"/"+f+"/pressure_"+g;
        if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading pressure_"<<i<<std::endl;LOG::filecout(ss.str());
            FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible_fluid_containers(i)->pressure);}
        filename=output_directory+"/"+f+"/temperature_"+g;
        if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading temperature_"<<i<<std::endl;LOG::filecout(ss.str());
            FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible_fluid_containers(i)->temperature_container.temperature);}
        if(fluids_parameters.water){
            filename=output_directory+"/"+f+"/chimera_levelset_"+g;
            if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading chimera_levelset_"<<i<<std::endl;LOG::filecout(ss.str());
                FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible_fluid_containers(i)->particle_levelset_evolution.particle_levelset.levelset);}
            filename=output_directory+"/"+f+"/positive_particles_"+g;
            if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading positive_particles_"<<i<<std::endl;LOG::filecout(ss.str());
                FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible_fluid_containers(i)->particle_levelset_evolution.particle_levelset.positive_particles);}
            filename=output_directory+"/"+f+"/negative_particles_"+g;
            if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading negative_particles_"<<i<<std::endl;LOG::filecout(ss.str());
                FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible_fluid_containers(i)->particle_levelset_evolution.particle_levelset.negative_particles);}
            filename=output_directory+"/"+f+"/removed_positive_particles_"+g;
            if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading removed_positive_particles_"<<i<<std::endl;LOG::filecout(ss.str());
                FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible_fluid_containers(i)->particle_levelset_evolution.particle_levelset.removed_positive_particles);}
            filename=output_directory+"/"+f+"/removed_negative_particles_"+g;
            if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading removed_negative_particles_"<<i<<std::endl;LOG::filecout(ss.str());
                FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible_fluid_containers(i)->particle_levelset_evolution.particle_levelset.removed_negative_particles);}
            filename=output_directory+"/"+f+"/last_unique_particle_id_"+g;
            if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading last_unique_particle_id_"<<i<<std::endl;LOG::filecout(ss.str());
                FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible_fluid_containers(i)->particle_levelset_evolution.particle_levelset.last_unique_particle_id);}
        }
    }

    //SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>* evolution;
    //evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID> * >(solids_evolution);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Write_Output_Files(const int frame) const
{
    LOG::SCOPE scope(STRING_UTILITIES::string_sprintf("writing example output files %d",chimera_grid->rank));
    FILE_UTILITIES::Create_Directory(output_directory);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    Write_Frame_Title(frame);

    //solids_parameters.write_deformable_body=false; //temporary workaround

    if(!chimera_grid->use_mpi || (chimera_grid->use_mpi && chimera_grid->rank==0))
        solid_body_collection.Write(stream_type,output_directory,frame,first_frame,solids_parameters.write_static_variables_every_frame,
            solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies,solids_parameters.write_deformable_body,solids_parameters.write_from_every_process,
            solids_parameters.triangle_collision_parameters.output_interaction_pairs);
    //if(NEWMARK_EVOLUTION<TV>* newmark=dynamic_cast<NEWMARK_EVOLUTION<TV>*>(solids_evolution))
    //  newmark->Write_Position_Update_Projection_Data(stream_type,output_directory+"/"+f+"/");

    SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID>* evolution;
    evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA<T_GRID> * >(solids_evolution);

/////////////////////////////////////////////////////////////////////////////////////////////////DEBUG ERROR
    /*for(int i=1;i<=incompressible_fluid_containers.Size();i++)
        for(CELL_ITERATOR iterator(rigid_grid_collection.Rigid_Grid(i).grid);iterator.Valid();iterator.Next())
        {
            T time=(T)frame/24;
            T PI=3.14159265,bump_radius=0.3;
            TV bump_center;bump_center(1)=0.5*cos(PI/2+time);bump_center(2)=0.5*sin(PI/2+time);
            TV location=incompressible_fluid_containers(i)->rigid_grid.Frame()*iterator.Location();
            T dist_to_center=(location-bump_center).Magnitude();
            T analytic_solution=dist_to_center-bump_radius;
            incompressible_fluid_containers(i)->temperature_container.temperature(iterator.Cell_Index())=abs(incompressible_fluid_containers(i)->density_container.density(iterator.Cell_Index())-analytic_solution);
        }*/
////////////////////////////////////////////////////////////////////////////////////////////////

    for(int i=1;i<=incompressible_fluid_containers.Size();i++){
        int grid_index=chimera_grid->local_to_global_grid_index_map(i);
        std::string g=chimera_grid->use_mpi?STRING_UTILITIES::string_sprintf("_%d",grid_index):STRING_UTILITIES::string_sprintf("_%d",i);
        if(frame==first_frame)
            FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common"+"/chimera_grid"+g,rigid_grid_collection.Rigid_Grid(i).grid);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/chimera_grid"+g,rigid_grid_collection.Rigid_Grid(i).grid);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/rigid_grid_frame"+g,rigid_grid_collection.Rigid_Grid(i).Frame());
        if(fluids_parameters.smoke) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density"+g,incompressible_fluid_containers(i)->density_container.density);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/valid_cell_points_current"+g,incompressible_fluid_containers(i)->density_container.valid_mask_current);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/valid_cell_points_next"+g,incompressible_fluid_containers(i)->density_container.valid_mask_next);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/face_valid_mask"+g,incompressible_fluid_containers(i)->face_velocities_valid_mask);
        if(solve_heat_equation_on_faces) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities"+g,face_temperatures(i));
        else FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities"+g,incompressible_fluid_containers(i)->face_velocities);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_N"+g,incompressible_fluid_containers(i)->psi_N);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D"+g,incompressible_fluid_containers(i)->psi_D);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure"+g,incompressible_fluid_containers(i)->pressure);
        if(solve_heat_equation) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/temperature"+g,incompressible_fluid_containers(i)->temperature_container.temperature);

        //FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/vorticity"+g,evolution->vorticity(grid_index));
    }
    if(!chimera_grid->use_mpi){
        if(write_voronoi)
            FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/voronoi",evolution->laplace_grid.voronoi);}
    else{
        if(frame==first_frame) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common"+"/side_neighbor_ranks",chimera_grid->side_neighbor_ranks_chimera);
        std::string g=STRING_UTILITIES::string_sprintf("_%d",chimera_grid->rank+1);
        if(write_voronoi)
            FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/voronoi"+g,evolution->laplace_grid.voronoi);}
    //particle levelset
    if(fluids_parameters.water) for(int i=1;i<=incompressible_fluid_containers.Size();i++){
        PARTICLE_LEVELSET_UNIFORM<T_GRID>& particle_levelset=incompressible_fluid_containers(i)->particle_levelset_evolution.particle_levelset;
        std::string g=chimera_grid->use_mpi?STRING_UTILITIES::string_sprintf("_%d",chimera_grid->local_to_global_grid_index_map(i)):STRING_UTILITIES::string_sprintf("_%d",i);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/chimera_levelset"+g,particle_levelset.levelset);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/positive_particles"+g,particle_levelset.positive_particles);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/negative_particles"+g,particle_levelset.negative_particles);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/removed_positive_particles"+g,particle_levelset.removed_positive_particles);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/removed_negative_particles"+g,particle_levelset.removed_negative_particles);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/last_unique_particle_id"+g,particle_levelset.last_unique_particle_id);}
}
//#####################################################################
// Function Log_Parameters 
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>::
Log_Parameters() const
{       
    LOG::SCOPE scope("SOLIDS_FLUIDS_EXAMPLE_CHIMERA parameters");
    BASE::Log_Parameters();
    fluids_parameters.Log_Parameters();
}
//#####################################################################
#define PROTECT(...) __VA_ARGS__
#define INSTANTIATION_HELPER_T_d_SHAPE(T,d,SHAPE) \
    template void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<GRID<VECTOR<T,d> > >::Get_Source_Velocities(const SHAPE&,const MATRIX<T,d+1>&,const VECTOR<T,d>&); \
    template void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<GRID<VECTOR<T,d> > >::Get_Source_Velocities(const SHAPE&,const MATRIX<T,d+1>&,const VECTOR<T,d>&,const T_FACE_ARRAYS_BOOL&); \
    template void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<GRID<VECTOR<T,d> > >::Adjust_Phi_With_Source(const SHAPE&,const MATRIX<T,d+1>&); \
    template void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<GRID<VECTOR<T,d> > >::Adjust_Phi_With_Source(const SHAPE&,const int,const MATRIX<T,d+1>&); \
    template void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<GRID<VECTOR<T,d> > >::Get_Source_Reseed_Mask(const SHAPE& source,const MATRIX<T,d+1>&,T_ARRAYS_BOOL*&,const bool); \
    template void SOLIDS_FLUIDS_EXAMPLE_CHIMERA<GRID<VECTOR<T,d> > >::Adjust_Density_And_Temperature_With_Sources(const SHAPE&,const MATRIX<T,d+1>&,const T,const T);
#define INSTANTIATION_HELPER_SHAPE(...) INSTANTIATION_HELPER_T_d_SHAPE(PROTECT(__VA_ARGS__::VECTOR_T::SCALAR),PROTECT(__VA_ARGS__::VECTOR_T::m),PROTECT(__VA_ARGS__))
#define INSTANTIATION_HELPER(T) \
    template class SOLIDS_FLUIDS_EXAMPLE_CHIMERA<GRID<VECTOR<T,1> > >;   \
    template class SOLIDS_FLUIDS_EXAMPLE_CHIMERA<GRID<VECTOR<T,2> > >;    \
    template class SOLIDS_FLUIDS_EXAMPLE_CHIMERA<GRID<VECTOR<T,3> > >;   \
    INSTANTIATION_HELPER_SHAPE(BOX<VECTOR<T,2> >) \
    INSTANTIATION_HELPER_SHAPE(BOX<VECTOR<T,3> >) \
    INSTANTIATION_HELPER_SHAPE(SPHERE<VECTOR<T,2> >) \
    INSTANTIATION_HELPER_SHAPE(SPHERE<VECTOR<T,3> >) \
    INSTANTIATION_HELPER_SHAPE(CYLINDER<T>)
INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
