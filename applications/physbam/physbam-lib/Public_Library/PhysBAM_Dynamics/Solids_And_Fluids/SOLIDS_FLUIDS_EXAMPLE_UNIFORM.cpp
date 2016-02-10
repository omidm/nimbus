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
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
SOLIDS_FLUIDS_EXAMPLE_UNIFORM(const STREAM_TYPE stream_type,const int number_of_regions,const typename FLUIDS_PARAMETERS<T_GRID>::TYPE type,const int array_collection_type)
    :SOLIDS_FLUIDS_EXAMPLE<TV>(stream_type,array_collection_type),fluids_parameters(number_of_regions,type,incompressible_fluid_container),resolution(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
~SOLIDS_FLUIDS_EXAMPLE_UNIFORM()
{
    if(dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(solids_evolution)) fluids_parameters.projection=0;
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Register_Options()
{
    BASE::Register_Options();
    this->parse_args->Add_Option_Argument("-skip_debug_data","turn off file io for debug data");
    this->parse_args->Add_Option_Argument("-use_fmm_extrapolation","use fast marching to extrapolate compressible flow data into solid state.");
    this->parse_args->Add_Integer_Argument("-resolution",1);
}
//#####################################################################
// Function Parse_Options
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Parse_Options()
{
    BASE::Parse_Options();
    if(this->parse_args->Is_Value_Set("-skip_debug_data")) fluids_parameters.write_debug_data=false;
    if(this->parse_args->Is_Value_Set("-use_fmm_extrapolation")) fluids_parameters.euler_solid_fluid_coupling_utilities->use_fast_marching=true;
    resolution=this->parse_args->Get_Integer_Value("-resolution");
}
//#####################################################################
// Function Add_Volumetric_Body_To_Fluid_Simulation
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Add_Volumetric_Body_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,bool add_collision,bool add_coupling)
{
    RIGID_COLLISION_GEOMETRY<TV>* collision_geometry=new RIGID_COLLISION_GEOMETRY<TV>(rigid_body);
    if(add_collision) fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(collision_geometry,rigid_body.particle_index,true);
    if(add_coupling)
        if(SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>*>(solids_evolution))
            coupled_evolution->iterator_info.coupling_bodies.Append(collision_geometry);
    if(rigid_body.simplicial_object) rigid_body.simplicial_object->Initialize_Hierarchy();
}
//#####################################################################
// Function Add_Thin_Shell_To_Fluid_Simulation
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
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
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
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
// Function Clamp_Dt_Adaptively
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Clamp_Dt_Adaptively(T_GRID& grid,const T_ARRAYS_DIMENSION_SCALAR& rhs,const T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_BOOL& psi,T_ARRAYS_SCALAR& rho_dt,T_ARRAYS_SCALAR& e_dt,const T& dt,T& clamp_rho,T& clamp_e)
{
    if(fluids_parameters.compressible_adaptive_time_step) for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(psi(cell_index)){
            T clamp_rho_cell=clamp_rho*U(cell_index)(1);
            if(rhs(cell_index)(1)>0 && U(cell_index)(1)>=clamp_rho_cell) rho_dt(cell_index)=(U(cell_index)(1)-clamp_rho_cell)/rhs(cell_index)(1);
            else rho_dt(cell_index)=dt;

            T e=EULER<T_GRID>::e(U,cell_index);int d=T_GRID::dimension+2;
            T clamp_e_cell=clamp_e*e,tmp_dt=dt,momentum_flux_sqr=0,momentum_sqr=0,momentum_flux_dot_product=0;
            for(int axis=1;axis<=TV::dimension;axis++){
                momentum_flux_sqr+=rhs(cell_index)(axis+1)*rhs(cell_index)(axis+1);
                momentum_sqr+=U(cell_index)(axis+1)*U(cell_index)(axis+1);
                momentum_flux_dot_product+=U(cell_index)(axis+1)*rhs(cell_index)(axis+1);}

                T rho_np1=U(cell_index)(1)-rho_dt(cell_index)*rhs(cell_index)(1);
                T rho_np1_sqr=rho_np1*rho_np1;
                assert(rho_np1>0);
                T E_np1=U(cell_index)(d)-rho_dt(cell_index)*rhs(cell_index)(d);
                T mom_np1_sqr=momentum_sqr-2*rho_dt(cell_index)*momentum_flux_dot_product+rho_dt(cell_index)*rho_dt(cell_index)*momentum_flux_sqr;
                T e_np1=E_np1/rho_np1-(T).5*mom_np1_sqr/rho_np1_sqr;

                T a=2*rhs(cell_index)(1)*rhs(cell_index)(d)-momentum_flux_sqr-2*clamp_e_cell*rhs(cell_index)(1)*rhs(cell_index)(1);
                T c=2*U(cell_index)(d)*U(cell_index)(1)-2*clamp_e_cell*U(cell_index)(1)*U(cell_index)(1)-momentum_sqr;

                
                if(e_np1<clamp_e_cell){
                    T b_over_two=momentum_flux_dot_product-U(cell_index)(1)*rhs(cell_index)(d)-U(cell_index)(d)*rhs(cell_index)(1)+2*clamp_e_cell*U(cell_index)(1)*rhs(cell_index)(1);
                    T b_sqr_over_four=b_over_two*b_over_two, ac=a*c;

                    if(b_sqr_over_four>ac){
                        if((a<1e-16&&b_over_two<1e-16)||(a>-1e-16&&b_over_two>-1e-16)) tmp_dt=dt;
                        else if(a<1e-16||a>-1e16){T tmp_dt1=(-(T).5*c)/b_over_two;tmp_dt=(tmp_dt1>0)?tmp_dt1:dt;}
                        else if(b_over_two<1e-16||b_over_two>-1e-16){T tmp_dt1=abs(sqrt(-c/a));tmp_dt=(tmp_dt1>0)?tmp_dt1:dt;}
                        else{T tmp_dt1=(-1*b_over_two-sqrt(b_sqr_over_four-ac))/a;
                        assert(tmp_dt1!=0);
                        T tmp_dt2=c/(-1*b_over_two-sqrt(b_sqr_over_four-ac));
                        tmp_dt=(tmp_dt1>0)?(tmp_dt2>0)?min(tmp_dt1,tmp_dt2):tmp_dt1:(tmp_dt2>0)?tmp_dt2:dt;}}
                    else tmp_dt=dt;
                assert(tmp_dt>0);
                e_dt(cell_index)=min(dt,tmp_dt);}}}
}
//#####################################################################
// Function Initialize_MPI
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
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
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization()
{
    fluids_parameters.use_poisson=true;
    if(fluids_parameters.use_slip){
        SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>* coupled_evolution=new SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>(solids_parameters,solid_body_collection,fluids_parameters,solids_fluids_parameters,
            incompressible_fluid_container);
        delete solids_evolution;
        solids_evolution=coupled_evolution;
        fluids_parameters.Set_Projection(coupled_evolution);}
    else{
        delete solids_evolution;
        solids_evolution=new SOLID_FLUID_COUPLED_EVOLUTION<TV>(solids_parameters,solid_body_collection,fluids_parameters,incompressible_fluid_container,solids_fluids_parameters);}
    // TODO: set up anything that needs to be set up for solid_affects_fluid only case
}
//#####################################################################
// Function Initialize_Solid_Fluid_Coupling
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
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
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Initialize_Compressible_Incompressible_Coupling()
{
    if(fluids_parameters.compressible && fluids_parameters.number_of_regions){
        fluids_parameters.euler->euler_projection.Set_Incompressible_Coupling_Callbacks(fluids_parameters.compressible_incompressible_coupling_utilities);}
}
//#####################################################################
// Function Set_Ghost_Density_And_Temperature_Inside_Flame_Core
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
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
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
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
            //coupled_evolution->Set_Dirichlet_Boundary_Conditions(incompressible_fluid_container.face_velocities,time);
        }
        else{
            SOLID_FLUID_COUPLED_EVOLUTION<TV>* coupled_evolution=dynamic_cast<SOLID_FLUID_COUPLED_EVOLUTION<TV>*>(solids_evolution);
            coupled_evolution->Set_Dirichlet_Boundary_Conditions(time);}}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Get_Source_Velocities(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity)
{
    for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) if(source.Lazy_Inside(world_to_source.Homogeneous_Times(iterator.Location()))){
        int axis=iterator.Axis();fluids_parameters.incompressible->projection.elliptic_solver->psi_N.Component(axis)(iterator.Face_Index())=true;
        incompressible_fluid_container.face_velocities.Component(axis)(iterator.Face_Index())=constant_source_velocity[axis];}
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Get_Source_Velocities(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity,const T_FACE_ARRAYS_BOOL& invalid_mask)
{
    for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){const int axis=iterator.Axis();const TV_INT face_index=iterator.Face_Index();
        if(!invalid_mask(axis,face_index) && source.Lazy_Inside(world_to_source.Homogeneous_Times(iterator.Location()))){
            fluids_parameters.incompressible->projection.elliptic_solver->psi_N.Component(axis)(face_index)=true;
            incompressible_fluid_container.face_velocities.Component(axis)(face_index)=constant_source_velocity[axis];}}
}
//#####################################################################
// Function Adjust_Phi_With_Source
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
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
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
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
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
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
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Adjust_Density_And_Temperature_With_Sources(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const T source_density,const T source_temperature)
{
    for(CELL_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()) if(source.Lazy_Inside(world_to_source.Homogeneous_Times(iterator.Location()))){
        if(fluids_parameters.use_density) fluids_parameters.density_container.density(iterator.Cell_Index())=source_density;
        if(fluids_parameters.use_temperature) fluids_parameters.temperature_container.temperature(iterator.Cell_Index())=source_temperature;}
}
//#####################################################################
// Function Revalidate_Fluid_Scalars
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Revalidate_Fluid_Scalars()
{
    for(int i=1;i<=fluids_parameters.number_of_regions;i++){
        T_FAST_LEVELSET& levelset=fluids_parameters.particle_levelset_evolution->Levelset(i);
        T_FAST_LEVELSET_ADVECTION& levelset_advection=fluids_parameters.particle_levelset_evolution->Levelset_Advection(i);
        int sign=1;if(fluids_parameters.number_of_regions>=2&&fluids_parameters.dirichlet_regions(i))sign=-1;
        if(levelset_advection.nested_semi_lagrangian_collidable)
            levelset_advection.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(*fluids_parameters.grid,sign*fluids_parameters.collidable_phi_replacement_value,levelset.phi);}
    if(fluids_parameters.use_density){
        if(fluids_parameters.density_container.nested_semi_lagrangian_collidable)
            fluids_parameters.density_container.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(fluids_parameters.density_container.grid,fluids_parameters.density_container.ambient_density,
            fluids_parameters.density_container.density);}
    if(fluids_parameters.use_temperature){
        if(fluids_parameters.temperature_container.nested_semi_lagrangian_collidable)
            fluids_parameters.temperature_container.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(fluids_parameters.temperature_container.grid,
                fluids_parameters.temperature_container.ambient_temperature,fluids_parameters.temperature_container.temperature);}
}
//#####################################################################
// Function Revalidate_Phi_After_Modify_Levelset
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Revalidate_Phi_After_Modify_Levelset()
{
    for(int i=1;i<=fluids_parameters.number_of_regions;i++){
        T_FAST_LEVELSET& levelset=fluids_parameters.particle_levelset_evolution->Levelset(i);
        T_FAST_LEVELSET_ADVECTION& levelset_advection=fluids_parameters.particle_levelset_evolution->Levelset_Advection(i);
        int sign=1;if(fluids_parameters.number_of_regions>=2&&fluids_parameters.dirichlet_regions(i))sign=-1;
        if(levelset_advection.nested_semi_lagrangian_collidable){
            levelset_advection.nested_semi_lagrangian_collidable->cell_valid_points_current=levelset_advection.nested_semi_lagrangian_collidable->cell_valid_points_next;
            levelset_advection.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(*fluids_parameters.grid,sign*fluids_parameters.collidable_phi_replacement_value,levelset.phi);}}
}
//#####################################################################
// Function Revalidate_Fluid_Velocity
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Revalidate_Fluid_Velocity(T_FACE_ARRAYS_SCALAR& face_velocities)
{
    if(fluids_parameters.incompressible->nested_nested_semi_lagrangian_fire_multiphase_collidable) 
        fluids_parameters.incompressible->nested_nested_semi_lagrangian_fire_multiphase_collidable->Average_To_Invalidated_Face(*fluids_parameters.grid,face_velocities);
    if(fluids_parameters.incompressible->nested_semi_lagrangian_collidable) 
        fluids_parameters.incompressible->nested_semi_lagrangian_collidable->Average_To_Invalidated_Face(*fluids_parameters.grid,face_velocities);
    //if(fluids_parameters.incompressible->nested_semi_lagrangian_collidable_slip)
    //    fluids_parameters.incompressible->nested_semi_lagrangian_collidable_slip->Average_To_Invalidated_Face(*fluids_parameters.grid,face_velocities);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
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
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Get_Levelset_Velocity(const T_GRID& grid,LEVELSET_MULTIPLE_UNIFORM<T_GRID>& levelset_multiple,T_FACE_ARRAYS_SCALAR& V_levelset,const T time) const
{
    if(!fluids_parameters.use_reacting_flow) T_FACE_ARRAYS_SCALAR::Put(incompressible_fluid_container.face_velocities,V_levelset);
    else{
        const PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection=fluids_parameters.incompressible_multiphase->projection;
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face=iterator.Face_Index();
            TV_INT cell1,cell2;grid.Cells_Touching_Face(axis,face,cell1,cell2);
            int region1,other_region1,region2,other_region2;T phi1,other_phi1,phi2,other_phi2;
            levelset_multiple.Two_Minimum_Regions(cell1,region1,other_region1,phi1,other_phi1);
            levelset_multiple.Two_Minimum_Regions(cell2,region2,other_region2,phi2,other_phi2);
            // TODO: break out of this if phi1 and ph2 are far from the interface
            int face_region=region1;
            if(region1==region2){ // not close enough to the interface
                if(phi1>phi2){region2=other_region1;phi2=other_phi1;} // cell1 is close to the interface
                else{region1=other_region2;phi1=other_phi2;}} // cell2 is closer to the interface
            else if(phi1>phi2) face_region=region2; // see LEVELSET_MULTIPLE::Inside_Region_Face()
            if(projection.flame_speed_constants(region1,region2).z==0) V_levelset.Component(axis)(face)=incompressible_fluid_container.face_velocities.Component(axis)(face);
            else{
                int fuel_region=region1,product_region=region2;if(fluids_parameters.densities(fuel_region)<fluids_parameters.densities(product_region)){fuel_region=region2;product_region=region1;}
                V_levelset.Component(axis)(face)=projection.Face_Velocity_With_Ghost_Value_Multiphase(incompressible_fluid_container.face_velocities,axis,face,fuel_region,face_region)-
                    projection.Flame_Speed_Face_Multiphase(axis,face,fuel_region,product_region)*
                    (levelset_multiple.Phi(fuel_region,cell2)-levelset_multiple.Phi(fuel_region,cell1))*grid.one_over_dX[axis];}}}
    if(fluids_parameters.pseudo_dirichlet_regions.Number_True()>0){
        T_ARRAYS_SCALAR phi_for_pseudo_dirichlet_regions;T_GRID grid_temp(grid);T_FAST_LEVELSET levelset_for_pseudo_dirichlet_regions(grid_temp,phi_for_pseudo_dirichlet_regions);
        fluids_parameters.particle_levelset_evolution_multiple->particle_levelset_multiple.levelset_multiple.Get_Single_Levelset(fluids_parameters.pseudo_dirichlet_regions,levelset_for_pseudo_dirichlet_regions,false);
        fluids_parameters.incompressible->Extrapolate_Velocity_Across_Interface(V_levelset,phi_for_pseudo_dirichlet_regions);}
}
//#####################################################################
// Function Initialize_Swept_Occupied_Blocks_For_Advection
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Initialize_Swept_Occupied_Blocks_For_Advection(const T dt,const T time,T maximum_fluid_speed,const bool include_removed_particle_velocities)
{
    T_GRID& grid=*fluids_parameters.grid;
    T maximum_particle_speed=0;
    if(fluids_parameters.fire){
        if(fluids_parameters.number_of_regions==1){
            Get_Levelset_Velocity(*fluids_parameters.grid,fluids_parameters.particle_levelset_evolution->particle_levelset.levelset,
                fluids_parameters.particle_levelset_evolution->V,time);
            maximum_fluid_speed=max(maximum_fluid_speed,fluids_parameters.particle_levelset_evolution->V.Maxabs().Max());}
        else{
            Get_Levelset_Velocity(*fluids_parameters.grid,fluids_parameters.particle_levelset_evolution_multiple->Levelset_Multiple(),
                fluids_parameters.particle_levelset_evolution_multiple->V,time);
            maximum_fluid_speed=max(maximum_fluid_speed,fluids_parameters.particle_levelset_evolution_multiple->V.Maxabs().Max());}}
    if(include_removed_particle_velocities){
        for(int i=1;i<=fluids_parameters.number_of_regions;i++){
            PARTICLE_LEVELSET_UNIFORM<T_GRID>& particle_levelset=fluids_parameters.particle_levelset_evolution->Particle_Levelset(i);
            if(particle_levelset.use_removed_negative_particles) for(CELL_ITERATOR iterator(particle_levelset.levelset.grid);iterator.Valid();iterator.Next()){
                PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles=particle_levelset.removed_negative_particles(iterator.Cell_Index());
                if(particles) maximum_particle_speed=max(maximum_particle_speed,ARRAYS_COMPUTATIONS::Maximum_Magnitude(particles->V));}
            if(particle_levelset.use_removed_positive_particles) for(CELL_ITERATOR iterator(particle_levelset.levelset.grid);iterator.Valid();iterator.Next()){
                PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles=particle_levelset.removed_positive_particles(iterator.Cell_Index());
                if(particles) maximum_particle_speed=max(maximum_particle_speed,ARRAYS_COMPUTATIONS::Maximum_Magnitude(particles->V));}}}
    T max_particle_collision_distance=0;
    for(int i=1;i<=fluids_parameters.number_of_regions;i++)
        max_particle_collision_distance=max(max_particle_collision_distance,fluids_parameters.particle_levelset_evolution->Particle_Levelset(i).max_collision_distance_factor*grid.dX.Max());
    fluids_parameters.collision_bodies_affecting_fluid->Compute_Occupied_Blocks(true,
        dt*max(maximum_fluid_speed,maximum_particle_speed)+2*max_particle_collision_distance+(T).5*fluids_parameters.p_grid.dX.Max(),10);

    if(fluids_parameters.use_maccormack_semi_lagrangian_advection && fluids_parameters.use_maccormack_compute_mask){
        typedef typename T_GRID::VECTOR_INT TV_INT;
        fluids_parameters.maccormack_cell_mask.Resize(grid.Domain_Indices(fluids_parameters.number_of_ghost_cells),false,false);
        fluids_parameters.maccormack_face_mask.Resize(grid,fluids_parameters.number_of_ghost_cells);
        // don't use maccormack near domain boundary conditions
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
            fluids_parameters.maccormack_cell_mask(iterator.Cell_Index())=TV_INT::Componentwise_Min(iterator.Cell_Index()-grid.Domain_Indices().Minimum_Corner(),
                grid.Domain_Indices().Maximum_Corner()-iterator.Cell_Index()).Min()>fluids_parameters.cfl;
        VECTOR<T_GRID,T_GRID::dimension> grids;
        for(int i=1;i<=T_GRID::dimension;i++) grids[i]=grid.Get_Face_Grid(i);
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            int axis=iterator.Axis();TV_INT index=iterator.Face_Index();
            fluids_parameters.maccormack_face_mask(axis,index)
                =TV_INT::Componentwise_Min(index-grids[axis].Domain_Indices().Minimum_Corner(),grids[axis].Domain_Indices().Maximum_Corner()-index).Min()>fluids_parameters.cfl;}
        if(fluids_parameters.mpi_grid){ // if mpi check turned off domain boundary regions to make sure they are global domain boundaries
            RANGE<TV> global_domain=fluids_parameters.mpi_grid->global_grid.domain;T double_cfl_dx=(T)fluids_parameters.cfl*fluids_parameters.p_grid.dX.Max();
            for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
                if(!fluids_parameters.maccormack_cell_mask(iterator.Cell_Index()) && global_domain.Inside(iterator.Location(),(T)double_cfl_dx))
                    fluids_parameters.maccormack_cell_mask(iterator.Cell_Index())=true;
            for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                int axis=iterator.Axis();TV_INT index=iterator.Face_Index();
                if(!fluids_parameters.maccormack_face_mask(axis,index) && global_domain.Inside(iterator.Location(),(T)double_cfl_dx))
                    fluids_parameters.maccormack_face_mask(axis,index)=true;}}
        // don't use maccormack near the interface
        if(fluids_parameters.bandwidth_without_maccormack_near_interface){
            T dx_times_band=grid.Maximum_Edge_Length()*fluids_parameters.bandwidth_without_maccormack_near_interface;
            for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next())if(fluids_parameters.particle_levelset_evolution->phi(iterator.Cell_Index())>-dx_times_band){
                fluids_parameters.maccormack_cell_mask(iterator.Cell_Index())=false;
                for(int axis=1;axis<=T_GRID::dimension;axis++){
                    fluids_parameters.maccormack_face_mask(axis,iterator.First_Face_Index(axis))=false;fluids_parameters.maccormack_face_mask(axis,iterator.Second_Face_Index(axis))=false;}}}
        // turn off maccormack near objects
        for(typename T_GRID::NODE_ITERATOR node_iterator(grid,2);node_iterator.Valid();node_iterator.Next()) 
            if(fluids_parameters.collision_bodies_affecting_fluid->occupied_blocks(node_iterator.Node_Index())){
                TV_INT block_index=node_iterator.Node_Index();BLOCK_UNIFORM<T_GRID> block(grid,block_index);
                for(int cell_index=0;cell_index<T_GRID::number_of_cells_per_block;cell_index++) fluids_parameters.maccormack_cell_mask(block.Cell(cell_index))=false;
                for(int axis=1;axis<=T_GRID::dimension;axis++) for(int face=0;face<T_GRID::number_of_incident_faces_per_block;face++)
                    fluids_parameters.maccormack_face_mask(axis,block.Incident_Face(axis,face))=false;}}
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Delete_Particles_Inside_Objects(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time)
{}
//#####################################################################
// Function Read_Output_Files_Fluids
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Read_Output_Files_Fluids(const int frame)
{
    fluids_parameters.Read_Output_Files(stream_type,output_directory,frame);
    incompressible_fluid_container.Read_Output_Files(stream_type,output_directory,frame);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    if(fluids_parameters.smoke||fluids_parameters.fire||fluids_parameters.water){
        if(fluids_parameters.solid_affects_fluid && fluids_parameters.fluid_affects_solid){std::string filename;
            /*
            if(fluids_parameters.smoke && fluids_parameters.use_density && fluids_parameters.semi_lagrangian_collidable_density.use_valid_mask){
                filename=output_directory+"/density_valid_mask."+f;
                if(FILE_UTILITIES::File_Exists(filename)){
                    LOG::cout<<"Reading "<<filename<<std::endl;
                    FILE_UTILITIES::Read_From_File(stream_type,filename,fluids_parameters.semi_lagrangian_collidable_density.valid_points_current);}}
            if(fluids_parameters.smoke && fluids_parameters.use_temperature && fluids_parameters.semi_lagrangian_collidable_temperature.use_valid_mask){
                filename=output_directory+"/temperature_valid_mask."+f;
                if(FILE_UTILITIES::File_Exists(filename)){
                    LOG::cout<<"Reading "<<filename<<std::endl;
                    FILE_UTILITIES::Read_From_File(stream_type,filename,fluids_parameters.semi_lagrangian_collidable_temperature.valid_points_current);}}
            if(fluids_parameters.semi_lagrangian_collidable_velocity.use_valid_mask){
                filename=output_directory+"/velocity_valid_mask."+f;
                if(FILE_UTILITIES::File_Exists(filename)){
                    LOG::cout<<"Reading "<<filename<<std::endl;
                    FILE_UTILITIES::Read_From_File(stream_type,filename,fluids_parameters.semi_lagrangian_collidable_velocity.valid_points_current);}}
            if((fluids_parameters.water || fluids_parameters.fire) && fluids_parameters.semi_lagrangian_collidable_phi.use_valid_mask){
                filename=output_directory+"/phi_valid_mask."+f;
                if(FILE_UTILITIES::File_Exists(filename)){
                    LOG::cout<<"Reading "<<filename<<std::endl;
                    FILE_UTILITIES::Read_From_File(stream_type,filename,fluids_parameters.semi_lagrangian_collidable_phi.valid_points_current);}}*/}}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Write_Output_Files(const int frame) const
{
    FILE_UTILITIES::Create_Directory(output_directory);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    Write_Frame_Title(frame);
    solid_body_collection.Write(stream_type,output_directory,frame,first_frame,solids_parameters.write_static_variables_every_frame,
        solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies,solids_parameters.write_deformable_body,solids_parameters.write_from_every_process,
        solids_parameters.triangle_collision_parameters.output_interaction_pairs);
    if(NEWMARK_EVOLUTION<TV>* newmark=dynamic_cast<NEWMARK_EVOLUTION<TV>*>(solids_evolution))
        newmark->Write_Position_Update_Projection_Data(stream_type,output_directory+"/"+f+"/");
    
    fluids_parameters.Write_Output_Files(stream_type,output_directory,first_frame,frame);
    if(fluids_parameters.incompressible) incompressible_fluid_container.Write_Output_Files(stream_type,output_directory,frame);
    if(fluids_parameters.smoke||fluids_parameters.fire||fluids_parameters.water){
        if(fluids_parameters.solid_affects_fluid && fluids_parameters.fluid_affects_solid){
            if(fluids_parameters.write_debug_data){
                /*
                FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/thin_shells_grid_visibility."+f,fluids_parameters.collision_bodies_affecting_fluid->cell_neighbors_visible,
                    fluids_parameters.collision_bodies_affecting_fluid->face_corners_visible_from_face_center);
                if(fluids_parameters.smoke && fluids_parameters.use_density && fluids_parameters.semi_lagrangian_collidable_density.use_valid_mask)
                    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/density_valid_mask."+f,fluids_parameters.semi_lagrangian_collidable_density.valid_points_current);
                if(fluids_parameters.smoke && fluids_parameters.use_temperature && fluids_parameters.semi_lagrangian_collidable_temperature.use_valid_mask)
                    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/temperature_valid_mask."+f,fluids_parameters.semi_lagrangian_collidable_temperature.valid_points_current);
                if(fluids_parameters.semi_lagrangian_collidable_velocity.use_valid_mask) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/velocity_valid_mask."+f,
                    fluids_parameters.semi_lagrangian_collidable_velocity.valid_points_current);
                if((fluids_parameters.water || fluids_parameters.fire) && fluids_parameters.semi_lagrangian_collidable_phi.use_valid_mask)
                    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/phi_valid_mask."+f,fluids_parameters.semi_lagrangian_collidable_phi.valid_points_current);*/}}}
}
//#####################################################################
// Function Get_Psi_D_Inside_Solids
//#####################################################################
template<class T_GRID> bool SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Get_Psi_D_Inside_Solids(T_ARRAYS_BOOL& psi_D)
{
    if(fluids_parameters.euler_solid_fluid_coupling_utilities && fluids_parameters.euler_solid_fluid_coupling_utilities->thinshell){
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(fluids_parameters.euler->grid,psi_D.Domain_Indices());iterator.Valid();iterator.Next()){const TV_INT& cell_index=iterator.Cell_Index();
            if(fluids_parameters.euler_solid_fluid_coupling_utilities->psi_np1.Valid_Index(cell_index))
                psi_D(cell_index)=!fluids_parameters.euler_solid_fluid_coupling_utilities->psi_np1(cell_index); else psi_D(cell_index)=true;}
        return true;}
    return false;
}
//#####################################################################
// Function Log_Parameters 
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>::
Log_Parameters() const
{       
    LOG::SCOPE scope("SOLIDS_FLUIDS_EXAMPLE_UNIFORM parameters");
    BASE::Log_Parameters();
    fluids_parameters.Log_Parameters();
}
//#####################################################################
#define PROTECT(...) __VA_ARGS__
#define INSTANTIATION_HELPER_T_d_SHAPE(T,d,SHAPE) \
    template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T,d> > >::Get_Source_Velocities(const SHAPE&,const MATRIX<T,d+1>&,const VECTOR<T,d>&); \
    template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T,d> > >::Get_Source_Velocities(const SHAPE&,const MATRIX<T,d+1>&,const VECTOR<T,d>&,const T_FACE_ARRAYS_BOOL&); \
    template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T,d> > >::Adjust_Phi_With_Source(const SHAPE&,const MATRIX<T,d+1>&); \
    template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T,d> > >::Adjust_Phi_With_Source(const SHAPE&,const int,const MATRIX<T,d+1>&); \
    template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T,d> > >::Get_Source_Reseed_Mask(const SHAPE& source,const MATRIX<T,d+1>&,T_ARRAYS_BOOL*&,const bool); \
    template void SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T,d> > >::Adjust_Density_And_Temperature_With_Sources(const SHAPE&,const MATRIX<T,d+1>&,const T,const T);
#define INSTANTIATION_HELPER_SHAPE(...) INSTANTIATION_HELPER_T_d_SHAPE(PROTECT(__VA_ARGS__::VECTOR_T::SCALAR),PROTECT(__VA_ARGS__::VECTOR_T::m),PROTECT(__VA_ARGS__))
#define INSTANTIATION_HELPER(T) \
    template class SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T,1> > >;   \
    template class SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T,2> > >;    \
    template class SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T,3> > >;   \
    INSTANTIATION_HELPER_SHAPE(BOX<VECTOR<T,2> >) \
    INSTANTIATION_HELPER_SHAPE(BOX<VECTOR<T,3> >) \
    INSTANTIATION_HELPER_SHAPE(SPHERE<VECTOR<T,2> >) \
    INSTANTIATION_HELPER_SHAPE(SPHERE<VECTOR<T,3> >) \
    INSTANTIATION_HELPER_SHAPE(CYLINDER<T>)
INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
