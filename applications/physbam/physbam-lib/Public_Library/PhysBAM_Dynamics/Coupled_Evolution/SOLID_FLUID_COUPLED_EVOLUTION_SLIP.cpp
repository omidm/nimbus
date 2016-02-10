//#####################################################################
// Copyright 2008-2009, Elliot English, Jon Gretarsson, Nipun Kwatra, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Krylov_Solvers/LANCZOS_ITERATION.h>
#include <PhysBAM_Tools/Krylov_Solvers/SYMMQMR.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Math_Tools/Inverse.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_ND.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_POINT_SIMPLICES_1D_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_SEGMENTED_CURVE_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_TRIANGULATED_SURFACE_INTERSECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SURFACE_TENSION_FORCE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Fracture/FRACTURE_PATTERN.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGIDS_NEWMARK_COLLISION_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_PROJECTION_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Grid_Based_Fields/INCOMPRESSIBLE_FLUID_CONTAINER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/VISCOSITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/LEVELSET_VISCOSITY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COLLISION_AWARE_INDEX_MAP.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COMPRESSIBLE_BOUNDARY_CONDITION_WALLS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COUPLED_SYSTEM_VECTOR.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/EXAMPLE_BOUNDARY_CONDITION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION_CUT.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION_COLLECTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION_COLLISIONS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION_SOURCES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_SOLID_INTERPOLATION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SURFACE_TENSION_BOUNDARY_CONDITION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_CALLBACKS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_CALLBACKS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_PARAMETERS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
SOLID_FLUID_COUPLED_EVOLUTION_SLIP(SOLIDS_PARAMETERS<TV>& solids_parameters_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,
    FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters_input,SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters_input,INCOMPRESSIBLE_FLUID_CONTAINER<GRID<TV> >& incompressible_fluid_container_input)
    :NEWMARK_EVOLUTION<TV>(solids_parameters_input,solid_body_collection_input),
    PROJECTION_DYNAMICS_UNIFORM<T_GRID>(*fluids_parameters_input.grid,fluids_parameters_input.fire,false,false,fluids_parameters_input.use_poisson),
    collision_bodies(*fluids_parameters_input.collision_bodies_affecting_fluid),
    fluids_parameters(fluids_parameters_input),solids_fluids_parameters(solids_fluids_parameters_input),incompressible_fluid_container(incompressible_fluid_container_input),
    iterator_info(*new UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>(collision_bodies)),
    boundary_condition_collection(*new IMPLICIT_BOUNDARY_CONDITION_COLLECTION<TV>(fluids_parameters.callbacks,!Simulate_Compressible_Fluids(),!Simulate_Compressible_Fluids(),
            fluids_parameters_input.use_psi_R,fluids_parameters.use_second_order_pressure,fluids_parameters.periodic_boundary)),
    grid(*fluids_parameters.grid),fracture_pattern(new FRACTURE_PATTERN<T>),number_of_coupling_faces_cached(0),time_cached((T)-1),cached_coupling_face_data(false),
    run_self_tests(false),print_poisson_matrix(false),print_index_map(false),print_rhs(false),print_each_matrix(false),output_iterators(false),
    use_viscous_forces(false),two_phase(false),use_full_ic(false)
{
    disable_thinshell=false;
    Initialize_Grid_Arrays();

    if(fluids_parameters.mpi_grid){
        VECTOR<VECTOR<bool,2>,TV::dimension> mpi_boundaries;mpi_boundaries.Fill(VECTOR<bool,2>(false,false));
        for(int axis=1;axis<=TV::dimension;++axis) for(int axis_side=1;axis_side<=2;++axis_side)
            if(fluids_parameters.mpi_grid->side_neighbor_ranks(2*(axis-1)+axis_side)>=0) mpi_boundaries(axis)(axis_side)=true;
        boundary_condition_collection.Set_Mpi_Boundaries(mpi_boundaries);}

    boundary_condition_collection.Add_Boundary_Condition(new IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<TV>(collision_bodies.collision_geometry_collection,false,
        *fluids_parameters.callbacks));
    boundary_condition_collection.Add_Boundary_Condition(new IMPLICIT_BOUNDARY_CONDITION_SOURCES<TV>(*fluids_parameters.callbacks));
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
~SOLID_FLUID_COUPLED_EVOLUTION_SLIP()
{
    delete &boundary_condition_collection;
    delete fracture_pattern;
    delete &iterator_info;
}
//#####################################################################
// Function Initialize_Grid_Arrays
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Initialize_Grid_Arrays()
{
    fluids_face_velocities.Resize(grid);
    pressure.Resize(grid.Domain_Indices(1));
    density.Resize(grid.Domain_Indices(1));
    if(Simulate_Compressible_Fluids()) centered_velocity.Resize(grid.Domain_Indices(1));
    one_over_rho_c_squared.Resize(grid.Domain_Indices(1));
    p_advected_over_rho_c_squared_dt.Resize(grid.Domain_Indices());
    p_advected.Resize(grid.Domain_Indices());
    solved_faces.Resize(grid);
}
//#####################################################################
// Function Backward_Euler_Step_Velocity_Helper
//#####################################################################
// assumes all solids_forces are linear in velocity, with a symmetric positive definite Jacobian.
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Backward_Euler_Step_Velocity_Helper(const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update)
{
    Solve(incompressible_fluid_container.face_velocities,dt,current_velocity_time,current_position_time,velocity_update,false);
    if(!velocity_update && solids_fluids_parameters.mpi_solid_fluid && solids_fluids_parameters.mpi_solid_fluid->Fluid_Node())
        Process_Collisions(dt,current_position_time,false);
}
//#####################################################################
// Function Process_Collisions
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Process_Collisions(const T dt,const T time,const bool advance_rigid_bodies)
{
    ARRAY<int> all_added_bodies(0),remove_body_particle_ids(0);
    ARRAY<COLLISION_GEOMETRY_ID> removed_bodies(0);
    if(solids_fluids_parameters.use_fluid_rigid_fracture){
        for(COLLISION_GEOMETRY_ID p(1);p<=fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Size();++p){
            ARRAY<int> added_bodies(0);
            if(RIGID_COLLISION_GEOMETRY<TV>* rigid_body_wrapper=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.bodies(p))){
                int rigid_particle_index=rigid_body_wrapper->rigid_geometry.particle_index;
                if(!solid_body_collection.rigid_body_collection.Is_Active(rigid_particle_index)){
                    LOG::filecout("THIS DOESN'T MAKE SENSE\n");
                    DEBUG_UTILITIES::Debug_Breakpoint();
                    continue;}
                RIGID_BODY<TV>& body=solid_body_collection.rigid_body_collection.Rigid_Body(rigid_particle_index);
                if(pressure_impulses_twist.Valid_Index(rigid_particle_index) && pressure_impulses_twist(rigid_particle_index).linear.Magnitude() > body.fracture_threshold){ // Valid_Index because we might be trying to fracture objects added in this substep.
                    std::stringstream ss;
                    ss<<"Maximum magnitude threshold exceeded: "<<pressure_impulses_twist(rigid_particle_index).linear.Magnitude()<<std::endl;
                    LOG::filecout(ss.str());

                    TV object_space_collision_location;
                    {
                        RAY<TV> intersection_ray(body.Object_Space_Point(body.X()),-body.Object_Space_Point(pressure_impulses_twist(rigid_particle_index).linear));
                        body.simplicial_object->Initialize_Hierarchy();body.simplicial_object->Refresh_Auxiliary_Structures();
                        if(!INTERSECTION::Intersects(intersection_ray,*body.simplicial_object,(T)1e-6)){LOG::filecout("Unable to find a point of impact\n");continue;}
                        object_space_collision_location=intersection_ray.Point(intersection_ray.t_max); // already in object space
                    }

                    if(!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Solid_Node())
                        rigids_evolution_callbacks.Begin_Fracture(rigid_particle_index);
                    fracture_pattern->Intersect_With_Rigid_Body(body,body.World_Space_Point(object_space_collision_location),added_bodies,solids_parameters.rigid_body_collision_parameters.allow_refracturing,
                                                                solids_parameters.rigid_body_collision_parameters.use_fracture_particle_optimization,true);

                    std::stringstream ss1;
                    ss1<<"Added "<<added_bodies.m<<" to the simulation"<<std::endl;
                    LOG::filecout(ss1.str());
                    removed_bodies.Append(solid_body_collection.rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(rigid_particle_index));
                    solid_body_collection.rigid_body_collection.rigid_body_particle.Remove_Body(rigid_particle_index);
                    for(int i=1;i<=added_bodies.m;++i){
                        COLLISION_GEOMETRY<TV>* collision_body=&(*solid_body_collection.rigid_body_collection.rigid_geometry_collection.collision_body_list)(
                                solid_body_collection.rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(added_bodies(i)));
                        fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(collision_body,added_bodies(i),false);
                        collision_body->Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE,time+dt);}

                    if(!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Solid_Node())
                        rigids_evolution_callbacks.End_Fracture(rigid_particle_index,added_bodies);
                    for(int i=1;i<=added_bodies.m;++i){
                        COLLISION_GEOMETRY<TV>* collision_body=&(*solid_body_collection.rigid_body_collection.rigid_geometry_collection.collision_body_list)(
                                solid_body_collection.rigid_body_collection.rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(added_bodies(i)));
                        collision_body->Save_State(COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE,time);}

                    all_added_bodies.Append_Elements(added_bodies);}}}

        if(all_added_bodies.m){
            for(int i=1;i<=removed_bodies.m;++i)
                fluids_parameters.collision_bodies_affecting_fluid->Remove_Body(removed_bodies(i));

            solid_body_collection.rigid_body_collection.rigid_geometry_collection.Destroy_Unreferenced_Geometry();
            solid_body_collection.rigid_body_collection.Update_Simulated_Particles();
            rigid_body_collisions->Initialize_Data_Structures(true);
            for(int i=1;i<=all_added_bodies.m;i++){
                rigid_body_collisions->added_bodies(1).Append(all_added_bodies(i));
                if(!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Solid_Node())
                    rigid_body_collisions->collision_callbacks.Euler_Step_Position_With_New_Velocity(all_added_bodies(i),dt,time);}
            if(solids_fluids_parameters.mpi_solid_fluid) solids_fluids_parameters.mpi_solid_fluid->Exchange_Solid_Positions_And_Velocities(solid_body_collection);
            rigid_body_collisions->spatial_partition->Reinitialize();}}

    if(!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Solid_Node()) BASE::Process_Collisions(dt,time,advance_rigid_bodies);
}
//#####################################################################
// Function Apply_Pressure
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Apply_Pressure(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,bool scale_by_dt)
{
}
//#####################################################################
// Function Setup_Solids
//#####################################################################
template<class TV> BACKWARD_EULER_SYSTEM<TV>* SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Setup_Solids(const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update,const bool leakproof_solve)
{
    PHYSBAM_ASSERT(!solids_parameters.use_trapezoidal_rule_for_velocities); // two-way coupling does not work with trapezoidal rule
    ARRAY<int> empty_list;
    MPI_SOLIDS<TV>* mpi_solids=solid_body_collection.deformable_body_collection.mpi_solids;
    GENERALIZED_VELOCITY<TV> B(B_full,solid_body_collection.deformable_body_collection.dynamic_particles,rigid_B_full,
        solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles,empty_list);
    //bool solids=Simulate_Solids() && (!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Solid_Node());

    if(!leakproof_solve){
        //if(!solids) return 0; // do nothing
        BACKWARD_EULER_SYSTEM<TV>* system=new BACKWARD_EULER_SYSTEM<TV>(*this,solid_body_collection,dt,current_velocity_time,current_position_time,
            &solid_body_collection.rigid_body_collection.articulated_rigid_body,(velocity_update && solids_parameters.enforce_repulsions_in_cg)?repulsions:0,mpi_solids,velocity_update);
        Prepare_Backward_Euler_System(*system,dt,current_velocity_time,current_position_time,velocity_update);
        if(solids_fluids_parameters.mpi_solid_fluid) solids_fluids_parameters.mpi_solid_fluid->Exchange_Solid_Positions_And_Velocities(solid_body_collection); // TODO: only need to exchange velocities
        return system;}

    // use F so we don't change the actual solid velocities
    DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>* collisions=0;
    for(COLLISION_GEOMETRY_ID i(1);i<=fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.bodies.m;i++)
        if(fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Is_Active(i)){
            collisions=dynamic_cast<DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>*>(fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.bodies(i));if(collisions) break;}

    if(solids_fluids_parameters.mpi_solid_fluid && solids_fluids_parameters.mpi_solid_fluid->Fluid_Node()) return 0;

    // TODO: make sure that V.V.Size() is in fact indexed only over the dynamic particles
    int state1=COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE;int state2=COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE;
    if(collisions) for(int i=1;i<=B.V.Size();i++)
        B.V(i)=collisions->Pointwise_Node_Pseudo_Velocity(B.V.indices(i),state1,state2);
    for(int i=1;i<=B.rigid_V.Size();i++){
        RIGID_COLLISION_GEOMETRY<TV>* collision_geometry=
            dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Get_Collision_Geometry(B.rigid_V.indices(i)));
        assert(collision_geometry);
        T one_over_dt=1/(collision_geometry->saved_states(state2).y-collision_geometry->saved_states(state1).y);
        FRAME<TV> a=collision_geometry->saved_states(state1).x,b=collision_geometry->saved_states(state2).x;
        B.rigid_V(i)=one_over_dt*TWIST<TV>(b.t-a.t,(b.r*a.r.Inverse()).Rotation_Vector());}
    return 0;
}
//#####################################################################
// Function Setup_Fluids
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Setup_Fluids(T_FACE_ARRAYS_SCALAR& incompressible_face_velocities,const T current_position_time,const T dt,const bool leakproof_solve)
{
    // TODO(kwatra): should we call with both incompressible and compressible face velocities to make sure psi_N faces are set correctly.
    if(solids_fluids_parameters.mpi_solid_fluid && solids_fluids_parameters.mpi_solid_fluid->Solid_Node()) return;
    const T one_over_cell_size=grid.One_Over_Cell_Size();
    boundary_condition_collection.periodic_boundary=fluids_parameters.periodic_boundary;
    boundary_condition_collection.Compute(grid,pressure,fluids_face_velocities,current_position_time);
    EULER_PROJECTION_UNIFORM<T_GRID>& euler_projection=fluids_parameters.euler->euler_projection;
    // TODO: in case of both compressible and incompressible fluids, use a levelset to decide which faces to fill with what.
    one_over_rho_c_squared.Fill(0);p_advected_over_rho_c_squared_dt.Fill(0);p_advected=pressure;
    if(Simulate_Incompressible_Fluids()){
        fluids_face_velocities=incompressible_face_velocities;
        if(fluids_parameters.use_density){
            fluids_parameters.density_container.Get_Ghost_Density(dt,current_position_time,1,density);}
        else if(!two_phase) density.Fill(fluids_parameters.density);
        else for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next())
            density(it.index)=boundary_condition_collection.psi_D(it.index)?fluids_parameters.outside_density:fluids_parameters.density;}
    if(Simulate_Compressible_Fluids()){
        euler_projection.Compute_Density_Weighted_Face_Velocities(dt,current_position_time,boundary_condition_collection.psi_N);
        fluids_face_velocities=euler_projection.face_velocities; // TODO(kwatra): fill only compressible faces
        euler_projection.Get_Ghost_Density(dt,current_position_time,1,density);
        euler_projection.Get_Ghost_Centered_Velocity(dt,current_position_time,1,centered_velocity);
        one_over_rho_c_squared=euler_projection.one_over_rho_c_squared;
        const T one_over_dt=(T)1/dt;
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            p_advected_over_rho_c_squared_dt(cell_index)=euler_projection.p_advected(cell_index)*one_over_rho_c_squared(cell_index)*one_over_dt;
            p_advected(cell_index)=euler_projection.p_advected(cell_index);}
    
        for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            euler_projection.density_scaling(cell_index) = fluids_parameters.euler_solid_fluid_coupling_utilities->cell_volumes_np1(cell_index) * one_over_cell_size;
            density(cell_index) *= euler_projection.density_scaling(cell_index);}}
}
//#####################################################################
// Function Solve
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Solve(T_FACE_ARRAYS_SCALAR& incompressible_face_velocities,const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update,const bool leakproof_solve)
{
    static int solve_id=0;solve_id++;
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particle;
    EULER_PROJECTION_UNIFORM<T_GRID>& euler_projection=fluids_parameters.euler->euler_projection;
    bool solids=Simulate_Solids();
    bool fluids=Simulate_Fluids();

    ARRAY<int> empty_list;

    if(solids_fluids_parameters.mpi_solid_fluid){
        PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
        V_save.Resize(particles.array_collection->Size(),false,false);
        rigid_velocity_save.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);}
    GENERALIZED_VELOCITY<TV> V_n(V_save,rigid_velocity_save,solid_body_collection);
    if(solids_fluids_parameters.mpi_solid_fluid) solids_fluids_parameters.mpi_solid_fluid->Distribute_Lists_From_Solid_Node(V_n);

    ARRAY<TWIST<TV> > twist;twist.Resize(rigid_body_particles.V.Size());
    for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(rigid_body_particles.V(i),rigid_body_particles.angular_velocity(i));

    SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>* coupled_system=new SYMMETRIC_POSITIVE_DEFINITE_COUPLING_SYSTEM<TV>(fluids_parameters.use_preconditioner_for_slip_system,
        0,boundary_condition_collection,iterator_info,solid_body_collection,incompressible_fluid_container,density,centered_velocity,one_over_rho_c_squared,leakproof_solve,false,
        fluids_parameters.use_coupled_implicit_viscosity);
    coupled_system->index_map.two_phase=two_phase;
    fluids_parameters.callbacks->Substitute_Coupling_Matrices(*coupled_system,dt,current_velocity_time,current_position_time,velocity_update,leakproof_solve);
    BACKWARD_EULER_SYSTEM<TV>* solid_system=Setup_Solids(dt,current_velocity_time,current_position_time,velocity_update,leakproof_solve);
    coupled_system->solid_system=solid_system;
    GENERALIZED_VELOCITY<TV> V(particles.V,twist,solid_body_collection),B(B_full,rigid_B_full,solid_body_collection),F(F_full,rigid_F_full,solid_body_collection);
    if(!leakproof_solve)
        if(solids_fluids_parameters.mpi_solid_fluid)
            solids_fluids_parameters.mpi_solid_fluid->Exchange_Solid_Positions_And_Velocities(solid_body_collection); // TODO: only need to exchange velocities
    if(fluids){
        Setup_Fluids(incompressible_face_velocities,current_position_time,dt,leakproof_solve);}

    if(coupled_system->fluid_to_solid_interpolation)
        if(FLUID_TO_SOLID_INTERPOLATION_CUT<TV>* temp_interpolate=dynamic_cast<FLUID_TO_SOLID_INTERPOLATION_CUT<TV>*>(coupled_system->fluid_to_solid_interpolation)){
            temp_interpolate->Setup_Before_Compute(coupled_system->index_map.boundary_condition_collection.psi_D,coupled_system->index_map.boundary_condition_collection.psi_N);
            if(SURFACE_TENSION_FORCE<VECTOR<T,2> >* force=solid_body_collection.deformable_body_collection.template Find_Force<SURFACE_TENSION_FORCE<VECTOR<T,2> >*>()) force->Dump_Curvatures();}

    coupled_system->debug_velocity=&incompressible_face_velocities;
    coupled_system->debug_generalized_velocity=&V;

    if(solids_fluids_parameters.mpi_solid_fluid) coupled_system->Set_MPI(*solids_fluids_parameters.mpi_solid_fluid,*fluids_parameters.mpi_grid);

    divergence_scaling=1;//sqrt(fluids_parameters.density*fluids_parameters.grid->Cell_Size())/fluids_parameters.grid->Face_Size(1);

    coupled_system->run_self_tests=run_self_tests;
    coupled_system->print_poisson_matrix=print_poisson_matrix;
    coupled_system->print_index_map=print_index_map;
    coupled_system->print_matrix=print_matrix;
    coupled_system->print_rhs=print_rhs;
    coupled_system->print_each_matrix=print_each_matrix;
    coupled_system->surface_tension_coefficient=fluids_parameters.surface_tension;
    coupled_system->use_full_ic=use_full_ic;
    coupled_system->Compute(solids_fluids_parameters.mpi_solid_fluid?1:0,dt,current_velocity_time,boundary_condition_collection.psi_N,disable_thinshell,Simulate_Compressible_Fluids(),
        fluids_face_velocities,fluids_parameters.viscosity,fluids_parameters.particle_levelset_evolution && fluids_parameters.second_order_cut_cell_method,
        fluids_parameters.particle_levelset_evolution?&fluids_parameters.particle_levelset_evolution->Levelset(1):0);

    coupled_system->Resize_Coupled_System_Vector(coupled_f);
    coupled_system->Resize_Coupled_System_Vector(coupled_r);
    coupled_system->Resize_Coupled_System_Vector(coupled_s);
    coupled_system->Resize_Coupled_System_Vector(coupled_ar);
    coupled_system->Resize_Coupled_System_Vector(coupled_z);
    if(solids_fluids_parameters.mpi_solid_fluid) solids_fluids_parameters.mpi_solid_fluid->Distribute_Lists_From_Solid_Node(F);
    if(solids_fluids_parameters.mpi_solid_fluid) solids_fluids_parameters.mpi_solid_fluid->Distribute_Lists_From_Solid_Node(B);
    // TODO(kwatra): switch to using B instead of V once use_post_cg_constraints work
    coupled_system->Set_Up_RHS(coupled_x,coupled_b,B,fluids_face_velocities,p_advected_over_rho_c_squared_dt,p_advected,pressure);
    //if(Simulate_Compressible_Fluids()) Set_Cached_Psi_N_And_Coupled_Face_Data(coupled_system->index_map,coupled_system->solid_interpolation,current_position_time); // dont cache. As velocities can change in the position update from collisions. Could cache solid_interpolation matrix.

    // fill psi_N,psi_D for viewing purposes
    if(Simulate_Compressible_Fluids() && (!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Fluid_Node())){
        euler_projection.elliptic_solver->psi_N=boundary_condition_collection.psi_N;
        coupled_system->Set_Coupling_Faces(0,euler_projection.elliptic_solver->psi_N);
        euler_projection.elliptic_solver->psi_D=boundary_condition_collection.psi_D;
        Warn_For_Exposed_Dirichlet_Cell(euler_projection.elliptic_solver->psi_D,euler_projection.elliptic_solver->psi_N);}

    T fluid_tolerance=(T)1e-6*divergence_scaling*fluids_parameters.grid->Cell_Size();
    if(!fluid_tolerance) fluid_tolerance=solids_parameters.implicit_solve_parameters.cg_tolerance;
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Before Cholesky factorization",0,1);

    int max_iterations=solids_parameters.implicit_solve_parameters.cg_iterations;
    static CONJUGATE_RESIDUAL<T> cr;
    static SYMMQMR<T> symmqmr;
    static CONJUGATE_GRADIENT<T> cg;
    KRYLOV_SOLVER<T>* solver=0;
    const char* solver_name=0;
    if(solids_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_cg){solver=&cg;solver_name="CG";}
    else if(solids_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_cr){solver=&cr;solver_name="CONJUGATE_RESIDUAL";}
    else if(solids_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_symmqmr){solver=&symmqmr;solver_name="SYMMQMR";}
    LOG::SCOPE solve_scope(solver_name);

    solver->nullspace_tolerance=(T)-1;
    solver->restart_iterations=solids_parameters.implicit_solve_parameters.cg_restart_iterations;
    solver->print_diagnostics=solid_body_collection.print_diagnostics;
    solver->print_residuals=solid_body_collection.print_residuals;
    solver->iterations_used=&solid_body_collection.iterations_used_diagnostic;
    solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;

    if(print_each_matrix) OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("b-%i.txt",solve_id).c_str()).Write("b",coupled_b);    // TODO: this isn't actually valid for all the solver types
    if(!solver->Solve(*coupled_system,coupled_x,coupled_b,coupled_f,coupled_s,coupled_r,coupled_ar,coupled_z,fluid_tolerance,1,max_iterations))
        PHYSBAM_DEBUG_WRITE_SUBSTEP("FAILED CONVERGENCE",0,1);
    if(print_each_matrix) OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("x-%i.txt",solve_id).c_str()).Write("x",coupled_x);
    std::stringstream ss;
    ss<<"Residual L_inf norm="<<coupled_system->Residual_Linf_Norm(coupled_x,coupled_b)<<std::endl;
    LOG::filecout(ss.str());

    LOG::Stop_Time();
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After coupled solve (sf coupled evolution)",0,1);

    bool want_fluid=fluids && (!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Fluid_Node()) && Simulate_Incompressible_Fluids();
    coupled_system->Apply_Velocity_Update(coupled_x,incompressible_face_velocities,pressure,V,F,solids && !leakproof_solve,want_fluid);

    if(solids && !leakproof_solve){
        pressure_impulses=coupled_system->pressure_impulses;
        pressure_impulses_twist=coupled_system->pressure_impulses_twist;}
        // TODO: only need to exchange velocities
    if(solids_fluids_parameters.mpi_solid_fluid) solids_fluids_parameters.mpi_solid_fluid->Exchange_Solid_Positions_And_Velocities(solid_body_collection);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After solid velocity update (sf coupled evolution)",0,1);

    if(fluids && (!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Fluid_Node()) && Simulate_Compressible_Fluids()){
        ARRAY<T,COUPLING_CONSTRAINT_ID> coupled_faces_solid_interpolated_velocity_n;
        ARRAY<T,COUPLING_CONSTRAINT_ID> coupled_faces_solid_interpolated_velocity_np1;
        coupled_system->Interpolate_Solid_Velocity_To_Coupled_Faces(V_n,coupled_faces_solid_interpolated_velocity_n);
        coupled_system->Interpolate_Solid_Velocity_To_Coupled_Faces(V,coupled_faces_solid_interpolated_velocity_np1);
        if(velocity_update || leakproof_solve){
            coupled_system->Get_Pressure(coupled_x,pressure);
            T_FACE_ARRAYS_SCALAR p_face;
            euler_projection.Get_Ghost_Pressures(dt,current_position_time,boundary_condition_collection.psi_D,boundary_condition_collection.psi_N,pressure,pressure);
            euler_projection.Get_Pressure_At_Faces(dt,current_position_time,pressure,p_face);
            coupled_system->Zero_Coupling_Faces_Values(p_face);
            euler_projection.Apply_Pressure(pressure,p_face,fluids_face_velocities,boundary_condition_collection.psi_D,boundary_condition_collection.psi_N,dt,current_position_time);
            euler_projection.p=pressure;// only for writing to file and viewing in debugger
            PHYSBAM_DEBUG_WRITE_SUBSTEP("After apply pressure at non-coupled faces (sf coupled evolution)",0,1);
            coupled_system->Apply_Lambda_To_Euler_State(coupled_x,coupled_faces_solid_interpolated_velocity_n,coupled_faces_solid_interpolated_velocity_np1,
                                                        fluids_parameters.euler_solid_fluid_coupling_utilities->cell_volumes_np1,fluids_parameters.euler->U);
            fluids_parameters.euler->Invalidate_Ghost_Cells();
            PHYSBAM_DEBUG_WRITE_SUBSTEP("After apply lambda to euler state (sf coupled evolution)",0,1);}}
    coupled_system->Get_Pressure(coupled_x,pressure);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After complete fluid state update (sf coupled evolution)",0,1);

    if(!leakproof_solve && (!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Solid_Node())){
        Finish_Backward_Euler_Step(*coupled_system,dt,current_position_time,velocity_update);}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("unscaled final pressures (sf coupled evolution)",0,1);

    coupled_system->Mark_Valid_Faces(solved_faces);
    for(int i=1;i<=twist.Size();i++){rigid_body_particles.V(i)=twist(i).linear;rigid_body_particles.angular_velocity(i)=twist(i).angular;}
    delete coupled_system;coupled_system=0;
    delete solid_system;solid_system=0;
}
//#####################################################################
// Function Output_Iterators
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Output_Iterators(const STREAM_TYPE stream_type,const char* output_directory,int frame) const
{
    // This isn't so great - it means the iterator might actually be different from the one used.
    LOG::filecout("OUTPUT ITERATORS\n");
    GEOMETRY_PARTICLES<TV> particles;
    particles.array_collection->template Add_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);
    ARRAY_VIEW<VECTOR<T,3> >& color_attribute=*particles.array_collection->template Get_Array<VECTOR<T,3> >(ATTRIBUTE_ID_COLOR);

    for(UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<TV> iterator(iterator_info,1);iterator.Valid();iterator.Next()){
        int p=particles.array_collection->Add_Element();
        particles.X(p)=iterator.Location();
        color_attribute(p)=VECTOR<T,3>(0,(T).5,0);}

    for(UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<TV> iterator(iterator_info,0,!two_phase);iterator.Valid();iterator.Next_Fluid()){
        int p=particles.array_collection->Add_Element();
        particles.X(p)=iterator.Location();
        color_attribute(p)=VECTOR<T,3>((T).5,(T).5,1);}

    FILE_UTILITIES::Create_Directory(STRING_UTILITIES::string_sprintf("%s/%i",output_directory,frame));
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%i/collision_iterators",output_directory,frame),particles);
}
//#####################################################################
// Function Make_Divergence_Free
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Make_Divergence_Free(T_FACE_ARRAYS_SCALAR& incompressible_face_velocities,const T dt,const T time)
{
    Solve(incompressible_face_velocities,dt,time,time,false,true);
}
//#####################################################################
// Function Apply_Second_Order_Cut_Cell_Method
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Apply_Second_Order_Cut_Cell_Method(const T_ARRAYS_INT& cell_index_to_divergence_matrix_index,const T_FACE_ARRAYS_INT& face_index_to_matrix_index,VECTOR_ND<T>& b)
{
/*
    // modify div_active to work for second order cut cell
    // assume only one color for the moment
    POISSON_COLLIDABLE_UNIFORM<T_GRID>& poisson=*poisson;
    T_FACE_ARRAYS_BOOL& psi_N=poisson->psi_N;
    T_ARRAYS_BOOL& psi_D=poisson->psi_D;
    // TODO: this will not work for multiphase, obviously
    assert(fluids_parameters.number_of_regions==1); // this should not be called for gas or multiphase
    T_LEVELSET* levelset=&fluids_parameters.particle_levelset_evolution->Levelset(1);
    for(FACE_ITERATOR iterator(poisson->grid);iterator.Valid();iterator.Next()){
        TV_INT face_index=iterator.Face_Index();int axis=iterator.Axis();
        TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
        if(!psi_N(axis,face_index) && LEVELSET_UTILITIES<T>::Interface(levelset->phi(first_cell_index),levelset->phi(second_cell_index)) && Cell_To_Cell_Visible(axis,first_cell_index,second_cell_index)){
            T theta=LEVELSET_UTILITIES<T>::Theta(levelset->phi(first_cell_index),levelset->phi(second_cell_index));
            if(psi_D(first_cell_index) && !psi_D(second_cell_index)){ // interface is to the negative side of second cell
                int matrix_index=cell_index_to_divergence_matrix_index(second_cell_index);
                // the lines below eliminate the usual contribution to the matrix and rhs
                int face_matrix_index=face_index_to_matrix_index(2,axis,face_index);
                T div_element=div_active(matrix_index,face_matrix_index);
                div_element/=sqrt(max((1-theta),poisson->second_order_cut_cell_threshold));
                if(div_precondition.Element_Present(matrix_index,face_matrix_index))
                    div_precondition.Set_Element(matrix_index,face_matrix_index,div_element);
                div.Set_Element(matrix_index,face_matrix_index,div_element);
                assert(!fluids_parameters.surface_tension);
                // must add this correctly before using surface tension
                //b(matrix_index)+=div_element*u_interface.Component(axis)(face_index);
            }
            else if(psi_D(second_cell_index) && !psi_D(first_cell_index)){ // interface is to the positive side of first cell
                int matrix_index=cell_index_to_divergence_matrix_index(first_cell_index);
                // the lines below eliminate the usual contribution to the matrix and rhs
                int face_matrix_index=face_index_to_matrix_index(1,axis,face_index);
                T div_element=div_active(matrix_index,face_matrix_index);
                div_element/=sqrt(max(theta,poisson->second_order_cut_cell_threshold));
                if(div_precondition.Element_Present(matrix_index,face_matrix_index))
                    div_precondition.Set_Element(matrix_index,face_matrix_index,div_element);
                div.Set_Element(matrix_index,face_matrix_index,div_element);
                assert(!fluids_parameters.surface_tension);
                // must add this correctly before using surface tension
                //b(matrix_index)+=div_element*u_interface.Component(axis)(face_index);
            }
        }
    }
*/
}
//#####################################################################
// Function Simulate_Fluids
//#####################################################################
template<class TV> bool SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Simulate_Incompressible_Fluids() const
{
    return (solids_fluids_parameters.mpi_solid_fluid || fluids_parameters.simulate) &&
        (fluids_parameters.smoke || fluids_parameters.fire || fluids_parameters.water || fluids_parameters.sph);
}
//#####################################################################
// Function Simulate_Fluids
//#####################################################################
template<class TV> bool SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Simulate_Compressible_Fluids() const
{
    return (solids_fluids_parameters.mpi_solid_fluid || fluids_parameters.simulate) && (fluids_parameters.compressible);
}
//#####################################################################
// Function Simulate_Fluids
//#####################################################################
template<class TV> bool SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Simulate_Fluids() const
{
    return (solids_fluids_parameters.mpi_solid_fluid || fluids_parameters.simulate) &&
        (fluids_parameters.smoke || fluids_parameters.fire || fluids_parameters.water || fluids_parameters.sph || fluids_parameters.compressible);
}
//#####################################################################
// Function Simulate_Solids
//#####################################################################
template<class TV> bool SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Simulate_Solids() const
{
    return (solids_fluids_parameters.mpi_solid_fluid || solid_body_collection.simulate);
}
//#####################################################################
// Function Apply_Viscosity
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Apply_Viscosity(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    if(fluids_parameters.use_coupled_implicit_viscosity) return;

    if(fluids_parameters.use_levelset_viscosity && fluids_parameters.viscosity && fluids_parameters.implicit_viscosity){
        LEVELSET_VISCOSITY_UNIFORM<TV> levelset_viscosity(fluids_parameters.callbacks,grid,dt,fluids_parameters.density,fluids_parameters.viscosity);
        levelset_viscosity.periodic_boundary=fluids_parameters.periodic_boundary;
        levelset_viscosity.print_matrix=fluids_parameters.print_viscosity_matrix;
        PHYSBAM_DEBUG_WRITE_SUBSTEP("before viscosity",0,2);
        levelset_viscosity.Apply_Viscosity(face_velocities,false,true,false);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after viscosity",0,2);
        return;}

    // TODO(kwatra): not handled for compressible flow
    if(fluids_parameters.viscosity && fluids_parameters.implicit_viscosity){
        boundary_condition_collection.Compute(grid,p,face_velocities,time);
        // fill these for viscosity
        poisson->psi_D=boundary_condition_collection.psi_D;
        poisson->psi_N=boundary_condition_collection.psi_N;
        if(boundary_condition_collection.use_psi_R){
            fluids_parameters.callbacks->Get_Reflection_Conditions(boundary_condition_collection.psi_R,time);
            poisson->psi_R=boundary_condition_collection.psi_R;}

        T_ARRAYS_SCALAR variable_viscosity;
        VISCOSITY<T_GRID> viscosity_helper(*poisson,variable_viscosity,fluids_parameters.density,fluids_parameters.viscosity,fluids_parameters.implicit_viscosity,false,false,1000,fluids_parameters.use_psi_R);

        std::stringstream ss;
        ss<<"KE before viscosity: "<<std::endl;
        TV KE=TV();
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
            KE[iterator.Axis()]+=sqr(face_velocities(iterator.Axis(),iterator.Face_Index()));
        for(int axis=1;axis<=TV::dimension;axis++)
            ss<<"axis "<<axis<<": "<<KE[axis]<<std::endl;
        viscosity_helper.Add_Implicit_Forces_Before_Projection(grid,face_velocities,face_velocities,dt,time);
        ss<<"KE after viscosity: "<<std::endl;
        KE=TV();
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next())
            KE[iterator.Axis()]+=sqr(face_velocities(iterator.Axis(),iterator.Face_Index()));
        for(int axis=1;axis<=TV::dimension;axis++)
            ss<<"axis "<<axis<<": "<<KE[axis]<<std::endl;
        LOG::filecout(ss.str());
    }
}
//#####################################################################
// Function Warn_For_Exposed_Dirichlet_Cell
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Warn_For_Exposed_Dirichlet_Cell(const T_ARRAYS_BOOL& psi_D,const T_FACE_ARRAYS_BOOL& psi_N)
{
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(!psi_D(cell_index)) continue;
        for(int axis=1;axis<=TV::dimension;axis++){
            TV_INT first_cell_index=iterator.Cell_Neighbor(2*axis-1),second_cell_index=iterator.Cell_Neighbor(2*axis);
            TV_INT first_face_index=iterator.First_Face_Index(axis),second_face_index=iterator.Second_Face_Index(axis);
            if((grid.Inside_Domain(first_cell_index) && !psi_D(first_cell_index) && !psi_N.Component(axis)(first_face_index)) || 
                (grid.Inside_Domain(second_cell_index) && !psi_D(second_cell_index) && !psi_N.Component(axis)(second_face_index))){
                std::stringstream ss;
                ss<<"Warning: Exposed dirichlet cell "<<cell_index<<std::endl;LOG::filecout(ss.str());}}}
}
//#####################################################################
// Function Set_Cached_Psi_N_And_Coupled_Face_Data
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Set_Cached_Psi_N_And_Coupled_Face_Data(const COLLISION_AWARE_INDEX_MAP<TV>& index_map,const MATRIX_SOLID_INTERPOLATION<TV>& solid_interpolation,
    const T time)
{
    time_cached=time;
    cached_coupling_face_data=true;
    Get_Coupled_Faces_And_Interpolated_Solid_Velocities(index_map,solid_interpolation,boundary_condition_collection.psi_N,cached_psi_N,coupling_face_velocities_cached);
    number_of_coupling_faces_cached=solid_interpolation.Number_Of_Constraints();
    indexed_faces_cached=index_map.indexed_faces;
}
//#####################################################################
// Function Fill_Coupled_Face_Data
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Fill_Coupled_Face_Data(const COUPLING_CONSTRAINT_ID number_of_coupling_faces,const ARRAY<FACE_INDEX<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::TV::dimension> >& indexed_faces,
    const ARRAY<T,COUPLING_CONSTRAINT_ID>& coupling_face_data,T_FACE_ARRAYS_SCALAR& face_data)
{
    for(COUPLING_CONSTRAINT_ID i(1);i<=number_of_coupling_faces;i++){
        FACE_INDEX<TV::dimension> face_index=indexed_faces(Value(i));
        face_data(face_index)=coupling_face_data(i);}
}
//#####################################################################
// Function Get_Coupled_Faces_And_Interpolated_Solid_Velocities
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Get_Coupled_Faces_And_Interpolated_Solid_Velocities(const COLLISION_AWARE_INDEX_MAP<TV>& index_map,
    const MATRIX_SOLID_INTERPOLATION<TV>& solid_interpolation,const T_FACE_ARRAYS_BOOL& psi_N_domain,T_FACE_ARRAYS_BOOL& psi_N,
    ARRAY<T,COUPLING_CONSTRAINT_ID>& coupling_face_velocities)
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particle;
    ARRAY<TWIST<TV> > twist;twist.Resize(rigid_body_particles.V.Size());
    for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(rigid_body_particles.V(i),rigid_body_particles.angular_velocity(i));

    psi_N=psi_N_domain;

    // Set coupling faces as Neumann
    for(UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<TV> iterator(index_map.iterator_info);iterator.Valid();iterator.Next()){
        psi_N(iterator.axis,iterator.index)=true;}

    // Project solid velocities at coupling faces. (JV_s)
    GENERALIZED_VELOCITY<TV> V_S(particles.V,twist,solid_body_collection);
    coupling_face_velocities.Resize(solid_interpolation.Number_Of_Constraints());
    solid_interpolation.Times(V_S,coupling_face_velocities);
    for(int i=1;i<=twist.Size();i++){rigid_body_particles.V(i)=twist(i).linear;rigid_body_particles.angular_velocity(i)=twist(i).angular;}
}
//#####################################################################
// Function Get_Coupled_Faces_And_Interpolated_Solid_Velocities
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Get_Coupled_Faces_And_Interpolated_Solid_Velocities(const T time,T_FACE_ARRAYS_BOOL& psi_N,T_FACE_ARRAYS_SCALAR& face_velocities)
{
    if(cached_coupling_face_data){
        assert(time==time_cached);
        psi_N=cached_psi_N;
        Fill_Coupled_Face_Data(number_of_coupling_faces_cached,indexed_faces_cached,coupling_face_velocities_cached,face_velocities);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after Get_Coupled_Faces_And_Velocities using cached data (sf coupled evolution)",0,2);}
    else{
        boundary_condition_collection.Compute(grid,pressure,face_velocities,time);

        collision_bodies.outside_fluid=&boundary_condition_collection.psi_D;
        
        COLLISION_AWARE_INDEX_MAP<TV> index_map(iterator_info,boundary_condition_collection);
        iterator_info.Initialize_Collision_Aware_Face_Iterator(boundary_condition_collection.psi_D,boundary_condition_collection.psi_N,0,disable_thinshell);
        index_map.Construct_Indices(solids_fluids_parameters.mpi_solid_fluid?1:0);
        MATRIX_SOLID_INTERPOLATION<TV> solid_interpolation(iterator_info);
        solid_interpolation.Compute(0);

        ARRAY<T,COUPLING_CONSTRAINT_ID> coupling_face_velocities;
        Get_Coupled_Faces_And_Interpolated_Solid_Velocities(index_map,solid_interpolation,boundary_condition_collection.psi_N,psi_N,coupling_face_velocities);
        Fill_Coupled_Face_Data(solid_interpolation.Number_Of_Constraints(),index_map.indexed_faces,coupling_face_velocities,face_velocities);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after Get_Coupled_Faces_And_Velocities (sf coupled evolution)",0,2);}
}
//#####################################################################
// Function Setup_Boundary_Condition_Collection
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::
Setup_Boundary_Condition_Collection()
{
    VECTOR<VECTOR<bool,2>,TV::dimension> mpi_boundaries;mpi_boundaries.Fill(VECTOR<bool,2>(false,false));
    if(fluids_parameters.mpi_grid) for(int axis=1;axis<=TV::dimension;++axis) for(int axis_side=1;axis_side<=2;++axis_side)
        if(fluids_parameters.mpi_grid->side_neighbor_ranks(2*(axis-1)+axis_side)>=0) mpi_boundaries(axis)(axis_side)=true;

    if(Simulate_Incompressible_Fluids()){
        if(fluids_parameters.particle_levelset_evolution){
            boundary_condition_collection.Add_Boundary_Condition(new INCOMPRESSIBLE_BOUNDARY_CONDITION_FREE_SURFACE<TV>(fluids_parameters.particle_levelset_evolution->Levelset(1).phi));
            if(fluids_parameters.surface_tension && false)
                boundary_condition_collection.Add_Boundary_Condition(
                    new SURFACE_TENSION_BOUNDARY_CONDITION<TV>(fluids_parameters.particle_levelset_evolution->Levelset(1),fluids_parameters.surface_tension));}
        boundary_condition_collection.Add_Boundary_Condition(new INCOMPRESSIBLE_BOUNDARY_CONDITION_WALLS<TV>(fluids_parameters.domain_walls,mpi_boundaries));}
    if(Simulate_Compressible_Fluids()){
        boundary_condition_collection.Add_Boundary_Condition(new COMPRESSIBLE_BOUNDARY_CONDITION_WALLS<TV>(fluids_parameters.domain_walls,mpi_boundaries,fluids_parameters));}
    boundary_condition_collection.Add_Boundary_Condition(new EXAMPLE_BOUNDARY_CONDITION<TV>(fluids_parameters.callbacks));
}
//#####################################################################
template class SOLID_FLUID_COUPLED_EVOLUTION_SLIP<VECTOR<float,1> >;
template class SOLID_FLUID_COUPLED_EVOLUTION_SLIP<VECTOR<float,2> >;
template class SOLID_FLUID_COUPLED_EVOLUTION_SLIP<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SOLID_FLUID_COUPLED_EVOLUTION_SLIP<VECTOR<double,1> >;
template class SOLID_FLUID_COUPLED_EVOLUTION_SLIP<VECTOR<double,2> >;
template class SOLID_FLUID_COUPLED_EVOLUTION_SLIP<VECTOR<double,3> >;
#endif
