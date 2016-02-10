//#####################################################################
// Copyright 2007-2009, Jon Gretarsson, Nipun Kwatra, Michael Lentine, Avi Robinson-Mosher, Craig Schroeder, Andrew Selle, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM hose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Krylov_Solvers/LANCZOS_ITERATION.h>
#include <PhysBAM_Tools/Krylov_Solvers/SYMMQMR.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Math_Tools/Inverse.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/IMPLICIT_OBJECT_INTERSECTOR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BINDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_1D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_PROJECTION_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Grid_Based_Fields/INCOMPRESSIBLE_FLUID_CONTAINER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_COMPRESSIBLE_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_SYSTEM.h>
#include <PhysBAM_Dynamics/Incompressible_Flows/INCOMPRESSIBLE_MULTIPHASE_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Parallel_Computation/CONJUGATE_RESIDUAL_SPARSE_MPI.h>
#include <PhysBAM_Dynamics/Parallel_Computation/FLUID_SYSTEM_MPI.h>
#include <PhysBAM_Dynamics/Parallel_Computation/MPI_SOLID_FLUID.h>
#include <PhysBAM_Dynamics/Parallel_Computation/SOLID_SYSTEM_MPI.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_CALLBACKS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLID_FLUID_COUPLED_EVOLUTION<TV>::
SOLID_FLUID_COUPLED_EVOLUTION(SOLIDS_PARAMETERS<TV>& solids_parameters_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >& fluids_parameters_input,
    INCOMPRESSIBLE_FLUID_CONTAINER<GRID<TV> >& incompressible_fluid_container_input,SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters_input)
    :NEWMARK_EVOLUTION<TV>(solids_parameters_input,solid_body_collection_input),rigid_body_count(0),print_matrix_rhs_and_solution(false),
    collision_bodies(*fluids_parameters_input.collision_bodies_affecting_fluid),dual_cell_weights(*fluids_parameters_input.grid),
    rigid_body_dual_cell_weights(*fluids_parameters_input.grid),dual_cell_fluid_volume(*fluids_parameters_input.grid),
    dual_cell_contains_solid(*fluids_parameters_input.grid,1),fluids_parameters(fluids_parameters_input),incompressible_fluid_container(incompressible_fluid_container_input),
    solids_fluids_parameters(solids_fluids_parameters_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLID_FLUID_COUPLED_EVOLUTION<TV>::
~SOLID_FLUID_COUPLED_EVOLUTION()
{
    dual_cell_weights.Delete_Pointers_And_Clean_Memory();
    rigid_body_dual_cell_weights.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Get_Density_At_Face
//#####################################################################
template<class TV> typename TV::SCALAR SOLID_FLUID_COUPLED_EVOLUTION<TV>::
Get_Density_At_Face(const int axis,const TV_INT& face_index)
{
    T density=0;
    if(fluids_parameters.compressible){
        TV_INT first_cell_index,second_cell_index;
        GRID<TV>::Cells_Touching_Face(axis,face_index,first_cell_index,second_cell_index);
        // TODO: use U_ghost instead
        if(fluids_parameters.euler->U.Valid_Index(first_cell_index) && fluids_parameters.euler->U.Valid_Index(second_cell_index))
            density=(fluids_parameters.euler->U(first_cell_index)(1)+fluids_parameters.euler->U(second_cell_index)(1))*(T).5;
        else if(fluids_parameters.euler->U.Valid_Index(first_cell_index))
            density=fluids_parameters.euler->U(first_cell_index)(1);
        else density=fluids_parameters.euler->U(second_cell_index)(1);}
    else density=fluids_parameters.density;
    return density;
}
//#####################################################################
// Function Backward_Euler_Step_Velocity_Helper
//#####################################################################
// assumes all solids_forces are linear in velocity, with a symmetric positive definite Jacobian.
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION<TV>::
Backward_Euler_Step_Velocity_Helper(const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update)
{
    PHYSBAM_ASSERT(!solids_parameters.use_trapezoidal_rule_for_velocities); // two-way coupling does not work with trapezoidal rule
    bool solids=Simulate_Solids() && (!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Solid_Node());
    bool fluids=Simulate_Fluids() && (!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Fluid_Node());
    ARRAY<int,VECTOR<int,1> > filled_region_cell_count;

    if(fluids_parameters.compressible && !fluids_parameters.euler->timesplit){
        //TODO(kwatra): move this logic to driver
        BASE::Backward_Euler_Step_Velocity_Helper(dt,current_velocity_time,current_position_time,velocity_update);
        return;}

    if(!fluids_parameters.fluid_affects_solid){
        if(fluids){
            POISSON_COLLIDABLE_UNIFORM<GRID<TV> >& poisson=*Get_Poisson();
            Compute_W(current_position_time);

            if(fluids_parameters.compressible){
                poisson.beta_given_on_faces=true;
                fluids_parameters.euler->euler_projection.Fill_Face_Weights_For_Projection(dt,current_position_time,poisson.beta_face);}// TODO(jontg): which time do we use?
            else{poisson.beta_given_on_faces=true;poisson.beta_face.Fill(1/fluids_parameters.density);}
            poisson.tolerance=Get_Grid().Cell_Size()*(T)1e-4;}
        if(solids) BASE::Backward_Euler_Step_Velocity_Helper(dt,current_velocity_time,current_position_time,velocity_update);
        return;}
    int colors=0;
    // Doing super fragments for s-f coupling requires connectivity through fluid.

    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particle;
    MPI_SOLIDS<TV>* mpi_solids=solid_body_collection.deformable_body_collection.mpi_solids;

    if(solids){
        F_full.Resize(particles.array_collection->Size(),false,false);rigid_F_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
        R_full.Resize(particles.array_collection->Size(),false,false);rigid_R_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
        S_full.Resize(particles.array_collection->Size(),false,false);rigid_S_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
        B_full.Resize(particles.array_collection->Size(),false,false);rigid_B_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
        ar_full.Resize(particles.array_collection->Size(),false,false);rigid_ar_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
        z_full.Resize(particles.array_collection->Size(),false,false);rigid_z_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
        zaq_full.Resize(particles.array_collection->Size(),false,false);rigid_zaq_full.Resize(rigid_body_particles.array_collection->Size(),false,false);}
    else{
        if(fluids && solids_fluids_parameters.mpi_solid_fluid){ // Gather the fluid terms for the RHS of the solid here
            F_full.Resize(particles.array_collection->Size(),false,false);rigid_F_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
            R_full.Resize(particles.array_collection->Size(),false,false);rigid_R_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
            B_full.Resize(particles.array_collection->Size(),false,false);rigid_B_full.Resize(rigid_body_particles.array_collection->Size(),false,false);}}
    ARRAY<TWIST<TV> > twist;twist.Resize(rigid_body_particles.V.Size());
    for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(rigid_body_particles.V(i),rigid_body_particles.angular_velocity(i));
    GENERALIZED_VELOCITY<TV> V(particles.V,twist,solid_body_collection),F(F_full,rigid_F_full,solid_body_collection),
        R(R_full,rigid_R_full,solid_body_collection),S(S_full,rigid_S_full,solid_body_collection),B(B_full,rigid_B_full,solid_body_collection),
        ar(ar_full,rigid_ar_full,solid_body_collection),z(z_full,rigid_z_full,solid_body_collection),zaq(zaq_full,rigid_zaq_full,solid_body_collection);
    GENERALIZED_MASS<TV> mass(solid_body_collection);

    T_ARRAYS_INT cell_index_to_matrix_index;
    ARRAY<INTERVAL<int> > interior_regions;
    int number_of_regions=0;
    ARRAY<SPARSE_MATRIX_FLAT_NXN<T> > A_array;
    KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > > kb_array;
    ARRAY<VECTOR_ND<T> >& b_array=kb_array.v;

    ARRAY<int> rigid_body_particles_to_dynamic_rigid_body_particles_map(rigid_body_particles.array_collection->Size());
    rigid_body_particles_to_dynamic_rigid_body_particles_map.Subset(solid_body_collection.rigid_body_collection.simulated_rigid_body_particles)=IDENTITY_ARRAY<int>(solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m);

    BACKWARD_EULER_SYSTEM<TV>* solid_system=0;
    if(solids){ // Set up solids RHS (B)
        rigid_body_collection.Update_Angular_Velocity(solid_body_collection.rigid_body_collection.simulated_rigid_body_particles); // make sure omega = I^{-1} L

        // compute solid momentum before
        /*TV dimensionwise_solid_momentum;
        for(int i=1;i<=V.V.Size();i++){
            dimensionwise_solid_momentum+=mass.mass(i)*V.V(i);
        }
        LOG::cout<<"Deformable solid momentum before explicit forces: "<<dimensionwise_solid_momentum<<std::endl;*/

        INDIRECT_ARRAY<ARRAY_VIEW<TV>,ARRAY<int>&> BV_subset=B.V.array.Subset(solid_body_collection.deformable_body_collection.dynamic_particles);ARRAYS_COMPUTATIONS::Fill(BV_subset,TV());
        ARRAYS_COMPUTATIONS::Fill(B.rigid_V.array,TWIST<TV>());
        solid_body_collection.example_forces_and_velocities->Add_External_Forces(B.V.array,current_velocity_time+dt);
        solid_body_collection.example_forces_and_velocities->Add_External_Forces(B.rigid_V.array,current_velocity_time+dt);
        solid_body_collection.Add_Velocity_Independent_Forces(B.V.array,B.rigid_V.array,current_velocity_time+dt); // this is a nop for binding forces
        solid_body_collection.deformable_body_collection.binding_list.Distribute_Force_To_Parents(B.V.array,B.rigid_V.array);
        solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Distribute_Force_To_Parents(rigid_B_full);
        if(solid_body_collection.deformable_body_collection.soft_bindings.Need_Bindings_Mapped()){
            solid_body_collection.deformable_body_collection.soft_bindings.Map_Forces_From_Parents(B.V.array,B.rigid_V.array);
            solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(B.V.array);
            for(int k=1;k<=solid_body_collection.solids_forces.m;k++){
                if(dynamic_cast<BINDING_SPRINGS<TV>*>(solid_body_collection.solids_forces(k)))
                    solid_body_collection.solids_forces(k)->Add_Force_Differential(particles.X,B.V.array,current_velocity_time+dt);}
            solid_body_collection.deformable_body_collection.binding_list.Distribute_Force_To_Parents(B.V.array);}

        // TODO: Make sure this runs for the solid node of an MPI coupled sim
        Initialize_World_Space_Masses();
        solid_system=new BACKWARD_EULER_SYSTEM<TV>(*this,solid_body_collection,dt,current_velocity_time,current_position_time,
            &solid_body_collection.rigid_body_collection.articulated_rigid_body,(velocity_update && solids_parameters.enforce_repulsions_in_cg)?repulsions:0,mpi_solids,velocity_update);
        
        if(solid_system->arb->constrain_pd_directions){
            for(int i=1;i<=B.rigid_V.Size();i++) B.rigid_V(i)=mass.world_space_rigid_mass_inverse(i)*B.rigid_V(i);
            solid_system->arb->Poststabilization_Projection(B.rigid_V.array,true);
            for(int i=1;i<=B.rigid_V.Size();i++) B.rigid_V(i)=mass.world_space_rigid_mass(i)*B.rigid_V(i);}

        if(Simulate_Fluids()){
            for(int i=1;i<=V.V.Size();i++) B.V(i)=mass.mass(i)*V.V(i)+dt*B.V(i);
            for(int i=1;i<=V.rigid_V.Size();i++) B.rigid_V(i)=mass.world_space_rigid_mass(i)*V.rigid_V(i)+dt*B.rigid_V(i);}
        else{
            for(int i=1;i<=V.V.Size();i++) B.V(i)=V.V(i)+dt*mass.one_over_mass(i)*B.V(i);
            for(int i=1;i<=V.rigid_V.Size();i++) B.rigid_V(i)=V.rigid_V(i)+dt*(mass.world_space_rigid_mass_inverse(i)*B.rigid_V(i));}

        if(solids_fluids_parameters.mpi_solid_fluid) solids_fluids_parameters.mpi_solid_fluid->Exchange_Solid_Positions_And_Velocities(solid_body_collection); // TODO: only need to exchange velocities
        }
    if(fluids){
        POISSON_COLLIDABLE_UNIFORM<GRID<TV> >& poisson=*Get_Poisson();
        T_FACE_ARRAYS_SCALAR& face_velocities=Get_Face_Velocities();
        const GRID<TV>& grid=Get_Grid();

        if(fluids_parameters.compressible){
            poisson.beta_given_on_faces=true;
            fluids_parameters.euler->euler_projection.Fill_Face_Weights_For_Projection(dt,current_position_time,poisson.beta_face);}// TODO(jontg): which time do we use?
        else{poisson.beta_given_on_faces=true;poisson.beta_face.Fill(1/fluids_parameters.density);}

        poisson.tolerance=grid.Cell_Size()*(T)1e-4;

        if(solids_fluids_parameters.mpi_solid_fluid) solids_fluids_parameters.mpi_solid_fluid->Exchange_Solid_Positions_And_Velocities(solid_body_collection); // TODO: only need to exchange velocities

        // TODO(jontg): compressible
        Compute_W(current_position_time); // set Neumann conditions on all dual cells which have object faces
        
        /*TV dimensionwise_fluid_momentum;
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            const int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            if(dual_cell_fluid_volume(axis,face_index)) dimensionwise_fluid_momentum(axis)+=Get_Density_At_Face(axis,face_index)*dual_cell_fluid_volume(axis,face_index)*face_velocities(axis,face_index);}
        LOG::cout<<"Dimensionwise fluid momentum after explicit forces: "<<dimensionwise_fluid_momentum<<std::endl;*/

        Set_Neumann_and_Dirichlet_Boundary_Conditions(current_position_time);


        // For energy calculation, get solids' V* velocity and project it to fluid-solid faces using the weights just computed
        if(fluids_parameters.compressible && velocity_update){
            for(int i=1;i<=V.V.Size();i++) V.V(i)=mass.one_over_mass(i)*B.V(i);
            for(int i=1;i<=V.rigid_V.Size();i++) V.rigid_V(i)=mass.world_space_rigid_mass_inverse(i)*B.rigid_V(i);
            solid_projected_face_velocities_star.Resize(Get_Grid(),true,false);
            Apply_Solid_Boundary_Conditions(current_velocity_time,false,solid_projected_face_velocities_star);}

        poisson.Find_Solution_Regions();
        poisson.Compute_beta_And_Add_Jumps_To_b(dt,current_position_time);
        
        // Set up fluids RHS (poisson->f)
        T_FACE_ARRAYS_SCALAR copy_face_velocities(face_velocities);
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            const int axis=iterator.Axis();const TV_INT face_index=iterator.Face_Index();
            if(dual_cell_contains_solid(axis,face_index)) face_velocities(axis,face_index)=(T)0;} // These faces should not be pulled to the RHS, so zero them out
        Transfer_Momentum_And_Set_Boundary_Conditions(current_position_time,&B);

        if(!fluids_parameters.compressible && fluids_parameters.second_order_cut_cell_method) poisson.Set_Up_Second_Order_Cut_Cell_Method();
        if(fluids_parameters.compressible){
            fluids_parameters.euler->euler_projection.Compute_Density_Weighted_Face_Velocities(dt,current_position_time,fluids_parameters.euler->euler_projection.elliptic_solver->psi_N);
            fluids_parameters.euler->euler_projection.Compute_Right_Hand_Side(face_velocities,dt,current_position_time);}
        else fluids_parameters.incompressible->projection.Compute_Divergence(T_FACE_LOOKUP(face_velocities),&poisson);

        // restore face velocities for mixed faces
        T_FACE_ARRAYS_SCALAR::Copy(copy_face_velocities,face_velocities);

        // Set up fluids matrix (A_array) and RHS (b_array from modifying poisson->f for boundary conditions in psi_D)
        number_of_regions=poisson.number_of_regions;
        matrix_index_to_cell_index_array.Resize(number_of_regions);cell_index_to_matrix_index.Resize(grid.Domain_Indices(1));
        ARRAY<int,VECTOR<int,1> > filled_region_cell_count(-1,number_of_regions);
        A_array.Resize(number_of_regions);b_array.Resize(number_of_regions);
        for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) filled_region_cell_count(poisson.filled_region_colors(iterator.Cell_Index()))++;
        for(int color=1;color<=number_of_regions;color++) if(poisson.filled_region_touches_dirichlet(color)||poisson.solve_neumann_regions){
            matrix_index_to_cell_index_array(color).Resize(filled_region_cell_count(color));}
        filled_region_cell_count.Fill(0); // reusing this array in order to make the indirection arrays

        if(!poisson.mpi_grid) poisson.Compute_Matrix_Indices(filled_region_cell_count,matrix_index_to_cell_index_array,cell_index_to_matrix_index);
        else poisson.laplace_mpi->Find_Matrix_Indices(filled_region_cell_count,cell_index_to_matrix_index,matrix_index_to_cell_index_array);

        PHYSBAM_DEBUG_WRITE_SUBSTEP("before poisson setup for coupled solve (sf coupled evolution)",0,1);

        T_ARRAYS_SCALAR& p=Get_Pressure();
        if(dt && !fluids_parameters.compressible) p*=dt; // RESCALE PRESSURE FOR A BETTER INITIAL GUESS!
        if(dt && fluids_parameters.compressible) fluids_parameters.euler->euler_projection.Scale_Pressure_By_Dt(dt);
 
        poisson.Set_Dt(dt);
        RANGE<TV_INT> domain(grid.Domain_Indices(1));
        poisson.Find_A(domain,A_array,b_array,filled_region_cell_count,cell_index_to_matrix_index);

        //TODO(kwatra): Why not just use number_of_regions
        colors=filled_region_cell_count.domain.max_corner.x;

        interior_regions.Resize(colors);
        for(int color=1;color<=colors;color++) if(filled_region_cell_count(color)>0){
            if(!poisson.laplace_mpi || color>poisson.laplace_mpi->partitions.m) interior_regions(color)=INTERVAL<int>(1,filled_region_cell_count(color));
            else interior_regions(color)=poisson.laplace_mpi->partitions(color).interior_indices;}}

    // Compute coupling matrices, modified solid masses. Fix angular part of rigid RHS(B) for modified com.
    if(Simulate_Fluids()){
        // interior_regions and cell_index_to_matrix_index are dummy here for solids mpi
        Compute_Coupling_Terms_Deformable(cell_index_to_matrix_index,interior_regions,number_of_regions);
        Compute_Coupling_Terms_Rigid(cell_index_to_matrix_index,interior_regions,number_of_regions);

        if(fluids){  // scale coupling matrices by V
            const GRID<TV>& grid=Get_Grid();
            J_deformable*=-grid.Cell_Size();J_rigid*=-grid.Cell_Size();
            if(solids_fluids_parameters.mpi_solid_fluid){  // in MPI, fluid nodes don't need these...
                ARRAYS_COMPUTATIONS::Fill(B.V,TV());ARRAYS_COMPUTATIONS::Fill(B.rigid_V,TWIST<TV>());}}
        if(solids) for(int i=1;i<=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m;i++){  // Update angular momentum for modified center of mass
            // TODO(kwatra): do we need to this for compressible case, as there is no mass lumping there?
            RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(solid_body_collection.rigid_body_collection.simulated_rigid_body_particles(i));
            B.rigid_V(i).angular+=TV::Cross_Product((rigid_body_updated_center_of_mass(i)-rigid_body.X()),B.rigid_V(i).linear);}}

    // Add fluids component to solids RHS(B)
    if(fluids){
        Add_Nondynamic_Solids_To_Right_Hand_Side(b_array,interior_regions,colors);
        const GRID<TV>& grid=Get_Grid();
        b_array*=-grid.Cell_Size();A_array*=-grid.Cell_Size();

        if(!fluids_parameters.compressible){
            PHYSBAM_DEBUG_WRITE_SUBSTEP("using face velocities to set up solid RHS (sf coupled evolution)",0,1);
            T_FACE_ARRAYS_SCALAR& face_velocities=Get_Face_Velocities();
            for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                const int axis=iterator.Axis();const TV_INT face_index=iterator.Face_Index();TV axis_vector=TV::Axis_Vector(axis);
                const TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
                TV momentum_chunk=axis_vector*dual_cell_fluid_volume(axis,face_index)*Get_Density_At_Face(axis,face_index)*face_velocities(axis,face_index);
                if(dual_cell_weights(axis,face_index) && dual_cell_fluid_volume(axis,face_index)){
                    FACE_WEIGHT_ELEMENTS& weights=*dual_cell_weights(axis,face_index);
                    for(int i=1;i<=weights.m;i++) B.V(weights(i).x)+=weights(i).y*momentum_chunk;}
                    if(rigid_body_dual_cell_weights(axis,face_index) && dual_cell_fluid_volume(axis,face_index)){
                        FACE_WEIGHT_ELEMENTS& weights=*rigid_body_dual_cell_weights(axis,face_index);
                        for(int i=1;i<=weights.m;i++){
                            int rigid_body_id=weights(i).x;T weight=weights(i).y;
                            int rigid_body_index=rigid_body_particles_to_dynamic_rigid_body_particles_map(rigid_body_id);
                            if(rigid_body_index!=0){
                                TV R_chunk=iterator.Location()-rigid_body_updated_center_of_mass(rigid_body_index);
                                B.rigid_V(rigid_body_index).linear+=weight*momentum_chunk;
                                B.rigid_V(rigid_body_index).angular+=TV::Cross_Product(R_chunk,weight*momentum_chunk);}}}}}}

    // System is now set. Solve it.
    bool preconditioned=true;
    SOLID_SYSTEM_MPI<TV>* solid_system_mpi=0;
    FLUID_SYSTEM_MPI<TV>* fluid_system_mpi=0;
    ARRAY<int> coupled_deformable_particle_indices;
    ARRAY<ARRAY<int> > coupled_deformable_particle_indices_array;

    if(solids_fluids_parameters.mpi_solid_fluid){
        if(fluids){
            RANGE<TV> grid_domain=fluids_parameters.grid->domain;
            grid_domain.Change_Size(fluids_parameters.grid->dX*(T).5);
            ARRAY<int> boundary_particles(particles.array_collection->Size());ARRAYS_COMPUTATIONS::Fill(boundary_particles,0);
            GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >& collision_bodies_affecting_fluid=*fluids_parameters.collision_bodies_affecting_fluid;
            for(COLLISION_GEOMETRY_ID i(1);i<=collision_bodies_affecting_fluid.collision_geometry_collection.bodies.m;i++)
                if(collision_bodies_affecting_fluid.collision_geometry_collection.Is_Active(i)){
                    COLLISION_GEOMETRY<TV>& body=*collision_bodies_affecting_fluid.collision_geometry_collection.bodies(i);
                    if(!body.active) continue;
                    if(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>* deformable=dynamic_cast<DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>*>(&body))
                        for(int simplex=1;simplex<=deformable->object.mesh.elements.m;simplex++){
                            RANGE<TV> simplex_bounding_box(particles.X(deformable->object.mesh.elements(simplex)(1)));
                            for(int node=2;node<=deformable->object.mesh.dimension;node++) simplex_bounding_box.Enlarge_To_Include_Point(particles.X(deformable->object.mesh.elements(simplex)(node)));
                            if(grid_domain.Lazy_Intersection(simplex_bounding_box)){
                                INDIRECT_ARRAY<ARRAY<int>,VECTOR<int,TV::m>&> subset=boundary_particles.Subset(deformable->object.mesh.elements(simplex)); // TODO(jontg): Check this
                                ARRAYS_COMPUTATIONS::Fill(subset,1);}}}
            for(int i=1;i<=boundary_particles.m;i++) if(boundary_particles(i)) coupled_deformable_particle_indices.Append(i);
            solids_fluids_parameters.mpi_solid_fluid->Exchange_Coupled_Deformable_Particle_List(&coupled_deformable_particle_indices,0);
            fluid_system_mpi=new FLUID_SYSTEM_MPI<TV>(J_deformable,J_rigid,A_array,interior_regions,
                solids_parameters.implicit_solve_parameters.cg_tolerance/Get_Poisson()->tolerance,solids_fluids_parameters.mpi_solid_fluid,F,R,coupled_deformable_particle_indices,preconditioned);
            fluid_system_mpi->Send_Generalized_Velocity_To_Solid(B);}
        else{
            coupled_deformable_particle_indices_array.Resize(solids_fluids_parameters.mpi_solid_fluid->fluid_ranks.n);
            solids_fluids_parameters.mpi_solid_fluid->Exchange_Coupled_Deformable_Particle_List(0,&coupled_deformable_particle_indices_array);
            solid_system_mpi=new SOLID_SYSTEM_MPI<TV>(*solid_system,nodal_fluid_mass,rigid_body_fluid_mass,rigid_body_fluid_inertia,solids_fluids_parameters.mpi_solid_fluid,
                coupled_deformable_particle_indices_array,*this,V.rigid_V.Size());
            // get fluid contribution to the solid RHS
            solid_system_mpi->Get_Generalized_Velocity_From_Fluid(B);}}

    if(solids && Simulate_Fluids()){B.V*=(T)-1;B.rigid_V*=(T)-1;}

    KRYLOV_VECTOR_WRAPPER<T,ARRAY<VECTOR_ND<T> > > x_array,p_array,ap_array,ar_array,r_array,z_array,zaq_array;
    x_array.v.Resize(A_array.m);p_array.v.Resize(A_array.m);ap_array.v.Resize(A_array.m);ar_array.v.Resize(A_array.m);
    r_array.v.Resize(A_array.m);z_array.v.Resize(A_array.m);zaq_array.v.Resize(A_array.m);
    if(fluids){
        for(int i=1;i<=A_array.m;i++){int n=A_array(i).n;
            x_array.v(i).Resize(n);p_array.v(i).Resize(n);ap_array.v(i).Resize(n);
            ar_array.v(i).Resize(n);r_array.v(i).Resize(n);z_array.v(i).Resize(n);zaq_array.v(i).Resize(n);}

        T_ARRAYS_SCALAR& p=Get_Pressure();
        for(int i=1;i<=A_array.m;i++){
            const SPARSE_MATRIX_FLAT_NXN<T>& A=A_array(i);
            for(int j=1;j<=A.n;j++) x_array.v(i)(j)=p(matrix_index_to_cell_index_array(i)(j));}}

    if(solids_fluids_parameters.mpi_solid_fluid){
        LOG::Time("conjugate residual parallel");
        if(fluids){
            POISSON_COLLIDABLE_UNIFORM<GRID<TV> >& poisson=*Get_Poisson();
            solids_fluids_parameters.mpi_solid_fluid->Parallel_Solve_Fluid_Part(*fluid_system_mpi,x_array,kb_array,p_array,ap_array,ar_array,r_array,z_array,zaq_array,
                1,solids_parameters.implicit_solve_parameters.cg_iterations,solids_parameters.implicit_solve_parameters.cg_tolerance,true,poisson.laplace_mpi->communicators,&poisson.laplace_mpi->partitions);}
        else{
            // TODO: unify with single proc case
            for(int i=1;i<=V.V.Size();i++) B.V(i)=solid_system_mpi->one_over_modified_mass(i)*B.V(i);
            for(int i=1;i<=V.rigid_V.Size();i++){B.rigid_V(i).linear=solid_system_mpi->modified_world_space_rigid_mass_inverse(i)*B.rigid_V(i).linear;
                B.rigid_V(i).angular=solid_system_mpi->modified_world_space_rigid_inertia_tensor_inverse(i)*B.rigid_V(i).angular;}
            V.V=B.V;V.rigid_V=B.rigid_V;
            V.V*=(T)-1;V.rigid_V*=(T)-1;
            solid_system->Set_Global_Boundary_Conditions(V,X_save,rigid_X_save,rigid_rotation_save,rigid_velocity_save,rigid_angular_momentum_save,V_save,
                solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);
            solids_fluids_parameters.mpi_solid_fluid->Parallel_Solve_Solid_Part(*solid_system_mpi,V,B,F,S,ar,R,z,zaq,1,solids_parameters.implicit_solve_parameters.cg_iterations,solids_parameters.implicit_solve_parameters.cg_tolerance);}
        LOG::Stop_Time();}
    else{
        T fluid_tolerance=0;
        if(fluids_parameters.compressible || fluids_parameters.incompressible) fluid_tolerance=Get_Poisson()->tolerance;

        SOLID_FLUID_SYSTEM<TV,SPARSE_MATRIX_FLAT_NXN<T> > solid_fluid_system(*solid_system,J_deformable,J_rigid,nodal_fluid_mass,rigid_body_fluid_mass,rigid_body_fluid_inertia,
            fluid_tolerance,solids_parameters.implicit_solve_parameters.cg_tolerance,A_array);
        if(preconditioned){solid_fluid_system.use_preconditioner=true;solid_fluid_system.preconditioner_commutes_with_projection=false;}

        if(fluids){
            // correct RHS for solids
            for(int i=1;i<=V.V.Size();i++) B.V(i)=solid_fluid_system.one_over_modified_mass(i)*B.V(i);
            for(int i=1;i<=V.rigid_V.Size();i++){B.rigid_V(i).linear=solid_fluid_system.modified_world_space_rigid_mass_inverse(i)*B.rigid_V(i).linear;
                B.rigid_V(i).angular=solid_fluid_system.modified_world_space_rigid_inertia_tensor_inverse(i)*B.rigid_V(i).angular;}}

        V.V=B.V;V.rigid_V=B.rigid_V;
        V.V*=(T)-1;V.rigid_V*=(T)-1;
        if(solids)
            solid_system->Set_Global_Boundary_Conditions(V,X_save,rigid_X_save,rigid_rotation_save,rigid_velocity_save,rigid_angular_momentum_save,V_save,
                solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);

        if(fluids){
            static CONJUGATE_RESIDUAL<T> cr;
            static SYMMQMR<T> symmqmr;
            KRYLOV_SOLVER<T>* solver=0;
            const char* solver_name=0;
            if(solids_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_cg){
                PHYSBAM_WARNING("Cannot use CG for the solid fluid coupling equations; using CONJUGATE_RESIDUAL instead.");solver=&cr;solver_name="CONJUGATE_RESIDUAL";}
            else if(solids_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_cr){solver=&cr;solver_name="CONJUGATE_RESIDUAL";}
            else if(solids_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_symmqmr){solver=&symmqmr;solver_name="SYMMQMR";}
            LOG::Time(solver_name);

            solver->nullspace_tolerance=(T)0;
            solver->restart_iterations=solids_parameters.implicit_solve_parameters.cg_restart_iterations;
            solver->print_diagnostics=solid_body_collection.print_diagnostics;
            solver->print_residuals=solid_body_collection.print_residuals;
            solver->iterations_used=&solid_body_collection.iterations_used_diagnostic;

            solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;

            // compute the preconditioner
            POISSON_COLLIDABLE_UNIFORM<GRID<TV> >& poisson=*Get_Poisson();
            for(int i=1;i<=A_array.m;i++)
                A_array(i).Construct_Incomplete_Cholesky_Factorization(poisson.pcg.modified_incomplete_cholesky,poisson.pcg.modified_incomplete_cholesky_coefficient,
                    poisson.pcg.preconditioner_zero_tolerance,poisson.pcg.preconditioner_zero_replacement); // check to see if the blocks can be preconditioned even though the whole

            ar_full.Resize(particles.array_collection->Size(),false,false);rigid_ar_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
            z_full.Resize(particles.array_collection->Size(),false,false);rigid_z_full.Resize(rigid_body_particles.array_collection->Size(),false,false);
            GENERALIZED_VELOCITY<TV> ar_V(ar_full,rigid_ar_full,solid_body_collection),z_V(z_full,rigid_z_full,solid_body_collection);
            PRESSURE_VELOCITY_VECTOR<TV> V_coupled(V,x_array.v),F_coupled(F,p_array.v),R_coupled(R,r_array.v),S_coupled(S,ap_array.v),B_coupled(B,b_array),ar_coupled(ar_V,ar_array.v),
                z_coupled(z_V,z_array.v);
            LOG::Time(solver_name);
            static int solve_id=0;
            solve_id++;
            if(print_matrix_rhs_and_solution) OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("bo-%i.txt",solve_id).c_str()).Write("bo",B_coupled);
            if(!solver->Solve(solid_fluid_system,V_coupled,B_coupled,F_coupled,S_coupled,ar_coupled,R_coupled,z_coupled,solids_parameters.implicit_solve_parameters.cg_tolerance,1,solids_parameters.implicit_solve_parameters.cg_iterations))
                PHYSBAM_DEBUG_WRITE_SUBSTEP("FAILED CONVERGENCE",0,1);
            if(print_matrix_rhs_and_solution){
                OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("xo-%i.txt",solve_id).c_str()).Write("xo",V_coupled);
                OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("matrixo-%i.txt",solve_id).c_str()).Write("MO",solid_fluid_system,S_coupled,F_coupled);}}
        else{
            static CONJUGATE_GRADIENT<T> cg;
            static CONJUGATE_RESIDUAL<T> cr;
            static SYMMQMR<T> symmqmr;
            KRYLOV_SOLVER<T>* solver=0;
            const char* solver_name=0;
            if(solids_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_cg){solver=&cg;solver_name="CG";}
            else if(solids_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_cr){solver=&cr;solver_name="CONJUGATE_RESIDUAL";}
            else if(solids_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_symmqmr){solver=&symmqmr;solver_name="SYMMQMR";}

            solver->nullspace_tolerance=(T)0;
            solver->restart_iterations=100;
            solver->print_diagnostics=solid_body_collection.print_diagnostics;
            solver->print_residuals=solid_body_collection.print_residuals;
            solver->iterations_used=&solid_body_collection.iterations_used_diagnostic;

            LOG::Time(solver_name);
            if(!solver->Solve(*solid_system,V,B,F,S,ar,R,z,solids_parameters.implicit_solve_parameters.cg_tolerance,1,solids_parameters.implicit_solve_parameters.cg_iterations) && solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure)
                throw std::runtime_error("Backward Euler Failed");}
        LOG::Stop_Time();

        rigid_deformable_collisions->rigid_body_collisions.Remove_Contact_Joints();}

    if(fluids){
        T_ARRAYS_SCALAR& p=Get_Pressure();
        // Copy pressures back into pressure array
        for(int i=1;i<=x_array.v.m;i++) for(int j=1;j<=x_array.v(i).n;j++) p(matrix_index_to_cell_index_array(i)(j))=x_array.v(i)(j);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("unscaled final pressures (sf coupled evolution)",0,1);

        // scale pressure back to get a real pressure
        if(dt && !fluids_parameters.compressible) p*=1/dt;
        if(dt && fluids_parameters.compressible) fluids_parameters.euler->euler_projection.Unscale_Pressure_By_Dt(dt);}
    if(solids && Simulate_Fluids()){
        for(int i=1;i<=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m;i++){
            RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(solid_body_collection.rigid_body_collection.simulated_rigid_body_particles(i));
            V.rigid_V(i).linear-=TV::Cross_Product(V.rigid_V(i).angular,rigid_body_updated_center_of_mass(i)-rigid_body.X());}}

    if(solids && velocity_update && solids_parameters.use_post_cg_constraints){ // save final force for friction calculation
        solid_system->Force(V,F);
        if(Simulate_Fluids()){ // TODO(jontg): Ask Craig if what is done here makes any sense...
            for(int i=1;i<=F.V.Size();i++) F.V(i)*=dt*mass.one_over_mass(i);} // save final force for friction calculation
        else{
            for(int i=1;i<=F.V.Size();i++) V.V(i)=B.V(i)+dt*mass.one_over_mass(i)*F.V(i);
            for(int i=1;i<=F.rigid_V.Size();i++) V.rigid_V(i)=B.rigid_V(i)+mass.world_space_rigid_mass_inverse(i)*F.rigid_V(i)*dt;}}

    if(solids){
        solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities();
        rigid_body_collection.Update_Angular_Momentum(solid_body_collection.rigid_body_collection.simulated_rigid_body_particles);} // make sure L = I omega


    if(solids_fluids_parameters.mpi_solid_fluid){
        if(solids) solid_system_mpi->Send_Generalized_Velocity_To_Fluid(V);
        else fluid_system_mpi->Get_Generalized_Velocity_From_Solid(V);}

    if(!solids)
        rigid_body_collection.Update_Angular_Momentum(solid_body_collection.rigid_body_collection.simulated_rigid_body_particles); // make sure L = I omega
    
    for(int i=1;i<=twist.Size();i++){rigid_body_particles.V(i)=twist(i).linear;rigid_body_particles.angular_velocity(i)=twist(i).angular;}

    delete solid_system_mpi;delete fluid_system_mpi;delete solid_system;
}
//#####################################################################
// Function Transfer_Momentum_And_Set_Boundary_Conditions
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION<TV>::
Transfer_Momentum_And_Set_Boundary_Conditions(const T time,GENERALIZED_VELOCITY<TV>* B)
{
    T_FACE_ARRAYS_SCALAR& face_velocities=Get_Face_Velocities();
    POISSON_COLLIDABLE_UNIFORM<GRID<TV> >& poisson=*Get_Poisson();
    const GRID<TV>& grid=Get_Grid();
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    ARRAY<int> deformable_simplices;
    ARRAY<PAIR<int,int> > rigid_simplices;
    ARRAY<T> deformable_simplex_forces;
    ARRAY<T> rigid_simplex_jet_velocities;
    T_FACE_ARRAYS_BOOL modified_face(grid);modified_face.Fill(false);
    T_FACE_ARRAYS_SCALAR modified_face_weight(grid);modified_face_weight.Fill((T)0);
    T_FACE_ARRAYS_SCALAR modified_face_velocities(grid);modified_face_velocities.Fill((T)0);
    // TODO: rewind velocities if we modify them here, otherwise may be doing extra contribution
    TV orientation;
    if(solids_evolution_callbacks->Get_Solid_Source_Velocities(deformable_simplices,deformable_simplex_forces,rigid_simplices,rigid_simplex_jet_velocities,orientation,time)){
        // let's say it gives a box of a certain size, set velocities on that? orientation needs to be correct, too.
        // TODO: add F to B appropriately
        ARRAY<T_THIN_SHELL_SIMPLEX> velocity_lines;
        ARRAY<TV> jet_velocities;
        for(int rigid_simplex=1;rigid_simplex<=rigid_simplices.m;rigid_simplex++){
            int rigid_body_id=rigid_simplices(rigid_simplex).x;
            int rigid_simplex_index=rigid_simplices(rigid_simplex).y;
            RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(rigid_body_id);

            T_THIN_SHELL_SIMPLEX velocity_line(rigid_body.simplicial_object->Get_Element(rigid_simplex_index));
            for(int i=1;i<=TV::dimension;i++)
                velocity_line.X(i)=rigid_body.Frame()*velocity_line.X(i);
            velocity_lines.Append(velocity_line);
            jet_velocities.Append(rigid_simplex_jet_velocities(rigid_simplex)*orientation);}

        //TODO (mlentine): This only works with the first collision body?
        DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>* collisions=0;
        for(COLLISION_GEOMETRY_ID i(1);i<=fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.bodies.m;i++)
            if(fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Is_Active(i)){
                collisions=dynamic_cast<DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>*>(fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.bodies(i));if(collisions) break;}
        for(int deformable_simplex=1;deformable_simplex<=deformable_simplices.m;deformable_simplex++){
            T_THIN_SHELL_SIMPLEX velocity_line(collisions->object.Get_Element(deformable_simplices(deformable_simplex)));velocity_lines.Append(velocity_line);
            jet_velocities.Append(deformable_simplex_forces(deformable_simplex)*orientation);}

        assert(velocity_lines.m==jet_velocities.m);
        for(int i=1;i<=velocity_lines.m;i++){
            T_THIN_SHELL_SIMPLEX velocity_line=velocity_lines(i);
            TV jet_velocity=jet_velocities(i);
            RANGE<TV> box=velocity_line.Bounding_Box();
            RANGE<TV_INT> bounding_grid_cells(grid.Clamp_To_Cell(box.min_corner,1),grid.Clamp_To_Cell(box.max_corner,1));
            bounding_grid_cells.Change_Size(1);
            // wonder if this works
            TV total_length_by_dimension;
            for(int axis=1;axis<=TV::dimension;axis++){
                for(FACE_ITERATOR iterator(grid,bounding_grid_cells,axis);iterator.Valid();iterator.Next()){
                    RANGE<TV> dual_cell=iterator.Dual_Cell();
                    const TV_INT& face_index=iterator.Face_Index();
                    if(((!rigid_body_dual_cell_weights.Component(axis).Valid_Index(face_index) || !rigid_body_dual_cell_weights(axis,face_index)) && 
                        (!dual_cell_weights.Component(axis).Valid_Index(face_index) || !dual_cell_weights(axis,face_index))) || !poisson.psi_N(axis,face_index)) continue;
                    ARRAY<T_THIN_SHELL_SIMPLEX> clipped_simplices;
                    velocity_line.Clip_To_Box(dual_cell,clipped_simplices);
                    T total_length=0;
                    for(int i=1;i<=clipped_simplices.m;i++)
                        total_length+=clipped_simplices(i).Size();
                    if(total_length){
                        if(!modified_face(axis,face_index)){
                            modified_face(axis,face_index)=true;
                            modified_face_velocities(axis,face_index)=(T)0;}
                        //dual_cell_contains_solid(axis,face_index)=false;
                        modified_face_velocities(axis,face_index)+=jet_velocity(axis)*total_length;
                        modified_face_weight(axis,face_index)+=total_length;}}}}}

    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        const int axis=iterator.Axis();
        const TV_INT& face_index=iterator.Face_Index();
        if(modified_face(axis,face_index)){
            modified_face_velocities(axis,face_index)/=modified_face_weight(axis,face_index);
            face_velocities(axis,face_index)+=modified_face_velocities(axis,face_index);
            dual_cell_fluid_volume(axis,face_index)=(T)0;}}
}
//#####################################################################
// Function Set_Neumann_and_Dirichlet_Boundary_Conditions
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION<TV>::
Set_Neumann_and_Dirichlet_Boundary_Conditions(const T time)
{
    POISSON_COLLIDABLE_UNIFORM<GRID<TV> >& poisson=*Get_Poisson();
    T_FACE_ARRAYS_SCALAR& face_velocities=Get_Face_Velocities();
    T_FACE_ARRAYS_BOOL::Copy(dual_cell_contains_solid,poisson.psi_N);poisson.psi_D.Fill(false);
    fluids_parameters.Set_Domain_Boundary_Conditions(poisson,face_velocities,time);

    // We aren't actually sourcing here, since these velocities get reset later
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before sourcing (sf coupled evolution)",0,1);  // TODO(jontg): What's the right thing to do here in the one-way coupled case?
    if(fluids_parameters.solids_override_source_velocities) fluids_parameters.callbacks->Get_Source_Velocities_Masked(time,dual_cell_contains_solid);
    else fluids_parameters.callbacks->Get_Source_Velocities(face_velocities,poisson.psi_N,time);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after sourcing (sf coupled evolution)",0,1);

    fluids_parameters.callbacks->Set_Dirichlet_Boundary_Conditions(time);

    // Don't use Neumann faces for coupling if they have dirichlet conditions on both sides - prevent nasty mass lumping bug
    for(FACE_ITERATOR iterator(Get_Grid(),0,GRID<TV>::INTERIOR_REGION);iterator.Valid();iterator.Next()){
        const int axis=iterator.Axis();const TV_INT face_index=iterator.Face_Index();
        if(poisson.psi_N(axis,face_index) && ((poisson.psi_D(iterator.First_Cell_Index()) && poisson.psi_D(iterator.Second_Cell_Index()))
    // Also, don't couple to thin-shell faces, as euler is having difficulties dealing with them.
           || (fluids_parameters.euler && !poisson.psi_D(iterator.First_Cell_Index()) && !poisson.psi_D(iterator.Second_Cell_Index())))){
            poisson.psi_N(axis,face_index)=false;dual_cell_contains_solid(axis,face_index)=false;}}
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION<TV>::
Set_Dirichlet_Boundary_Conditions(const T time)
{
    POISSON_COLLIDABLE_UNIFORM<GRID<TV> >& poisson=*Get_Poisson();
    T_ARRAYS_SCALAR& p=Get_Pressure();
    const GRID<TV>& grid=Get_Grid();

    // Set all cells inside a solid to be dirichlet
    COLLISION_GEOMETRY_ID inside_body_id;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();TV location=iterator.Location();
        if(fluids_parameters.collision_bodies_affecting_fluid->Inside_Any_Body(location,inside_body_id)){
            poisson.psi_D(cell_index)=true;p(cell_index)=(T)-100;}} // These dirichlet pressures should never be used.

    if(fluids_parameters.mpi_grid) fluids_parameters.mpi_grid->Exchange_Boundary_Cell_Data(poisson.psi_D,1);
    if(!fluids_parameters.solve_neumann_regions){ // Turn off cells which should not be solved for...
        poisson.Find_Solution_Regions();
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT cell_index=iterator.Cell_Index();int color=poisson.filled_region_colors(cell_index);
            if(!poisson.psi_D(cell_index) && (color==-1 || !poisson.filled_region_touches_dirichlet(color))){
                poisson.psi_D(cell_index)=true;p(cell_index)=-100;}}}

    if(fluids_parameters.mpi_grid) fluids_parameters.mpi_grid->Exchange_Boundary_Cell_Data(poisson.psi_D,1);
    if(!fluids_parameters.compressible) poisson.beta_plus=poisson.beta_minus=(T)1/fluids_parameters.density;
    
    // TODO(jontg): I'm pretty sure this call is extraneous...
    poisson.Set_Relative_Tolerance((T)1e-7);
    if(!fluids_parameters.compressible){poisson.beta_face.Fill((T)1/fluids_parameters.density);poisson.beta_given_on_faces=true;}
}
//#####################################################################
// Function Compute_W
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION<TV>::
Compute_W(const T current_position_time)
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particle;
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    T_FACE_ARRAYS_SCALAR& face_velocities=Get_Face_Velocities();
    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >& collision_bodies_affecting_fluid=*fluids_parameters.collision_bodies_affecting_fluid;
    const GRID<TV>& grid=Get_Grid();
    rigid_body_count=0;kinematic_rigid_bodies.Remove_All();

    T_ARRAYS_STRUCTURE_SIMPLEX_LIST structure_simplex_list(grid.Domain_Indices(1));
    for(COLLISION_GEOMETRY_ID i(1);i<=collision_bodies_affecting_fluid.collision_geometry_collection.bodies.m;i++)
        if(collision_bodies_affecting_fluid.collision_geometry_collection.Is_Active(i)){
            COLLISION_GEOMETRY<TV>& body=*collision_bodies_affecting_fluid.collision_geometry_collection.bodies(i);
            if(!body.active) continue;
            T_THIN_SHELL_MESH* mesh=0;
            if(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>* deformable=dynamic_cast<DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>*>(&body))
                mesh=&deformable->object.mesh;
            else if(RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(&body)){
                ++rigid_body_count;
                const RIGID_GEOMETRY<TV>& rigid_geometry=rigid_collision_geometry->rigid_geometry;
                if(dynamic_cast<const RIGID_BODY<TV>&>(rigid_geometry).Has_Infinite_Inertia() || !fluids_parameters.fluid_affects_solid){
                    kinematic_rigid_bodies.Append(rigid_geometry.particle_index);}
                const T_THIN_SHELL* thin_shell=rigid_geometry.simplicial_object;
                if(!thin_shell) PHYSBAM_NOT_IMPLEMENTED();
                mesh=&thin_shell->mesh;}
            else PHYSBAM_FATAL_ERROR();
            for(int e=1;e<=mesh->elements.m;e++){
                const RANGE<TV>& box=body.World_Space_Simplex(e).Bounding_Box();
                RANGE<TV_INT> bounding_grid_cells(grid.Clamp_To_Cell(box.min_corner,1),grid.Clamp_To_Cell(box.max_corner,1));
                for(CELL_ITERATOR iterator(grid,bounding_grid_cells);iterator.Valid();iterator.Next())
                    structure_simplex_list(iterator.Cell_Index()).Append(PAIR<COLLISION_GEOMETRY_ID,int>(i,e));}}

    // =====================================================================================
    // DEFORMABLE SETUP
    // =====================================================================================
    // populates total_nodal_volume and dual_cell_weights
    // TODO: optimize this by first clipping all simplices against all the face dual cell face planes and storing barycentric weights for the clipped simplices
    // NOTE: assume non-embedded collision surface for now
    ARRAY<int> particles_to_dynamic_particles_map(particles.array_collection->Size());
    particles_to_dynamic_particles_map.Subset(solid_body_collection.deformable_body_collection.dynamic_particles)=
        IDENTITY_ARRAY<int>(solid_body_collection.deformable_body_collection.dynamic_particles.m);

    // =====================================================================================
    // FLUID SECTION
    // =====================================================================================
    bool using_levelset=fluids_parameters.particle_levelset_evolution!=0;

    // TODO: this rasterization will not function for objects which are not in the collision body list!!! static bodies are included in this
    collision_bodies.Rasterize_Objects();

    dual_cell_contains_solid.Fill(false);
    collision_bodies.Compute_Psi_N_Zero_Velocity(dual_cell_contains_solid,0);
    POISSON_COLLIDABLE_UNIFORM<GRID<TV> >* poisson=Get_Poisson();
    fluids_parameters.Set_Domain_Boundary_Conditions(*poisson,face_velocities,current_position_time);
    const T one_over_number_nodes_thin_shell=(T)1/TV::dimension;
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){const int axis=iterator.Axis();const TV_INT face_index=iterator.Face_Index();
        const TV axis_vector=TV::Axis_Vector(axis);const TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
        const RANGE<TV> dual_cell=iterator.Dual_Cell().Thickened(grid.Minimum_Edge_Length()*(T)1e-2);
        if(!using_levelset || Negative(grid,axis,face_index,fluids_parameters.particle_levelset_evolution->particle_levelset.levelset.phi)){
            // compute solid volume in dual cell
            RANGE<TV> clamped_dual_cell=RANGE<TV>::Intersect(dual_cell,grid.domain);

            dual_cell_fluid_volume(axis,face_index)=clamped_dual_cell.Size();
            dual_cell_fluid_volume(axis,face_index)=max((T)0,dual_cell_fluid_volume(axis,face_index));}  // TODO(jontg): This is no longer necessary
        else dual_cell_fluid_volume(axis,face_index)=(T)0;

        if(dual_cell_weights(axis,face_index)) dual_cell_weights(axis,face_index)->Remove_All();
        if(rigid_body_dual_cell_weights(axis,face_index)) rigid_body_dual_cell_weights(axis,face_index)->Remove_All();

        ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> > dual_cell_structure_simplex_list=structure_simplex_list(first_cell_index);
        dual_cell_structure_simplex_list.Append_Unique_Elements(structure_simplex_list(second_cell_index));

        if(!dual_cell_fluid_volume(axis,face_index) || !dual_cell_structure_simplex_list.m) continue; // TODO: alternatively, if fluid volume is sufficiently small, clamp it

        // this needs to be surface simplices
        ARRAY<T_THIN_SHELL_SIMPLEX> clipped_simplices_thin_shell;

        T overall_weight=(T)0;
        if(!dual_cell_weights(axis,face_index)) dual_cell_weights(axis,face_index)=new FACE_WEIGHT_ELEMENTS();
        if(!rigid_body_dual_cell_weights(axis,face_index)) rigid_body_dual_cell_weights(axis,face_index)=new FACE_WEIGHT_ELEMENTS();
        FACE_WEIGHT_ELEMENTS& deformable_face_weights=*dual_cell_weights(axis,face_index);
        FACE_WEIGHT_ELEMENTS& rigid_dual_cell_weights=*rigid_body_dual_cell_weights(axis,face_index);

        for(int s=1;s<=dual_cell_structure_simplex_list.m;s++){
            const COLLISION_GEOMETRY_ID object_id=dual_cell_structure_simplex_list(s).x;
            const int object_simplex_index=dual_cell_structure_simplex_list(s).y;
            
            COLLISION_GEOMETRY<TV>& body=*collision_bodies_affecting_fluid.collision_geometry_collection.bodies(object_id);
            if(object_simplex_index>0){
                const T_THIN_SHELL_SIMPLEX& simplex=body.World_Space_Simplex(object_simplex_index);
                simplex.Clip_To_Box(dual_cell,clipped_simplices_thin_shell); // TODO: don't clip against every plane twice
                if(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>* deformable=dynamic_cast<DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>*>(&body)){
                    TV total_weight;
                    for(int i=1;i<=clipped_simplices_thin_shell.m;i++){
                        const T_THIN_SHELL_SIMPLEX& clipped_simplex=clipped_simplices_thin_shell(i);
                        total_weight+=simplex.Sum_Barycentric_Coordinates(clipped_simplex)*clipped_simplex.Size()*one_over_number_nodes_thin_shell;}
                    const T_THIN_SHELL_ELEMENT& simplex_indices=deformable->object.mesh.elements(object_simplex_index);
                    T sum_total_weight=total_weight.Sum();
                    overall_weight+=sum_total_weight;
                    if(sum_total_weight){
                        for(int p=1;p<=simplex_indices.m;p++){
                            int dynamic_particle_index=particles_to_dynamic_particles_map(simplex_indices(p)); // TODO: make this work for embedded particles
                            if(dynamic_particle_index) deformable_face_weights.Append(PAIR<int,T>(dynamic_particle_index,total_weight(p)));}}}
                else if(RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry=dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>*>(&body)){
                    T total_weight=(T)0;
                    for(int i=1;i<=clipped_simplices_thin_shell.m;i++) total_weight+=clipped_simplices_thin_shell(i).Size();
                    overall_weight+=total_weight;
                    if(total_weight){
                        const RIGID_GEOMETRY<TV>& rigid_geometry=rigid_collision_geometry->rigid_geometry;
                        rigid_dual_cell_weights.Append(PAIR<int,T>(rigid_geometry.particle_index,total_weight));}}
                else PHYSBAM_FATAL_ERROR();}
            else{ // TODO: non-simplicial volumetric rigid body
                PHYSBAM_NOT_IMPLEMENTED();
            }}

        // normalize weights for this dual cell
        if(dual_cell_contains_solid(axis,face_index) && overall_weight){
            rigid_dual_cell_weights.template Project<T,&PAIR<int,T>::y>()*=(T)1/overall_weight;
            deformable_face_weights.template Project<T,&PAIR<int,T>::y>()*=(T)1/overall_weight;}
        else{
            rigid_dual_cell_weights.Remove_All();deformable_face_weights.Remove_All();
            dual_cell_contains_solid(axis,face_index)=false;
        }

        if(using_levelset && fluids_parameters.fluid_affects_solid){
            ARRAY<TWIST<TV> > twist;twist.Resize(rigid_body_particles.V.Size());
            for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(rigid_body_particles.V(i),rigid_body_particles.angular_velocity(i));
            GENERALIZED_VELOCITY<TV> V(particles.V,twist,solid_body_collection);
            // check whether this face can see fluid
            // cast a ray left and right
            // TODO: must exchange solid velocities to do this! for parallel.  Alt, just do solid explicit part on all procs.
            T_ARRAYS_SCALAR& phi=fluids_parameters.particle_levelset_evolution->particle_levelset.levelset.phi;
            ARRAY<COLLISION_GEOMETRY_ID> objects;collision_bodies.objects_in_cell.Get_Objects_For_Cells(first_cell_index,second_cell_index,collision_bodies.collision_geometry_collection.bodies.m,objects);if(!objects.m) continue;
            COLLISION_GEOMETRY_ID body_id;bool occluded=true;
            RAY<TV> ray(iterator.First_Cell_Center(),TV::Axis_Vector(axis),true);ray.t_max=(T).5*grid.dX[axis];ray.semi_infinite=false;
            if(phi(first_cell_index)<=0 && !collision_bodies.Intersection_With_Any_Simplicial_Object(ray,body_id,&objects)) occluded=false;
            else{
                ray.Initialize(iterator.Second_Cell_Center(),-TV::Axis_Vector(axis),true);ray.t_max=(T).5*grid.dX[axis];ray.semi_infinite=false;
                if(phi(second_cell_index)<=0 && !collision_bodies.Intersection_With_Any_Simplicial_Object(ray,body_id,&objects)) occluded=false;}
            if(occluded){
                T& velocity=face_velocities(axis,face_index);
                velocity=(T)0;
                if(dual_cell_weights(axis,face_index)){
                    FACE_WEIGHT_ELEMENTS& face_weights=*dual_cell_weights(axis,face_index);
                    for(int i=1;i<=face_weights.m;i++) velocity+=face_weights(i).y*V.V(face_weights(i).x)(axis);}
                if(rigid_body_dual_cell_weights(axis,face_index)){
                    FACE_WEIGHT_ELEMENTS& face_weights=*rigid_body_dual_cell_weights(axis,face_index);
                    for(int i=1;i<=face_weights.m;i++) velocity+=face_weights(i).y*solid_body_collection.rigid_body_collection.Rigid_Body(face_weights(i).x).Pointwise_Object_Velocity(iterator.Location())(axis);}}    
            for(int i=1;i<=twist.Size();i++){rigid_body_particles.V(i)=twist(i).linear;rigid_body_particles.angular_velocity(i)=twist(i).angular;}}}
}
//#####################################################################
// Function Compute_Coupling_Terms_Deformable
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION<TV>::
Compute_Coupling_Terms_Deformable(const T_ARRAYS_INT& cell_index_to_matrix_index,const ARRAY<INTERVAL<int> >& interior_regions,const int colors)
{
    ARRAY<ARRAY<int> > row_counts(colors);
    for(int i=1;i<=colors;i++) row_counts(i).Resize(solid_body_collection.deformable_body_collection.dynamic_particles.Size()*TV::dimension); // TODO: fix // TODO: is this right?

    // have face weights.  Want to set up W^T * \nabla.  dual cells correspond to rows of W.
    // stencil for grad p looks like (p_i+1 - p_i)/dx
    // W^T * grad p 

    J_deformable.Resize(colors);

    nodal_fluid_mass.Resize(solid_body_collection.deformable_body_collection.dynamic_particles.m);ARRAYS_COMPUTATIONS::Fill(nodal_fluid_mass,T_DIAGONAL_MATRIX());
    if(!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Fluid_Node()){
        POISSON_COLLIDABLE_UNIFORM<GRID<TV> >* poisson=Get_Poisson();
        const GRID<TV>& grid=Get_Grid();
        RANGE<TV_INT> domain_indices=grid.Domain_Indices();
        // estimate the number of non-zero entries per row - DEFORMABLE PART
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            if(!dual_cell_weights(axis,face_index)) continue;
            FACE_WEIGHT_ELEMENTS& dual_cell_weight=*dual_cell_weights(axis,face_index);
            TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
            int left_column_index=-1,right_column_index=-1;
            int left_column_color=poisson->filled_region_colors(first_cell_index),right_column_color=poisson->filled_region_colors(second_cell_index);
            if(left_column_color>0 && first_cell_index(axis)>=domain_indices.min_corner(axis))
                left_column_index=cell_index_to_matrix_index(first_cell_index);
            if(right_column_color>0 && second_cell_index(axis)<=domain_indices.max_corner(axis))
                right_column_index=cell_index_to_matrix_index(second_cell_index);
            for(int i=1;i<=dual_cell_weight.m;i++){
                int j_row=TV::dimension*(dual_cell_weight(i).x-1)+axis;
                if(left_column_index>0) row_counts(left_column_color)(j_row)++;
                if(right_column_index>0) row_counts(right_column_color)(j_row)++;
                if(!fluids_parameters.compressible){
                    const T dual_cell_fluid_mass=Get_Density_At_Face(axis,face_index)*dual_cell_fluid_volume(axis,face_index);
                    nodal_fluid_mass(dual_cell_weight(i).x)(axis)+=dual_cell_fluid_mass*dual_cell_weight(i).y;}}}

        // assume we've done everything in terms of dynamic particles
        for(int i=1;i<=colors;i++){
            J_deformable(i).n=interior_regions(i).Size()+1;
            J_deformable(i).Set_Row_Lengths(row_counts(i));
            PROJECTED_ARRAY<ARRAY<SPARSE_MATRIX_ENTRY<T> >,T_PROJECTED_INDEX> J_deformableiAj=J_deformable(i).A.template Project<int,&SPARSE_MATRIX_ENTRY<T>::j>();
            ARRAYS_COMPUTATIONS::Fill(J_deformableiAj,0);
            PROJECTED_ARRAY<ARRAY<SPARSE_MATRIX_ENTRY<T> >,T_PROJECTED_VALUE> J_deformableiAa=J_deformable(i).A.template Project<T,&SPARSE_MATRIX_ENTRY<T>::a>();
            ARRAYS_COMPUTATIONS::Fill(J_deformableiAa,(T)0);}}

    if(!fluids_parameters.compressible && solids_fluids_parameters.mpi_solid_fluid){ // TODO see if we can change this to a single direction send
        ARRAY<T> serial_nodal_fluid_mass(nodal_fluid_mass.m*TV::m),global_nodal_fluid_mass(nodal_fluid_mass.m*TV::m);
        for(int node=1;node<=nodal_fluid_mass.m;node++) for(int axis=1;axis<=TV::m;axis++) serial_nodal_fluid_mass((node-1)*TV::m+axis)=nodal_fluid_mass(node)(axis);
        solids_fluids_parameters.mpi_solid_fluid->Reduce_Add(serial_nodal_fluid_mass,global_nodal_fluid_mass);
        for(int node=1;node<=nodal_fluid_mass.m;node++) for(int axis=1;axis<=TV::m;axis++) nodal_fluid_mass(node)(axis)=global_nodal_fluid_mass((node-1)*TV::m+axis);}

    // populate the entries of J
    if(!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Fluid_Node()){
        POISSON_COLLIDABLE_UNIFORM<GRID<TV> >* poisson=Get_Poisson();
        const GRID<TV>& grid=Get_Grid();
        RANGE<TV_INT> domain_indices=grid.Domain_Indices();
        TV one_over_dx=grid.one_over_dX;
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            if(!dual_cell_weights(axis,face_index)) continue;
            FACE_WEIGHT_ELEMENTS& dual_cell_weight=*dual_cell_weights(axis,face_index);
            TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
            T left_column_weight=(T)0,right_column_weight=(T)0;
            int left_column_index=-1,right_column_index=-1;
            int left_column_color=poisson->filled_region_colors(first_cell_index),right_column_color=poisson->filled_region_colors(second_cell_index);
            if(left_column_color>0 && first_cell_index(axis)>=domain_indices.min_corner(axis)){
                left_column_weight=(T)-one_over_dx(axis);left_column_index=cell_index_to_matrix_index(first_cell_index);}
            if(right_column_color>0 && second_cell_index(axis)<=domain_indices.max_corner(axis)){
                right_column_weight=(T)one_over_dx(axis);right_column_index=cell_index_to_matrix_index(second_cell_index);}
            for(int i=1;i<=dual_cell_weight.m;i++){
                int j_row=TV::dimension*(dual_cell_weight(i).x-1)+axis;T weight=dual_cell_weight(i).y;
                if(left_column_index>0) J_deformable(left_column_color).Add_Element(j_row,left_column_index,weight*left_column_weight);
                if(right_column_index>0) J_deformable(right_column_color).Add_Element(j_row,right_column_index,weight*right_column_weight);}}
        for(int i=1;i<=colors;i++){SPARSE_MATRIX_FLAT_MXN<T> temp(J_deformable(i));temp.Compress(J_deformable(i));}}
}
//#####################################################################
// Function Compute_Coupling_Terms_Rigid
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION<TV>::
Compute_Coupling_Terms_Rigid(const T_ARRAYS_INT& cell_index_to_matrix_index,const ARRAY<INTERVAL<int> >& interior_regions,const int colors)
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particle;
    ARRAY<ARRAY<int> > row_counts(colors),kinematic_row_counts(colors);
    ARRAY<int> rigid_body_particles_to_dynamic_rigid_body_particles_map(rigid_body_particles.array_collection->Size());
    rigid_body_particles_to_dynamic_rigid_body_particles_map.Subset(kinematic_rigid_bodies)=IDENTITY_ARRAY<int>(kinematic_rigid_bodies.m);

    for(int i=1;i<=colors;i++) kinematic_row_counts(i).Resize(rigid_body_count*rows_per_rigid_body);
    J_rigid_kinematic.Resize(colors);

    if(fluids_parameters.fluid_affects_solid){
        rigid_body_particles_to_dynamic_rigid_body_particles_map.Subset(solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles)=IDENTITY_ARRAY<int>(solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m);
        for(int i=1;i<=colors;i++) row_counts(i).Resize(solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m*rows_per_rigid_body); // TODO: fix
        rigid_body_fluid_mass.Resize(solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m);ARRAYS_COMPUTATIONS::Fill(rigid_body_fluid_mass,T_DIAGONAL_MATRIX());
        rigid_body_updated_center_of_mass.Resize(solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m);ARRAYS_COMPUTATIONS::Fill(rigid_body_updated_center_of_mass,TV());
        J_rigid.Resize(colors);}

    if(!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Fluid_Node()){
        POISSON_COLLIDABLE_UNIFORM<GRID<TV> >* poisson=Get_Poisson();
        const GRID<TV>& grid=Get_Grid();
        // estimate the number of non-zero entries per row - RIGID PART
        // think of this as a separate block of the matrix - call it J_rigid
        RANGE<TV_INT> domain_indices=grid.Domain_Indices();
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            if(!rigid_body_dual_cell_weights(axis,face_index)) continue;
            FACE_WEIGHT_ELEMENTS& dual_cell_weight=*rigid_body_dual_cell_weights(axis,face_index);
            TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
            int left_column_index=-1,right_column_index=-1;
            int left_column_color,right_column_color;
            if((left_column_color=poisson->filled_region_colors(first_cell_index))>0 && first_cell_index(axis)>=domain_indices.min_corner(axis))
                left_column_index=cell_index_to_matrix_index(first_cell_index);
            if((right_column_color=poisson->filled_region_colors(second_cell_index))>0 && second_cell_index(axis)<=domain_indices.max_corner(axis))
                right_column_index=cell_index_to_matrix_index(second_cell_index);
            // each one goes into four rows in 3d (with nonzero weight)
            for(int i=1;i<=dual_cell_weight.m;i++){
                int rigid_body_id=dual_cell_weight(i).x;
                RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(rigid_body_id);
                const int base_row=rows_per_rigid_body*(rigid_body_particles_to_dynamic_rigid_body_particles_map(rigid_body_id)-1);
                if(!fluids_parameters.fluid_affects_solid || rigid_body.Has_Infinite_Inertia()){
                    if(left_column_index>0){
                        kinematic_row_counts(left_column_color)(base_row+axis)++;
                        for(int j=1;j<=T_SPIN::dimension;j++) kinematic_row_counts(left_column_color)(base_row+TV::dimension+j)++;}
                    if(right_column_index>0){
                        kinematic_row_counts(right_column_color)(base_row+axis)++;
                        for(int j=1;j<=T_SPIN::dimension;j++) kinematic_row_counts(right_column_color)(base_row+TV::dimension+j)++;}}
                else{
                    const T weight=dual_cell_weight(i).y;
                    const int rigid_body_index=rigid_body_particles_to_dynamic_rigid_body_particles_map(rigid_body_id);
                    const TV location=iterator.Location();
                    if(!fluids_parameters.compressible){
                        const T dual_cell_fluid_mass_on_body=Get_Density_At_Face(axis,face_index)*dual_cell_fluid_volume(axis,face_index)*weight;
                        rigid_body_fluid_mass(rigid_body_index)(axis)+=dual_cell_fluid_mass_on_body;
                        rigid_body_updated_center_of_mass(rigid_body_index)(axis)+=location(axis)*dual_cell_fluid_mass_on_body;}
                    if(left_column_index>0){
                        row_counts(left_column_color)(base_row+axis)++;kinematic_row_counts(left_column_color)(base_row+axis)++;
                        for(int j=1;j<=T_SPIN::dimension;j++){
                            row_counts(left_column_color)(base_row+TV::dimension+j)++;kinematic_row_counts(left_column_color)(base_row+axis)++;}}
                    if(right_column_index>0){
                        row_counts(right_column_color)(base_row+axis)++;kinematic_row_counts(right_column_color)(base_row+axis)++;
                        for(int j=1;j<=T_SPIN::dimension;j++){
                            row_counts(right_column_color)(base_row+TV::dimension+j)++;kinematic_row_counts(right_column_color)(base_row+axis)++;}}}}}

        for(int i=1;i<=J_rigid_kinematic.m;i++){
            J_rigid_kinematic(i).n=interior_regions(i).Size()+1;J_rigid_kinematic(i).Set_Row_Lengths(kinematic_row_counts(i));
            PROJECTED_ARRAY<ARRAY<SPARSE_MATRIX_ENTRY<T> >,T_PROJECTED_INDEX> J_rigidiAj=J_rigid_kinematic(i).A.template Project<int,&SPARSE_MATRIX_ENTRY<T>::j>();ARRAYS_COMPUTATIONS::Fill(J_rigidiAj,0);}
        for(int i=1;i<=J_rigid.m;i++){
            J_rigid(i).n=interior_regions(i).Size()+1;J_rigid(i).Set_Row_Lengths(row_counts(i));
            PROJECTED_ARRAY<ARRAY<SPARSE_MATRIX_ENTRY<T> >,T_PROJECTED_VALUE> J_rigidiAa=J_rigid(i).A.template Project<T,&SPARSE_MATRIX_ENTRY<T>::a>();
            ARRAYS_COMPUTATIONS::Fill(J_rigidiAa,(T)0);}}

    if(!fluids_parameters.compressible && solids_fluids_parameters.mpi_solid_fluid){ // TODO see if we can change this to a single direction send
        ARRAY<T> serial_rigid_body_updated_center_of_mass(rigid_body_updated_center_of_mass.m*TV::m),global_rigid_body_updated_center_of_mass(rigid_body_updated_center_of_mass.m*TV::m);
        ARRAY<T> serial_rigid_body_fluid_mass(rigid_body_fluid_mass.m*TV::m),global_rigid_body_fluid_mass(rigid_body_fluid_mass.m*TV::m);
        for(int rigid_body=1;rigid_body<=rigid_body_updated_center_of_mass.m;rigid_body++)
            for(int axis=1;axis<=TV::m;axis++){
                serial_rigid_body_updated_center_of_mass((rigid_body-1)*TV::m+axis)=rigid_body_updated_center_of_mass(rigid_body)(axis);
                serial_rigid_body_fluid_mass((rigid_body-1)*TV::m+axis)=rigid_body_fluid_mass(rigid_body)(axis);}
        solids_fluids_parameters.mpi_solid_fluid->Reduce_Add(serial_rigid_body_updated_center_of_mass,global_rigid_body_updated_center_of_mass);
        solids_fluids_parameters.mpi_solid_fluid->Reduce_Add(serial_rigid_body_fluid_mass,global_rigid_body_fluid_mass);
        for(int rigid_body=1;rigid_body<=rigid_body_updated_center_of_mass.m;rigid_body++)
            for(int axis=1;axis<=TV::m;axis++){
                rigid_body_updated_center_of_mass(rigid_body)(axis)=global_rigid_body_updated_center_of_mass((rigid_body-1)*TV::m+axis);
                rigid_body_fluid_mass(rigid_body)(axis)=global_rigid_body_fluid_mass((rigid_body-1)*TV::m+axis);}}

    if(fluids_parameters.fluid_affects_solid){
        rigid_body_fluid_inertia.Resize(solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m);
        // get updated center of mass for each rigid body, compute modified inertia tensor<
        for(int i=1;i<=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){
            RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i));
            if(fluids_parameters.compressible) rigid_body_updated_center_of_mass(i)=rigid_body.X();
            else rigid_body_updated_center_of_mass(i)=(rigid_body_fluid_mass(i)+rigid_body.Mass()).Solve_Linear_System(
                    rigid_body_updated_center_of_mass(i)+rigid_body.X()*rigid_body.Mass());}
        ARRAYS_COMPUTATIONS::Fill(rigid_body_fluid_inertia,T_INERTIA_TENSOR());}

    // populate the entries of J_rigid
    if(!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Fluid_Node()){
        POISSON_COLLIDABLE_UNIFORM<GRID<TV> >* poisson=Get_Poisson();
        const GRID<TV>& grid=Get_Grid();
        TV one_over_dx=grid.one_over_dX;
        RANGE<TV_INT> domain_indices=grid.Domain_Indices();
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            if(!rigid_body_dual_cell_weights(axis,face_index)) continue;
            const TV axis_vector=TV::Axis_Vector(axis);
            FACE_WEIGHT_ELEMENTS& dual_cell_weight=*rigid_body_dual_cell_weights(axis,face_index);
            
            TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
            T left_column_weight=(T)0,right_column_weight=(T)0;
            int left_column_index=-1,right_column_index=-1;
            int left_column_color,right_column_color;
            if((left_column_color=poisson->filled_region_colors(first_cell_index))>0 && first_cell_index(axis)>=domain_indices.min_corner(axis)){
                left_column_weight=(T)-one_over_dx(axis);left_column_index=cell_index_to_matrix_index(first_cell_index);}
            if((right_column_color=poisson->filled_region_colors(second_cell_index))>0 && second_cell_index(axis)<=domain_indices.max_corner(axis)){
                right_column_weight=(T)one_over_dx(axis);right_column_index=cell_index_to_matrix_index(second_cell_index);}
            // each one goes into four rows in 3d (with nonzero weight)
            for(int i=1;i<=dual_cell_weight.m;i++){
                int rigid_body_id=dual_cell_weight(i).x;T weight=dual_cell_weight(i).y;
                RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(rigid_body_id);
                TV center_of_mass=((!fluids_parameters.fluid_affects_solid || rigid_body.Has_Infinite_Inertia()) ? rigid_body.X() : rigid_body_updated_center_of_mass(rigid_body_particles_to_dynamic_rigid_body_particles_map(rigid_body_id)));
                TV R=iterator.Location()-center_of_mass; // vector from center of mass to dual cell
                // get inertia contribution, fluid mass contribution
                // inertia
                // TODO: minus outer product of ith column by ith row of r* - form is simple
                T_SPIN cross_product=TV::Cross_Product(R,axis_vector);
                
                const int base_row=rows_per_rigid_body*(rigid_body_particles_to_dynamic_rigid_body_particles_map(rigid_body_id)-1);
                if(!fluids_parameters.fluid_affects_solid || rigid_body.Has_Infinite_Inertia()){
                    if(left_column_index>0){
                        J_rigid_kinematic(left_column_color).Add_Element(base_row+axis,left_column_index,weight*left_column_weight);
                        for(int j=1;j<=T_SPIN::dimension;j++) J_rigid_kinematic(left_column_color).Add_Element(base_row+TV::dimension+j,left_column_index,weight*left_column_weight*cross_product(j));}
                    if(right_column_index>0){
                        J_rigid_kinematic(right_column_color).Add_Element(base_row+axis,right_column_index,weight*right_column_weight);
                        for(int j=1;j<=T_SPIN::dimension;j++) J_rigid_kinematic(right_column_color).Add_Element(base_row+TV::dimension+j,right_column_index,weight*right_column_weight*cross_product(j));}}
                else{
                    if(!fluids_parameters.compressible){
                        const T dual_cell_fluid_mass_on_body=Get_Density_At_Face(axis,face_index)*dual_cell_fluid_volume(axis,face_index)*weight;
                        rigid_body_fluid_inertia(rigid_body_particles_to_dynamic_rigid_body_particles_map(rigid_body_id))+=dual_cell_fluid_mass_on_body*T_INERTIA_TENSOR::Outer_Product(cross_product);}
                    if(left_column_index>0){
                        J_rigid(left_column_color).Add_Element(base_row+axis,left_column_index,weight*left_column_weight);
                        for(int j=1;j<=T_SPIN::dimension;j++) J_rigid(left_column_color).Add_Element(base_row+TV::dimension+j,left_column_index,weight*left_column_weight*cross_product(j));}
                    if(right_column_index>0){
                        J_rigid(right_column_color).Add_Element(base_row+axis,right_column_index,weight*right_column_weight);
                        for(int j=1;j<=T_SPIN::dimension;j++) J_rigid(right_column_color).Add_Element(base_row+TV::dimension+j,right_column_index,weight*right_column_weight*cross_product(j));}}}}}

    if(!fluids_parameters.compressible && solids_fluids_parameters.mpi_solid_fluid){ // TODO see if we can change this to a single direction send
        int inertia_tensor_size=T_INERTIA_TENSOR::m*T_INERTIA_TENSOR::n;
        ARRAY<T> serial_rigid_body_fluid_inertia(rigid_body_fluid_inertia.m*inertia_tensor_size),global_rigid_body_fluid_inertia(rigid_body_fluid_inertia.m*inertia_tensor_size);
        for(int rigid_body=1;rigid_body<=rigid_body_fluid_inertia.m;rigid_body++)
            for(int row=1;row<=T_INERTIA_TENSOR::m;row++)
                for(int col=1;col<=T_INERTIA_TENSOR::n;col++)
                    serial_rigid_body_fluid_inertia((rigid_body-1)*inertia_tensor_size+(row-1)*T_INERTIA_TENSOR::m+col)=rigid_body_fluid_inertia(rigid_body)(row,col);
        solids_fluids_parameters.mpi_solid_fluid->Reduce_Add(serial_rigid_body_fluid_inertia,global_rigid_body_fluid_inertia);
        for(int rigid_body=1;rigid_body<=rigid_body_fluid_inertia.m;rigid_body++)
            for(int row=1;row<=T_INERTIA_TENSOR::m;row++)
                for(int col=1;col<=T_INERTIA_TENSOR::n;col++)
                    rigid_body_fluid_inertia(rigid_body)(row,col)=global_rigid_body_fluid_inertia((rigid_body-1)*inertia_tensor_size+(row-1)*T_INERTIA_TENSOR::m+col);}
    if(fluids_parameters.fluid_affects_solid){
        for(int i=1;i<=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){
            RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i));
            rigid_body_fluid_inertia(i)+=rigid_body.World_Space_Inertia_Tensor(rigid_body_updated_center_of_mass(i));}}

    if(!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Fluid_Node()){
        // Go through and eliminate empty elements
        for(int i=1;i<=J_rigid.m;i++){SPARSE_MATRIX_FLAT_MXN<T> temp(J_rigid(i));temp.Compress(J_rigid(i));}
        for(int i=1;i<=J_rigid_kinematic.m;i++){SPARSE_MATRIX_FLAT_MXN<T> temp(J_rigid_kinematic(i));temp.Compress(J_rigid_kinematic(i));}}
}
//#####################################################################
// Function Apply_Pressure
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION<TV>::
Apply_Pressure(const T dt,const T time)
{
    EULER_PROJECTION_UNIFORM<GRID<TV> >& euler_projection=fluids_parameters.euler->euler_projection;
    T_FACE_ARRAYS_SCALAR& face_velocities=Get_Face_Velocities();
    //TODO(kwatra): move this logic to driver
    if(fluids_parameters.compressible && !fluids_parameters.euler->timesplit) return;
    Apply_Solid_Boundary_Conditions(time,false,face_velocities);
    if(fluids_parameters.compressible){
        Average_Solid_Projected_Face_Velocities_For_Energy_Update(solid_projected_face_velocities_star,face_velocities,face_velocities);
        T_ARRAYS_SCALAR p_ghost;
        euler_projection.Get_Ghost_Pressures(dt,time,euler_projection.elliptic_solver->psi_D,euler_projection.elliptic_solver->psi_N,euler_projection.p,p_ghost);
        T_FACE_ARRAYS_SCALAR p_face;
        euler_projection.Get_Pressure_At_Faces(dt,time,p_ghost,p_face);
        fluids_parameters.euler_solid_fluid_coupling_utilities->Project_Fluid_Pressure_At_Neumann_Faces(p_ghost,p_face); // Bp
        euler_projection.Apply_Pressure(p_ghost,p_face,face_velocities,euler_projection.elliptic_solver->psi_D,euler_projection.elliptic_solver->psi_N,dt,time);}
    else fluids_parameters.incompressible->projection.Apply_Pressure(face_velocities,dt,time,true);
}
//#####################################################################
// Function Add_Nondynamic_Solids_To_Right_Hand_Side
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION<TV>::
Add_Nondynamic_Solids_To_Right_Hand_Side(ARRAY<VECTOR_ND<T> >& right_hand_side,const ARRAY<INTERVAL<int> >& interior_regions,const int colors)
{
    POISSON_COLLIDABLE_UNIFORM<GRID<TV> >* poisson=Get_Poisson();
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particle;
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    ARRAY<TWIST<TV> > twist;twist.Resize(rigid_body_particles.V.Size());
    for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(rigid_body_particles.V(i),rigid_body_particles.angular_velocity(i));
    GENERALIZED_VELOCITY<TV> V(particles.V,twist,solid_body_collection);
    for(int i=1;i<=colors;++i){
        if(poisson->filled_region_touches_dirichlet(i)||poisson->solve_neumann_regions){
            VECTOR_ND<T> right_hand_side_i;right_hand_side_i.Set_Subvector_View(right_hand_side(i),interior_regions(i));
            SOLID_FLUID_SYSTEM<TV,SPARSE_MATRIX_FLAT_NXN<T> >::Add_J_Rigid_Transpose_Times_Velocity(J_rigid_kinematic(i),V,right_hand_side_i);}}
    for(int i=1;i<=twist.Size();i++){rigid_body_particles.V(i)=twist(i).linear;rigid_body_particles.angular_velocity(i)=twist(i).angular;}
}
//#####################################################################
// Function Average_Solid_Projected_Face_Velocities_For_Energy_Update
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION<TV>::
Average_Solid_Projected_Face_Velocities_For_Energy_Update(const T_FACE_ARRAYS_SCALAR& solid_projected_face_velocities_star,const T_FACE_ARRAYS_SCALAR& solid_projected_face_velocities_np1,T_FACE_ARRAYS_SCALAR& face_velocities)
{
    for(FACE_ITERATOR iterator(Get_Grid());iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
        if(dual_cell_contains_solid(axis,face_index)){
            face_velocities.Component(axis)(face_index)=(solid_projected_face_velocities_star.Component(axis)(face_index)+solid_projected_face_velocities_np1.Component(axis)(face_index))*(T).5;}}
}
//#####################################################################
// Function Apply_Solid_Boundary_Conditions
//#####################################################################
template<class TV> void SOLID_FLUID_COUPLED_EVOLUTION<TV>::
Apply_Solid_Boundary_Conditions(const T time,const bool use_pseudo_velocities,T_FACE_ARRAYS_SCALAR& face_velocities)
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particle;
    POISSON_COLLIDABLE_UNIFORM<GRID<TV> >& poisson=*Get_Poisson();

    // TODO: it is possible that we will end up with strange face weight behavior in cells on domain boundaries, depending on how we're setting up our Ws.
    // if that happens, look here.

    // use effective velocities if available
    if(use_pseudo_velocities){
        int state1=COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_OLD_STATE;int state2=COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE;
        // TODO: fix once DEFORMABLE_OBJECT_FLUID_COLLISIONS::saved_states is moved
        DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>* collisions=0;
        for(COLLISION_GEOMETRY_ID i(1);i<=fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.bodies.m;i++)
            if(fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Is_Active(i)){
                collisions=dynamic_cast<DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>*>(fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.bodies(i));if(collisions) break;}
        for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            const int axis=iterator.Axis();const TV_INT face_index=iterator.Face_Index();
            if(dual_cell_contains_solid(axis,face_index)){
                poisson.psi_N(axis,face_index)=true;
                T& velocity=face_velocities(axis,face_index);
                velocity=(T)0;
                if(dual_cell_weights(axis,face_index) && dual_cell_weights(axis,face_index)->m){
                    assert(collisions);
                    FACE_WEIGHT_ELEMENTS& face_weights=*dual_cell_weights(axis,face_index);
                    for(int i=1;i<=face_weights.m;i++) velocity+=face_weights(i).y*collisions->Pointwise_Node_Pseudo_Velocity(solid_body_collection.deformable_body_collection.dynamic_particles(face_weights(i).x),state1,state2)(axis);}
                if(rigid_body_dual_cell_weights(axis,face_index)){
                    FACE_WEIGHT_ELEMENTS& face_weights=*rigid_body_dual_cell_weights(axis,face_index);
                    for(int i=1;i<=face_weights.m;i++){
                        COLLISION_GEOMETRY<TV>& collision_geometry=*fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Get_Collision_Geometry(face_weights(i).x);
                        velocity+=face_weights(i).y*collision_geometry.Pointwise_Object_Pseudo_Velocity(0,iterator.Location(),state1,state2)(axis);}}}}}
    else{
        for(FACE_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){
            const int axis=iterator.Axis();const TV_INT face_index=iterator.Face_Index();
            if(dual_cell_contains_solid(axis,face_index)){
                poisson.psi_N(axis,face_index)=true;
                T& velocity=face_velocities(axis,face_index);
                velocity=(T)0;
                if(dual_cell_weights(axis,face_index)){
                    FACE_WEIGHT_ELEMENTS& face_weights=*dual_cell_weights(axis,face_index);
                    ARRAY<TWIST<TV> > twist;twist.Resize(rigid_body_particles.V.Size());
                    for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(rigid_body_particles.V(i),rigid_body_particles.angular_velocity(i));
                    GENERALIZED_VELOCITY<TV> V(particles.V,twist,solid_body_collection);
                    for(int i=1;i<=face_weights.m;i++) velocity+=face_weights(i).y*V.V(face_weights(i).x)(axis);
                    for(int i=1;i<=twist.Size();i++){rigid_body_particles.V(i)=twist(i).linear;rigid_body_particles.angular_velocity(i)=twist(i).angular;}}
                if(rigid_body_dual_cell_weights(axis,face_index)){
                    FACE_WEIGHT_ELEMENTS& face_weights=*rigid_body_dual_cell_weights(axis,face_index);
                    for(int i=1;i<=face_weights.m;i++) velocity+=face_weights(i).y*solid_body_collection.rigid_body_collection.Rigid_Body(face_weights(i).x).Pointwise_Object_Velocity(iterator.Location())(axis);}}}}

    // add after
    Transfer_Momentum_And_Set_Boundary_Conditions(time,0);
}
//#####################################################################
template class SOLID_FLUID_COUPLED_EVOLUTION<VECTOR<float,1> >;
template class SOLID_FLUID_COUPLED_EVOLUTION<VECTOR<float,2> >;
template class SOLID_FLUID_COUPLED_EVOLUTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SOLID_FLUID_COUPLED_EVOLUTION<VECTOR<double,1> >;
template class SOLID_FLUID_COUPLED_EVOLUTION<VECTOR<double,2> >;
template class SOLID_FLUID_COUPLED_EVOLUTION<VECTOR<double,3> >;
#endif
