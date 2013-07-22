//#####################################################################
// Copyright 2006-2007, Nipun Kwatra, Frank Losasso, Nick Rasmussen, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPH_EVOLUTION_UNIFORM
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UTILITIES.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Incompressible_Flows/SPH_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SPH_CALLBACKS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> SPH_EVOLUTION_UNIFORM<T_GRID>::
SPH_EVOLUTION_UNIFORM(T_GRID& grid_input,INCOMPRESSIBLE_UNIFORM<T_GRID>& incompressible_input,FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters_input,
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>* particle_levelset_evolution_input)
    :grid(grid_input.Get_MAC_Grid()),callbacks(0),incompressible(incompressible_input),fluids_parameters(fluids_parameters_input),particle_levelset_evolution(particle_levelset_evolution_input),
    particle_radius(1.5),target_particles_per_unit_volume((T)1000),ballistic_particles_per_unit_volume(0),use_analytic_divergence(false),use_analytic_divergence_for_expansion_only(false),
    adjust_cell_weights_on_neumann_boundaries(false),enforce_density_near_interface(true),particle_targeting_time(1),divergence_expansion_multiplier(0),grid_to_particle_slip_multiplier(1),
    neumann_boundary_slip_multiplier((T)0),attraction_strength(0),attraction_restlength_cells((T).2),attraction_forward_move_multiplier((T)1),attraction_skip_number(0),
    stable_attraction_integration(true),uniform_attraction_force(false),flip_ratio((T)1),use_variable_density_solve(false),use_two_way_coupling(false),
    convert_particles_to_fluid(false),maximum_phi_refinement_depth(1),radius(0),radius_plus_half_dx_squared(0),one_over_radius_squared(FLT_MAX),radius_vector(TV()),target_particles_per_cell(0),
    target_minus_ballistic_particles_per_cell(0),one_over_target_particles_per_cell(0),one_over_target_minus_ballistic_particles_per_cell(0),ballistic_particles_per_cell(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> SPH_EVOLUTION_UNIFORM<T_GRID>::
~SPH_EVOLUTION_UNIFORM()
{}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T_GRID> void SPH_EVOLUTION_UNIFORM<T_GRID>::
Euler_Step(const T dt,const T time)
{
    std::stringstream ss;
    ss<<"Number of SPH particles: "<<sph_particles.array_collection->Size()<<std::endl;
    LOG::filecout(ss.str());
    for(int i=sph_particles.array_collection->Size();i>=1;i--){
        callbacks->Adjust_SPH_Particle_For_Domain_Boundaries(sph_particles,i,sph_particles.V(i),dt,time);
        if(!callbacks->Adjust_SPH_Particle_For_Objects(sph_particles,i,sph_particles.V(i),dt,time)) sph_particles.array_collection->Delete_Element(i);
        else sph_particles.X(i)+=dt*sph_particles.V(i);}
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR SPH_EVOLUTION_UNIFORM<T_GRID>::
CFL() const
{
    T max_particle_velocity=0;
    if(particle_levelset_evolution) for(NODE_ITERATOR node_iterator(grid);node_iterator.Valid();node_iterator.Next()){
        TV_INT block=node_iterator.Node_Index();
        if(particle_levelset_evolution->particle_levelset.removed_negative_particles(block)){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*particle_levelset_evolution->particle_levelset.removed_negative_particles(block);
            for(int p=1;p<=particles.array_collection->Size();p++) max_particle_velocity=max(max_particle_velocity,particles.V(p).Max_Abs());}}
    else for(int i=1;i<=sph_particles.array_collection->Size();i++)max_particle_velocity=max(max_particle_velocity,sph_particles.V(i).Max_Abs());
    if(max_particle_velocity) return 3*grid.Minimum_Edge_Length()/max_particle_velocity;
    else return 1e8;
}
//#####################################################################
// Function Copy_Particle_Attributes_From_Array
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_PARTICLES> void SPH_EVOLUTION_UNIFORM<T_GRID>::
Copy_Particle_Attributes_From_Array(T_ARRAYS_PARTICLES& particles)
{
    int particle_index=0;sph_particles.array_collection->Delete_All_Elements();
    for(NODE_ITERATOR iterator(grid,2);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();
        if(particles(block)){sph_particles.array_collection->Add_Elements(particles(block)->array_collection->Size());
            for(int p=1;p<=particles(block)->array_collection->Size();p++){particle_index++;sph_particles.X(particle_index)=particles(block)->X(p);sph_particles.V(particle_index)=particles(block)->V(p);}}}
}
//#####################################################################
// Function Copy_Particle_Attributes_To_Array
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_PARTICLES> void SPH_EVOLUTION_UNIFORM<T_GRID>::
Copy_Particle_Attributes_To_Array(T_ARRAYS_PARTICLES& particles) const
{
    int particle_index=0;
    for(NODE_ITERATOR iterator(grid,2);iterator.Valid();iterator.Next()){TV_INT block=iterator.Node_Index();
        if(particles(block)) for(int p=1;p<=particles(block)->array_collection->Size();p++){particle_index++;
            particles(block)->X(p)=sph_particles.X(particle_index);particles(block)->V(p)=sph_particles.V(particle_index);}}
}
//#####################################################################
// Function Make_Incompressible
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_PARTICLES> void SPH_EVOLUTION_UNIFORM<T_GRID>::
Make_Incompressible(T_ARRAYS_PARTICLES& particles,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    Copy_Particle_Attributes_From_Array(particles);
    Make_Incompressible(face_velocities,dt,time);Copy_Particle_Attributes_To_Array(particles);
}
//#####################################################################
// Function Make_Incompressible
//#####################################################################
template<class T_GRID> void SPH_EVOLUTION_UNIFORM<T_GRID>::
Make_Incompressible(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    LOG::SCOPE scope("SPH_INCOMPRESSIBLE","sph incompressible (time=%f,dt=%f)",time,dt);
    Set_Up_For_Projection(face_velocities,time);
    LOG::Time("make divergence free");
    incompressible.projection.Make_Divergence_Free(face_velocities,dt,time);
    Postprocess_Particles(face_velocities,dt,time);
}
//#####################################################################
// Function Calculate_SPH_Constants
//#####################################################################
template<class T_GRID> void SPH_EVOLUTION_UNIFORM<T_GRID>::
Calculate_SPH_Constants()
{
    radius=particle_radius*grid.Minimum_Edge_Length();
    radius_plus_half_dx_squared=sqr(radius+(T)root_three/(T)2*grid.Maximum_Edge_Length());
    one_over_radius_squared=1/sqr(radius);
    radius_vector=(radius+grid.Minimum_Edge_Length())*TV::All_Ones_Vector();
    target_particles_per_cell=target_particles_per_unit_volume*grid.Cell_Size();
    target_minus_ballistic_particles_per_cell=(target_particles_per_unit_volume-ballistic_particles_per_unit_volume)*grid.Cell_Size();
    one_over_target_particles_per_cell=1/target_particles_per_cell;
    one_over_target_minus_ballistic_particles_per_cell=1/target_minus_ballistic_particles_per_cell;
    ballistic_particles_per_cell=ballistic_particles_per_unit_volume*grid.Cell_Size();
}
//#####################################################################
// Function Set_Up_For_Projection
//#####################################################################
template<class T_GRID> void SPH_EVOLUTION_UNIFORM<T_GRID>::
Set_Up_For_Projection(T_FACE_ARRAYS_SCALAR& face_velocities,const T time)
{
    Calculate_SPH_Constants();

    LOG::Time("initializing particles in cell");
    particles_in_cell.Resize(grid.Domain_Indices(2));one_over_total_particle_cell_weight.Resize(sph_particles.array_collection->Size());one_over_total_particle_face_weight.Resize(sph_particles.array_collection->Size());

    for(int p=1;p<=sph_particles.array_collection->Size();p++){
        RANGE<TV_INT> particle_cells(grid.Cell(sph_particles.X(p)-radius_vector,2),grid.Cell(sph_particles.X(p)+radius_vector,2));
        for(CELL_ITERATOR iterator(grid,particle_cells);iterator.Valid();iterator.Next()){
            TV_INT cell=iterator.Cell_Index();
            TV X_minus_Xp=grid.Center(cell)-sph_particles.X(p);T distance_squared=X_minus_Xp.Magnitude_Squared();
            T weight=1-distance_squared*one_over_radius_squared;
            if(weight>0)one_over_total_particle_cell_weight(p)+=weight;
            if(distance_squared<radius_plus_half_dx_squared && particles_in_cell.Valid_Index(cell)){
                if(!particles_in_cell(cell).m)particles_in_cell(cell).Preallocate(15);
                particles_in_cell(cell).Append(p);}}
        for(int axis=1;axis<=T_GRID::dimension;axis++)
            for(FACE_ITERATOR iterator(grid,particle_cells+RANGE<TV_INT>(TV_INT(),TV_INT::Axis_Vector(axis)),axis);iterator.Valid();iterator.Next()){
                TV_INT face=iterator.Face_Index();
                TV X_minus_Xp=grid.Face(axis,face)-sph_particles.X(p);T distance_squared=X_minus_Xp.Magnitude_Squared();
                T weight=1-distance_squared*one_over_radius_squared;
                if(weight>0)one_over_total_particle_face_weight(p)+=weight;}}
    for(int p=1;p<=sph_particles.array_collection->Size();p++){
        if(one_over_total_particle_cell_weight(p))one_over_total_particle_cell_weight(p)=1/one_over_total_particle_cell_weight(p);
        if(one_over_total_particle_face_weight(p))one_over_total_particle_face_weight(p)=1/one_over_total_particle_face_weight(p);}

    LOG::Time("initializing");
    PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection=incompressible.projection;
    T_FACE_ARRAYS_SCALAR preset_velocities=face_velocities;
    particle_velocities.Resize(grid.Domain_Indices(3),false,false);particle_velocities.Fill(0);
    valid_particle_face_velocities.Resize(grid.Domain_Indices(3),false,false);valid_particle_face_velocities.Fill(false);
    cell_weight=T_ARRAYS_SCALAR(grid.Domain_Indices(1));face_weight=T_FACE_ARRAYS_SCALAR(grid,1);
    
    LOG::Time("rasterize velocities to grid");
    Rasterize_Velocities_To_Grid(particle_velocities,cell_weight,face_weight);
    for(FACE_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){
        int axis=iterator.Axis();TV_INT face=iterator.Face_Index();
        if(face_weight(axis,face)) particle_velocities(axis,face)/=face_weight(axis,face);}

    if(use_two_way_coupling){
        LOG::Time("adding levelset weight to cells"); // TODO: Determine the correct value for interface cells when using 2nd order cut cell
        for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) if(!projection.elliptic_solver->psi_D(iterator.Cell_Index()))
            cell_weight(iterator.Cell_Index())=target_particles_per_cell;}

    if(fluids_parameters.mpi_grid){ // TODO: Assert that received cell weights are within tolerance of those generated locally
        LOG::Time("syncing sph cell weights across MPI boundaries"); 
        fluids_parameters.mpi_grid->Exchange_Boundary_Cell_Data(cell_weight,1);}    

    LOG::Time("doing something with the density");
    callbacks->Do_Something_With_Density(grid,cell_weight);

    LOG::Time("copying fluid velocities");
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        int axis=iterator.Axis();TV_INT face=iterator.Face_Index();TV_INT cell_1=iterator.First_Cell_Index(),cell_2=iterator.Second_Cell_Index();
        if(projection.elliptic_solver->psi_D(cell_1)&&projection.elliptic_solver->psi_D(cell_2)) face_velocities(axis,face)=particle_velocities(axis,face);}

    if(use_variable_density_solve){projection.poisson->Set_Variable_beta();projection.poisson->variable_beta.Fill(0);}

    LOG::Time("setting divergence variables");
    projection.Use_Non_Zero_Divergence(true);projection.Use_Divergence_Multiplier(true);
    projection.divergence_multiplier.Fill(0);projection.divergence.Fill(0);
    T_ARRAYS_BOOL cells_valid(grid.Domain_Indices(1));cells_valid.Fill(true); 
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(projection.elliptic_solver->psi_D(cell)) for(int axis=1;axis<=T_GRID::dimension;axis++) for(int k=0;k<=1;k++){
            TV_INT face(cell);face[axis]+=k;if(!face_weight(axis,face)) cells_valid(cell)=false;}}
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(cells_valid(cell) && cell_weight(cell)>ballistic_particles_per_cell) projection.elliptic_solver->psi_D(cell)=false;
        if(projection.divergence.Valid_Index(cell)) Set_Divergence_And_Multiplier(cell,cells_valid,time);
        if(cell_weight(cell) && !projection.elliptic_solver->psi_D(cell) && use_variable_density_solve)
            projection.poisson->variable_beta(cell)=1/(cell_weight(cell)*one_over_target_particles_per_cell);}

    if(use_variable_density_solve){LOG::Time("setting betas for poisson solve");projection.poisson->Find_Variable_beta();}

    LOG::Time("restoring object and water boundaries");
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        int axis=iterator.Axis();TV_INT face=iterator.Face_Index();
        if(projection.elliptic_solver->psi_N(axis,face)) face_velocities(axis,face)=preset_velocities(axis,face);}
    for(int axis=1;axis<=T_GRID::dimension;axis++)for(int axis_side=1;axis_side<=2;axis_side++){int side=2*(axis-1)+axis_side;
        if(!fluids_parameters.mpi_grid || !fluids_parameters.mpi_grid->Neighbor(axis,axis_side)){
            TV_INT interior_cell_offset=axis_side==1?TV_INT():-TV_INT::Axis_Vector(axis);
            for(FACE_ITERATOR iterator(grid,1,T_GRID::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Face_Index()+interior_cell_offset;
                projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;}}}
}
//#####################################################################
// Function Set_Divergence_And_Multiplier
//#####################################################################
template<class T_GRID> void SPH_EVOLUTION_UNIFORM<T_GRID>::
Set_Divergence_And_Multiplier(const TV_INT cell,const T_ARRAYS_BOOL& cells_valid,const T time)
{
    Calculate_SPH_Constants();

    PROJECTION_DYNAMICS_UNIFORM<T_GRID>& projection=incompressible.projection;
    if(projection.elliptic_solver->psi_D(cell)){projection.divergence_multiplier(cell)=1;projection.divergence(cell)=0;return;}
    
    if(!use_analytic_divergence){
        T normalized_cell_weight=max((T)0,cell_weight(cell)-ballistic_particles_per_cell)/target_minus_ballistic_particles_per_cell;
        projection.divergence_multiplier(cell)=min((T)1,normalized_cell_weight);
        projection.divergence(cell)=max((T)0,normalized_cell_weight-1)*divergence_expansion_multiplier;}
    else{
        T target_density_factor=callbacks->Target_Density_Factor(grid.Center(cell),time);bool near_interface=false;
        if(!enforce_density_near_interface){
            RANGE<TV_INT> neighbor_cells(grid.Clamp_To_Cell(grid.Center(cell)-radius_vector,1),grid.Clamp_To_Cell(grid.Center(cell)+radius_vector,1));
            for(CELL_ITERATOR iterator_neighbor(grid,neighbor_cells);iterator_neighbor.Valid();iterator_neighbor.Next()){
                if(!cells_valid(iterator_neighbor.Cell_Index()) || cell_weight(iterator_neighbor.Cell_Index())<=ballistic_particles_per_cell){near_interface=true;break;}}}
        if(!near_interface){
            projection.divergence_multiplier(cell)=1; 
            if(cell_weight(cell)){ // TODO: This divergence may be missing an advection term
                if(use_analytic_divergence_for_expansion_only)
                    projection.divergence(cell)=max((T)0,(T)log(cell_weight(cell)*one_over_target_minus_ballistic_particles_per_cell/target_density_factor)/particle_targeting_time);
                else projection.divergence(cell)=log(cell_weight(cell)*one_over_target_minus_ballistic_particles_per_cell/target_density_factor)/particle_targeting_time;}}}
}
//#####################################################################
// Function Postprocess_Particles
//#####################################################################
template<class T_GRID> void SPH_EVOLUTION_UNIFORM<T_GRID>::
Postprocess_Particles(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    Calculate_SPH_Constants();

    LOG::Time("calculate particle deltas");
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face=iterator.Face_Index();
        particle_velocities(axis,face)-=face_velocities(axis,face);}
    ARRAY<TV> delta_velocity(sph_particles.array_collection->Size());ARRAY<TV> delta_weight(sph_particles.array_collection->Size());
    Calculate_Particle_Deltas(particle_velocities,delta_velocity,delta_weight);

    if(neumann_boundary_slip_multiplier){ // TODO:  Probably can yank this if we use adjust_cell_weight_for
        LOG::Time("neumann slip multiplier");
        T_ARRAYS_INT cell_faces_blocked(grid.Domain_Indices(1));
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            if(incompressible.projection.elliptic_solver->psi_N(iterator.Axis(),iterator.Face_Index())){
                cell_faces_blocked(iterator.First_Cell_Index())++;cell_faces_blocked(iterator.Second_Cell_Index())++;}}
        T multiplier_over_number_of_faces_per_cell=neumann_boundary_slip_multiplier/(T)T_GRID::number_of_faces_per_cell;
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            cell_weight(iterator.Cell_Index())*=1+cell_faces_blocked(iterator.Cell_Index())*multiplier_over_number_of_faces_per_cell;}}

    LOG::Time("applying deltas");
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T> interpolation;
    for(int p=1;p<=sph_particles.array_collection->Size();p++){
        T delta_slip=grid_to_particle_slip_multiplier*max((T)0,interpolation.Clamped_To_Array_Cell(grid,cell_weight,sph_particles.X(p))-ballistic_particles_per_cell);
        delta_slip=min((T)1,delta_slip/target_minus_ballistic_particles_per_cell);
        for(int axis=1;axis<=T_GRID::dimension;axis++) if(delta_weight(p)[axis]){
        sph_particles.V(p)[axis]+=delta_slip*delta_velocity(p)[axis]/delta_weight(p)[axis];}}

    if(flip_ratio!=1){
        LOG::Time(STRING_UTILITIES::string_sprintf("applying FLIP to PIC with ratio %f",flip_ratio));
        int ghost_cells=3;
        T_FACE_ARRAYS_SCALAR face_velocities_ghost(grid,ghost_cells);
        incompressible.boundary->Fill_Ghost_Cells_Face(grid,face_velocities,face_velocities_ghost,time,ghost_cells);
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face=iterator.Face_Index();
            if(incompressible.projection.elliptic_solver->psi_N(axis,face)){int valid_neighbors=0;T face_velocity=0;
                for(int i=1;i<=T_GRID::number_of_one_ring_neighbors_per_cell;i++){
                    TV_INT face_neighbor=T_GRID::One_Ring_Neighbor(face,i); // TODO: Check that this is correct for one-way coupling
                    if(face_weight(axis,face_neighbor)&&!incompressible.projection.elliptic_solver->psi_N(axis,face_neighbor)){
                        valid_neighbors++;face_velocity+=face_velocities_ghost(axis,face_neighbor);}}
                if(valid_neighbors) face_velocities_ghost(axis,face)=face_velocity/(T)valid_neighbors;}}
        const RANGE<TV>& domain=grid.domain;
        for(int p=1;p<=sph_particles.array_collection->Size();p++)if(domain.Lazy_Inside(sph_particles.X(p))){
            sph_particles.V(p)=flip_ratio*sph_particles.V(p)+(1-flip_ratio)*interpolation.Clamped_To_Array_Face(grid,face_velocities_ghost,sph_particles.X(p));}}

    if(attraction_strength){
        int skip_number=1+attraction_skip_number;
        LOG::Time(STRING_UTILITIES::string_sprintf("attracting particles with skip %d",skip_number));
        T attraction_restlength=attraction_restlength_cells*grid.Minimum_Edge_Length();
        T halfway_between_restlength_and_radius=(T).5*(attraction_restlength+radius);
        const RANGE<TV>& domain=grid.domain;
        ARRAY<TV> forces(sph_particles.array_collection->Size());
        int skip_offset=random.Get_Uniform_Integer(1,skip_number);
        for(int p=skip_offset;p<=sph_particles.array_collection->Size();p+=skip_number){
            TV& X=sph_particles.X(p);if(!domain.Inside(X,(T)1e-8)) continue;
            TV_INT cell=grid.Clamp_To_Cell(sph_particles.X(p));
            T particles_in_cell_ratio=0;
            if(!uniform_attraction_force)particles_in_cell_ratio=interpolation.Clamped_To_Array_Cell(grid,cell_weight,sph_particles.X(p))*one_over_target_particles_per_cell;
            if(particles_in_cell_ratio>=1) continue;
            for(int i=1;i<=particles_in_cell(cell).m;i++){
                int k=particles_in_cell(cell)(i);if(k==p)continue;
                T strength=attraction_strength*(1-particles_in_cell_ratio);
                TV X_minus_Xp=X-sph_particles.X(k);T distance=X_minus_Xp.Magnitude();X_minus_Xp.Normalize();T x=distance-attraction_restlength;T force_magnitude=0;
                // TODO: make force computed with relative deformation for non-constant restlength
                if(distance<halfway_between_restlength_and_radius) force_magnitude=strength*x;
                else if(distance<radius) force_magnitude=strength*(radius-distance);
                else force_magnitude=0;
                TV half_force=(T).5*force_magnitude*X_minus_Xp;
                forces(p)+=-half_force;forces(k)+=half_force;}}
        T dt_squared=sqr(dt);
        if(attraction_forward_move_multiplier==1){
            if(stable_attraction_integration) for(int p=1;p<=sph_particles.array_collection->Size();p++)sph_particles.X(p)+=forces(p)*dt_squared;
            else for(int p=1;p<=sph_particles.array_collection->Size();p++)sph_particles.V(p)+=forces(p)*dt;}
        else if(stable_attraction_integration) for(int p=1;p<=sph_particles.array_collection->Size();p++){
            TV dx=forces(p)*dt_squared;
            TV velocity_normalized=sph_particles.V(p).Normalized();
            TV velocity_component=TV::Dot_Product(dx,velocity_normalized)*velocity_normalized;
            TV tangent_component=dx-velocity_component;
            sph_particles.X(p)+=velocity_component*attraction_forward_move_multiplier+tangent_component;}
        else for(int p=1;p<=sph_particles.array_collection->Size();p++){
            TV dv=forces(p)*dt;
            TV velocity_normalized=sph_particles.V(p).Normalized();
            TV velocity_component=TV::Dot_Product(dv,velocity_normalized)*velocity_normalized;
            TV tangent_component=dv-velocity_component;
            sph_particles.V(p)+=velocity_component*attraction_forward_move_multiplier+tangent_component;}}
    particles_in_cell.Clean_Memory();one_over_total_particle_cell_weight.Resize(0);one_over_total_particle_face_weight.Resize(0);
}
//#####################################################################
// Function Rasterize_Velocities_To_Grid
//#####################################################################
template<class T_GRID> void SPH_EVOLUTION_UNIFORM<T_GRID>::
Rasterize_Velocities_To_Grid(T_FACE_ARRAYS_SCALAR& velocities,T_ARRAYS_SCALAR& cell_weight,T_FACE_ARRAYS_SCALAR& face_weight)
{
    Calculate_SPH_Constants();

    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){
        TV_INT cell=iterator.Cell_Index();TV X=grid.Center(cell);
        for(int k=1;k<=particles_in_cell(cell).m;k++){
            int p=particles_in_cell(cell)(k);
            TV X_minus_Xp=X-sph_particles.X(p);T distance_squared=X_minus_Xp.Magnitude_Squared();
            T weight=(1-distance_squared*one_over_radius_squared)*one_over_total_particle_cell_weight(p); 
            int weight_multiplier=1;if(adjust_cell_weights_on_neumann_boundaries) for(int axis=1;axis<=T_GRID::dimension;axis++){ // TODO: Make this more accurate
                TV_INT face1=grid.First_Face_Index_In_Cell(axis,cell),face2=grid.Second_Face_Index_In_Cell(axis,cell);
                if(incompressible.projection.elliptic_solver->psi_N(axis,face1)) if(sph_particles.X(p)(axis)>X(axis)+.5*grid.dX(axis)) weight_multiplier++;
                if(incompressible.projection.elliptic_solver->psi_N(axis,face2)) if(sph_particles.X(p)(axis)<X(axis)-.5*grid.dX(axis)) weight_multiplier++;}
            if(weight>0)cell_weight(cell)+=weight_multiplier*weight;}}
    for(FACE_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){
        TV_INT face=iterator.Face_Index();int axis=iterator.Axis();TV X=grid.Face(axis,face);TV_INT cell=face;
        for(int k=1;k<=particles_in_cell(cell).m;k++){
            int p=particles_in_cell(cell)(k);
            TV X_minus_Xp=X-sph_particles.X(p);T distance_squared=X_minus_Xp.Magnitude_Squared();
            T weight=(1-distance_squared*one_over_radius_squared)*one_over_total_particle_face_weight(p);
            if(weight>0){
                valid_particle_face_velocities(axis,face)=true;
                velocities(axis,face)+=weight*sph_particles.V(p)[axis];
                face_weight(axis,face)+=weight;}}}
}
//#####################################################################
// Function Calculate_Particle_Deltas
//#####################################################################
template<class T_GRID> void SPH_EVOLUTION_UNIFORM<T_GRID>::
Calculate_Particle_Deltas(const T_FACE_ARRAYS_SCALAR& minus_face_delta,ARRAY<TV>& delta_velocity,ARRAY<TV>& delta_weight)
{
    Calculate_SPH_Constants();

    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT face=iterator.Face_Index();int axis=iterator.Axis();TV X=grid.Face(axis,face);TV_INT cell=face;grid.Clamp(cell);
        for(int k=1;k<=particles_in_cell(cell).m;k++){
            int p=particles_in_cell(cell)(k);
            TV X_minus_Xp=X-sph_particles.X(p);T distance_squared=X_minus_Xp.Magnitude_Squared();
            T weight=1-distance_squared*one_over_radius_squared;
            if(weight>0)if(minus_face_delta.Component(axis).Valid_Index(face)){
                delta_velocity(p)[axis]+=weight*-minus_face_delta(axis,face);
                delta_weight(p)[axis]+=weight;}}}
}
//#####################################################################
// Function Create_Fluid_From_Particles
//#####################################################################
template<class T_GRID> void SPH_EVOLUTION_UNIFORM<T_GRID>::
Modify_Levelset_And_Particles_To_Create_Fluid(const T time,T_FACE_ARRAYS_SCALAR* face_velocities)
{
    Calculate_SPH_Constants();

    LOG::Time("modifying levelset");
    PARTICLE_LEVELSET_UNIFORM<T_GRID>& particle_levelset=particle_levelset_evolution->particle_levelset;
    particle_levelset.Modify_Levelset_Using_Escaped_Particles(face_velocities);

    LOG::Time("creating fluid from particles");
    T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES& removed_negative_particles=particle_levelset.removed_negative_particles;
    T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES& removed_positive_particles=particle_levelset.removed_positive_particles;
    T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& negative_particles=particle_levelset.negative_particles;
    T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& positive_particles=particle_levelset.positive_particles;
    T_ARRAYS_ARRAY_PAIR_TV_INT_INT removed_negative_particles_affecting_cell(grid.Domain_Indices());
    T_ARRAYS_ARRAY_SCALAR one_over_total_removed_negative_particle_cell_weight(grid.Block_Indices(1));
    T_ARRAYS_SCALAR removed_negative_particle_cell_weight(grid.Domain_Indices());
    
    for(NODE_ITERATOR node_iterator(grid,1);node_iterator.Valid();node_iterator.Next()){
        TV_INT block=node_iterator.Node_Index();
        if(removed_negative_particles(block)){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*removed_negative_particles(block);
            one_over_total_removed_negative_particle_cell_weight(block).Resize(particles.array_collection->Size());
            for(int p=1;p<=particles.array_collection->Size();p++){
                PAIR<TV_INT,int> particle_record(block,p);
                TV_INT lower=grid.Clamp_To_Cell(particles.X(p)-radius_vector),upper=grid.Clamp_To_Cell(particles.X(p)+radius_vector);
                for(CELL_ITERATOR cell_iterator(grid,RANGE<TV_INT>(lower,upper));cell_iterator.Valid();cell_iterator.Next()){
                    TV_INT cell=cell_iterator.Cell_Index();
                    TV X_minus_Xp=grid.Center(cell)-particles.X(p);T distance_squared=X_minus_Xp.Magnitude_Squared();
                    T weight=1-distance_squared*one_over_radius_squared;
                    if(weight>0)one_over_total_removed_negative_particle_cell_weight(block)(p)+=weight;
                    if(distance_squared<radius_plus_half_dx_squared){
                        if(!removed_negative_particles_affecting_cell(cell).m)removed_negative_particles_affecting_cell(cell).Preallocate(15);
                        removed_negative_particles_affecting_cell(cell).Append(particle_record);}}}}}

    for(NODE_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){
        TV_INT block=iterator.Node_Index(); if(removed_negative_particles(block)){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*removed_negative_particles(block);
            for(int p=1;p<=particles.array_collection->Size();p++) if(one_over_total_removed_negative_particle_cell_weight(block)(p))
                one_over_total_removed_negative_particle_cell_weight(block)(p)=1/one_over_total_removed_negative_particle_cell_weight(block)(p);}}

    T_ARRAYS_BOOL converting_cells(grid.Domain_Indices(1)),converting_cells_neighborhood(grid.Domain_Indices(1));

    bool created_fluid_from_particles=false;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell=iterator.Cell_Index();TV X=grid.Center(cell);
        for(int k=1;k<=removed_negative_particles_affecting_cell(cell).m;k++){
            PAIR<TV_INT,int>& particle_record=removed_negative_particles_affecting_cell(cell)(k);
            TV X_minus_Xp=X-removed_negative_particles(particle_record.x)->X(particle_record.y);T distance_squared=X_minus_Xp.Magnitude_Squared();
            T weight=(1-distance_squared*one_over_radius_squared)*one_over_total_removed_negative_particle_cell_weight(particle_record.x)(particle_record.y);
            if(weight>0)removed_negative_particle_cell_weight(cell)+=weight;}
        if(removed_negative_particle_cell_weight(cell)>=target_particles_per_cell){created_fluid_from_particles=true;converting_cells(cell)=true;}}

    if(created_fluid_from_particles){
        ARRAYS_UTILITIES<T_GRID,T>::Make_Ghost_Mask_From_Active_Mask(grid.Get_Regular_Grid_At_MAC_Positions(),converting_cells,converting_cells_neighborhood,2,1);
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();converting_cells(cell)=converting_cells(cell)||converting_cells_neighborhood(cell);}
        
        T blending_particle_radius=(T).5*grid.Minimum_Edge_Length();
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){ 
            TV_INT cell=iterator.Cell_Index();
            if(converting_cells(cell)){ // TODO: Use smooth density function to determine if a particle contributes to the levelset
                for(int p=1;p<=removed_negative_particles_affecting_cell(cell).m;p++){
                    PAIR<TV_INT,int>& particle_record=removed_negative_particles_affecting_cell(cell)(p);
                    TV Xp=removed_negative_particles(particle_record.x)->X(particle_record.y);
                    T radius_plus_phi=blending_particle_radius+particle_levelset_evolution->phi(cell);
                    if(radius_plus_phi > 0){TV X=grid.Center(cell);
                        T distance_squared=(X-Xp).Magnitude_Squared();
                        if(distance_squared < sqr(radius_plus_phi)) particle_levelset_evolution->phi(cell)=sqrt(distance_squared)-blending_particle_radius;}}}}}
    
    LOG::Time("reinitializing levelset");
    particle_levelset_evolution->Make_Signed_Distance();
    
    if(created_fluid_from_particles){
        LOG::Time("reclassifying sph particles");
        for(NODE_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){
            TV_INT block=iterator.Node_Index();
            if(positive_particles(block)){
                for(int p=1;p<=positive_particles(block)->array_collection->Size();p++) 
                    if(particle_levelset_evolution->particle_levelset.levelset.Phi(positive_particles(block)->X(p))<0){
                        positive_particles(block)->array_collection->Add_To_Deletion_List(p);
                        if(!removed_positive_particles(block)) removed_positive_particles(block)=particle_levelset.template_removed_particles.Clone();
                        removed_positive_particles(block)->array_collection->Append(*positive_particles(block)->array_collection,p);}
                positive_particles(block)->array_collection->Delete_Elements_On_Deletion_List();}
            if(removed_negative_particles(block)){
                for(int p=1;p<=removed_negative_particles(block)->array_collection->Size();p++)
                    if(particle_levelset_evolution->particle_levelset.levelset.Phi(removed_negative_particles(block)->X(p))<0){
                        removed_negative_particles(block)->array_collection->Add_To_Deletion_List(p);
                        if(!negative_particles(block)) negative_particles(block)=particle_levelset.template_particles.Clone();
                        negative_particles(block)->array_collection->Append(*removed_negative_particles(block)->array_collection,p);}
                removed_negative_particles(block)->array_collection->Delete_Elements_On_Deletion_List();}}}

    LOG::Time("reseeding around new fluid");
    ARRAYS_UTILITIES<T_GRID,T>::Make_Ghost_Mask_From_Active_Mask(grid.Get_Regular_Grid_At_MAC_Positions(),converting_cells,converting_cells_neighborhood,2,1);
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();converting_cells(cell)=converting_cells(cell)||converting_cells_neighborhood(cell);}
    particle_levelset_evolution->particle_levelset.Reseed_Particles(time,&converting_cells);
        
    LOG::Time("modifying levelset");        
    particle_levelset_evolution->particle_levelset.Modify_Levelset_Using_Escaped_Particles(face_velocities);
    LOG::Time("adjusting particle radii");
    particle_levelset_evolution->particle_levelset.Adjust_Particle_Radii();
    LOG::Stop_Time();
}
//#####################################################################
// Function Move_Particles_Off_Grid_Boundaries
//#####################################################################
template<class T_GRID> template<class T_ARRAYS_PARTICLES> void SPH_EVOLUTION_UNIFORM<T_GRID>::
Move_Particles_Off_Grid_Boundaries(T_ARRAYS_PARTICLES& particles,const T tolerance) const
{
    assert(fluids_parameters.mpi_grid);
    for(int axis=1;axis<=T_GRID::dimension;axis++)for(int axis_side=1;axis_side<=2;axis_side++){int side=2*(axis-1)+axis_side;
        if(fluids_parameters.mpi_grid->Neighbor(axis,axis_side)){
            for(NODE_ITERATOR iterator(grid,0,T_GRID::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){
                TV_INT block=iterator.Node_Index();
                if(particles(block)) for(int p=1;p<=particles(block)->array_collection->Size();p++)
                    if(axis_side==1) particles(block)->X(p)[axis]=max(grid.domain.Minimum_Corner()[axis]+tolerance,particles(block)->X(p)[axis]);
                    else particles(block)->X(p)[axis]=min(grid.domain.Maximum_Corner()[axis]-tolerance,particles(block)->X(p)[axis]);}}}
}
//#####################################################################
#define P(...) __VA_ARGS__
#define INSTANTIATION_HELPER(T_GRID) \
    template class SPH_EVOLUTION_UNIFORM<T_GRID >; \
    template void SPH_EVOLUTION_UNIFORM<T_GRID >::Copy_Particle_Attributes_To_Array(GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::REBIND<PARTICLE_LEVELSET_REMOVED_PARTICLES<T_GRID::VECTOR_T>*>::TYPE&) const; \
    template void SPH_EVOLUTION_UNIFORM<T_GRID >::Copy_Particle_Attributes_From_Array(GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::REBIND<PARTICLE_LEVELSET_REMOVED_PARTICLES<T_GRID::VECTOR_T>*>::TYPE&); \
    template void SPH_EVOLUTION_UNIFORM<T_GRID >::Make_Incompressible(GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::REBIND<PARTICLE_LEVELSET_REMOVED_PARTICLES<T_GRID::VECTOR_T>*>::TYPE&,GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS&,const T_GRID::SCALAR,const T_GRID::SCALAR); \
    template void SPH_EVOLUTION_UNIFORM<T_GRID >::Move_Particles_Off_Grid_Boundaries(GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::REBIND<PARTICLE_LEVELSET_REMOVED_PARTICLES<T_GRID::VECTOR_T>*>::TYPE&,const T_GRID::SCALAR) const
INSTANTIATION_HELPER(P(GRID<VECTOR<float,1> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<float,2> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<float,3> >));
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(P(GRID<VECTOR<double,1> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<double,2> >));
INSTANTIATION_HELPER(P(GRID<VECTOR<double,3> >));
#endif
