#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//####################################################################/
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUIDS_PARAMETERS_DYADIC
//#####################################################################
#include <PhysBAM_Tools/Fourier_Transforms_Calculations/TURBULENCE_2D.h>
#include <PhysBAM_Tools/Fourier_Transforms_Calculations/TURBULENCE_3D.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic_Advection/ADVECTION_SEMI_LAGRANGIAN_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_OCTREE_HELPER.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_QUADTREE_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Collisions/GRID_BASED_COLLISION_GEOMETRY_DYADIC.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Dyadic_Level_Sets/READ_WRITE_LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Dyadic_Level_Sets/READ_WRITE_LEVELSET_QUADTREE.h>
#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_PHI_WATER.h>
#include <PhysBAM_Dynamics/Interpolation/FACE_LOOKUP_FIRE_DYADIC.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_CALLBACKS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_DYADIC.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> FLUIDS_PARAMETERS_DYADIC<T_GRID>::
FLUIDS_PARAMETERS_DYADIC(const typename FLUIDS_PARAMETERS<T_GRID>::TYPE type)
    :FLUIDS_PARAMETERS<T_GRID>(type),particle_levelset_evolution(*grid),incompressible(*grid,fire)
{
    if(smoke){
        density_refinement_thresholds.Resize(3);
        density_refinement_thresholds(1)=VECTOR<T,2>((T).3,(T).4);density_refinement_thresholds(2)=VECTOR<T,2>((T).2,(T).8);
        density_refinement_thresholds(3)=VECTOR<T,2>((T).1,2);
        temperature_container.Set_Velocity(&incompressible.projection.face_velocities);
        density_container.Set_Velocity(&incompressible.projection.face_velocities);}
    else if(fire){
        temperature_refinement_thresholds.Resize(4);
        temperature_refinement_thresholds(1)=VECTOR<T,2>(925,975);temperature_refinement_thresholds(2)=VECTOR<T,2>(800,1100);
        temperature_refinement_thresholds(3)=VECTOR<T,2>(600,1300);temperature_refinement_thresholds(4)=VECTOR<T,2>(400,4000);}
    else if(water){
        phi_boundary=&phi_boundary_water; // override default
        phi_boundary_water.Set_Velocity_Pointer(incompressible.projection.face_velocities);}
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> FLUIDS_PARAMETERS_DYADIC<T_GRID>::
~FLUIDS_PARAMETERS_DYADIC()
{}
//#####################################################################
// Function Use_Fluid_Coupling_Defaults
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Use_Fluid_Coupling_Defaults()
{
    solid_affects_fluid=true;fluid_affects_solid=true;
    density_container.Use_Semi_Lagrangian_Collidable_Advection(*collision_bodies_affecting_fluid,incompressible.valid_mask);
    temperature_container.Use_Semi_Lagrangian_Collidable_Advection(*collision_bodies_affecting_fluid,incompressible.valid_mask);
    particle_levelset_evolution.levelset_advection.Use_Semi_Lagrangian_Collidable_Advection(*collision_bodies_affecting_fluid,collidable_phi_replacement_value,incompressible.valid_mask);
    if(!fire) incompressible.Use_Semi_Lagrangian_Collidable_Advection(*collision_bodies_affecting_fluid);
    else PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Use_Fluid_Slip_Coupling_Defaults
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Use_Fluid_Slip_Coupling_Defaults()
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Use_No_Fluid_Coupling_Defaults
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Use_No_Fluid_Coupling_Defaults()
{
    density_container.Set_Custom_Advection(semi_lagrangian);
    temperature_container.Set_Custom_Advection(semi_lagrangian);
    particle_levelset_evolution.levelset_advection.Set_Custom_Advection(semi_lagrangian);
    if(!fire) incompressible.Set_Custom_Advection(semi_lagrangian);
    else PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Initialize_Grids
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Initialize_Grids()
{
    if(!grid->uniform_grid.Is_Isotropic((T)1e-5)) PHYSBAM_FATAL_ERROR("Dyadic cells must be square/cube");
    callbacks->Initialize_Fluids_Grids();
}    
//#####################################################################
// Function Get_Levelset_Velocity
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Get_Levelset_Velocity(const T_GRID& grid,T_LEVELSET& levelset,ARRAY<T>& face_velocities_levelset,const T time) const
{
    if(water) face_velocities_levelset=incompressible.projection.face_velocities; 
    else if(fire) for(int i=1;i<=grid.number_of_nodes;i++){
        const PROJECTION_DYADIC<T_GRID>& projection=incompressible.projection;
        ARRAY<T> flame_speed_face(grid.number_of_faces);T_LINEAR_INTERPOLATION_DYADIC_HELPER::Interpolate_From_Cells_To_Faces(grid,incompressible.flame_speed,flame_speed_face);
        ARRAY<TV> normals_face(grid.number_of_faces);T_LINEAR_INTERPOLATION_DYADIC_HELPER::Interpolate_From_Cells_To_Faces(grid,*levelset.normals,normals_face);
        for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
            TV normal_at_face=normals_face(iterator.Face_Index()).Normalized();
            face_velocities_levelset(iterator.Face_Index())=projection.Face_Velocity_With_Ghost_Value(projection.face_velocities,iterator.Face_Index(),-1)-
                flame_speed_face(iterator.Face_Index())*normal_at_face[iterator.Axis()];}}
}
//#####################################################################
// Function Adjust_Particle_For_Domain_Boundaries
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time)
{
    // remove this for speed - don't call the function for the other particles
    if(particle_type == PARTICLE_LEVELSET_POSITIVE || particle_type == PARTICLE_LEVELSET_REMOVED_POSITIVE) return;

    TV& X=particles.X(index);TV X_new=X+dt*V;
    T max_collision_distance=particle_levelset_evolution.particle_levelset.Particle_Collision_Distance(particles.quantized_collision_distance(index));
    TV min_corner=grid->Domain().Minimum_Corner(),max_corner=grid->Domain().Maximum_Corner();
    for(int axis=1;axis<=T_GRID::dimension;axis++){
        if(domain_walls[axis][1] && X_new[axis] < min_corner[axis]+max_collision_distance){
            T collision_distance=X[axis]-min_corner[axis];if(collision_distance > max_collision_distance) collision_distance=X_new[axis]-min_corner[axis];
            collision_distance=max((T)0,collision_distance);
            X_new[axis]+=max((T)0,min_corner[axis]-X_new[axis]+collision_distance);
            V[axis]=max((T)0,V[axis]);X=X_new-dt*V;}
        if(domain_walls[axis][2] && X_new[axis] > max_corner[axis]-max_collision_distance){
            T collision_distance=max_corner[axis]-X[axis];if(collision_distance > max_collision_distance) collision_distance=max_corner[axis]-X_new[axis];
            collision_distance=max((T)0,collision_distance);
            X_new[axis]-=max((T)0,X_new[axis]-max_corner[axis]+collision_distance);
            V[axis]=min((T)0,V[axis]);X=X_new-dt*V;}}
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Delete_Particles_Inside_Objects(const T time)
{
    Delete_Particles_Inside_Objects(particle_levelset_evolution.particle_levelset.positive_particles,PARTICLE_LEVELSET_POSITIVE,time);
    Delete_Particles_Inside_Objects(particle_levelset_evolution.particle_levelset.negative_particles,PARTICLE_LEVELSET_NEGATIVE,time);
    if(particle_levelset_evolution.particle_levelset.use_removed_positive_particles)
        Delete_Particles_Inside_Objects(particle_levelset_evolution.particle_levelset.removed_positive_particles,PARTICLE_LEVELSET_REMOVED_POSITIVE,time);
    if(particle_levelset_evolution.particle_levelset.use_removed_negative_particles)
        Delete_Particles_Inside_Objects(particle_levelset_evolution.particle_levelset.removed_negative_particles,PARTICLE_LEVELSET_REMOVED_NEGATIVE,time);
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
template<class T_GRID> template<class T_PARTICLES> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Delete_Particles_Inside_Objects(ARRAY<T_PARTICLES*>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time)
{
    ARRAY<T_CELL*>& cell_pointer_from_index=grid->Cell_Pointer_From_Index();
    for(int base_cell=1;base_cell<=grid->number_of_cells;base_cell++)if(particles(base_cell)){
        BLOCK_DYADIC<T_GRID> block(*grid,cell_pointer_from_index(base_cell));
        if(collision_bodies_affecting_fluid->Occupied_Block(block)) Delete_Particles_Inside_Objects(block,*particles(base_cell),particle_type,time);}
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Delete_Particles_Inside_Objects(const BLOCK_DYADIC<T_GRID>& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time)
{
    //T contour_value=particle_type==PARTICLE_LEVELSET_NEGATIVE || particle_type==PARTICLE_LEVELSET_REMOVED_NEGATIVE?-grid.Minimum_Edge_Length():0;
    COLLISION_GEOMETRY_ID body_id;int aggregate_id;
    for(int k=particles.array_collection->Size();k>=1;k--)if(collision_bodies_affecting_fluid->Inside_Any_Simplex_Of_Any_Body(block,particles.X(k),body_id,aggregate_id)) particles.array_collection->Delete_Element(k);
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Update_Fluid_Parameters(const T dt,const T time)
{
    incompressible.Set_Gravity(gravity,gravity_direction);
    incompressible.Set_Body_Force(use_body_force);
    if(use_body_force) callbacks->Get_Body_Force(incompressible.force,dt,time);
    incompressible.projection.Set_Density(density);
    if(fire) PHYSBAM_NOT_IMPLEMENTED();//incompressible.projection.Set_Density_Fuel(density_fuel);
    if(fire){
        PHYSBAM_NOT_IMPLEMENTED();
        /*incompressible.projection.Set_Normal_Flame_Speed(normal_flame_speed);
        if(curvature_flame_speed) incompressible.projection.Set_Curvature_Flame_Speed(curvature_flame_speed);*/}
    incompressible.Set_Viscosity(viscosity);
    if(fire) PHYSBAM_NOT_IMPLEMENTED();//incompressible.Set_Viscosity_Fuel(viscosity_fuel);
    if(!fire){
        //incompressible.Set_Variable_Viscosity(variable_viscosity);
        incompressible.Use_Explicit_Part_Of_Implicit_Viscosity(use_explicit_part_of_implicit_viscosity);
        if(implicit_viscosity && implicit_viscosity_iterations) incompressible.Set_Maximum_Implicit_Viscosity_Iterations(implicit_viscosity_iterations); 
        /*if(variable_viscosity) callbacks->Get_Variable_Viscosity(incompressible.variable_viscosity,time);*/}
    if(use_vorticity_confinement) incompressible.Set_Vorticity_Confinement(confinement_parameter);
    if(fire && use_vorticity_confinement_fuel) PHYSBAM_NOT_IMPLEMENTED();//incompressible.Set_Vorticity_Confinement_Fuel(confinement_parameter_fuel);
    if(!fire){
        incompressible.Use_Variable_Vorticity_Confinement(use_variable_vorticity_confinement);
        if(use_variable_vorticity_confinement) callbacks->Get_Variable_Vorticity_Confinement(incompressible.variable_vorticity_confinement,time);  
        incompressible.projection.Use_Non_Zero_Divergence(use_non_zero_divergence);
        if(use_non_zero_divergence) callbacks->Get_Divergence(incompressible.projection.divergence,dt,time);}
    incompressible.projection.elliptic_solver->Solve_Neumann_Regions(solve_neumann_regions);
}
//#####################################################################
// Function Set_Domain_Boundary_Conditions
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Set_Domain_Boundary_Conditions(const T time)
{
    PROJECTION_DYADIC<T_GRID>& projection=incompressible.projection;
    ARRAY<bool> &psi_D=projection.elliptic_solver->psi_D,&psi_N=projection.elliptic_solver->psi_N;

    for(int axis=1;axis<=T_GRID::dimension;axis++)for(int axis_side=1;axis_side<=2;axis_side++){
        int side=2*(axis-1)+axis_side,interior_cell=2-axis_side,exterior_cell=axis_side-1;
        if(domain_walls[axis][axis_side])
            for(FACE_ITERATOR iterator(*grid,grid->Map_Individual_Side_Boundary_Faces(side));iterator.Valid();iterator.Next()){
                if(!water || particle_levelset_evolution.phi(iterator.Cell(interior_cell)->Cell())<=0){
                    int face=iterator.Face_Index();psi_N(face)=true;projection.face_velocities(face)=0;}
                else{int cell=iterator.Cell(exterior_cell)->Cell();
                    psi_D(cell)=true;projection.p(cell)=0;}}
        else
            for(FACE_ITERATOR iterator(*grid,grid->Map_Individual_Side_Boundary_Faces(side));iterator.Valid();iterator.Next()){
                int cell=iterator.Cell(exterior_cell)->Cell();
                psi_D(cell)=true;projection.p(cell)=0;}}
}
//#####################################################################
// Function Get_Body_Force
//#####################################################################
// flames and smoke have a buoyancy force
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Get_Body_Force(ARRAY<TV>& force,const T dt,const T time)
{
    if(smoke){
        ARRAYS_COMPUTATIONS::Fill(force,TV());
        for(int i=1;i<=grid->number_of_nodes;i++){ // y-direction forces only
            T rho_atm=rho_bottom+(rho_top-rho_bottom)*(grid->uniform_grid.Axis_X(i,2)-grid->uniform_grid.domain.min_corner.y)/(grid->uniform_grid.domain.max_corner.y-grid->uniform_grid.domain.min_corner.y);
            if(density_container.density(i)>density_buoyancy_threshold){
                T density_difference=density_container.density(i)-rho_atm,temperature_difference=temperature_container.temperature(i)-temperature_container.ambient_temperature;
                if(density_difference>0||temperature_difference>0)
                    force(i).y=temperature_buoyancy_constant*temperature_difference-density_buoyancy_constant*density_difference;}}}
    else if(fire){
        ARRAYS_COMPUTATIONS::Fill(force,TV());
        for(int i=1;i<=grid->number_of_nodes;i++){ // y-direction forces only
            if(particle_levelset_evolution.particle_levelset.levelset.phi(i)<0) force(i).y=temperature_buoyancy_constant*(temperature_fuel-temperature_container.ambient_temperature); // fuel
            else force(i).y=temperature_buoyancy_constant*(temperature_container.temperature(i)-temperature_container.ambient_temperature);}} // products 
}
//#####################################################################
// Function Get_Neumann_And_Dirichlet_Boundary_Conditions
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Get_Neumann_And_Dirichlet_Boundary_Conditions(const T dt,const T time)
{
    LAPLACE_COLLIDABLE_DYADIC<T_GRID>& elliptic_solver=*incompressible.projection.elliptic_solver;
    POISSON_COLLIDABLE_DYADIC<T_GRID>* poisson=incompressible.projection.poisson;
    ARRAYS_COMPUTATIONS::Fill(elliptic_solver.psi_N,false);ARRAYS_COMPUTATIONS::Fill(elliptic_solver.psi_D,false);
    if(poisson) ARRAYS_COMPUTATIONS::Fill(poisson->beta_face,(T)1/density);
    Set_Domain_Boundary_Conditions(time);
    callbacks->Get_Source_Velocities(time);
    callbacks->Get_Object_Velocities(incompressible.projection,incompressible.projection.face_velocities,dt,time);
    callbacks->Set_Dirichlet_Boundary_Conditions(time);
}
//#####################################################################
// Function Blend_In_External_Velocity
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Blend_In_External_Velocity(const T dt,const T time)
{  
    ARRAY<T>& face_velocities=incompressible.projection.face_velocities;
    if(use_external_velocity){
        ARRAY<TV> V_blend(grid->number_of_cells);ARRAY<T> blend(grid->number_of_cells);callbacks->Get_External_Velocity(V_blend,blend,time);
        for(FACE_ITERATOR iterator(*grid,grid->Map_Regular_Faces());iterator.Valid();iterator.Next()){int axis=iterator.Axis()+1;T& face_velocity=face_velocities(iterator.Face_Index());
            T b=(T).5*(blend(iterator.First_Cell_Index())+blend(iterator.Second_Cell_Index()));b=(T)1-pow(max((T)1-b,(T)0),dt);
            face_velocity=(1-b)*face_velocity+b*(T).5*(V_blend(iterator.First_Cell_Index())[axis]+V_blend(iterator.Second_Cell_Index())[axis]);}}
    if(kolmogorov){
        while(time > turbulence.time_end) turbulence.Advance_Turbulence();
        T b=(T)1-pow(max((T)1-kolmogorov,(T)0),dt),fraction=(time-turbulence.time_start)/(turbulence.time_end-turbulence.time_start);
        for(FACE_ITERATOR iterator(*grid,grid->Map_Regular_Faces());iterator.Valid();iterator.Next()){int axis=iterator.Axis()+1;T& face_velocity=face_velocities(iterator.Face_Index());
            face_velocity=(1-b)*face_velocity+b*turbulence.Turbulent_Face_Velocity(axis,iterator.Location(),fraction);}}
}
//#####################################################################
// Function Move_Grid
//#####################################################################        
template<class T> static void
Move_Grid_Helper(FLUIDS_PARAMETERS_DYADIC<QUADTREE_GRID<T> >& fluids_parameters)
{
    typedef VECTOR<T,2> TV;
    QUADTREE_GRID<T>& grid=*fluids_parameters.grid;
    PARTICLE_LEVELSET_EVOLUTION_DYADIC<QUADTREE_GRID<T> >& particle_levelset_evolution=fluids_parameters.particle_levelset_evolution;

    bool shift_left(false),shift_right(false),shift_bottom(false),shift_top(false);
    T x_min_plus_dx(grid.uniform_grid.domain.min_corner.x+6*grid.uniform_grid.dX.x),x_max_minus_dx(grid.uniform_grid.domain.max_corner.x-6*grid.uniform_grid.dX.x),
        y_min_plus_dy(grid.uniform_grid.domain.min_corner.y+6*grid.uniform_grid.dX.y),y_max_minus_dy(grid.uniform_grid.domain.max_corner.y-6*grid.uniform_grid.dX.y);
    for(int i=1;i<=grid.number_of_nodes;i++)if(particle_levelset_evolution.phi(i)<0){
        VECTOR<T,2> node_location=grid.Node_Location(i);
        if(node_location.x<x_min_plus_dx)shift_left=true;if(node_location.x>x_max_minus_dx)shift_right=true;
        if(node_location.y<y_min_plus_dy)shift_bottom=true;if(node_location.y>y_max_minus_dy)shift_top=true;}
    if(shift_left && shift_right){shift_left=false;shift_right=false;}
    if(shift_bottom && shift_top){shift_bottom=false;shift_top=false;}
    
    RANGE<TV> domain=grid.uniform_grid.domain;
    if(shift_left){
        domain.min_corner.x-=2*grid.uniform_grid.dX.x;domain.max_corner.x-=2*grid.uniform_grid.dX.x;grid.Move_Contents_Right_Two();
        grid.uniform_grid.Initialize(grid.uniform_grid.counts,domain);
        grid.Tree_Topology_Changed();}
    if(shift_right){
        domain.min_corner.x+=2*grid.uniform_grid.dX.x;domain.max_corner.x+=2*grid.uniform_grid.dX.x;grid.Move_Contents_Left_Two();
        grid.uniform_grid.Initialize(grid.uniform_grid.counts,domain);
        grid.Tree_Topology_Changed();}
    if(shift_bottom){
        domain.min_corner.y-=2*grid.uniform_grid.dX.y;domain.max_corner.y-=2*grid.uniform_grid.dX.y;grid.Move_Contents_Up_Two();
        grid.uniform_grid.Initialize(grid.uniform_grid.counts,domain);
        grid.Tree_Topology_Changed();}
    if(shift_top){
        domain.min_corner.y+=2*grid.uniform_grid.dX.y;domain.max_corner.y+=2*grid.uniform_grid.dX.y;grid.Move_Contents_Down_Two();
        grid.uniform_grid.Initialize(grid.uniform_grid.counts,domain);
        grid.Tree_Topology_Changed();}
}
template<class T> static void
Move_Grid_Helper(FLUIDS_PARAMETERS_DYADIC<OCTREE_GRID<T> >& fluids_parameters)
{   
    typedef VECTOR<T,3> TV;
    OCTREE_GRID<T>& grid=*fluids_parameters.grid;
    PARTICLE_LEVELSET_EVOLUTION_DYADIC<OCTREE_GRID<T> >& particle_levelset_evolution=fluids_parameters.particle_levelset_evolution;

    bool shift_left(false),shift_right(false),shift_bottom(false),shift_top(false),shift_front(false),shift_back(false);
    T x_min_plus_dx(grid.uniform_grid.domain.min_corner.x+6*grid.uniform_grid.dX.x),x_max_minus_dx(grid.uniform_grid.domain.max_corner.x-6*grid.uniform_grid.dX.x),
        y_min_plus_dy(grid.uniform_grid.domain.min_corner.y+6*grid.uniform_grid.dX.y),y_max_minus_dy(grid.uniform_grid.domain.max_corner.y-6*grid.uniform_grid.dX.y),
        z_min_plus_dz(grid.uniform_grid.domain.min_corner.z+6*grid.uniform_grid.dX.z),z_max_minus_dz(grid.uniform_grid.domain.max_corner.z-6*grid.uniform_grid.dX.z);
    for(int i=1;i<=grid.number_of_nodes;i++)if(particle_levelset_evolution.phi(i)<0){
        VECTOR<T,3> node_location=grid.Node_Location(i);
        if(node_location.x<x_min_plus_dx)shift_left=true;if(node_location.x>x_max_minus_dx)shift_right=true;
        if(node_location.y<y_min_plus_dy)shift_bottom=true;if(node_location.y>y_max_minus_dy)shift_top=true;
        if(node_location.z<z_min_plus_dz)shift_front=true;if(node_location.z>z_max_minus_dz)shift_back=true;}
    if(shift_left && shift_right){shift_left=false;shift_right=false;}
    if(shift_bottom && shift_top){shift_bottom=false;shift_top=false;}
    if(shift_front && shift_back){shift_front=false;shift_back=false;}

    RANGE<TV> domain=grid.uniform_grid.domain;
    if(shift_left){
        domain.min_corner.x-=2*grid.uniform_grid.dX.x;domain.max_corner.x-=2*grid.uniform_grid.dX.x;grid.Move_Contents_Right_Two();
        grid.uniform_grid.Initialize(grid.uniform_grid.counts,domain);
        grid.Tree_Topology_Changed();}
    if(shift_right){
        domain.min_corner.x+=2*grid.uniform_grid.dX.x;domain.max_corner.x+=2*grid.uniform_grid.dX.x;grid.Move_Contents_Left_Two();
        grid.uniform_grid.Initialize(grid.uniform_grid.counts,domain);
        grid.Tree_Topology_Changed();}
    if(shift_bottom){
        domain.min_corner.y-=2*grid.uniform_grid.dX.y;domain.max_corner.y-=2*grid.uniform_grid.dX.y;grid.Move_Contents_Up_Two();
        grid.uniform_grid.Initialize(grid.uniform_grid.counts,domain);
        grid.Tree_Topology_Changed();}
    if(shift_top){
        domain.min_corner.y+=2*grid.uniform_grid.dX.y;domain.max_corner.y+=2*grid.uniform_grid.dX.y;grid.Move_Contents_Down_Two();
        grid.uniform_grid.Initialize(grid.uniform_grid.counts,domain);
        grid.Tree_Topology_Changed();}
    if(shift_front){
        domain.min_corner.z-=2*grid.uniform_grid.dX.z;domain.max_corner.z-=2*grid.uniform_grid.dX.z;grid.Move_Contents_Backward_Two();
        grid.uniform_grid.Initialize(grid.uniform_grid.counts,domain);
        grid.Tree_Topology_Changed();}
    if(shift_back){
        domain.min_corner.z+=2*grid.uniform_grid.dX.z;domain.max_corner.z+=2*grid.uniform_grid.dX.z;grid.Move_Contents_Forward_Two();
        grid.uniform_grid.Initialize(grid.uniform_grid.counts,domain);
        grid.Tree_Topology_Changed();}
}
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Move_Grid(const T time)
{
    assert(!smoke && !fire);
    
    if(move_grid_explicitly) callbacks->Move_Grid_Explicitly(time); // otherwise move automatically
    else Move_Grid_Helper(*this);
}
//#####################################################################
// Function Specify_Refinement_With_Thresholds_Helper
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Specify_Refinement_For_Cell_With_Thresholds_Helper(const T_CELL& cell,ARRAY<typename T_CELL::REFINE_ACTION>& refine_action,const ARRAY<VECTOR<T,2> >& refinement_thresholds,const ARRAY<T>& refinement_variable,const T dt,const T time)
{
    if(cell.Depth_Of_This_Cell()<maximum_tree_depth){
        int refinement_criteria_to_use=min(maximum_tree_depth-cell.Depth_Of_This_Cell(),refinement_thresholds.m);
        if(refinement_thresholds(refinement_criteria_to_use)[1]<refinement_variable(cell.Cell())&&refinement_variable(cell.Cell())<refinement_thresholds(refinement_criteria_to_use)[2]){
            refine_action(cell.Cell())=T_CELL::REFINE;return;}}
    if(refine_action(cell.Cell())<=T_CELL::COARSEN){
        int refinement_criteria_to_use=min(maximum_tree_depth-cell.Depth_Of_This_Cell()+1,refinement_thresholds.m);
        if((T).9*refinement_thresholds(refinement_criteria_to_use)[1]>refinement_variable(cell.Cell())||refinement_variable(cell.Cell())>(T)1.1*refinement_thresholds(refinement_criteria_to_use)[2]){
            refine_action(cell.Cell())=T_CELL::COARSEN;return;}}
    refine_action(cell.Cell())=max(refine_action(cell.Cell()),T_CELL::NONE); // this cell is at the right level. Don't allow coarsening
}
//#####################################################################
// Function Specify_Refinement_With_Levelset_Helper
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Specify_Refinement_For_Cell_With_Levelset_Helper(const T_CELL& cell,ARRAY<typename T_CELL::REFINE_ACTION>& refine_action,const ARRAY<T>& phi_array,const T band_width_times_small_dx,const T dt,const T time)
{
    T phi=phi_array(cell.Cell());
    T dx=cell.Diagonal_Length()+band_width_times_small_dx;
    if(!(phi<-dx||phi>dx)){refine_action(cell.Cell())=T_CELL::REFINE;return;}
    if(refine_action(cell.Cell())<=T_CELL::COARSEN){
        T dx=(T)2.01*cell.Diagonal_Length()+band_width_times_small_dx;
        if(phi<-dx||phi>dx){refine_action(cell.Cell())=T_CELL::COARSEN;return;}}
    refine_action(cell.Cell())=max(refine_action(cell.Cell()),T_CELL::NONE); // this cell is at the right level. Don't allow coarsening
}
//#####################################################################
// Function Specify_Refinement_With_Object_Helper
//#####################################################################
template<class T_GRID> template<class OBJECT> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Specify_Refinement_For_Cell_With_Object_Helper(const T_CELL& cell,ARRAY<typename T_CELL::REFINE_ACTION>& refine_action,const OBJECT& object,const TRANSFORMATION_MATRIX& world_to_source,const T band_width_times_small_dx,const T dt,const T time)
{
    T phi=object.Signed_Distance(world_to_source.Homogeneous_Times(cell.Center()));
    T dx=cell.Diagonal_Length()+band_width_times_small_dx;
    if(particle_levelset_evolution.phi(cell.Cell())<dx&&!(phi<-dx||phi>dx)){refine_action(cell.Cell())=T_CELL::REFINE;return;}
    if(refine_action(cell.Cell())<=T_CELL::COARSEN){
        T dx=(T)2.01*cell.Diagonal_Length()+band_width_times_small_dx;
        if(particle_levelset_evolution.phi(cell.Cell())>dx||phi<-dx||phi>dx){refine_action(cell.Cell())=T_CELL::COARSEN;return;}}
    refine_action(cell.Cell())=max(refine_action(cell.Cell()),T_CELL::NONE); // this cell is at the right level. Don't allow coarsening
}
//#####################################################################
// Function Specify_Refinement
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Specify_Refinement(ARRAY<typename T_CELL::REFINE_ACTION>& refine_action,const ARRAY<T_CELL*>* cells_to_specify_for,const T dt,const T time)
{
    ARRAY<T_CELL*>& cell_pointer_from_index=grid->Cell_Pointer_From_Index();
    ARRAY<T>& phi=particle_levelset_evolution.particle_levelset.levelset.phi;
    if(fire){
        T refinement_distance=(T).5*levelset_refinement_bandwidth*grid->Minimum_Edge_Length();
        if(cells_to_specify_for){for(int i=1;i<=cells_to_specify_for->m;i++){
            Specify_Refinement_For_Cell_With_Thresholds_Helper(*(*cells_to_specify_for)(i),refine_action,temperature_refinement_thresholds,temperature_container.temperature,dt,time);
            Specify_Refinement_For_Cell_With_Levelset_Helper(*(*cells_to_specify_for)(i),refine_action,phi,refinement_distance,dt,time);}}
        else{for(int i=1;i<=grid->number_of_cells;i++){
            Specify_Refinement_For_Cell_With_Thresholds_Helper(*cell_pointer_from_index(i),refine_action,temperature_refinement_thresholds,temperature_container.temperature,dt,time);
            Specify_Refinement_For_Cell_With_Levelset_Helper(*cell_pointer_from_index(i),refine_action,phi,refinement_distance,dt,time);}}}
    else if(water){
        T refinement_distance=(T).5*levelset_refinement_bandwidth*grid->Minimum_Edge_Length();
        if(cells_to_specify_for){for(int i=1;i<=cells_to_specify_for->m;i++)Specify_Refinement_For_Cell_With_Levelset_Helper(*(*cells_to_specify_for)(i),refine_action,phi,refinement_distance,dt,time);}
        else{for(int i=1;i<=grid->number_of_cells;i++)Specify_Refinement_For_Cell_With_Levelset_Helper(*cell_pointer_from_index(i),refine_action,phi,refinement_distance,dt,time);}}
    else if(smoke){
        if(cells_to_specify_for){for(int i=1;i<=cells_to_specify_for->m;i++)
            Specify_Refinement_For_Cell_With_Thresholds_Helper(*(*cells_to_specify_for)(i),refine_action,density_refinement_thresholds,density_container.density,dt,time);}
        else{for(int i=1;i<=grid->number_of_cells;i++)
            Specify_Refinement_For_Cell_With_Thresholds_Helper(*cell_pointer_from_index(i),refine_action,density_refinement_thresholds,density_container.density,dt,time);}}
    else PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Setup_Initial_Refinement
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Setup_Initial_Refinement(const T initial_time)
{
    ARRAY<T>& phi=particle_levelset_evolution.particle_levelset.levelset.phi;
    grid->Node_Iterator_Data();grid->Face_Iterator_Data();
    if(smoke || fire) density_container.Initialize_Array();
    if(fire) temperature_container.Initialize_Array();
    Initialize_Density_And_Temperature(initial_time);
    while(true){
        if(fire || water){
            T_LEVELSET levelset(*grid,phi);
            phi.Resize(grid->number_of_cells);
            callbacks->Initialize_Phi();
            /*levelset.Fast_Marching_Method();*/} // TODO: this will fail because the nodal fast march does not have a band to use yet....
        bool refined_any_cells=false;int number_of_cells=grid->number_of_cells;ARRAY<T_CELL*>& cell_pointer_from_index=grid->Cell_Pointer_From_Index();
        ARRAY<typename T_CELL::REFINE_ACTION> refine_action(grid->number_of_cells);Specify_Refinement(refine_action,0,0,initial_time);
        for(int i=1;i<=number_of_cells;i++)if(refine_action(i)==T_CELL::REFINE){
            T_CELL* cell=cell_pointer_from_index(i);
            if(grid->Domain().Lazy_Inside(cell->Center())&&cell->Depth_Of_This_Cell()<grid->maximum_depth&&!cell->Has_Children()){
                refined_any_cells=true;cell->Create_Children(grid->number_of_cells,0,grid->number_of_nodes,0,grid->number_of_faces,0,grid);}}
        if(!refined_any_cells)break;
        grid->Tree_Topology_Changed();
        grid->Node_Iterator_Data();grid->Face_Iterator_Data();
        callbacks->Topology_Changed();
        if(smoke || fire) density_container.Initialize_Array();
        if(fire) temperature_container.Initialize_Array();
        Initialize_Density_And_Temperature(initial_time);
    }
}
//#####################################################################
// Function Total_Number_Of_Particles
//#####################################################################
template<class T_GRID> template<class T_PARTICLES> int FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Total_Number_Of_Particles(const ARRAY<T_PARTICLES*>& particles) const
{
    int total=0;
    for(int i=1;i<=particles.m;i++) if(particles(i)) total+=particles(i)->array_collection->Size();
    return total;
}
//#####################################################################
// Function Write_Particles
//#####################################################################
template<class T_GRID> template<class T_PARTICLES> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Write_Particles(const STREAM_TYPE stream_type,const ARRAY<T_PARTICLES*>& particles,const std::string& output_directory,const std::string& prefix,const int frame) const
{
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,prefix.c_str()),particles);
    LOG::cout << "Writing " << Total_Number_Of_Particles(particles) << " " << prefix << std::endl;
}
//#####################################################################
// Function Read_Particles
//#####################################################################
template<class T_GRID> template<class T_PARTICLES> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Read_Particles(const STREAM_TYPE stream_type,ARRAY<T_PARTICLES*>& particles,const std::string& output_directory,const std::string& prefix,const int frame)
{
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,prefix.c_str()),particles);
    LOG::cout << "Reading " << Total_Number_Of_Particles(particles) << " " << prefix << std::endl;
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Read_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame),prefix=output_directory+"/"+grid->name;
    if(smoke || fire || water){
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"_grid."+f,*grid);
        // scalar fields
        if(smoke||fire){
            if(use_density) FILE_UTILITIES::Read_From_File(stream_type,prefix+"_density."+f,density_container.density);
            if(use_temperature) FILE_UTILITIES::Read_From_File(stream_type,prefix+"_temperature."+f,temperature_container.temperature);}
        // velocities
        if(write_velocity)if(fire)FILE_UTILITIES::Read_From_File(stream_type,prefix+"_face_velocities."+f,incompressible.projection.face_velocities);
        // particle levelset 
        if(write_levelset && (water||fire))
            FILE_UTILITIES::Read_From_File(stream_type,prefix+"_levelset."+f,particle_levelset_evolution.particle_levelset.levelset.phi);
        if(water){
            if(write_particles){
                Read_Particles(stream_type,particle_levelset_evolution.particle_levelset.positive_particles,output_directory,"positive_particles",frame);
                Read_Particles(stream_type,particle_levelset_evolution.particle_levelset.negative_particles,output_directory,"negative_particles",frame);}
            if(write_removed_positive_particles)
                Read_Particles(stream_type,particle_levelset_evolution.particle_levelset.removed_positive_particles,output_directory,"removed_positive_particles",frame);
            if(write_removed_negative_particles)
                Read_Particles(stream_type,particle_levelset_evolution.particle_levelset.removed_negative_particles,output_directory,"removed_negative_particles",frame);
            if(store_particle_ids) FILE_UTILITIES::Read_From_Text_File(output_directory+"/last_unique_particle_id."+f,particle_levelset_evolution.particle_levelset.last_unique_particle_id);}
        // for more accurate restarts
        std::string filename;
        filename=prefix+"_pressure."+f;
        if(FILE_UTILITIES::File_Exists(filename)){LOG::cout << "Reading dyadic_pressure " << filename << std::endl;
            FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible.projection.p);}

        // read next bodies states for fluid collision bodies
        if(solid_affects_fluid && fluid_affects_solid && collision_bodies_affecting_fluid->collision_geometry_collection.bodies.m){
            try{
                std::istream* input_raw(FILE_UTILITIES::Safe_Open_Input(output_directory+"/fluid_collision_bodies_new_state."+f));TYPED_ISTREAM input(*input_raw,stream_type);
                collision_bodies_affecting_fluid->Read_State(input,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);delete input_raw;}
            catch(FILESYSTEM_ERROR&){LOG::cerr<<"PHYSBAM_WARNING: Could not read fluid collision bodies new state, restart may be incorrect!"<<std::endl;}}}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T_GRID> void FLUIDS_PARAMETERS_DYADIC<T_GRID>::
Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int first_frame,const int frame) const
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame),prefix=output_directory+"/"+grid->name;
    if(smoke || fire || water){
        FILE_UTILITIES::Write_To_File(stream_type,prefix+"_grid."+f,*grid);
        // scalar fields
        if(smoke||fire){
            if(use_density) FILE_UTILITIES::Write_To_File(stream_type,prefix+"_density."+f,density_container.density);
            if(use_temperature) FILE_UTILITIES::Write_To_File(stream_type,prefix+"_temperature."+f,temperature_container.temperature);}
        // velocities
        if(write_velocity && frame%restart_data_write_rate==0)if(fire)FILE_UTILITIES::Write_To_File(stream_type,prefix+"_face_velocities."+f,incompressible.projection.face_velocities);
        // particle levelset 
        if(write_levelset && (water||fire)){
            assert(!write_ghost_values);
            FILE_UTILITIES::Write_To_File(stream_type,prefix+"_levelset."+f,particle_levelset_evolution.particle_levelset.levelset.phi);}
        if(water){
            if(write_particles && frame%restart_data_write_rate==0){
                Write_Particles(stream_type,particle_levelset_evolution.particle_levelset.positive_particles,output_directory,"positive_particles",frame);
                Write_Particles(stream_type,particle_levelset_evolution.particle_levelset.negative_particles,output_directory,"negative_particles",frame);}
            if(write_removed_positive_particles)
                Write_Particles(stream_type,particle_levelset_evolution.particle_levelset.removed_positive_particles,output_directory,"removed_positive_particles",frame);
            if(write_removed_negative_particles)
                Write_Particles(stream_type,particle_levelset_evolution.particle_levelset.removed_negative_particles,output_directory,"removed_negative_particles",frame);
            if(store_particle_ids) FILE_UTILITIES::Write_To_Text_File(output_directory+"/last_unique_particle_id."+f,particle_levelset_evolution.particle_levelset.last_unique_particle_id);}
        // write next bodies states for fluid collision bodies
        if(solid_affects_fluid && fluid_affects_solid && collision_bodies_affecting_fluid->collision_geometry_collection.bodies.m && frame%restart_data_write_rate==0){
            std::ostream* output_raw(FILE_UTILITIES::Safe_Open_Output(output_directory+"/fluid_collision_bodies_new_state."+f));TYPED_OSTREAM output(*output_raw,stream_type);
            collision_bodies_affecting_fluid->Write_State(output,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);delete output_raw;}
        if(write_debug_data || (write_restart_data && frame%restart_data_write_rate==0)){
            FILE_UTILITIES::Write_To_File(stream_type,prefix+"_pressure."+f,incompressible.projection.p);
            FILE_UTILITIES::Write_To_File(stream_type,prefix+"_face_velocities."+f,incompressible.projection.face_velocities);}
        // debugging
        if(write_debug_data){
            FILE_UTILITIES::Write_To_File(stream_type,prefix+"_psi_N."+f,incompressible.projection.elliptic_solver->psi_N);
            FILE_UTILITIES::Write_To_File(stream_type,prefix+"_psi_D."+f,incompressible.projection.elliptic_solver->psi_D);
            FILE_UTILITIES::Write_To_File(stream_type,prefix+"_colors."+f,incompressible.projection.elliptic_solver->filled_region_colors);}}
}
//#####################################################################
#define COMMA ,
#define INSTANTIATION_HELPER_SHAPE(SHAPE) \
    template void FLUIDS_PARAMETERS_DYADIC<DYADIC_GRID_POLICY<SHAPE::VECTOR_T>::DYADIC_GRID>::Specify_Refinement_For_Cell_With_Object_Helper(const T_CELL&,ARRAY<T_CELL::REFINE_ACTION>&, \
        const SHAPE&,const TRANSFORMATION_MATRIX&,const SHAPE::VECTOR_T::SCALAR,const SHAPE::VECTOR_T::SCALAR,const SHAPE::VECTOR_T::SCALAR);
#define INSTANTIATION_HELPER(T) \
    template class FLUIDS_PARAMETERS_DYADIC<OCTREE_GRID<T> >; \
    template class FLUIDS_PARAMETERS_DYADIC<QUADTREE_GRID<T> >; \
    INSTANTIATION_HELPER_SHAPE(BOX<VECTOR<T COMMA 2> >) \
    INSTANTIATION_HELPER_SHAPE(BOX<VECTOR<T COMMA 3> >) \
    INSTANTIATION_HELPER_SHAPE(CYLINDER<T>) \
    INSTANTIATION_HELPER_SHAPE(IMPLICIT_OBJECT<VECTOR<T COMMA 2> >) \
    INSTANTIATION_HELPER_SHAPE(IMPLICIT_OBJECT<VECTOR<T COMMA 3> >)
INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
#endif
