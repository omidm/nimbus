#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_TRANSFER_ITERATOR.h>
#include <PhysBAM_Tools/Grids_RLE_Advection/ADVECTION_SEMI_LAGRANGIAN_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Boundaries/BOUNDARY_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_RLE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_FRAME.h>
#include <PhysBAM_Tools/Read_Write/Random_Numbers/READ_WRITE_RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Utilities/Find_Type.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Grids_RLE_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_CELL_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Advection_Collidable/ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_FACE_RLE.h>
#include <PhysBAM_Geometry/Read_Write/Grids_RLE_Level_Sets/READ_WRITE_LEVELSET_RLE.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/FLUID_COLLISION_BODY_INACCURATE_UNION.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_RLE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
SOLIDS_FLUIDS_EXAMPLE_RLE(const STREAM_TYPE stream_type)
    :SOLIDS_FLUIDS_EXAMPLE<TV>(stream_type),FLUIDS_PARAMETERS<T_GRID>(FLUIDS_PARAMETERS<T_GRID>::WATER),fluids_parameters(*this),incompressible(*fluids_parameters.grid),
    particle_levelset(*fluids_parameters.grid),levelset_advection(&particle_levelset.levelset),mpi_grid(0),
    refine_all_water(false),use_phi_for_vertical_refinement(true),vertical_refinement_depth(0),refine_near_some_objects(true),
    refine_near_negative_particles(true),refine_near_removed_negative_particles(false),enforce_refinement_slope(false),use_incompressible_cfl(true),
    clamp_long_velocities_to_short_velocities(false)
{
    fluids_parameters.phi_boundary=&fluids_parameters.phi_boundary_water; // override default, TODO: make RLE aware
    // TODO: fluids_parameters.phi_boundary_water.Set_Velocity_Pointer(incompressible.V);
    fluids_parameters.implicit_viscosity=false;
    fluids_parameters.second_order_cut_cell_method=true;
    fluids_parameters.enforce_divergence_free_extrapolation=false;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
~SOLIDS_FLUIDS_EXAMPLE_RLE()
{
    delete mpi_grid;
}
//#####################################################################
// Function Use_Fluid_Coupling_Defaults
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Use_Fluid_Coupling_Defaults()
{
    // TODO(jontg): This function belongs in FLUIDS_PARAMETERS, like th DYADIC and UNIFORM cases...
    fluids_parameters.solid_affects_fluid=true;fluids_parameters.fluid_affects_solid=true;
    levelset_advection.Use_Semi_Lagrangian_Collidable_Advection(*fluids_parameters.collision_bodies_affecting_fluid,fluids_parameters.collidable_phi_replacement_value,
        incompressible.valid_mask);
    incompressible.Use_Semi_Lagrangian_Collidable_Advection(*fluids_parameters.collision_bodies_affecting_fluid);
}
//#####################################################################
// Function Use_Fluid_Slip_Coupling_Defaults
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Use_Fluid_Slip_Coupling_Defaults()
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Use_No_Fluid_Coupling_Defaults
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Use_No_Fluid_Coupling_Defaults()
{
    levelset_advection.Set_Custom_Advection(fluids_parameters.semi_lagrangian);
    incompressible.Set_Custom_Advection(fluids_parameters.semi_lagrangian);
}
//#####################################################################
// Function Initialize_Ground
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Initialize_Ground()
{   
    T_GRID& grid=*fluids_parameters.grid;
    ground_j.Resize(grid.horizontal_grid.Get_MAC_Grid().Domain_Indices(grid.number_of_ghost_cells+1));
    for(HORIZONTAL_CELL_ITERATOR iterator(grid.horizontal_grid,grid.number_of_ghost_cells+1);iterator.Valid();iterator.Next()){TV_HORIZONTAL_INT cell=iterator.Cell_Index();
        T ground=Initial_Ground(iterator.Location());
        ground_j(cell)=(int)floor((ground-grid.uniform_grid.domain.min_corner.y)*grid.uniform_grid.one_over_dX.y+(T)1.5);}
}
//#####################################################################
// Function Initialize_Grid
//#####################################################################
template<class T_GRID> class HELPER_IMPLICIT_OBJECT:public RLE_INITIAL_PHI_HELPER<typename T_GRID::SCALAR, T_GRID::VECTOR_T::dimension>
{
public:
    SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>* example;
    ~HELPER_IMPLICIT_OBJECT() {};
    typename T_GRID::SCALAR operator()(const typename T_GRID::VECTOR_T& X)const{return max(example->Initial_Phi(X),-example->Initial_Phi_Object(X));}
};
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Initialize_Grid()
{
    T_GRID& grid=*fluids_parameters.grid;
    HELPER_IMPLICIT_OBJECT<T_GRID> helper;helper.example=this;
    grid.Initialize(helper,&ground_j);
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Initialize_Phi()
{
    T_GRID& grid=*fluids_parameters.grid;
    ARRAY<T>& phi=particle_levelset.phi;
    for(CELL_ITERATOR cell(grid,1);cell;cell++){int c=cell.Cell();
        phi(c)=Initial_Phi(cell.Center());if(cell.Long()) phi(c+1)=phi(c);}
}
//#####################################################################
// Function Initialize_MPI
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Initialize_MPI()
{
    T_GRID& grid=*fluids_parameters.grid;
    VECTOR<VECTOR<bool,2>,T_GRID::dimension-1> horizontal_domain_walls=fluids_parameters.domain_walls.Remove_Index(2);
    mpi_grid->Initialize(horizontal_domain_walls);
    fluids_parameters.domain_walls(1)=horizontal_domain_walls(1);
    if(T_GRID::dimension==3) fluids_parameters.domain_walls(3)=horizontal_domain_walls(2);
    grid.mpi_grid=mpi_grid;
    fluids_parameters.phi_boundary=new BOUNDARY_MPI<T_GRID>(mpi_grid,*fluids_parameters.phi_boundary,grid.number_of_ghost_cells);
    fluids_parameters.fluid_boundary=new BOUNDARY_MPI<T_GRID>(mpi_grid,*fluids_parameters.fluid_boundary,grid.number_of_ghost_cells); 
    particle_levelset.mpi_grid=mpi_grid;
    incompressible.mpi_grid=mpi_grid;
    incompressible.projection.laplace.mpi_grid=mpi_grid;
    if(!restart && fluids_parameters.store_particle_ids) particle_levelset.last_unique_particle_id=mpi_grid->rank*30000000;
}
//#####################################################################
// Function Update_Fluid_Parameters
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Update_Fluid_Parameters(const T dt,const T time)
{
    if(fluids_parameters.fire || fluids_parameters.viscosity || fluids_parameters.variable_viscosity || fluids_parameters.use_vorticity_confinement) PHYSBAM_FATAL_ERROR();
    incompressible.Set_Gravity(fluids_parameters.gravity,fluids_parameters.gravity_direction);
    incompressible.Set_Body_Force(fluids_parameters.use_body_force);
    if(fluids_parameters.use_body_force) Get_Body_Force(incompressible.force,dt,time);
    incompressible.projection.Set_Density(fluids_parameters.density);
    incompressible.projection.Use_Non_Zero_Divergence(fluids_parameters.use_non_zero_divergence);
    if(fluids_parameters.use_non_zero_divergence) Get_Divergence(incompressible.projection.divergence,dt,time);
    incompressible.projection.laplace.Solve_Neumann_Regions(fluids_parameters.solve_neumann_regions);
}
//#####################################################################
// Function Adjust_Particle_For_Domain_Boundaries
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time)
{
    // remove this for speed - don't call the function for the other particles
    if(particle_type == PARTICLE_LEVELSET_POSITIVE || particle_type == PARTICLE_LEVELSET_REMOVED_POSITIVE) return;
    
    TV& X=particles.X(index);TV X_new=X+dt*V;
    T max_collision_distance=particle_levelset.Particle_Collision_Distance(particles.quantized_collision_distance(index));
    T_GRID::template Horizontal_Face_Loop<Adjust_Particle_For_Domain_Boundaries_Helper>(fluids_parameters,X,X_new,V,max_collision_distance,dt);
}
template<class T_GRID> template<class T_FACE> inline void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Adjust_Particle_For_Domain_Boundaries_Helper::Apply(const FLUIDS_PARAMETERS<T_GRID>& fluids_parameters,TV& X,TV& X_new,TV& V,const T max_collision_distance,const T dt)
{
    const RANGE<TV>& domain=fluids_parameters.grid->domain;
    TV domain_min=domain.Minimum_Corner(),domain_max=domain.Maximum_Corner();
    VECTOR<VECTOR<bool,2>,T_GRID::dimension-1> horizontal_walls=fluids_parameters.domain_walls.Remove_Index(2);
    int axis=T_FACE::Axis(),wall_axis_index=T_FACE::Horizontal_Axis();
    if(horizontal_walls[wall_axis_index][1] && X_new[axis]<domain_min[axis]+max_collision_distance){
        T collision_distance=X[axis]-domain_min[axis];if(collision_distance > max_collision_distance) collision_distance=X_new[axis]-domain_min[axis];
        collision_distance=max((T)0,collision_distance);
        X_new[axis]+=max((T)0,domain_min[axis]-X_new[axis]+collision_distance);
        V[axis]=max((T)0,V[axis]);X=X_new-dt*V;}
    if(horizontal_walls[wall_axis_index][2] && X_new[axis]>domain_max[axis]-max_collision_distance){
        T collision_distance=domain_max[axis]-X[axis];if(collision_distance > max_collision_distance) collision_distance=domain_max[axis]-X_new[axis];
        collision_distance=max((T)0,collision_distance);
        X_new[axis]-=max((T)0,X_new[axis]-domain_max[axis]+collision_distance);
        V[axis]=min((T)0,V[axis]);X=X_new-dt*V;}
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Delete_Particles_Inside_Objects(const T time)
{
    Delete_Particles_Inside_Objects(particle_levelset.positive_particles,PARTICLE_LEVELSET_POSITIVE,time);
    Delete_Particles_Inside_Objects(particle_levelset.negative_particles,PARTICLE_LEVELSET_NEGATIVE,time);
    if(particle_levelset.use_removed_positive_particles)
        Delete_Particles_Inside_Objects(particle_levelset.removed_positive_particles,PARTICLE_LEVELSET_REMOVED_POSITIVE,time);
    if(particle_levelset.use_removed_negative_particles)
        Delete_Particles_Inside_Objects(particle_levelset.removed_negative_particles,PARTICLE_LEVELSET_REMOVED_NEGATIVE,time);
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
template<class T_GRID> template<class T_PARTICLES> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Delete_Particles_Inside_Objects(ARRAY<T_PARTICLES*>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time)
{
    if(!fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.bodies.m) return;
    T_GRID& grid=*fluids_parameters.grid;
    for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();
        if(particles(b) && fluids_parameters.collision_bodies_affecting_fluid->Occupied_Block(block)) Delete_Particles_Inside_Objects(block,*particles(b),particle_type,time);}
}
//#####################################################################
// Function Delete_Particles_Inside_Objects
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Delete_Particles_Inside_Objects(const BLOCK_ITERATOR& block,PARTICLE_LEVELSET_PARTICLES<TV>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T time)
{
    COLLISION_GEOMETRY_ID body_id;int aggregate_id;
    for(int k=particles.array_collection->Size();k>=1;k--)
        if(fluids_parameters.collision_bodies_affecting_fluid->Inside_Any_Simplex_Of_Any_Body(block,particles.X(k),body_id,aggregate_id))
            particles.array_collection->Delete_Element(k);
}
//#####################################################################
// Function Get_Neumann_And_Dirichlet_Boundary_Conditions
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Get_Neumann_And_Dirichlet_Boundary_Conditions(const T dt,const T time)
{
    T_GRID& grid=*fluids_parameters.grid;
    LAPLACE_COLLIDABLE_RLE<T_GRID>& laplace=incompressible.projection.laplace;
    laplace.psi_N.Fill(false);laplace.psi_D.Fill(false);
    Set_Domain_Boundary_Conditions(time);
    Get_Source_Velocities(time);
    Get_Object_Velocities(dt,time);
    Get_Ground_Velocities(dt,time);
    if(clamp_long_velocities_to_short_velocities) T_GRID::template Face_Loop<Clamp_Long_Velocities_To_Short_Velocities>(grid,incompressible.V);
    Set_Dirichlet_Boundary_Conditions(time);
}
//#####################################################################
// Function Set_Domain_Boundary_Conditions
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Set_Domain_Boundary_Conditions(const T time)
{   
    T pressure=0;
    T_GRID& grid=*fluids_parameters.grid;
    ARRAY<bool>& psi_D=incompressible.projection.laplace.psi_D;
    ARRAY<T_BOX_HORIZONTAL_INT> regions;grid.Find_Ghost_Regions(regions,CELL_ITERATOR::Sentinels(),1,false);
    for(int r=1;r<=regions.m;r++)for(CELL_ITERATOR cell(grid,regions(r));cell;cell++){int c=cell.Cell();
        psi_D(c)=true;incompressible.projection.p(c)=pressure;if(cell.Long()){psi_D(c+1)=true;incompressible.projection.p(c+1)=pressure;}}
    T_GRID::template Horizontal_Face_Loop<Set_Domain_Boundary_Conditions_Helper>(*this);
}
template<class T_GRID> template<class T_FACE> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Set_Domain_Boundary_Conditions_Helper::Apply(SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>& example)
{
    FLUIDS_PARAMETERS<T_GRID>& fluids_parameters=example.fluids_parameters;
    T_GRID& grid=*fluids_parameters.grid;
    ARRAY<T>& V=example.incompressible.V;
    ARRAY<bool>& psi_N=example.incompressible.projection.laplace.psi_N;
    VECTOR<VECTOR<bool,2>,T_GRID::dimension-1> horizontal_walls=fluids_parameters.domain_walls.Remove_Index(2);
    ARRAY<T_BOX_HORIZONTAL_INT> cell_regions;grid.Find_Ghost_Regions(cell_regions,CELL_ITERATOR::Sentinels(),1,false);
    ARRAY<T_BOX_HORIZONTAL_INT> face_regions;grid.Find_Ghost_Regions(face_regions,T_FACE::Sentinels()+T_FACE::Boundary_Sentinels(false),1,false);
    int horizontal_axis_index=T_FACE::Horizontal_Axis();
    int region_index=2*T_FACE::Horizontal_Axis()-1;
    for(int side=0;side<2;side++)if(horizontal_walls[horizontal_axis_index][side+1]){
        for(T_FACE face(grid,face_regions(region_index+side));face;face++)if(example.particle_levelset.phi(face.Cell(1-side))<=0){int f=face.Face();
            psi_N(f)=true;V(f)=0;if(face.Long()){psi_N(f+1)=true;V(f+1)=0;}}
        for(CELL_ITERATOR cell(grid,cell_regions(region_index+side));cell;cell++)if(cell.Long()) psi_N(cell.Face_Y()+1)=true;}
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Set_Dirichlet_Boundary_Conditions(const T time)
{
    incompressible.Set_Dirichlet_Boundary_Conditions(particle_levelset.phi);
}
//#####################################################################
// Function Revalidate_Fluid_Scalars
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Revalidate_Fluid_Scalars()
{
    if(!levelset_advection.nested_semi_lagrangian_collidable) return;
    T_GRID& grid=*fluids_parameters.grid;
    ARRAY<T>& phi=particle_levelset.phi;
    levelset_advection.nested_semi_lagrangian_collidable->Average_To_Invalidated_Cells(grid,fluids_parameters.collidable_phi_replacement_value,phi);
    // make sure revalidation didn't screw up values outside the band
    T positive_threshold=(T).999*grid.positive_bandwidth,negative_threshold=(T)-.999*grid.negative_bandwidth;
    for(int c=1;c<=grid.number_of_cells;c++){
        if(phi(c)<=negative_threshold) phi(c)=-grid.negative_bandwidth;
        else if(phi(c)>=positive_threshold) phi(c)=grid.positive_bandwidth;}
}
//#####################################################################
// Function Revalidate_Phi_After_Modify_Levelset
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Revalidate_Phi_After_Modify_Levelset()
{
    if(!levelset_advection.nested_semi_lagrangian_collidable) return;
    levelset_advection.nested_semi_lagrangian_collidable->cell_valid_points_current=levelset_advection.nested_semi_lagrangian_collidable->cell_valid_points_next;
    Revalidate_Fluid_Scalars();
}
//#####################################################################
// Function Revalidate_Fluid_Velocity
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Revalidate_Fluid_Velocity()
{
    if(incompressible.nested_semi_lagrangian_collidable)
        incompressible.nested_semi_lagrangian_collidable->Average_To_Invalidated_Face(*fluids_parameters.grid,incompressible.V);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Get_Object_Velocities(const T dt,const T time)
{
    fluids_parameters.collision_bodies_affecting_fluid->Compute_Psi_N(incompressible.projection.laplace.psi_N,&incompressible.V);
}
//#####################################################################
// Function Get_Ground_Velocities
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Get_Ground_Velocities(const T dt,const T time)
{
    T_GRID::template Horizontal_Face_Loop<Get_Ground_Velocities_Horizontal>(*this);
    Get_Ground_Velocities_Vertical();
}
//#####################################################################
// Function Get_Ground_Velocities_Horizontal
//#####################################################################
template<class T_GRID> template<class T_FACE> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Get_Ground_Velocities_Horizontal::Apply(SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>& example)
{
    T_GRID& grid=*example.fluids_parameters.grid;
    ARRAY<bool>& psi_N=example.incompressible.projection.laplace.psi_N;ARRAY<T>& V=example.incompressible.V;T_ARRAYS_HORIZONTAL_INT& ground_j=example.ground_j;
    for(T_FACE face(grid,0);face;face++)
        if(face.cell1.j<ground_j(face.cell1.Horizontal_Index()) || face.cell2.j<ground_j(face.cell2.Horizontal_Index())){int f=face.Face();
            psi_N(f)=true;V(f)=0;
            if(face.Long()){psi_N(f+1)=true;V(f+1)=0;}}
}
//#####################################################################
// Function Get_Ground_Velocities_Vertical
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Get_Ground_Velocities_Vertical()
{
    typedef typename T_GRID::ARRAYS_HORIZONTAL T_ARRAYS_HORIZONTAL_T;
    T_GRID& grid=*fluids_parameters.grid;
    ARRAY<bool>& psi_N=incompressible.projection.laplace.psi_N;ARRAY<T>& V=incompressible.V;

    // get ground velocity field
    T_ARRAYS_HORIZONTAL_T v(grid.horizontal_grid.Get_MAC_Grid().Domain_Indices(),false);
    if(use_deep_water) deep_water.Get_Vertical_Velocity(v);
    else v.Fill(0);

    for(FACE_Y_ITERATOR face(grid,0,true);face;face++){TV_HORIZONTAL_INT column=face.cell2.Horizontal_Index();
        int ground=ground_j(column);
        if(face.cell2.j<=ground){int f=face.Face();
            psi_N(f)=true;V(f)=v(column);
            if(face.Upper_Cell_Long() && face.cell2.jmax()<=ground){psi_N(f+1)=true;V(f+1)=v(column);}}}
}
//#####################################################################
// Function Clamp_Long_Velocities_To_Short_Velocities
//#####################################################################
template<class T_GRID> template<class T_FACE> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Clamp_Long_Velocities_To_Short_Velocities::Apply(const T_GRID& grid,ARRAY<T>& V)
{
    T max_short_V=0;
    for(T_FACE face(grid,0,false);face;face++){int f=face.Face();if(face.Both_Cells_Short()) max_short_V=max(max_short_V,(T)abs(V(f)));}
    LOG::cout<<"maximum short velocity = "<<max_short_V<<std::endl;
    for(T_FACE face(grid,0,false);face;face++){int f=face.Face();
        if(abs(V(f))>max_short_V) V(f)=sign(V(f))*max_short_V;
        if(face.Long() && abs(V(f+1))>max_short_V) V(f+1)=sign(V(f+1))*max_short_V;}
}
//#####################################################################
// Function Initialize_Swept_Occupied_Cell_For_Advection
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Initialize_Swept_Occupied_Blocks_For_Advection(const T dt,const T time,const bool include_removed_particle_velocities)
{
    T_GRID& grid=*fluids_parameters.grid;
    T maximum_fluid_speed=incompressible.V.Maxabs(),maximum_particle_speed=0;
    if(include_removed_particle_velocities){
        for(int i=1;i<=particle_levelset.removed_negative_particles.m;i++){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles=particle_levelset.removed_negative_particles(i);
            if(particles) maximum_particle_speed=max(maximum_particle_speed,particles->V.Maximum_Magnitude());}
        for(int i=1;i<=particle_levelset.removed_positive_particles.m;i++){
            PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles=particle_levelset.removed_positive_particles(i);
            if(particles) maximum_particle_speed=max(maximum_particle_speed,particles->V.Maximum_Magnitude());}}
    LOG::cout<<"maximum_fluid_speed = "<<maximum_fluid_speed<<", maximum_particle_speed = "<<maximum_particle_speed<<std::endl;
    T max_particle_collision_distance=particle_levelset.max_collision_distance_factor*grid.Minimum_Edge_Length();
    fluids_parameters.collision_bodies_affecting_fluid->Compute_Occupied_Blocks(true,fluids_parameters.cfl*dt*max(maximum_fluid_speed,maximum_particle_speed)+
        2*max_particle_collision_distance+(T).5*grid.Minimum_Edge_Length(),10);
//    fluids_parameters.collision_bodies_affecting_fluid->Compute_Occupied_Cells(fluids_parameters.grid,occupied_cell_for_advection,true,
//    2*dt*max(maximum_fluid_speed,maximum_particle_speed)+2*max_collision_distance+fluids_parameters.grid.Minimum_Edge_Length(),10);
}
//#####################################################################
// Function Get_Cell_Should_Be_Long
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Get_Cell_Should_Be_Long(ARRAY<bool>& cell_should_be_long,const T time) const
{
    const T_GRID& grid=*fluids_parameters.grid;
    const ARRAY<T>& phi=particle_levelset.phi;
    if(refine_all_water)
        for(CELL_ITERATOR cell(grid,grid.number_of_ghost_cells);cell;cell++){int c=cell.Cell();
            cell_should_be_long(c)=phi(c)>=grid.positive_bandwidth || cell.j<ground_j(cell.Horizontal_Index())-1;}
    else for(int c=1;c<=grid.number_of_cells;c++)cell_should_be_long(c)=grid.positive_bandwidth<=phi(c) || phi(c)<=-grid.negative_bandwidth;
    T vertical_threshold=use_phi_for_vertical_refinement?-grid.negative_bandwidth:0;
    if(vertical_refinement_depth>-vertical_threshold){
        const ARRAY<VECTOR<bool,T_GRID::dimension> >& cell_neighbors_visible=fluids_parameters.collision_bodies_affecting_fluid->cell_neighbors_visible;
        int extra_band=(int)ceil((vertical_refinement_depth+vertical_threshold)/grid.uniform_grid.dX.Max());
        for(FACE_Y_ITERATOR face(grid,0,false);face;face++)if(face.Both_Cells_Short()){int c2=face.cell2.Cell(),c1=c2-1;
            if(phi(c1)<=vertical_threshold && phi(c2)>vertical_threshold){
                int first_cell=c1-min(extra_band-1,face.cell1.dj);
                for(int c=c1;c>=first_cell;c--){
                    if(!cell_neighbors_visible(c)(2) || phi(c)>vertical_threshold) break;
                    cell_should_be_long(c)=false;}}
            else if(phi(c1)>vertical_threshold && phi(c2)<=vertical_threshold){
                int last_cell=c2+min(extra_band-1,(face.cell2.run+1)->jmin-face.cell2.j-1);
                for(int c=c2;c<=last_cell;c++){
                    if(!cell_neighbors_visible(c-1)(2) || phi(c)>vertical_threshold) break;
                    cell_should_be_long(c)=false;}}}}
    if(refine_near_negative_particles){
        bool consider_removed_particles=refine_near_removed_negative_particles && particle_levelset.use_removed_negative_particles;
        for(BLOCK_ITERATOR block(grid,1);block;block++){int b=block.Block();
            if(particle_levelset.negative_particles(b) || (consider_removed_particles && particle_levelset.removed_negative_particles(b)))
                for(int i=0;i<T_GRID::number_of_cells_per_block;i++){
                    T phi_cell=phi(block.Cell(i));if(0<phi_cell && phi_cell<grid.positive_bandwidth) cell_should_be_long(block.Cell(i))=false;}}}
    if(enforce_refinement_slope){
        for(CELL_ITERATOR cell(grid,0);cell;cell++)if(min(cell.dj,cell.j_end-cell.j)>=1){int c=cell.Cell();
            if(particle_levelset.phi(c)<0 && cell.j>ground_j(cell.Horizontal_Index())+2) cell_should_be_long(c)=false;}}
    if(refine_near_some_objects) T_GRID::template Face_Loop<Get_Cell_Should_Be_Long_For_Objects>(*this,cell_should_be_long);
}
//#####################################################################
// Function Get_Cell_Should_Be_Long_For_Objects
//#####################################################################
template<class T_GRID> template<class T_FACE> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Get_Cell_Should_Be_Long_For_Objects::Apply(const SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>& example,ARRAY<bool>& cell_should_be_long)
{
    const T_GRID& grid=*example.fluids_parameters.grid;
    const GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>& body_list=*example.fluids_parameters.collision_bodies_affecting_fluid;
    if(!body_list.collision_geometry_collection.bodies.m) return;
    const ARRAY<T>& phi=example.particle_levelset.phi;
    for(T_FACE face(grid,1,false);face;face++)if(face.Both_Cells_Short()){
        int c1=face.cell1.Cell(),c2=face.cell2.Cell();
        if((!cell_should_be_long(c1) && !cell_should_be_long(c2)) || (phi(c1)>0 && phi(c2)>0)) continue;
        if(body_list.Refine_Due_To_Objects(c1,c2)) cell_should_be_long(c1)=cell_should_be_long(c2)=false;}
}
//#####################################################################
// Function Add_Volumetric_Body_To_Fluid_Simulation
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Add_Volumetric_Body_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,const bool affects_fluid)
{
    if(affects_fluid){
        RIGID_COLLISION_GEOMETRY<TV>* collision_geometry=new RIGID_COLLISION_GEOMETRY<TV>(rigid_body);
        fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(collision_geometry,rigid_body.particle_index,true);}
}
//#####################################################################
// Function Add_Thin_Shell_To_Fluid_Simulation
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Add_Thin_Shell_To_Fluid_Simulation(RIGID_BODY<TV>& rigid_body,const bool affects_fluid)
{
    rigid_body.thin_shell=true;
    if(affects_fluid){
        RIGID_COLLISION_GEOMETRY<TV>* collision_geometry=new RIGID_COLLISION_GEOMETRY<TV>(rigid_body);
        if(collision_geometry->collision_thickness<minimum_collision_thickness) collision_geometry->Set_Collision_Thickness(minimum_collision_thickness);
        fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(collision_geometry,rigid_body.particle_index,true);}
}
//#####################################################################
// Function Add_To_Fluid_Simulation
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Add_To_Fluid_Simulation(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>& deformable_collisions,const bool affects_fluid)
{
    if(deformable_collisions.collision_thickness<minimum_collision_thickness) deformable_collisions.Set_Collision_Thickness(minimum_collision_thickness);
    if(affects_fluid){
        fluids_parameters.collision_bodies_affecting_fluid->collision_geometry_collection.Add_Body(&deformable_collisions,0,false);
        deformable_collisions.Initialize_For_Thin_Shells_Fluid_Coupling();}
}
//#####################################################################
// Function Initialize_Solid_Fluid_Coupling
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Initialize_Solid_Fluid_Coupling()
{
    fluids_parameters.collision_bodies_affecting_fluid->Initialize_Grids();
}
//#####################################################################
// Function Get_Source_Velocities
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Get_Source_Velocities(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity)
{
    T_GRID::template Horizontal_Face_Loop<Get_Source_Velocities_Horizontal<GEOMETRY> >(incompressible.projection,source,world_to_source,constant_source_velocity);
    Get_Source_Velocities_Vertical<GEOMETRY>(source,world_to_source,constant_source_velocity);
}
//#####################################################################
// Function Get_Source_Velocities_Horizontal
//#####################################################################
template<class T_GRID> template<class GEOMETRY> template<class T_FACE> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Get_Source_Velocities_Horizontal<GEOMETRY>::Apply(PROJECTION_RLE<T_GRID>& projection,const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,
    const TV& constant_source_velocity)
{
    ARRAY<bool>& psi_N=projection.laplace.psi_N;ARRAY<T>& V=projection.V;T source_velocity=constant_source_velocity[T_FACE::Axis()];
    for(T_FACE face(projection.grid,2);face;face++)
        if(face.Both_Cells_Short()?source.Lazy_Inside(world_to_source.Homogeneous_Times(face.cell1.X())) || source.Lazy_Inside(world_to_source.Homogeneous_Times(face.cell2.X())):source.Lazy_Inside(world_to_source.Homogeneous_Times(face.X()))){
            int f=face.Face();psi_N(f)=true;V(f)=source_velocity;
            if(face.Long()){psi_N(f+1)=true;V(f+1)=source_velocity;}}
}
//#####################################################################
// Function Get_Source_Velocities_Vertical
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Get_Source_Velocities_Vertical(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const TV& constant_source_velocity)
{
    T_GRID& grid=*fluids_parameters.grid;
    ARRAY<bool>& psi_N=incompressible.projection.laplace.psi_N;ARRAY<T>& V=incompressible.V;T source_velocity=constant_source_velocity.y;
    for(FACE_Y_ITERATOR face(grid,2);face;face++)
        if(face.Both_Cells_Short()?source.Lazy_Inside(world_to_source.Homogeneous_Times(face.Cell_X(0))) || source.Lazy_Inside(world_to_source.Homogeneous_Times(face.Cell_X(1))):source.Lazy_Inside(world_to_source.Homogeneous_Times(face.X()))){
            int f=face.Face();psi_N(f)=true;V(f)=source_velocity;if(face.Upper_Cell_Long()){psi_N(f+1)=true;V(f+1)=source_velocity;}}
}
//#####################################################################
// Function Adjust_Phi_With_Source
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Adjust_Phi_With_Source(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,const bool only_near_existing_water)
{   
    T_GRID& grid=*fluids_parameters.grid;
    ARRAY<T>& phi=particle_levelset.phi;
    for(CELL_ITERATOR cell(grid,3);cell;cell++)if(cell.Short()){int c=cell.Cell();
        if(only_near_existing_water && phi(c)>=grid.positive_bandwidth) continue;
        phi(c)=min(phi(c),source.Signed_Distance(world_to_source.Homogeneous_Times(cell.X())));}
}
//#####################################################################
// Function Get_Source_Reseed_Mask
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Get_Source_Reseed_Mask(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,ARRAY<bool>*& cell_centered_mask,const bool reset_mask) const
{   
    const T_GRID& grid=*fluids_parameters.grid;
    if(reset_mask){if(cell_centered_mask) delete cell_centered_mask;cell_centered_mask=new ARRAY<bool>(grid.number_of_cells);}
    T padding=3*grid.uniform_grid.dX.Max();
    for(CELL_ITERATOR cell(grid,3);cell;cell++)if(cell.Short() && !source.Outside(world_to_source.Homogeneous_Times(cell.X()),padding)) (*cell_centered_mask)(cell.Cell())=true;
}
//#####################################################################
// Function Get_Source_Cell_Should_Be_Long
//#####################################################################
template<class T_GRID> template<class GEOMETRY> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Get_Source_Cell_Should_Be_Long(const GEOMETRY& source,const T_TRANSFORMATION_MATRIX& world_to_source,ARRAY<bool>& cell_should_be_long,const bool only_near_existing_water) const
{   
    const T_GRID& grid=*fluids_parameters.grid;
    T padding=3*grid.uniform_grid.dX.Max();
    for(CELL_ITERATOR cell(grid,3);cell;cell++)if(cell.Short() && abs(source.Signed_Distance(world_to_source.Homogeneous_Times(cell.X())))<=padding){int c=cell.Cell();
        if(only_near_existing_water && particle_levelset.phi(c)>=grid.positive_bandwidth) continue;
        cell_should_be_long(c)=false;}
}
//#####################################################################
// Function Total_Number_Of_Particles
//#####################################################################
template<class T_GRID> template<class T_PARTICLES> int SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Total_Number_Of_Particles(const ARRAY<T_PARTICLES*>& particles) const
{
    int total=0;
    for(int c=1;c<=particles.m;c++)if(particles(c)) total+=particles(c)->array_collection->Size();
    return total;
}
//#####################################################################
// Function Read_Particles
//#####################################################################
template<class T_GRID> template<class T_PARTICLES> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Read_Particles(const T_PARTICLES& template_particles,ARRAY<T_PARTICLES*>& particles,const std::string& prefix,const int frame)
{
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,prefix.c_str()),particles);
    if(typeid(template_particles)!=typeid(T_PARTICLES)) // swap in clones of template_particle for pure T_PARTICLESs
        for(int i=1;i<=particles.Size();i++) if(particles(i)){
            T_PARTICLES* replacement=template_particles.Clone();
            replacement->array_collection->Initialize(*particles(i)->array_collection);
            delete particles(i);
            particles(i)=replacement;}
    LOG::cout<<"Reading "<<Total_Number_Of_Particles(particles)<<" "<<prefix<<std::endl;
}
//#####################################################################
// Function Read_Particles
//#####################################################################
template<class T_GRID> template<class T_PARTICLES> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Read_Particles(const T_PARTICLES& template_particles,ARRAY<T_PARTICLES*>& particles,T_PARTICLES*& particles_in_long_cells,const std::string& prefix,const int frame)
{
    T_GRID& grid=*fluids_parameters.grid;
    Read_Particles(template_particles,particles,prefix,frame);
    if(grid.number_of_blocks>particles.m && particles.m<grid.number_of_blocks+1) PHYSBAM_FATAL_ERROR();
    if(particles.m==grid.number_of_blocks+1){delete particles_in_long_cells;particles_in_long_cells=particles(particles.m);}
    particles.Resize(grid.number_of_blocks);
}
//#####################################################################
// Function Write_Particles
//#####################################################################
template<class T_GRID> template<class T_PARTICLES> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Write_Particles(const ARRAY<T_PARTICLES*>& particles,const std::string& prefix,const int frame) const
{
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,prefix.c_str()),particles);
    LOG::cout<<"Writing "<<Total_Number_Of_Particles(particles)<<" "<<prefix<<std::endl;
}
//#####################################################################
// Function Write_Particles
//#####################################################################
template<class T_GRID> template<class T_PARTICLES> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Write_Particles(const ARRAY<T_PARTICLES*>& particles,const T_PARTICLES* particles_in_long_cells,const std::string& prefix,const int frame) const
{
    ARRAY<const T_PARTICLES*> all_particles(particles.m+1,false);
    for(int b=1;b<=particles.m;b++) all_particles(b)=particles(b);
    all_particles(all_particles.m)=particles_in_long_cells;
    Write_Particles(all_particles,prefix,frame);
}
//#####################################################################
// Function Read_Output_Files_Fluids
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Read_Output_Files_Fluids(const int frame)
{
    T_GRID& grid=*fluids_parameters.grid;
    GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>& collision_bodies_affecting_fluid=*fluids_parameters.collision_bodies_affecting_fluid;
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/rle_grid."+f,grid);
    // particle levelset
    if(fluids_parameters.write_levelset) FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/rle_levelset."+f,particle_levelset.levelset.phi);
    if(fluids_parameters.write_particles){
        Read_Particles(particle_levelset.template_particles,particle_levelset.positive_particles,"positive_particles",frame);
        Read_Particles(particle_levelset.template_particles,particle_levelset.negative_particles,"negative_particles",frame);}
    if(fluids_parameters.write_removed_positive_particles)
        Read_Particles(particle_levelset.template_removed_particles,particle_levelset.removed_positive_particles,particle_levelset.removed_positive_particles_in_long_cells,
            "removed_positive_particles",frame);
    if(fluids_parameters.write_removed_negative_particles)
        Read_Particles(particle_levelset.template_removed_particles,particle_levelset.removed_negative_particles,particle_levelset.removed_negative_particles_in_long_cells,
            "removed_negative_particles",frame);
    // deep water
    if(use_deep_water){
        if(deep_water.lambda) PHYSBAM_FATAL_ERROR();
        FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/deep_water_height."+f,deep_water.h);
        FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/deep_water_h_hats."+f,deep_water.h_hat1,deep_water.h_hat2);}
    // restarts
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/particle_levelset_random_seed."+f,particle_levelset.random);
    if(fluids_parameters.store_particle_ids) FILE_UTILITIES::Read_From_Text_File(output_directory+"/last_unique_particle_id."+f,particle_levelset.last_unique_particle_id);
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/rle_pressure."+f,incompressible.projection.p);
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/rle_face_velocities."+f,incompressible.projection.V);
    // read next bodies states for fluid collision bodies
    if(fluids_parameters.solid_affects_fluid && fluids_parameters.fluid_affects_solid && collision_bodies_affecting_fluid.collision_geometry_collection.bodies.m){
        try{
            std::istream* input_raw=FILE_UTILITIES::Safe_Open_Input(output_directory+"/fluid_collision_bodies_new_state."+f);TYPED_ISTREAM input(*input_raw,stream_type);
            collision_bodies_affecting_fluid.Read_State(input,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);delete input_raw;}
        catch(FILESYSTEM_ERROR&){LOG::cerr<<"PHYSBAM_WARNING: Could not read fluid collision bodies new state, restart may be incorrect!"<<std::endl;}}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T_GRID> void SOLIDS_FLUIDS_EXAMPLE_RLE<T_GRID>::
Write_Output_Files(const int frame) const
{
    const T_GRID& grid=*fluids_parameters.grid;
    const GRID_BASED_COLLISION_GEOMETRY_RLE<T_GRID>& collision_bodies_affecting_fluid=*fluids_parameters.collision_bodies_affecting_fluid;
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Create_Directory(output_directory);
    FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    Write_Frame_Title(frame);
    solids_parameters.Write_Output_Files(stream_type,output_directory,first_frame,frame);
    if(mpi_grid) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/global_uniform_grid",mpi_grid->global_uniform_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/rle_grid."+f,grid);
    // particle levelset
    if(fluids_parameters.write_levelset){
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/rle_levelset."+f,particle_levelset.levelset.phi);
        if(FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>* inaccurate_union=Find_Type<FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>*>(collision_bodies_affecting_fluid.collision_geometry_collection.bodies))
            FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/rle_object_levelset."+f,inaccurate_union->levelset.levelset.phi);}
    if(fluids_parameters.write_particles){
        Write_Particles(particle_levelset.positive_particles,"positive_particles",frame);
        Write_Particles(particle_levelset.negative_particles,"negative_particles",frame);}
    if(fluids_parameters.write_removed_positive_particles)
        Write_Particles(particle_levelset.removed_positive_particles,particle_levelset.removed_positive_particles_in_long_cells,"removed_positive_particles",frame);
    if(fluids_parameters.write_removed_negative_particles)
        Write_Particles(particle_levelset.removed_negative_particles,particle_levelset.removed_negative_particles_in_long_cells,"removed_negative_particles",frame);
    // deep water
    if(use_deep_water){
        if(deep_water.lambda) PHYSBAM_FATAL_ERROR();
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/deep_water_height."+f,deep_water.h);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/deep_water_h_hats."+f,deep_water.h_hat1,deep_water.h_hat2);}
    // restarts
    if(fluids_parameters.write_debug_data || (fluids_parameters.write_restart_data && frame%fluids_parameters.restart_data_write_rate==0)){
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/particle_levelset_random_seed."+f,particle_levelset.random);
        if(fluids_parameters.store_particle_ids) FILE_UTILITIES::Write_To_Text_File(output_directory+"/last_unique_particle_id."+f,particle_levelset.last_unique_particle_id);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/rle_pressure."+f,incompressible.projection.p);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/rle_face_velocities."+f,incompressible.projection.V);}
    // write next bodies states for fluid collision bodies
    if(fluids_parameters.solid_affects_fluid && fluids_parameters.fluid_affects_solid && collision_bodies_affecting_fluid.collision_geometry_collection.bodies.m && frame%fluids_parameters.restart_data_write_rate==0){
        std::ostream* output_raw=FILE_UTILITIES::Safe_Open_Output(output_directory+"/fluid_collision_bodies_new_state."+f);TYPED_OSTREAM output(*output_raw,stream_type);
        collision_bodies_affecting_fluid.Write_State(output,COLLISION_GEOMETRY<TV>::FLUID_COLLISION_GEOMETRY_NEW_STATE);delete output_raw;}
    // debugging
    if(fluids_parameters.write_debug_data){
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/rle_psi_D."+f,incompressible.projection.laplace.psi_D);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/rle_psi_N."+f,incompressible.projection.laplace.psi_N);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/rle_colors."+f,incompressible.projection.laplace.filled_region_colors);}
}
//#####################################################################
#define INSTANTIATION_HELPER_SHAPE(...) \
    template void SOLIDS_FLUIDS_EXAMPLE_RLE<RLE_GRID_POLICY<__VA_ARGS__::VECTOR_T>::RLE_GRID>::Adjust_Phi_With_Source(const __VA_ARGS__&,const T_TRANSFORMATION_MATRIX&,const bool); \
    template void SOLIDS_FLUIDS_EXAMPLE_RLE<RLE_GRID_POLICY<__VA_ARGS__::VECTOR_T>::RLE_GRID>::Get_Source_Velocities(const __VA_ARGS__&,const T_TRANSFORMATION_MATRIX&,const __VA_ARGS__::VECTOR_T&); \
    template void SOLIDS_FLUIDS_EXAMPLE_RLE<RLE_GRID_POLICY<__VA_ARGS__::VECTOR_T>::RLE_GRID>::Get_Source_Reseed_Mask(const __VA_ARGS__&,const T_TRANSFORMATION_MATRIX&,ARRAY<bool>*&,const bool) const; \
    template void SOLIDS_FLUIDS_EXAMPLE_RLE<RLE_GRID_POLICY<__VA_ARGS__::VECTOR_T>::RLE_GRID>::Get_Source_Cell_Should_Be_Long(const __VA_ARGS__&,const T_TRANSFORMATION_MATRIX&,ARRAY<bool>&,const bool) const;
#define INSTANTIATION_HELPER(T) \
    template class SOLIDS_FLUIDS_EXAMPLE_RLE<RLE_GRID_2D<T> >; \
    template class SOLIDS_FLUIDS_EXAMPLE_RLE<RLE_GRID_3D<T> >; \
    INSTANTIATION_HELPER_SHAPE(BOX<VECTOR<T,2> >) \
    INSTANTIATION_HELPER_SHAPE(BOX<VECTOR<T,3> >) \
    INSTANTIATION_HELPER_SHAPE(SPHERE<VECTOR<T,2> >) \
    INSTANTIATION_HELPER_SHAPE(SPHERE<VECTOR<T,3> >) \
    INSTANTIATION_HELPER_SHAPE(CYLINDER<T>)
INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
#endif 
