//#####################################################################
// Copyright 2008-2009, Elliot English, Nipun Kwatra, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA
//#####################################################################
#ifndef __SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA__
#define __SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS_BINARY_UNIFORM.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COUPLED_SYSTEM_VECTOR.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION_SLIP.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_CHIMERA.h>
#include <PhysBAM_Dynamics/Grids_Chimera/Grids_Chimera_Boundaries/BOUNDARY_CHIMERA.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology/TETRAHEDRON_MESH.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/CHIMERA_SYSTEM.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Dynamics/Grids_Chimera_PDE_Linear/Parallel_Computation/LAPLACE_CHIMERA_MPI.h>
#include <PhysBAM_Dynamics/Grids_Chimera_PDE_Linear/Parallel_Computation/LAPLACE_CHIMERA_GRID_MPI.h>
#include <PhysBAM_Dynamics/Grids_Chimera_Interpolation/AVERAGING_CHIMERA.h>
namespace PhysBAM{

template<class T_GRID> class SOLIDS_FLUIDS_DRIVER_CHIMERA;
template<class T> class FRACTURE_PATTERN;
template<class TV> class COLLISION_AWARE_INDEX_MAP;
template<class TV> class FLUIDS_PARAMETERS_UNIFORM;
template<class T_GRID> class POISSON_COLLIDABLE_UNIFORM;
template<class TV> class GENERALIZED_VELOCITY;
template<class TV> class GENERALIZED_MASS;
template<class TV> class MATRIX_SOLID_INTERPOLATION;
template<class TV> class SOLIDS_FLUIDS_PARAMETERS;
template<class TV> class SLIP_SYSTEM;
template<class TV> class BACKWARD_EULER_SYSTEM;
template<class T_GRID> class INCOMPRESSIBLE_FLUID_CONTAINER;
template<class T_GRID> class EULER_PROJECTION_UNIFORM;
template<class TV> class IMPLICIT_BOUNDARY_CONDITION_COLLECTION;
template<class TV> class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO;
template<class T> class SPARSE_MATRIX_FLAT_MXN;
template<class T> class SPARSE_MATRIX_FLAT_NXN;

template<class T_GRID>
class SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA:public NEWMARK_EVOLUTION<typename T_GRID::VECTOR_T>//public SOLID_FLUID_COUPLED_EVOLUTION_SLIP<typename T_GRID::VECTOR_T>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    typedef FACE_INDEX<TV::dimension> D_FACE_INDEX;
    typedef NEWMARK_EVOLUTION<TV> BASE;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_BASE::template REBIND<bool>::TYPE T_ARRAYS_BASE_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<ARRAY<int> >::TYPE T_ARRAYS_ARRAYS_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<int>::TYPE T_FACE_ARRAYS_INT;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_BOOL::template REBIND<VECTOR<bool,T_GRID::dimension> >::TYPE T_FACE_ARRAYS_BOOL_DIMENSION;
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;typedef typename LEVELSET_POLICY<T_GRID>::FAST_LEVELSET_T T_FAST_LEVELSET;
    typedef typename T_ARRAYS_SCALAR::template REBIND<PARTICLE_LEVELSET_PARTICLES<TV>*>::TYPE T_ARRAYS_PARTICLE_LEVELSET_PARTICLES;
    typedef typename T_ARRAYS_SCALAR::template REBIND<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>::TYPE T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES;
    typedef EXTRAPOLATION_UNIFORM<T_GRID,T> T_EXTRAPOLATION_SCALAR;
    typedef ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T, AVERAGING_UNIFORM<GRID<TV>, FACE_LOOKUP_UNIFORM<GRID<TV> > >,LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T,FACE_LOOKUP_UNIFORM<GRID<TV> > > > T_ADVECTION_UNIFORM;
    typedef ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_ALE<GRID<TV>,T,FACE_LOOKUP_UNIFORM<T_GRID>,FACE_LOOKUP_COLLIDABLE_UNIFORM<T_GRID> > T_ADVECTION_COLLIDABLE_ALE;
    typedef ADVECTION_WRAPPER_ALE<GRID<TV>,T,T_ADVECTION_UNIFORM,FACE_LOOKUP_UNIFORM<T_GRID> > T_ADVECTION_WRAPPER_ALE;
    typedef ADVECTION_WRAPPER_MACCORMACK_CHIMERA<GRID<TV>,T,T_ADVECTION_UNIFORM,FACE_LOOKUP_UNIFORM<T_GRID> > T_ADVECTION_MACCORMACK_CHIMERA;
    typedef PAIR<int,TV_INT> GRID_CELL_INDEX;
    typedef VECTOR<GRID_CELL_INDEX,2> VORONOI_FACE_INDICES;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::CONVEX_POLYGON CONVEX_POLYGON;
    typedef PAIR<VORONOI_FACE_INDICES,CONVEX_POLYGON> VORONOI_FACE;
    typedef LAPLACE_CHIMERA_MPI<T_GRID> T_LAPLACE;

public:
    SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>& example;

    LAPLACE_CHIMERA_GRID_MPI<T_GRID> laplace_grid;
    LAPLACE_CHIMERA_MPI<T_GRID> laplace_incompressible;
    LAPLACE_CHIMERA_MPI<T_GRID> laplace_diffusion;
    
    //coupled face solve
    ARRAY<LAPLACE_CHIMERA_MPI<T_GRID>*> laplace_diffusion_components;
    VECTOR_ND<T> face_temperature_packed,cell_temperature_packed;
    AVERAGING_CHIMERA<T_GRID> averaging;

    VECTOR<ARRAY<T_ARRAYS_SCALAR>,TV::dimension> cell_velocity_components,cell_velocity_components_updated;
    VECTOR<ARRAY<ARRAY<T> >,TV::dimension> cell_velocity_boundary_components,cell_velocity_boundary_components_updated;
    
    BOUNDARY_CHIMERA<T_GRID,T>* boundary;
    ARRAY<T_FACE_ARRAYS_BOOL> above_water_face_masks;
    
    SOLIDS_FLUIDS_DRIVER_CHIMERA<T_GRID>* driver;
    int n_local_grids;
    int n_global_grids;
    
    VECTOR_ND<T> rhs_incompressible;
    VECTOR_ND<T> x_incompressible;
    
    //SOLUTION VARIABLES
    VECTOR_ND<T> face_velocities_packed;
    ARRAY<T> voronoi_face_velocities;
    //ARRAY<T_ARRAYS_SCALAR> pressures;
    ARRAY<ARRAY<T> > pressures_boundary;

    SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA(SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>& example_input):
        BASE(example_input.solids_parameters,example_input.solid_body_collection),
            example(example_input),laplace_grid(*example.chimera_grid),laplace_incompressible(*example.chimera_grid,laplace_grid),laplace_diffusion(*example.chimera_grid,laplace_grid),averaging(laplace_grid)
    {
        //minimum_dual_cell_volume_fractions=1e-5;
        n_global_grids=example.chimera_grid->number_of_global_grids;
        n_local_grids=example.rigid_grid_collection.particles.array_collection->Size();

        vorticity.Resize(n_global_grids);

        example.solids_evolution=this;
        boundary=new BOUNDARY_CHIMERA<T_GRID>(example.chimera_grid,*(new BOUNDARY_UNIFORM<T_GRID,T>()));example.boundary_chimera=boundary;
        if(example.fluids_parameters.water) example.Initialize_Boundary_Phi_Water(&(example.incompressible_fluid_containers(1)->face_velocities));

        //example.advection_uniform.number_of_ghost_cells=example.fluids_parameters.number_of_ghost_cells;
        
        if(example.use_maccormack_advection_scalar)
            example.advection_maccormack_scalar=new T_ADVECTION_MACCORMACK_CHIMERA(example.advection_uniform,example,example.clamp_extrema_maccormack_scalar,example.enforce_second_order_maccormack_scalar);
        if(example.use_maccormack_advection_vector)
            example.advection_maccormack_vector=new T_ADVECTION_MACCORMACK_CHIMERA(example.advection_uniform,example,example.clamp_extrema_maccormack_vector,example.enforce_second_order_maccormack_vector);

        if(!example.simulate_solids || !example.use_collidable_advection){
            for(int grid_index=1;grid_index<=n_local_grids;grid_index++) 
                example.advection.Append(new T_ADVECTION_WRAPPER_ALE(example.advection_uniform,example.rigid_grid_collection.Rigid_Grid(grid_index),example.fluids_parameters.number_of_ghost_cells));}
        else{
            T default_cell_replacement_value=example.fluids_parameters.water?(T)example.fluids_parameters.collidable_phi_replacement_value:(T)0;
            for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
                if(example.fluids_parameters.smoke){
                    example.advection_collidable.Append(new T_ADVECTION_COLLIDABLE_ALE(example.rigid_grid_collection.Rigid_Grid(grid_index),example.incompressible_fluid_containers(grid_index)->collision_bodies_affecting_fluid,example.incompressible_fluid_containers(grid_index)->face_velocities_valid_mask,example.incompressible_fluid_containers(grid_index)->density_container.valid_mask_current,example.incompressible_fluid_containers(grid_index)->density_container.valid_mask_next,default_cell_replacement_value,false));
                    example.advection.Append(example.advection_collidable.Last());}
                else if(example.fluids_parameters.water){
                    example.advection_collidable.Append(new T_ADVECTION_COLLIDABLE_ALE(example.rigid_grid_collection.Rigid_Grid(grid_index),example.incompressible_fluid_containers(grid_index)->collision_bodies_affecting_fluid,example.incompressible_fluid_containers(grid_index)->face_velocities_valid_mask,example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.levelset.valid_mask_current,example.incompressible_fluid_containers(grid_index)->particle_levelset_evolution.particle_levelset.levelset.valid_mask_next,default_cell_replacement_value,false));
                    example.advection.Append(example.advection_collidable.Last());}}}

        for(int axis=1;axis<=TV::dimension;axis++)
            laplace_diffusion_components.Append(new LAPLACE_CHIMERA_MPI<T_GRID>(*example.chimera_grid,laplace_grid));
    };
    virtual ~SOLID_FLUID_COUPLED_EVOLUTION_CHIMERA(){
        for(int i=1;i<=laplace_diffusion_components.Size();i++)
            delete laplace_diffusion_components(i);
        laplace_diffusion_components.Remove_All();
        delete boundary;boundary=0;
        //should delete advection here as well
    };

    bool Simulate_Fluids() const
    {const FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters=example.fluids_parameters;SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters=example.solids_fluids_parameters;
    return (solids_fluids_parameters.mpi_solid_fluid || fluids_parameters.simulate) && (fluids_parameters.smoke || fluids_parameters.fire || fluids_parameters.water || fluids_parameters.sph || fluids_parameters.compressible);}

    bool Simulate_Solids() const
    {SOLIDS_PARAMETERS<TV>& solids_parameters=example.solids_parameters;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=example.solid_body_collection.deformable_body_collection;
    return (deformable_body_collection.simulate && deformable_body_collection.particles.array_collection->Size()) || (solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies && example.solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size());}
    
    void Advance_One_Time_Step(const T time,const T dt);
    T Calculate_Maximum_Allowable_dt(const T time);
    void Setup_Solids(const T time);
    void Setup_Fluids(const T time);
    void Advect_Fluid(const T time,const T dt);

    void Fill_Ghost_Regions_Before_Advection(ARRAY<T_FACE_ARRAYS_SCALAR>& face_velocities_ghost,ARRAY<T_ARRAYS_SCALAR>& scalar_ghost,const T time,const T dt);
    void Advect_Scalar_Fields(ARRAY<T_FACE_ARRAYS_SCALAR>& face_velocities_ghost,ARRAY<T_ARRAYS_SCALAR>& scalar_ghost,const T time,const T dt);
    void Advect_Vector_Fields(ARRAY<T_FACE_ARRAYS_SCALAR>& face_velocities_ghost,const T time,const T dt);
    void Revalidation();
    void Modify_Levelset_And_Particles_After_Advection(ARRAY<T_FACE_ARRAYS_SCALAR>& face_velocities_ghost,const T time,const T dt);
    void Coupling_Scalar_Fields(const T time,const T dt);
    void Postprocess_After_Advection(const T time,const T dt);

    void Coupling_Vector_Fields();
    void Revalidate_Fluid_Scalars();
    void Revalidate_Phi_After_Modify_Levelset();
    void Revalidate_Fluid_Velocity();
    void Adjust_Particles_After_Updating_Grid(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,RIGID_GRID<T_GRID>& rigid_grid,const T dt,const T time);
    void Adjust_Removed_Particles_After_Updating_Grid(T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES& particles,RIGID_GRID<T_GRID>& rigid_grid,const T dt,const T time);
    void Adjust_Face_Scalar_Field_After_Updating_Grid(T_FACE_ARRAYS_SCALAR& u, T_FACE_ARRAYS_SCALAR u_ghost, RIGID_GRID<T_GRID>& rigid_grid, const T dt, const T time);
    void Extrapolate_Velocity_Across_Interface(T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_velocities,const T_ARRAYS_SCALAR& phi_ghost,const T number_of_ghost_cells,const T band_width=3,const T damping=0,const TV& air_speed=TV(),const T_FACE_ARRAYS_BOOL_DIMENSION* face_neighbors_visible=0,const T_FACE_ARRAYS_BOOL* fixed_faces_input=0);
    void Integrate_Fluid_Non_Advection_Forces(T_FACE_ARRAYS_SCALAR& face_velocities,const RIGID_GRID<T_GRID>& rigid_grid,const T dt,const int substep);
    void Solid_Position_Update(const T time,const T dt);
    void Solid_Velocity_Update(const T time,const T dt);
    
    void Advance_Fluid_One_Time_Step_Implicit_Part(const T time,const T dt);

    T Phi(const int grid_index,const TV_INT& cell_index)
    {
        //ASSUME THAT GHOST PHI IS VALID
        //ONE GHOST CELL IS NEEDED FOR COMPUTING BOUNDARY FACE VOLUME FRACTIONS
        if(laplace_grid.Local_Grid(grid_index)) return example.incompressible_fluid_containers(laplace_grid.Local_Grid_Index(grid_index))->particle_levelset_evolution.phi(cell_index);
        else return laplace_incompressible.boundary_cell_phi(grid_index)(laplace_grid.boundary_cell_indices_to_linear_index(grid_index).Get(cell_index));
    }

    void Pack_And_Compute_Voronoi_Face_Values(LAPLACE_CHIMERA_MPI<T_GRID>& laplace,ARRAY<T_FACE_ARRAYS_SCALAR*>& face_values,VECTOR_ND<T>& packed_values,const T time,const T dt);
    void Build_Diffusion_System(LAPLACE_CHIMERA_MPI<T_GRID>& laplace,LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time,const T dt,ARRAY<T_ARRAYS_SCALAR*>& u,ARRAY<ARRAY<T> >& u_boundary,VECTOR_ND<T>& x,VECTOR_ND<T>& rhs,const bool compute_explicit_forces_only=false);
    void Build_Poisson_System(LAPLACE_CHIMERA_MPI<T_GRID>& laplace,LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time,const T dt,VECTOR_ND<T>* packed_values,VECTOR_ND<T>& x,VECTOR_ND<T>& rhs);
    void Solve_System(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& rhs,VECTOR_ND<T>& x,T tolerance,int iterations,SPARSE_MATRIX_PARTITION* partition=0); //COMMUNICATION //SOLVE THE SYSTEM
                                                                                                                                                       //AND UPDATE BOUNDARY CELLS OF
                                                                                                                                                       //NEIGHBORING GRIDS EACH ITERATION
    void Unpack_Cell_Values(LAPLACE_CHIMERA_MPI<T_GRID>& laplace,LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time,ARRAY<T_ARRAYS_SCALAR*>& cell_values,VECTOR_ND<T>* packed_values);
    void Set_Dirichlet_Boundary_Conditions(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time,ARRAY<T_ARRAYS_SCALAR*>& cell_values);
    void Pack_And_Exchange_Boundary_Cell_Values(ARRAY<T_ARRAYS_SCALAR*>& cell_values,ARRAY<ARRAY<T> >& boundary_cell_values);
    void Interpolate_Cell_Values_From_Delaunay_Simplices(ARRAY<T_ARRAYS_SCALAR*>& cell_values,ARRAY<ARRAY<T> >& boundary_cell_values);
    void Unpack_And_Interpolate_And_Inject_Cell_Values(LAPLACE_CHIMERA_MPI<T_GRID>& laplace,LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time,ARRAY<T_ARRAYS_SCALAR*>& cell_values,ARRAY<ARRAY<T> >& boundary_cell_values,VECTOR_ND<T>* packed_values);
    
    typedef MATRIX<T,(1+TV::dimension)*TV::dimension,(1+TV::dimension)*TV::dimension> D2_MATRIX;
    typedef VECTOR<T,(1+TV::dimension)*TV::dimension> TV2;
    void Add_Cell_Faces_Least_Squares_Terms(ARRAY<T_FACE_ARRAYS_SCALAR*>& face_velocities,const GRID_CELL_INDEX& grid_cell_index,const TV& x,D2_MATRIX& A,TV2& rhs,HASHTABLE<int>& added_voronoi_faces);
    TV Compute_Cell_Center_Velocity(ARRAY<T_FACE_ARRAYS_SCALAR*>& face_velocities,const int grid_index,const TV_INT& cell_index);
    void Update_Velocities(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,ARRAY<T_FACE_ARRAYS_SCALAR*>& face_velocities,ARRAY<T_ARRAYS_SCALAR*>& pressures,const T time,const T dt);
    TV Compute_Cell_Center_Velocity_Local(ARRAY<T_FACE_ARRAYS_SCALAR*>& face_velocities,const int grid_index,const TV_INT& cell_index);    

    void Set_Zero_Air_Velocities(); //SET AIR FACE VELOCITIES TO ZERO //DO NOT CALL THIS AFTER FILLED GHOST CELL FOR EXTRAPOLATION - MIGHT BLOW CFL
    void Extrapolate_Velocity_Across_Interface_On_All_Local_Grids(); //EXTRAPOLATE VELOCITY ACROSS INTERFACE //DONE AT THE END OF FRAME TO ENSURE CORRECT CFL CONDITION
    
    void Write_Substep(const std::string& title,const int substep,const int level);

    ARRAY<T_ARRAYS_SCALAR> vorticity;
    void Compute_Vorticity(const T time,const T dt);
//#####################################################################
};
}
#endif
