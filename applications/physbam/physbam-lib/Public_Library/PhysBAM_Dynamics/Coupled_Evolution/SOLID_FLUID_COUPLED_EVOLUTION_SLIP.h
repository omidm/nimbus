//#####################################################################
// Copyright 2008-2009, Elliot English, Nipun Kwatra, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_FLUID_COUPLED_EVOLUTION_SLIP
//#####################################################################
#ifndef __SOLID_FLUID_COUPLED_EVOLUTION_SLIP__
#define __SOLID_FLUID_COUPLED_EVOLUTION_SLIP__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS_BINARY_UNIFORM.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_DYNAMICS_UNIFORM.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COUPLED_SYSTEM_VECTOR.h>
namespace PhysBAM{

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

template<class TV_input>
class SOLID_FLUID_COUPLED_EVOLUTION_SLIP:public NEWMARK_EVOLUTION<TV_input>,public PROJECTION_DYNAMICS_UNIFORM<GRID<TV_input> >
{
    typedef TV_input TV;typedef typename TV::SCALAR T;typedef GRID<TV> T_GRID;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<int>::TYPE T_FACE_ARRAYS_INT;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;

protected:
    typedef NEWMARK_EVOLUTION<TV> BASE;
    typedef PROJECTION_DYNAMICS_UNIFORM<T_GRID> FLUID_BASE;
    using BASE::solid_body_collection;using BASE::solids_parameters;using BASE::F_full;using BASE::rigid_F_full;using BASE::R_full;using BASE::rigid_R_full;using BASE::S_full;
    using BASE::rigid_S_full;using BASE::B_full;using BASE::rigid_B_full;using BASE::repulsions;using BASE::rigids_evolution_callbacks;using BASE::rigid_body_collisions;
    using FLUID_BASE::p;using FLUID_BASE::poisson;using BASE::V_save;using BASE::rigid_velocity_save;

    COUPLED_SYSTEM_VECTOR<TV> coupled_x,coupled_f,coupled_r,coupled_b,coupled_s,coupled_ar,coupled_z;

    GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >& collision_bodies;
    FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters;
    SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters;
    INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container;
    ARRAY<TV> ar_full,z_full,zaq_full,leakproof_empty_V,temp_solid_full,second_temp_solid_full;
    ARRAY<TWIST<TV> > rigid_ar_full,rigid_z_full,rigid_zaq_full,rigid_temp_solid_full,rigid_second_temp_solid_full,leakproof_empty_rigid_V;
    
    T divergence_scaling;
    bool disable_thinshell;

public:
    UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& iterator_info;
    IMPLICIT_BOUNDARY_CONDITION_COLLECTION<TV>& boundary_condition_collection;
protected:

    T_GRID& grid;
    T_FACE_ARRAYS_SCALAR fluids_face_velocities; // Stores combined compressible-incompressible face velocities.
public:
    T_ARRAYS_SCALAR pressure;
    T_ARRAYS_SCALAR density;
    T_FACE_ARRAYS_BOOL solved_faces;
protected:
    T_ARRAYS_VECTOR centered_velocity;
    T_ARRAYS_SCALAR one_over_rho_c_squared;
    T_ARRAYS_SCALAR p_advected_over_rho_c_squared_dt;
    T_ARRAYS_SCALAR p_advected;

public:
    FRACTURE_PATTERN<T>* fracture_pattern;
private:
    ARRAY<TV> pressure_impulses;
    ARRAY<TWIST<TV> > pressure_impulses_twist;

    T_FACE_ARRAYS_BOOL cached_psi_N;
    ARRAY<T,COUPLING_CONSTRAINT_ID> coupling_face_velocities_cached;
    COUPLING_CONSTRAINT_ID number_of_coupling_faces_cached;
    ARRAY<FACE_INDEX<TV::dimension> > indexed_faces_cached;
    T time_cached;
    bool cached_coupling_face_data;
public:
    using BASE::print_matrix;
    bool run_self_tests;
    bool print_poisson_matrix;
    bool print_index_map;
    bool print_rhs;
    bool print_each_matrix;
    bool output_iterators;
    bool use_viscous_forces;
    bool two_phase;
    bool use_full_ic;

    SOLID_FLUID_COUPLED_EVOLUTION_SLIP(SOLIDS_PARAMETERS<TV>& solids_parameters_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,
        FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters_input,SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters_input,INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& incompressible_fluid_container_input);
    virtual ~SOLID_FLUID_COUPLED_EVOLUTION_SLIP();

    bool Cell_To_Cell_Visible(const int axis,const TV_INT& first_cell,const TV_INT& second_cell) const
    {assert(first_cell==second_cell+TV_INT::Axis_Vector(axis));return collision_bodies.cell_neighbors_visible(second_cell)(axis);}

//#####################################################################
    void Initialize_Grid_Arrays();
    bool Simulate_Incompressible_Fluids() const;
    bool Simulate_Compressible_Fluids() const;
    bool Simulate_Fluids() const;
    bool Simulate_Solids() const;
    void Backward_Euler_Step_Velocity_Helper(const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update) PHYSBAM_OVERRIDE;
    void Process_Collisions(const T dt,const T time,const bool advance_rigid_bodies) PHYSBAM_OVERRIDE;
    void Apply_Pressure(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,bool scale_by_dt=false);
    BACKWARD_EULER_SYSTEM<TV>* Setup_Solids(const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update,const bool leakproof_solve);
    void Setup_Fluids(T_FACE_ARRAYS_SCALAR& incompressible_face_velocities,const T current_position_time,const T dt,const bool leakproof_solve);
    void Solve(T_FACE_ARRAYS_SCALAR& incompressible_face_velocities,const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update,const bool leakproof_solve);
    void Make_Divergence_Free(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Apply_Second_Order_Cut_Cell_Method(const T_ARRAYS_INT& cell_index_to_divergence_matrix_index,const T_FACE_ARRAYS_INT& face_index_to_matrix_index,VECTOR_ND<T>& b);
    void Apply_Viscosity(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Setup_Boundary_Condition_Collection();
private:
    void Warn_For_Exposed_Dirichlet_Cell(const T_ARRAYS_BOOL& psi_D,const T_FACE_ARRAYS_BOOL& psi_N);
    void Set_Cached_Psi_N_And_Coupled_Face_Data(const COLLISION_AWARE_INDEX_MAP<TV>& index_map,
        const MATRIX_SOLID_INTERPOLATION<TV>& solid_interpolation,const T time);
    void Fill_Coupled_Face_Data(const COUPLING_CONSTRAINT_ID number_of_coupling_faces,const ARRAY<FACE_INDEX<SOLID_FLUID_COUPLED_EVOLUTION_SLIP<TV>::TV::dimension> >& indexed_faces,
        const ARRAY<T,COUPLING_CONSTRAINT_ID>& coupling_face_data,T_FACE_ARRAYS_SCALAR& face_data);
    void Get_Coupled_Faces_And_Interpolated_Solid_Velocities(const COLLISION_AWARE_INDEX_MAP<TV>& index_map,
        const MATRIX_SOLID_INTERPOLATION<TV>& solid_interpolation,const T_FACE_ARRAYS_BOOL& psi_N_domain,T_FACE_ARRAYS_BOOL& psi_N,
        ARRAY<T,COUPLING_CONSTRAINT_ID>& coupling_face_velocities);
public:
    void Get_Coupled_Faces_And_Interpolated_Solid_Velocities(const T time,T_FACE_ARRAYS_BOOL& psi_N,T_FACE_ARRAYS_SCALAR& face_velocities);
    void Output_Iterators(const STREAM_TYPE stream_type,const char* output_directory,int frame) const;
//#####################################################################
};
}
#endif
