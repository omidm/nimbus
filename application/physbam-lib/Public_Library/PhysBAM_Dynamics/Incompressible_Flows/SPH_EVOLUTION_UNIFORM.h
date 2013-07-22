//#####################################################################
// Copyright 2006, Nipun Kwatra, Frank Losasso, Nick Rasmussen, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPH_EVOLUTION_UNIFORM
//#####################################################################
#ifndef __SPH_EVOLUTION_UNIFORM__
#define __SPH_EVOLUTION_UNIFORM__

#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_FORWARD.h>
#include <PhysBAM_Dynamics/Particles/PARTICLES_FORWARD.h>
#include <PhysBAM_Dynamics/Particles/SPH_PARTICLES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID> class PARTICLE_LEVELSET_EVOLUTION_UNIFORM;

template<class T_GRID>
class SPH_EVOLUTION_UNIFORM:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    STATIC_ASSERT((IS_SAME<typename T_GRID::GRID_TAG,UNIFORM_TAG<TV> >::value));
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<ARRAY<int> >::TYPE T_ARRAYS_ARRAY_INT;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>::TYPE T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES;
    typedef typename T_ARRAYS_SCALAR::template REBIND<PARTICLE_LEVELSET_PARTICLES<TV>*>::TYPE T_ARRAYS_PARTICLE_LEVELSET_PARTICLES;
    typedef typename T_ARRAYS_SCALAR::template REBIND<ARRAY<PAIR<TV_INT,int> > >::TYPE T_ARRAYS_ARRAY_PAIR_TV_INT_INT;
    typedef typename T_ARRAYS_SCALAR::template REBIND<ARRAY<T> >::TYPE T_ARRAYS_ARRAY_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_INT>::TYPE T_ARRAYS_TV_INT;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
public:

    T_GRID grid;
    const SPH_CALLBACKS<T_GRID>* callbacks;
    INCOMPRESSIBLE_UNIFORM<T_GRID>& incompressible;
    FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>* particle_levelset_evolution;
    SPH_PARTICLES<TV> sph_particles;
    T_FACE_ARRAYS_BOOL valid_particle_face_velocities;

    T particle_radius;
    T target_particles_per_unit_volume;
    T ballistic_particles_per_unit_volume;
    bool use_analytic_divergence,use_analytic_divergence_for_expansion_only,adjust_cell_weights_on_neumann_boundaries,enforce_density_near_interface;
    T particle_targeting_time;
    T divergence_expansion_multiplier,grid_to_particle_slip_multiplier,neumann_boundary_slip_multiplier;
    T attraction_strength,attraction_restlength_cells,attraction_forward_move_multiplier;
    int attraction_skip_number;
    bool stable_attraction_integration,uniform_attraction_force;
    T flip_ratio;
    bool use_variable_density_solve;
    bool use_two_way_coupling;
    bool convert_particles_to_fluid;
    int maximum_phi_refinement_depth;
    T_ARRAYS_SCALAR cell_weight;
private:
    T radius,radius_plus_half_dx_squared,one_over_radius_squared;
    TV radius_vector;
    T target_particles_per_cell,target_minus_ballistic_particles_per_cell,one_over_target_particles_per_cell,one_over_target_minus_ballistic_particles_per_cell,ballistic_particles_per_cell;

    T_ARRAYS_ARRAY_INT particles_in_cell;
    T_FACE_ARRAYS_SCALAR face_weight;
    ARRAY<T> one_over_total_particle_cell_weight,one_over_total_particle_face_weight;
    T_FACE_ARRAYS_SCALAR particle_velocities;
    RANDOM_NUMBERS<T> random;
public:

    SPH_EVOLUTION_UNIFORM(T_GRID& grid_input,INCOMPRESSIBLE_UNIFORM<T_GRID>& incompressible_input,FLUIDS_PARAMETERS_UNIFORM<T_GRID>& fluids_parameters_input,
        PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>* particle_levelset_evolution=0);
    ~SPH_EVOLUTION_UNIFORM();

    void Initialize_Grids(const T_GRID& grid_input)
    {grid=grid_input.Get_MAC_Grid();}

    void Set_SPH_Callbacks(const SPH_CALLBACKS<T_GRID>& callbacks_input)
    {callbacks=&callbacks_input;}

//#####################################################################
    void Euler_Step(const T dt,const T time);
    T CFL() const;
    template<class T_ARRAYS_PARTICLES> void Copy_Particle_Attributes_From_Array(T_ARRAYS_PARTICLES& particles);
    template<class T_ARRAYS_PARTICLES> void Copy_Particle_Attributes_To_Array(T_ARRAYS_PARTICLES& particles) const;
    template<class T_ARRAYS_PARTICLES> void Make_Incompressible(T_ARRAYS_PARTICLES& particles,T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Make_Incompressible(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Calculate_SPH_Constants();
    void Set_Up_For_Projection(T_FACE_ARRAYS_SCALAR& face_velocities,const T time);
    void Set_Divergence_And_Multiplier(const TV_INT cell,const T_ARRAYS_BOOL& cells_valid,const T time);
    void Postprocess_Particles(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Rasterize_Velocities_To_Grid(T_FACE_ARRAYS_SCALAR& velocities,T_ARRAYS_SCALAR& cell_weight,T_FACE_ARRAYS_SCALAR& face_weight);
    void Calculate_Particle_Deltas(const T_FACE_ARRAYS_SCALAR& minus_face_delta,ARRAY<TV>& delta_velocity,ARRAY<TV>& delta_weight);
    void Modify_Levelset_And_Particles_To_Create_Fluid(const T time,T_FACE_ARRAYS_SCALAR* face_velocities);
    template<class T_ARRAYS_PARTICLES> void Move_Particles_Off_Grid_Boundaries(T_ARRAYS_PARTICLES& particles,const T tolerance) const;
//#####################################################################
};
}
#endif
