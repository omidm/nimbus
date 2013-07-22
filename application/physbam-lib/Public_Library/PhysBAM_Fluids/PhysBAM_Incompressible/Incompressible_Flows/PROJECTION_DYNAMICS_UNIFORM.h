//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTION_DYNAMICS_UNIFORM  
//#####################################################################
#ifndef __PROJECTION_DYNAMICS_UNIFORM__
#define __PROJECTION_DYNAMICS_UNIFORM__

#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/PROJECTION_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_DYNAMICS.h>
#include <PhysBAM_Dynamics/Interpolation/FIRE_INTERPOLATION_POLICY.h>
namespace PhysBAM{

template<class T_GRID> class DETONATION_SHOCK_DYNAMICS;

template<class T_GRID>
class PROJECTION_DYNAMICS_UNIFORM:public PROJECTION_COLLIDABLE_UNIFORM<T_GRID>,public PROJECTION_DYNAMICS<typename T_GRID::SCALAR>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_ARRAYS_SLIP T_FACE_ARRAYS_SLIP_SCALAR;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;typedef typename FIRE_INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP_FIRE_MULTIPHASE T_FACE_LOOKUP_FIRE_MULTIPHASE;
public:
    typedef PROJECTION_COLLIDABLE_UNIFORM<T_GRID> BASE;
    using BASE::use_non_zero_divergence;using BASE::p_grid;using BASE::poisson;using BASE::elliptic_solver;using BASE::laplace;using BASE::p;using BASE::collidable_solver;using BASE::use_divergence_multiplier;using BASE::divergence;
    using BASE::divergence_multiplier;using BASE::poisson_collidable;using BASE::laplace_collidable;
    using BASE::Use_Divergence_Multiplier;using BASE::Use_Non_Zero_Divergence;using BASE::Compute_Divergence;using BASE::density;
    using PROJECTION_DYNAMICS<T>::flame;using PROJECTION_DYNAMICS<T>::flame_speed_constants;
    
    bool use_flame_speed_multiplier;
    T_FACE_ARRAYS_SCALAR flame_speed_multiplier;
    DETONATION_SHOCK_DYNAMICS<T_GRID>* dsd;

protected:
    bool use_divergence_multiplier_save_for_sph,use_non_zero_divergence_save_for_sph;
    T_ARRAYS_SCALAR *p_save_for_sph,*divergence_save_for_sph,*divergence_multiplier_save_for_sph;
    T_FACE_ARRAYS_SCALAR *face_velocities_save_for_sph;
    LAPLACE_UNIFORM<T_GRID>* elliptic_solver_save_for_sph;
    LAPLACE_COLLIDABLE_UNIFORM<T_GRID>* laplace_save_for_sph;
    POISSON_COLLIDABLE_UNIFORM<T_GRID>* poisson_save_for_sph;
    LAPLACE_COLLIDABLE<T_GRID>* collidable_solver_save_for_sph;
public:

    PROJECTION_DYNAMICS_UNIFORM(const T_GRID& mac_grid,const bool flame_input=false,const bool multiphase=false,const bool use_variable_beta=false,const bool use_poisson=false,THREAD_QUEUE* thread_queue=0);
    PROJECTION_DYNAMICS_UNIFORM(const T_GRID& mac_grid,T_LEVELSET& levelset_input);
    virtual ~PROJECTION_DYNAMICS_UNIFORM();

    T Face_Velocity_With_Ghost_Value_Multiphase(const T_ARRAYS_BASE& face_velocities_ghost,const int axis,const TV_INT& face_index,const int current_region) const
    {assert(flame);return face_velocities_ghost(face_index)-Face_Jump_Multiphase(axis,face_index,current_region);}

    T Face_Velocity_With_Ghost_Value_Multiphase(const T_ARRAYS_BASE& face_velocities_ghost,const int axis,const TV_INT& face_index,const int current_region,const int face_region) const
    {assert(flame);return face_velocities_ghost(face_index)-Face_Jump_Multiphase(axis,face_index,current_region,face_region);}

    T Face_Velocity_With_Ghost_Value_Multiphase(const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,const int axis,const TV_INT& face_index,const int current_region) const
    {return Face_Velocity_With_Ghost_Value_Multiphase(face_velocities_ghost.Component(axis),axis,face_index,current_region);}

    T Face_Velocity_With_Ghost_Value_Multiphase(const T_FACE_ARRAYS_SCALAR& face_velocities_ghost,const int axis,const TV_INT& face_index,const int current_region,const int face_region) const
    {return Face_Velocity_With_Ghost_Value_Multiphase(face_velocities_ghost.Component(axis),axis,face_index,current_region,face_region);}

    T Face_Jump_Multiphase(const int axis,const TV_INT& face_index,const int current_region) const
    {int face_region=poisson_collidable->levelset_multiple->Inside_Region_Face(axis,face_index);return Face_Jump_Multiphase(axis,face_index,current_region,face_region);}

//#####################################################################
    virtual void Initialize_Grid(const T_GRID& mac_grid);
    void Initialize_Dsd(const LEVELSET_MULTIPLE_UNIFORM<T_GRID>& levelset_multiple,const ARRAY<bool>& is_fuel_region);
    void Initialize_Dsd(const T_LEVELSET& levelset,const ARRAY<bool>& fuel_region);
    virtual void Make_Divergence_Free(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Update_Potential_Energy(T_FACE_ARRAYS_SCALAR& face_velocities,T_FACE_ARRAYS_SCALAR& potential_energy,const T dt);
    void Update_Potential_Energy(T_FACE_ARRAYS_SCALAR& face_velocities,T_ARRAYS_SCALAR& potential_energy,const T dt);
    void Update_Potential_Energy(T_FACE_ARRAYS_SCALAR& face_velocities_new,T_FACE_ARRAYS_SCALAR& face_velocities_old,T_FACE_ARRAYS_SCALAR& potential_energy,const T dt);
    void Update_Potential_Energy(T_FACE_ARRAYS_SCALAR& face_velocities_new,T_FACE_ARRAYS_SCALAR& face_velocities_old,T_ARRAYS_SCALAR& potential_energy,const T dt);
    template<class T_POTENTIAL> void Update_Potential_Energy(T_FACE_ARRAYS_SCALAR& face_velocities_old,T_FACE_ARRAYS_SCALAR& face_velocities_new,T_POTENTIAL& potential_energy,T& allowed_energy_gained,const T dt);
    VECTOR<T,2> Compare_Potential_Energy(T_FACE_ARRAYS_SCALAR& face_velocities_old,T_FACE_ARRAYS_SCALAR& face_velocities_new,const T dt);
    void Compute_Divergence(const T_FACE_LOOKUP_FIRE_MULTIPHASE& face_lookup,LAPLACE_UNIFORM<T_GRID>* solver);
    void Compute_Divergence_For_Energy_Correction(const T_FACE_ARRAYS_SCALAR& face_velocities);
    void Apply_Pressure(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,bool scale_by_dt=false);
    void Set_Up_For_SPH(T_FACE_ARRAYS_SCALAR& face_velocities,const bool use_variable_density_solve=false,const bool use_one_way_coupling=false);
    void Restore_After_SPH(T_FACE_ARRAYS_SCALAR& face_velocities,const bool use_variable_density_solve=false,const bool use_one_way_coupling=false);
    void Update_Phi_And_Move_Velocity_Discontinuity(T_FACE_ARRAYS_SCALAR& face_velocities,LEVELSET_MULTIPLE_UNIFORM<T_GRID>& levelset_multiple,const T time,const bool update_phi_only=false);
    template<class FACE_LOOKUP> void Compute_Divergence(const FACE_LOOKUP &face_lookup,LAPLACE_UNIFORM<T_GRID>* solver);
    T Flame_Speed_Face_Multiphase(const int axis,const TV_INT& face_index,const int fuel_region,const int product_region) const;
    void Use_Flame_Speed_Multiplier(const bool use_flame_speed_multiplier_input=true);
    T Face_Jump_Multiphase(const int axis,const TV_INT& face_index,const int current_region,const int face_region) const;
//#####################################################################
};
}
#endif
