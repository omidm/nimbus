#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_EVOLUTION_DYADIC
//##################################################################### 
#ifndef __PARTICLE_LEVELSET_EVOLUTION_DYADIC__
#define __PARTICLE_LEVELSET_EVOLUTION_DYADIC__

#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_DYADIC.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_POLICY.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_DYADIC.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION.h>
namespace PhysBAM{

template<class T_ARRAY> class RUNGEKUTTA;

template<class T_GRID>
class PARTICLE_LEVELSET_EVOLUTION_DYADIC:public PARTICLE_LEVELSET_EVOLUTION<typename T_GRID::SCALAR>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::MAP_MESH MAP_MESH;
    typedef typename T_GRID::CELL CELL;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FLOOD_FILL FLOOD_FILL;
    typedef typename T_GRID::UNIFORM_GRID UNIFORM_GRID;typedef typename GRID_ARRAYS_POLICY<UNIFORM_GRID>::ARRAYS_SCALAR UNIFORM_ARRAYS;
    typedef typename UNIFORM_GRID::CELL_ITERATOR UNIFORM_CELL_ITERATOR;typedef typename UNIFORM_GRID::NODE_ITERATOR UNIFORM_NODE_ITERATOR;
    typedef typename LEVELSET_POLICY<T_GRID>::PARTICLE_LEVELSET T_PARTICLE_LEVELSET;
    typedef typename LEVELSET_ADVECTION_POLICY<T_GRID>::LEVELSET_ADVECTION T_LEVELSET_ADVECTION;
public:
    typedef PARTICLE_LEVELSET_EVOLUTION<T> BASE;
    using BASE::track_mass;using BASE::initial_mass;using BASE::runge_kutta_order_levelset;using BASE::runge_kutta_order_particles;using BASE::use_particle_levelset;
    using BASE::use_frozen_velocity;using BASE::cfl_number;using BASE::reseeding_frequency;using BASE::time;using BASE::use_fmm;using BASE::use_reinitialization;

    T_GRID& grid;
    ARRAY<T> phi;
    ARRAY<T> face_velocities,face_velocities2;
    T_PARTICLE_LEVELSET particle_levelset;
    RUNGEKUTTA<ARRAY<T> >* rungekutta_phi;

    T_LEVELSET_ADVECTION levelset_advection;

    PARTICLE_LEVELSET_EVOLUTION_DYADIC(T_GRID& grid_input);
    virtual ~PARTICLE_LEVELSET_EVOLUTION_DYADIC();

    void Use_Semi_Lagrangian_Advection()
    {levelset_advection.Use_Semi_Lagrangian_Advection();}

    void Track_Mass(const bool track_mass_input=true)
    {track_mass=track_mass_input;if(track_mass){initial_mass=particle_levelset.levelset.Approximate_Negative_Material();LOG::cout << "negative material = " << initial_mass << std::endl;}}

    void Initialize_Domain()
    {phi.Resize(grid.number_of_cells);face_velocities.Resize(grid.number_of_faces);
    if(runge_kutta_order_particles == 2 && !use_frozen_velocity) face_velocities2.Resize(grid.number_of_faces);
    particle_levelset.Initialize_Particle_Levelset_Grid_Values();}

    void Make_Signed_Distance()
    {assert(use_fmm);particle_levelset.levelset.Fast_Marching_Method();} // only FMM is supported

//#####################################################################
    void Initialize_Runge_Kutta();
    void Advance_To_Time(const T stopping_time,const bool verbose=true);
    T Time_Step(const T stopping_time,bool& limited_by_stopping_time);
    T CFL(const bool need_to_get_velocity=true);
    void Advance_One_Time_Step(const T dt);
    void Advance_Levelset(const T dt);
    void Advance_Particles(const T dt);
    T Advance_Particles(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T input_time);
    void Modify_Levelset_And_Particles();
    void Reseed_Particles(const T time,const int time_step=0,ARRAY<bool>* cell_centered_mask=0,const bool verbose=true);
    void Fill_Levelset_Ghost_Cells(const T time);
//#####################################################################
};
}
#endif
#endif
