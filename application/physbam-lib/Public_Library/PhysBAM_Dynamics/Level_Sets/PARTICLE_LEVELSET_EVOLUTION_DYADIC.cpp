#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic/BLOCK_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Dynamics/Level_Sets/FAST_LEVELSET_ADVECTION.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_DYADIC.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_DYADIC.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> PARTICLE_LEVELSET_EVOLUTION_DYADIC<T_GRID>::
PARTICLE_LEVELSET_EVOLUTION_DYADIC(T_GRID& grid_input)
    :grid(grid_input),particle_levelset(grid_input,phi),rungekutta_phi(0),levelset_advection(&particle_levelset.levelset)
{
    Use_Semi_Lagrangian_Advection();
    Track_Mass(false);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> PARTICLE_LEVELSET_EVOLUTION_DYADIC<T_GRID>::
~PARTICLE_LEVELSET_EVOLUTION_DYADIC()
{
    delete rungekutta_phi;
}
//#####################################################################
// Function Initialize_Runge_Kutta
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_DYADIC<T_GRID>::
Initialize_Runge_Kutta()
{
    if(runge_kutta_order_levelset > 1){
        delete rungekutta_phi;rungekutta_phi=new RUNGEKUTTA<ARRAY<T> >(phi);
        rungekutta_phi->Set_Order(runge_kutta_order_levelset);rungekutta_phi->Set_Time(time);
        rungekutta_phi->Set_Grid_And_Boundary_Condition(grid,*particle_levelset.levelset.boundary);}
}
//#####################################################################
// Function Time_Step
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR PARTICLE_LEVELSET_EVOLUTION_DYADIC<T_GRID>::
Time_Step(const T stopping_time,bool& limited_by_stopping_time)
{
    T dt=CFL();limited_by_stopping_time=false;
    if(time+dt >= stopping_time){dt=stopping_time-time;limited_by_stopping_time=true;}else if(time+2*dt >= stopping_time) dt=(T).51*(stopping_time-time);
    return dt;
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR PARTICLE_LEVELSET_EVOLUTION_DYADIC<T_GRID>::
CFL(bool need_to_get_velocity)
{
    if(need_to_get_velocity) particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,face_velocities,time);
    return cfl_number*particle_levelset.levelset.CFL(face_velocities);
}
//#####################################################################
// Function Advance_To_Time
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_DYADIC<T_GRID>::
Advance_To_Time(const T stopping_time,const bool verbose)
{
    int substep=0;bool done=false;
    while(!done){substep++;
    T dt=Time_Step(stopping_time,done);
    if(verbose) LOG::cout << "substep = " << substep << ", dt = " << dt << std::endl;
    Advance_One_Time_Step(dt);}
    if(track_mass){
        T mass=particle_levelset.levelset.Approximate_Negative_Material();
        LOG::cout << "negative area = " << mass << " - change = " << (mass-initial_mass)/initial_mass*100 << "%" << std::endl;}
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_DYADIC<T_GRID>::
Advance_One_Time_Step(const T dt)
{
    LOG::Time("advancing levelset");
    Advance_Levelset(dt);
    LOG::Time("advancing particles");
    Advance_Particles(dt);
    Modify_Levelset_And_Particles();
}
//#####################################################################
// Function Advance_Levelset
//#####################################################################
// only does an Euler step for the level set
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_DYADIC<T_GRID>::
Advance_Levelset(const T dt)
{
    if(runge_kutta_order_levelset == 1){
        particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,face_velocities,time);
        levelset_advection.Euler_Step(face_velocities,dt,time);time+=dt;}
    else{
        rungekutta_phi->Start(dt);
        for(int k=1;k<=runge_kutta_order_levelset;k++){
            if(k == 1 || !use_frozen_velocity) particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,face_velocities,time);
            levelset_advection.Euler_Step(face_velocities,dt,time);
            time=rungekutta_phi->Main();}}
}
//#####################################################################
// Function Advance_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_DYADIC<T_GRID>::
Advance_Particles(const T dt)
{
    if(use_particle_levelset){
        time-=dt; // to fix up time advancement due to Advance_Levelset()
        if(runge_kutta_order_particles == 1 || (runge_kutta_order_particles == 2 && use_frozen_velocity)){
            particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,face_velocities,time);
            particle_levelset.Euler_Step_Particles(face_velocities,dt,time,runge_kutta_order_particles==2);time+=dt;}
        else if(runge_kutta_order_particles == 2 || runge_kutta_order_particles == 3){
            T start_time=time;
            time=Advance_Particles(particle_levelset.positive_particles,PARTICLE_LEVELSET_POSITIVE,dt,start_time);
            time=Advance_Particles(particle_levelset.negative_particles,PARTICLE_LEVELSET_NEGATIVE,dt,start_time);
            particle_levelset.Euler_Step_Removed_Particles(dt,start_time,true);}} // can only Euler Step removed particles
}
//#####################################################################
// Function Advance_Particles
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR PARTICLE_LEVELSET_EVOLUTION_DYADIC<T_GRID>::
Advance_Particles(ARRAY<PARTICLE_LEVELSET_PARTICLES<TV>*>& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T input_time)
{
    T time=input_time;
    ARRAY<RUNGEKUTTA<ARRAY_VIEW<TV> >*> rungekutta_particles;rungekutta_particles.Resize(grid.number_of_cells);
    for(int i=1;i<=grid.number_of_cells;i++) if(particles(i)) rungekutta_particles(i)=RUNGEKUTTA<ARRAY_VIEW<TV> >::Create(particles(i)->X,runge_kutta_order_particles,time,dt);
    for(int k=1;k<=runge_kutta_order_particles;k++){
        if(k == 1 || !use_frozen_velocity) particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,face_velocities,time);
        particle_levelset.Euler_Step_Particles(face_velocities,particles,particle_type,dt,time,false,k==1);
        for(int i=1;i<=grid.number_of_cells;i++) if(particles(i)) time=rungekutta_particles(i)->Main();}
    particle_levelset.Update_Particle_Cells(particles);
    rungekutta_particles.Delete_Pointers_And_Clean_Memory();
    return time;
}
//#####################################################################
// Function Modify_Levelset_And_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_DYADIC<T_GRID>::
Modify_Levelset_And_Particles()
{
    if(use_particle_levelset){
        LOG::Time("modifying levelset");
        particle_levelset.Modify_Levelset_Using_Escaped_Particles();}
    LOG::Time("reinitializing levelset");
    Make_Signed_Distance();
    if(use_particle_levelset){
        LOG::Time("modifying levelset");        
        particle_levelset.Modify_Levelset_Using_Escaped_Particles();
        LOG::Time("compacting particles");
        particle_levelset.Compact_Particles_Into_Single_Particle_Bin();
        LOG::Time("adjusting particle radii");
        particle_levelset.Adjust_Particle_Radii();
        LOG::Stop_Time();}
}
//#####################################################################
// Function Reseed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_DYADIC<T_GRID>::
Reseed_Particles(const T time,const int time_step,ARRAY<bool>* cell_centered_mask,const bool verbose)
{
    if(use_particle_levelset && (!time_step || (reseeding_frequency && time_step%reseeding_frequency == 0))){
        int new_particles=particle_levelset.Reseed_Particles(time,cell_centered_mask);
        if(verbose) LOG::cout << "Reseeding... " << new_particles << " new particles" << std::endl;
        Initialize_Runge_Kutta();} // need to reset based on new number of particles
}
//#####################################################################
// Function Fill_Levelset_Ghost_Cells
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_DYADIC<T_GRID>::
Fill_Levelset_Ghost_Cells(const T time)
{
    particle_levelset.levelset.boundary->Fill_Ghost_Cells_Cell(grid,phi,particle_levelset.levelset.phi,time);
}
//#####################################################################
template class PARTICLE_LEVELSET_EVOLUTION_DYADIC<OCTREE_GRID<float> >;
template class PARTICLE_LEVELSET_EVOLUTION_DYADIC<QUADTREE_GRID<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PARTICLE_LEVELSET_EVOLUTION_DYADIC<OCTREE_GRID<double> >;
template class PARTICLE_LEVELSET_EVOLUTION_DYADIC<QUADTREE_GRID<double> >;
#endif
#endif
