//#####################################################################
// Copyright 2005, Ron Fedkiw, Eran Guendelman, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_MULTIPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM(const T_GRID& grid_input,const int number_of_ghost_cells_input)
    :PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>(grid_input,number_of_ghost_cells_input),particle_levelset_multiple(grid,phis,number_of_ghost_cells_input),levelset_advection_multiple(&particle_levelset_multiple.levelset_multiple)
{
    Use_Semi_Lagrangian_Advection();
    Track_Mass(false);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
~PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM()
{
    for(int i=1;i<=rungekutta_phis.m;i++)delete rungekutta_phis(i);
}
//#####################################################################
// Function Initialize_Runge_Kutta
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Initialize_Runge_Kutta()
{
    if(runge_kutta_order_levelset > 1) for(int i=1;i<=phis.m;i++){
        delete rungekutta_phis(i);rungekutta_phis(i)=new RUNGEKUTTA<T_ARRAYS_SCALAR>(phis(i));
        rungekutta_phis(i)->Set_Order(runge_kutta_order_levelset);rungekutta_phis(i)->Set_Time(time);
        rungekutta_phis(i)->Set_Grid_And_Boundary_Condition(grid,*particle_levelset_multiple.particle_levelsets(i)->levelset.boundary);}
}
//#####################################################################
// Function Time_Step
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Time_Step(const T stopping_time,bool& limited_by_stopping_time)
{
    T dt=CFL();limited_by_stopping_time=false;
    if(time+dt >= stopping_time){dt=stopping_time-time;limited_by_stopping_time=true;}else if(time+2*dt >= stopping_time) dt=(T).51*(stopping_time-time);
    return dt;
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
CFL(const bool need_to_get_velocity,const bool analytic_test)
{
    if(need_to_get_velocity) particle_levelset_multiple.levelset_multiple.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset_multiple.levelset_multiple,V,time);
    if(analytic_test){T max_time_step=particle_levelset_multiple.levelset_multiple.levelsets(1)->max_time_step;
        for(int i=2;i<=particle_levelset_multiple.levelset_multiple.levelsets.m;i++) max_time_step=min(max_time_step,particle_levelset_multiple.levelset_multiple.levelsets(i)->max_time_step);
        return cfl_number/max(V.Maxabs().Max(),1/max_time_step);}
    return cfl_number*particle_levelset_multiple.levelset_multiple.CFL(V);
}
//#####################################################################
// Function Advance_To_Time
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Advance_To_Time(T_FACE_ARRAYS_SCALAR* face_velocities,const T stopping_time,const bool verbose)
{
    int substep=0;bool done=false;
    while(!done){substep++;
        T dt=Time_Step(stopping_time,done);
        if(verbose) {std::stringstream ss;ss << "substep = " << substep << ", dt = " << dt << std::endl;LOG::filecout(ss.str());}
        Advance_One_Time_Step(face_velocities,dt);}
    if(track_mass)for(int i=1;i<=particle_levelset_multiple.particle_levelsets.m;i++){
        T mass=Levelset_Advection(i).Approximate_Negative_Material();
        std::stringstream ss;ss << "negative material(" << i << ") = " << mass << " - change = " << (mass-initial_mass(i))/initial_mass(i)*100 << "%" << std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Advance_One_Time_Step(T_FACE_ARRAYS_SCALAR* face_velocities,const T dt)
{
    LOG::Time("advancing levelset");
    Advance_Levelset(dt);
    LOG::Time("advancing particles");
    Advance_Particles(*face_velocities,dt);
    Modify_Levelset_And_Particles(face_velocities);
}
//#####################################################################
// Function Advance_Levelset
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Advance_Levelset(const T dt)
{
    if(runge_kutta_order_levelset == 1){
        particle_levelset_multiple.levelset_multiple.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset_multiple.levelset_multiple,V,time);
        levelset_advection_multiple.Euler_Step(V,dt,time,particle_levelset_multiple.number_of_ghost_cells);time+=dt;}
    else{
        for(int i=1;i<=rungekutta_phis.m;i++)rungekutta_phis(i)->Start(dt);
        for(int k=1;k<=runge_kutta_order_levelset;k++){
            if(k == 1 || !use_frozen_velocity) particle_levelset_multiple.levelset_multiple.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset_multiple.levelset_multiple,V,time);
            levelset_advection_multiple.Euler_Step(V,dt,time,particle_levelset_multiple.number_of_ghost_cells);
            for(int i=1;i<=rungekutta_phis.m;i++)time=rungekutta_phis(i)->Main();}}
}
//#####################################################################
// Function Advance_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Advance_Particles(const T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const bool analytic_test)
{
    if(analytic_test) PHYSBAM_NOT_IMPLEMENTED("analytic_test");
    if(use_particle_levelset){
        time-=dt; // to fix up time advancement due to Advance_Levelset()
        if(runge_kutta_order_particles == 1 || (runge_kutta_order_particles == 2 && use_frozen_velocity)){
            // TODO: still needed?
            //particle_levelset_multiple.levelset_multiple.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset_multiple.levelset_multiple,V,time);
            particle_levelset_multiple.Euler_Step_Particles(face_velocities,dt,time,runge_kutta_order_particles==2);time+=dt;}
        else if(runge_kutta_order_particles == 2 || runge_kutta_order_particles == 3){
            T start_time=time;
            for(int i=1;i<=particle_levelset_multiple.particle_levelsets.m;i++)
                time=Advance_Particles(particle_levelset_multiple.particle_levelsets(i)->negative_particles,PARTICLE_LEVELSET_NEGATIVE,dt,start_time);
            particle_levelset_multiple.Euler_Step_Removed_Particles(dt,start_time,true);}} // can only Euler Step removed particles
}
//#####################################################################
// Function Advance_Particles
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Advance_Particles(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T input_time)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Advance_Particles
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Advance_Particles(T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T input_time)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Modify_Levelset_And_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Modify_Levelset_And_Particles(T_FACE_ARRAYS_SCALAR* face_velocities)
{
    if(use_particle_levelset){
        LOG::Time("modifying levelset");
        particle_levelset_multiple.Modify_Levelset_Using_Escaped_Particles(face_velocities);}
    LOG::Time("reinitializing and projecting levelset");
    particle_levelset_multiple.levelset_multiple.Project_Levelset();
    Make_Signed_Distance();
    if(use_particle_levelset){
        LOG::Time("modifying levelset");
        particle_levelset_multiple.Modify_Levelset_Using_Escaped_Particles(face_velocities);
        LOG::Time("adjusting particle radii");
        particle_levelset_multiple.Adjust_Particle_Radii();
        LOG::Stop_Time();}
    particle_levelset_multiple.levelset_multiple.Project_Levelset();
}
//#####################################################################
// Function Modify_Levelset_And_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Modify_Levelset_And_Particles(T_FACE_ARRAYS_SCALAR* face_velocities,const T stopping_distance)
{
    if(use_particle_levelset){
        LOG::Time("modifying levelset");
        particle_levelset_multiple.Modify_Levelset_Using_Escaped_Particles(face_velocities);}
    LOG::Time("reinitializing and projecting levelset");
    particle_levelset_multiple.levelset_multiple.Project_Levelset();
    Make_Signed_Distance(stopping_distance);
    if(use_particle_levelset){
        LOG::Time("modifying levelset");
        particle_levelset_multiple.Modify_Levelset_Using_Escaped_Particles(face_velocities);
        LOG::Time("adjusting particle radii");
        particle_levelset_multiple.Adjust_Particle_Radii();
        LOG::Stop_Time();}
    particle_levelset_multiple.levelset_multiple.Project_Levelset();
}
//#####################################################################
// Function Reseed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Reseed_Particles(const T time,const int time_step,T_ARRAYS_BOOL* cell_centered_mask,const bool verbose)
{
    if((use_particle_levelset && cell_centered_mask) || (!time_step || (reseeding_frequency && time_step%reseeding_frequency == 0))){
        int new_particles=particle_levelset_multiple.Reseed_Particles(time,cell_centered_mask);
        if(verbose) {std::stringstream ss;ss << "Reseeding... " << new_particles << " new particles" << std::endl;LOG::filecout(ss.str());}
        Initialize_Runge_Kutta();} // need to reset based on new number of particles
}
//#####################################################################
// Function Fill_Levelset_Ghost_Cells
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Fill_Levelset_Ghost_Cells(const T time)
{
    for(int i=1;i<=phis.m;i++) particle_levelset_multiple.levelset_multiple.levelsets(i)->boundary->Fill_Ghost_Cells(grid,phis(i),phis(i),0,time,particle_levelset_multiple.number_of_ghost_cells);
}
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Use_Semi_Lagrangian_Advection()
{
    for(int i=1;i<=phis.m;i++) levelset_advection_multiple.levelset_advections(i).Use_Local_Semi_Lagrangian_Advection();
}
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Use_Hamilton_Jacobi_Weno_Advection()
{
    for(int i=1;i<=phis.m;i++) levelset_advection_multiple.levelset_advections(i).Use_Local_WENO_For_Advection();
}
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Use_Hamilton_Jacobi_Eno_Advection(const int order)
{
    assert(order >= 1 && order <= 3);
    for(int i=1;i<=phis.m;i++) levelset_advection_multiple.levelset_advections(i).Use_Local_ENO_For_Advection(order);
}
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Track_Mass(const bool track_mass_input)
{
    track_mass=track_mass_input;
    if(track_mass) for(int i=1;i<=phis.m;i++){
        initial_mass(i)=levelset_advection_multiple.levelset_advections(i).Approximate_Negative_Material();
        std::stringstream ss;ss << "negative material(" << i << ") = " << initial_mass(i) << std::endl;LOG::filecout(ss.str());}
}
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Initialize_Domain(const T_GRID& grid_input,const int number_of_regions)  // don't call up to the base class here because we don't need those variables initialized OVERRIDE PROBLEM
{
    assert(grid_input.Is_MAC_Grid());
    grid=grid_input;
    phis.Resize(number_of_regions);
    for(int i=1;i<=phis.m;i++) phis(i).Resize(grid.Domain_Indices(particle_levelset_multiple.number_of_ghost_cells));
    V.Resize(grid);
    particle_levelset_multiple.Initialize_Particle_Levelsets_And_Grid_Values(grid,phis,number_of_regions);
    //TODO: patkar - call this function somewhere in solids fluids driver after the advections are resized and set
    /*for(int i=1;i<=phis.m;i++) 
        if(levelset_advection_multiple.levelset_advections(i).semi_lagrangian_collidable)
            particle_levelset_multiple.particle_levelsets(i)->levelset.Initialize_Valid_Masks(grid);*/
    initial_mass.Resize(number_of_regions);
    rungekutta_phis.Resize(number_of_regions);
}
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Initialize_Domain(const T_GRID& grid_input)
{
    Initialize_Domain(grid_input,2);
}
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Make_Signed_Distance()
{
    if(use_fmm){
        ARRAY<int> local_advection_spatial_orders(phis.m);
        for(int i=1;i<=phis.m;i++)
            local_advection_spatial_orders(i)=levelset_advection_multiple.levelset_advections(i).local_advection_spatial_order;
        particle_levelset_multiple.levelset_multiple.Fast_Marching_Method(local_advection_spatial_orders);}
    else if(use_reinitialization) levelset_advection_multiple.Reinitialize();
}
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Make_Signed_Distance(const T stopping_distance)
{
    if(use_fmm){
        ARRAY<int> local_advection_spatial_orders(phis.m);
        for(int i=1;i<=phis.m;i++)
            local_advection_spatial_orders(i)=levelset_advection_multiple.levelset_advections(i).local_advection_spatial_order;
        particle_levelset_multiple.levelset_multiple.Fast_Marching_Method(stopping_distance,local_advection_spatial_orders);}
    else if(use_reinitialization) levelset_advection_multiple.Reinitialize();
}
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Set_Number_Particles_Per_Cell(const int number_particles_per_cell)
{
    for(int i=1;i<=phis.m;i++) particle_levelset_multiple.particle_levelsets(i)->Set_Number_Particles_Per_Cell(number_particles_per_cell);
}
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Set_Levelset_Callbacks(LEVELSET_CALLBACKS<T_GRID>& levelset_callbacks)
{
    particle_levelset_multiple.levelset_multiple.Set_Levelset_Callbacks(levelset_callbacks);
}
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Initialize_FMM_Initialization_Iterative_Solver(const bool refine_fmm_initialization_with_iterative_solver_input,const int fmm_initialization_iterations_input,
        const T fmm_initialization_iterative_tolerance_input,const T fmm_initialization_iterative_drift_fraction_input)
{
    for(int i=1;i<=phis.m;i++) particle_levelset_multiple.levelset_multiple.levelsets(i)->Initialize_FMM_Initialization_Iterative_Solver(refine_fmm_initialization_with_iterative_solver_input,
        fmm_initialization_iterations_input,fmm_initialization_iterative_tolerance_input,fmm_initialization_iterative_drift_fraction_input);
}
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Bias_Towards_Negative_Particles(const bool bias_towards_negative_particles)
{
    for(int i=1;i<=phis.m;i++) particle_levelset_multiple.particle_levelsets(i)->Bias_Towards_Negative_Particles(bias_towards_negative_particles);
}
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Set_Seed(const int seed)
{
    for(int i=1;i<=phis.m;i++) particle_levelset_multiple.particle_levelsets(i)->random.Set_Seed(seed);
}
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Seed_Particles(const T time)
{
    for(int i=1;i<=phis.m;i++) particle_levelset_multiple.particle_levelsets(i)->Seed_Particles(time);
}
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Delete_Particles_Outside_Grid()
{
    for(int i=1;i<=phis.m;i++) particle_levelset_multiple.particle_levelsets(i)->Delete_Particles_Outside_Grid();
}
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Set_CFL_Number(const T cfl_number_input)
{
    PARTICLE_LEVELSET_EVOLUTION<T>::Set_CFL_Number(cfl_number_input);
    for(int i=1;i<=phis.m;i++) particle_levelset_multiple.particle_levelsets(i)->cfl_number=cfl_number_input;
}
template<class T_GRID> PARTICLE_LEVELSET_MULTIPLE_UNIFORM<T_GRID>& PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Particle_Levelset_Multiple()
{
    return particle_levelset_multiple;
}
template<class T_GRID> LEVELSET_MULTIPLE_UNIFORM<T_GRID>& PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Levelset_Multiple()
{
    return particle_levelset_multiple.levelset_multiple;
}
template<class T_GRID> PARTICLE_LEVELSET_UNIFORM<T_GRID>& PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Particle_Levelset(const int i)
{
    return *particle_levelset_multiple.particle_levelsets(i);
}
template<class T_GRID> typename LEVELSET_POLICY<T_GRID>::FAST_LEVELSET_T& PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Levelset(const int i)
{
    return *particle_levelset_multiple.levelset_multiple.levelsets(i);
}
template<class T_GRID> typename LEVELSET_ADVECTION_POLICY<T_GRID>::FAST_LEVELSET_ADVECTION_T& PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<T_GRID>::
Levelset_Advection(const int i)
{
    return levelset_advection_multiple.levelset_advections(i);
}
//#####################################################################
template class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<GRID<VECTOR<float,1> > >;
template class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<GRID<VECTOR<float,2> > >;
template class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<GRID<VECTOR<double,1> > >;
template class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<GRID<VECTOR<double,2> > >;
template class PARTICLE_LEVELSET_EVOLUTION_MULTIPLE_UNIFORM<GRID<VECTOR<double,3> > >;
#endif
