//#####################################################################
// Copyright 2002-2009, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Avi Robinson-Mosher, Tamar Shinar, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_QUADTREE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/VOF_ADVECTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
PARTICLE_LEVELSET_EVOLUTION_UNIFORM(const T_GRID& grid_input,const int number_of_ghost_cells_input)
    :grid(grid_input),particle_levelset(grid,phi,number_of_ghost_cells_input),rungekutta_phi(0),levelset_advection(&particle_levelset.levelset)
{
    Use_Semi_Lagrangian_Advection();
    Track_Mass(false);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
~PARTICLE_LEVELSET_EVOLUTION_UNIFORM()
{
    delete rungekutta_phi;
}
//#####################################################################
// Function Initialize_Runge_Kutta
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Initialize_Runge_Kutta()
{
    if(runge_kutta_order_levelset > 1){
        delete rungekutta_phi;rungekutta_phi=new RUNGEKUTTA<T_ARRAYS_SCALAR>(phi);
        rungekutta_phi->Set_Order(runge_kutta_order_levelset);rungekutta_phi->Set_Time(time);
        rungekutta_phi->Set_Grid_And_Boundary_Condition(grid,*particle_levelset.levelset.boundary);}
}
//#####################################################################
// Function Time_Step
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Time_Step(const T stopping_time,bool& limited_by_stopping_time)
{
    T dt=CFL();limited_by_stopping_time=false;
    if(time+dt >= stopping_time){dt=stopping_time-time;limited_by_stopping_time=true;}else if(time+2*dt >= stopping_time) dt=(T).51*(stopping_time-time);
    return dt;
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
CFL(const bool need_to_get_velocity,const bool analytic_test)
{
    if(need_to_get_velocity) particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,V,time);
    if(analytic_test) return cfl_number/max((V.Maxabs()/grid.dX).Max(),1/particle_levelset.levelset.max_time_step);
    return cfl_number*particle_levelset.levelset.CFL(V);
}
//#####################################################################
// Function Advance_To_Time
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Advance_To_Time(T_FACE_ARRAYS_SCALAR* face_velocities,const T stopping_time,const bool verbose)
{
    int substep=0;bool done=false;
    while(!done){substep++;
        T dt=Time_Step(stopping_time,done);
        if(verbose) {std::stringstream ss;ss << "substep = " << substep << ", dt = " << dt << std::endl;LOG::filecout(ss.str());}
        Advance_One_Time_Step(face_velocities,dt);}
    if(track_mass){
        T mass=levelset_advection.Approximate_Negative_Material();
        std::stringstream ss;
        ss << "negative material = " << mass << " - change = " << (mass-initial_mass)/initial_mass*100 << "%" << std::endl;
        LOG::filecout(ss.str());}
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Advance_One_Time_Step(T_FACE_ARRAYS_SCALAR* face_velocities,const T dt)
{
    Advance_Levelset(dt);
    Advance_Particles(*face_velocities,dt);
    Modify_Levelset_And_Particles(face_velocities);
}
//#####################################################################
// Function Advance_Levelset
//#####################################################################
// only does an Euler step for the level set
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Advance_Levelset(const T dt)
{
    if(runge_kutta_order_levelset == 1){
        particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,V,time);
        levelset_advection.Euler_Step(V,dt,time,particle_levelset.number_of_ghost_cells);time+=dt;}
    else{
        if(!rungekutta_phi) Initialize_Runge_Kutta();
        rungekutta_phi->Start(dt);
        for(int k=1;k<=runge_kutta_order_levelset;k++){
            if(k == 1 || !use_frozen_velocity) particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,V,time);
            levelset_advection.Euler_Step(V,dt,time,particle_levelset.number_of_ghost_cells);
            time=rungekutta_phi->Main();}}
}
//#####################################################################
// Function Advance_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Advance_Particles(const T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const bool analytic_test)
{
    if(use_particle_levelset){
        time-=dt; // to fix up time advancement due to Advance_Levelset()
        if(runge_kutta_order_particles == 1 || (runge_kutta_order_particles == 2 && use_frozen_velocity)){
            // TODO: still needed?
            //particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,V,time);
            particle_levelset.Euler_Step_Particles(face_velocities,dt,time,runge_kutta_order_particles==2,true,analytic_test);time+=dt;}
        else if(runge_kutta_order_particles == 2 || runge_kutta_order_particles == 3){
            T start_time=time;
            time=Advance_Particles(particle_levelset.positive_particles,PARTICLE_LEVELSET_POSITIVE,dt,start_time);
            time=Advance_Particles(particle_levelset.negative_particles,PARTICLE_LEVELSET_NEGATIVE,dt,start_time);
            if(analytic_test){
                time=Advance_Particles(particle_levelset.removed_positive_particles,PARTICLE_LEVELSET_POSITIVE,dt,start_time);
                time=Advance_Particles(particle_levelset.removed_negative_particles,PARTICLE_LEVELSET_NEGATIVE,dt,start_time);}
            else particle_levelset.Euler_Step_Removed_Particles(dt,start_time,true);}} // can only Euler Step removed particles
}
//#####################################################################
// Function Advance_Particles
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Advance_Particles(T_ARRAYS_PARTICLE_LEVELSET_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T input_time)
{
    T_GRID mac_grid=grid.Get_MAC_Grid();
    T current_time=input_time;
    T_ARRAYS_RUNGEKUTTA rungekutta_particles;rungekutta_particles.Resize(mac_grid.Domain_Indices(1));
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index())){
        PARTICLE_LEVELSET_PARTICLES<TV>* particles_list=particles(iterator.Cell_Index());
        while(particles_list){
            rungekutta_particles(iterator.Cell_Index()).Append(RUNGEKUTTA<ARRAY_VIEW<TV> >::Create(particles_list->X,runge_kutta_order_particles,dt,current_time));
            rungekutta_particles(iterator.Cell_Index()).Last()->Start(dt);
            particles_list=particles_list->next;}}
    for(int k=1;k<=runge_kutta_order_particles;k++){
        if(k == 1 || !use_frozen_velocity) particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,V,current_time);
        particle_levelset.Euler_Step_Particles_Wrapper(V,particles,particle_type,dt,current_time,false,k==1,false);
        for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index()))
            for(int i=1;i<=rungekutta_particles(iterator.Cell_Index()).m;i++) current_time=rungekutta_particles(iterator.Cell_Index())(i)->Main();}
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index())){
        PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(iterator.Cell_Index());
        for(int k=1;k<=cell_particles.array_collection->Size();k++){
            TV velocity;
            particle_levelset.levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,k,velocity,particle_type,(T)1,current_time);
            cell_particles.X(k)+=velocity;}}
    particle_levelset.Update_Particle_Cells(particles);
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) rungekutta_particles(iterator.Cell_Index()).Delete_Pointers_And_Clean_Memory();
    return current_time;
}
//#####################################################################
// Function Advance_Particles
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Advance_Particles(T_ARRAYS_PARTICLE_LEVELSET_REMOVED_PARTICLES& particles,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T input_time)
{
    T_GRID mac_grid=grid.Get_MAC_Grid();
    T current_time=input_time;
    T_ARRAYS_RUNGEKUTTA rungekutta_particles;rungekutta_particles.Resize(mac_grid.Domain_Indices(1));
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index())){
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles_list=particles(iterator.Cell_Index());
        while(particles_list){
            rungekutta_particles(iterator.Cell_Index()).Append(RUNGEKUTTA<ARRAY_VIEW<TV> >::Create(particles(iterator.Cell_Index())->X,runge_kutta_order_particles,dt,current_time));
            rungekutta_particles(iterator.Cell_Index()).Last()->Start(dt);
            particles_list=particles_list->next;}}
    for(int k=1;k<=runge_kutta_order_particles;k++){
        if(k == 1 || !use_frozen_velocity) particle_levelset.levelset.levelset_callbacks->Get_Levelset_Velocity(grid,particle_levelset.levelset,V,current_time);
        particle_levelset.Euler_Step_Particles_Wrapper(V,particles,particle_type,dt,current_time,false,k==1,false);
        for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index())) 
            for(int i=1;i<=rungekutta_particles(iterator.Cell_Index()).m;i++) current_time=rungekutta_particles(iterator.Cell_Index())(i)->Main();}
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) if(particles(iterator.Cell_Index())){
        PARTICLE_LEVELSET_PARTICLES<TV>& cell_particles=*particles(iterator.Cell_Index());
        for(int k=1;k<=cell_particles.array_collection->Size();k++){
            TV velocity;
            particle_levelset.levelset.levelset_callbacks->Adjust_Particle_For_Domain_Boundaries(cell_particles,k,velocity,particle_type,(T)1,current_time);
            cell_particles.X(k)+=velocity;}}
    particle_levelset.Update_Particle_Cells(particles);
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) rungekutta_particles(iterator.Cell_Index()).Delete_Pointers_And_Clean_Memory();
    return input_time+dt; // there may be no removed particles
}
//#####################################################################
// Note: This function is added for correct execution of 3d water simulation
// with Nimbus. There are no MPI calls hidden in this function. However, this
// is not expected to give correct results with other PhysBAM simulations. The
// function is not tested for any case apart from 3d water simulation with
// Nimbus. Calls unrequired by the simple 3d water simulation are deleted for
// convenience.
// -- Chinmayee
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Modify_Levelset_And_Particles_Nimbus(T_FACE_ARRAYS_SCALAR* face_velocities,
                                     T_ARRAYS_SCALAR* phi_ghost)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Modify_Levelset_And_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Modify_Levelset_And_Particles(T_FACE_ARRAYS_SCALAR* face_velocities)
{
    // TODO: a call for creating particles from the geometry if necessary
    if(use_particle_levelset){
        particle_levelset.Modify_Levelset_Using_Escaped_Particles(face_velocities);
        if(particle_levelset.vof_advection){
            particle_levelset.vof_advection->Rasterize_Material_Postimages();
            particle_levelset.vof_advection->Adjust_Levelset_With_Material_Volumes();}}
    Make_Signed_Distance();
    if(use_particle_levelset){
        particle_levelset.Modify_Levelset_Using_Escaped_Particles(face_velocities);
        if(particle_levelset.vof_advection) particle_levelset.vof_advection->Adjust_Levelset_With_Material_Volumes();
        particle_levelset.Adjust_Particle_Radii();}
}
//#####################################################################
// Function Modify_Levelset_And_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Modify_Levelset_And_Particles(T_FACE_ARRAYS_SCALAR* face_velocities,const T stopping_distance)
{
    // TODO: a call for creating particles from the geometry if necessary
    if(use_particle_levelset){
        particle_levelset.Modify_Levelset_Using_Escaped_Particles(face_velocities);
        if(particle_levelset.vof_advection){
            particle_levelset.vof_advection->Rasterize_Material_Postimages();
            particle_levelset.vof_advection->Adjust_Levelset_With_Material_Volumes();}}
    Make_Signed_Distance(stopping_distance);
    if(use_particle_levelset){
        particle_levelset.Modify_Levelset_Using_Escaped_Particles(face_velocities);
        if(particle_levelset.vof_advection) particle_levelset.vof_advection->Adjust_Levelset_With_Material_Volumes();
        particle_levelset.Adjust_Particle_Radii();}
}
//#####################################################################
// Function Reseed_Particles
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Reseed_Particles(const T time,const int time_step,T_ARRAYS_BOOL* cell_centered_mask,const bool verbose)
{
    if((use_particle_levelset && cell_centered_mask) || (!time_step || (reseeding_frequency && time_step%reseeding_frequency == 0))){
        int new_particles=particle_levelset.Reseed_Particles(time,cell_centered_mask);
        if(verbose) {std::stringstream ss;ss << "Reseeding... " << new_particles << " new particles" << std::endl;LOG::filecout(ss.str());}
        Initialize_Runge_Kutta();} // need to reset based on new number of particles
}
//#####################################################################
// Function Apply_Mass_Conservation
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Apply_Mass_Conservation(const int number_of_regions,const T time,const T dt,T_FACE_ARRAYS_SCALAR& face_velocities)
{
    // TODO: determine conditions for when to perform this correction
    for(int i=1;i<=number_of_regions;i++){
        Particle_Levelset(i).vof_advection->Make_Approximately_Incompressible(face_velocities,dt,time);
        Particle_Levelset(i).vof_advection->Refine_Or_Coarsen_Geometry();}
}
//#####################################################################
// Function Reinitialize_Geometry
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Reinitialize_Geometry(const int number_of_regions)
{
    // precondition: rasterize must have been called
    for(int i=1;i<=number_of_regions;i++) Particle_Levelset(i).vof_advection->Create_Geometry();
}
//#####################################################################
// Function Perform_Conservative_Advection
//#####################################################################
template<class T_GRID> void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID>::
Perform_Conservative_Advection(const int number_of_regions,const T time,const T dt,T_FACE_ARRAYS_SCALAR& face_velocities)
{
    typename INTERPOLATION_POLICY<GRID<TV> >::FACE_LOOKUP V_lookup(face_velocities);
    for(int i=1;i<=number_of_regions;i++) Particle_Levelset(i).vof_advection->Perform_Conservative_Advection(V_lookup,dt,time);
}
//#####################################################################
template class PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<VECTOR<float,1> > >;
template class PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<VECTOR<float,2> > >;
template class PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<VECTOR<double,1> > >;
template class PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<VECTOR<double,2> > >;
template class PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<VECTOR<double,3> > >;
#endif

//#####################################################################
// Specialized implementation for 3d water simulation as mentioned before.
// -- Chinmayee
//#####################################################################
template <>
void PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<VECTOR<float, 3> > >::
Modify_Levelset_And_Particles_Nimbus(T_FACE_ARRAYS_SCALAR* face_velocities,
                                     T_ARRAYS_SCALAR* phi_ghost)
{
    std::cout << "#### CHINMAYEE: In specialized implementation for "
              << "modify levelset for 3d water simulation\n";
    typedef float T;
    typedef VECTOR<T, 3> TV;
    typedef VECTOR<int, 3> TV_INT;
    typedef GRID<VECTOR<float, 3> > T_GRID;
    if(use_particle_levelset) {
        particle_levelset.Modify_Levelset_Using_Escaped_Particles(face_velocities);
    }
    {
        // Make_Signed_Distance()
        if(use_fmm) {
            {
                // T_FAST_LEVELSET::Fast_Marching_Method()
                T_FAST_LEVELSET* ls = &particle_levelset.levelset;
                const int local_advection_spatial_order =
                    levelset_advection.local_advection_spatial_order;
                T time = 0;
                {
                    // T_FAST_LEVELSET::BASE::Fast_Marching_Method()
                    // That is, LEVELSET_3D::Fast_Marching_Method()
                    T stopping_distance = ls->half_band_width +
                                          ls->grid.dX.Max() *
                                          (1 + min(3, local_advection_spatial_order));
                    {
                        // LEVELSET_3D::Get_Signed_Distance_Using_FMM
                        const ARRAY<TV_INT>* seed_indices = NULL;
                        const bool add_seed_indices_for_ghost_cells = false;
                        const int ghost_cells = 7;
                        T_ARRAYS_SCALAR pg(ls->grid.Domain_Indices(ghost_cells));
                        ls->boundary->Fill_Ghost_Cells(ls->grid,
                                                       ls->phi,
                                                       pg,
                                                       0,
                                                       time,
                                                       ghost_cells);
                        FAST_MARCHING_METHOD_UNIFORM<T_GRID> fmm(*ls,
                                                                 ghost_cells,
                                                                 ls->thread_queue);
                        fmm.Fast_Marching_Method(pg,
                                                 stopping_distance,
                                                 seed_indices,
                                                 add_seed_indices_for_ghost_cells);
                        ARRAY<T, TV_INT>::Get(ls->phi, pg);
                        ls->boundary->Apply_Boundary_Condition(ls->grid, ls->phi, time);
                    }
                }
                // Unrequired??
                ls->boundary->Apply_Boundary_Condition(ls->grid, ls->phi, time);
            }
        }
        else if(use_reinitialization)
            levelset_advection.Reinitialize();
    }
    if(use_particle_levelset) {
        particle_levelset.Modify_Levelset_Using_Escaped_Particles(face_velocities);
        particle_levelset.Adjust_Particle_Radii();
    }
}
