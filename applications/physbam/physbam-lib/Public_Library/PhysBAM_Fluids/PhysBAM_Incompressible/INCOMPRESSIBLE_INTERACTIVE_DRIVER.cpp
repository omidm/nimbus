//#####################################################################
// Copyright 2009-2010, Mridul Aanjaneya, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_THREADED.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_INTERACTIVE_DRIVER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_EXAMPLE.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM_FORWARD.h>

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#ifdef WIN32
    #include <windows.h>
#endif
using namespace PhysBAM;
namespace{
template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
{
    ((INCOMPRESSIBLE_INTERACTIVE_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
}
};
//#####################################################################
// Initialize
//#####################################################################
template<class TV> INCOMPRESSIBLE_INTERACTIVE_DRIVER<TV>::
INCOMPRESSIBLE_INTERACTIVE_DRIVER(INCOMPRESSIBLE_EXAMPLE<TV>& example)
    :example(example),kinematic_evolution(example.rigid_geometry_collection,true),visualizer(0),rigid_component(0),mac_component(0),scalar_component(0)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> INCOMPRESSIBLE_INTERACTIVE_DRIVER<TV>::
~INCOMPRESSIBLE_INTERACTIVE_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void INCOMPRESSIBLE_INTERACTIVE_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void INCOMPRESSIBLE_INTERACTIVE_DRIVER<TV>::
Initialize()
{
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);

    // setup time
    if(example.restart) current_frame=example.restart;else current_frame=example.first_frame;
    output_number=current_frame;
    time=example.Time_At_Frame(current_frame);

    // initialize collision objects
    example.Initialize_Bodies();
    kinematic_evolution.Get_Current_Kinematic_Keyframes(1/example.frame_rate,time);
    kinematic_evolution.Set_External_Positions(example.rigid_geometry_collection.particles.X,example.rigid_geometry_collection.particles.rotation,time);
    kinematic_evolution.Set_External_Velocities(example.rigid_geometry_collection.particles.V,example.rigid_geometry_collection.particles.angular_velocity,time,time);

    // mpi
    if(example.mpi_grid) example.mpi_grid->Initialize(example.domain_boundary);
    ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>* advection_conservative=dynamic_cast<ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>*>(example.incompressible.advection);
    if(advection_conservative){
        for(int i=1;i<=TV::dimension;i++){advection_conservative->mpi_boundary(i)(1)=!example.domain_boundary(i)(1);advection_conservative->mpi_boundary(i)(2)=!example.domain_boundary(i)(2);}
        if(example.mpi_grid) example.mpi_grid->Initialize(advection_conservative->solid_walls);}
    example.incompressible.mpi_grid=example.mpi_grid;
    example.projection->elliptic_solver->mpi_grid=example.mpi_grid;
    if(example.mpi_grid) example.boundary=new BOUNDARY_MPI<GRID<TV> >(example.mpi_grid,example.boundary_scalar);
    else if(example.thread_queue) example.boundary=new BOUNDARY_THREADED<GRID<TV> >(*example.thread_queue,example.boundary_scalar);    
    else example.boundary=&example.boundary_scalar;
    example.incompressible.Set_Custom_Boundary(*example.boundary);
    
    // setup grids and velocities
    example.incompressible.Initialize_Grids(example.mac_grid);
    example.projection->Initialize_Grid(example.mac_grid);
    example.face_velocities.Resize(example.mac_grid);
    example.face_velocities_save.Resize(example.mac_grid);
    example.density.Resize(example.mac_grid.Domain_Indices());
    example.incompressible.kinetic_energy.Resize(example.mac_grid);
    example.Initialize_Fields();
    example.Initialize_Confinement();

    // setup laplace
    example.projection->elliptic_solver->Set_Relative_Tolerance((T)1e-8);
    example.projection->elliptic_solver->pcg.Set_Maximum_Iterations(40);
    example.projection->elliptic_solver->pcg.evolution_solver_type=krylov_solver_cg;
    example.projection->elliptic_solver->pcg.cg_restart_iterations=40;

    if(example.restart) example.Read_Output_Files(example.restart);
    
    // setup domain boundaries
    VECTOR<VECTOR<bool,2>,TV::dimension> constant_extrapolation;constant_extrapolation.Fill(VECTOR<bool,2>::Constant_Vector(true));
    example.incompressible.boundary->Set_Constant_Extrapolation(constant_extrapolation);
    example.boundary->Set_Constant_Extrapolation(constant_extrapolation);
    example.incompressible.Set_Custom_Boundary(*example.boundary);
    example.Set_Boundary_Conditions(time); // get so CFL is correct

    if(example.use_viscosity){
        PHYSBAM_FATAL_ERROR();
        example.incompressible.Use_Explicit_Part_Of_Implicit_Viscosity(false);    
        example.incompressible.Set_Viscosity(1);
        example.incompressible.Set_Variable_Viscosity(false);}

    example.Get_Scalar_Field_Sources(time);
    if(!example.restart) Write_Output_Files(example.first_frame);
    since_last_frame=std::clock();
}
//#####################################################################
// Scalar_Advance
//#####################################################################
template<class TV> void INCOMPRESSIBLE_INTERACTIVE_DRIVER<TV>::
Scalar_Advance(const T dt,const T time)
{
    example.Get_Scalar_Field_Sources(time);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before advection",0,1);
    ARRAY<T,TV_INT> density_ghost(example.mac_grid.Domain_Indices(2*example.number_of_ghost_cells+1));
    example.boundary->Set_Fixed_Boundary(true,0);
    example.boundary->Fill_Ghost_Cells(example.mac_grid,example.density,density_ghost,dt,time,2*example.number_of_ghost_cells+1);
    example.incompressible.advection->Update_Advection_Equation_Cell(example.mac_grid,example.density,density_ghost,example.face_velocities,*example.boundary,dt,time);
    example.boundary->Set_Fixed_Boundary(false);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after advection",0,1);
}
//#####################################################################
// Convect
//#####################################################################
template<class TV> void INCOMPRESSIBLE_INTERACTIVE_DRIVER<TV>::
Convect(const T dt,const T time)
{
    example.boundary->Set_Fixed_Boundary(true,0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before convection",0,1);
    example.incompressible.Advance_One_Time_Step_Convection(dt,time,example.face_velocities,example.face_velocities,example.number_of_ghost_cells);
    example.boundary->Set_Fixed_Boundary(false);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after convection",0,1);
}
//#####################################################################
// Add_Forces
//#####################################################################
template<class TV> void INCOMPRESSIBLE_INTERACTIVE_DRIVER<TV>::
Add_Forces(const T dt,const T time)
{
    if(example.incompressible.use_force){
        ARRAY<T,TV_INT> density_ghost(example.mac_grid.Domain_Indices(example.number_of_ghost_cells));
        example.boundary->Fill_Ghost_Cells(example.mac_grid,example.density,density_ghost,dt,time,example.number_of_ghost_cells);
        example.Get_Body_Force(example.incompressible.force,density_ghost,dt,time);}
    example.incompressible.Advance_One_Time_Step_Forces(example.face_velocities,dt,time,true,0,example.number_of_ghost_cells);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after forces",0,1);
    kinematic_evolution.Set_External_Velocities(example.rigid_geometry_collection.particles.V,example.rigid_geometry_collection.particles.angular_velocity,time+dt,time+dt);
}
//#####################################################################
// Project
//#####################################################################
template<class TV> void INCOMPRESSIBLE_INTERACTIVE_DRIVER<TV>::
Project(const T dt,const T time)
{
    example.Set_Boundary_Conditions(time+dt);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after boundary",0,1);
    example.projection->p*=dt; // rescale pressure for guess
    example.incompressible.Advance_One_Time_Step_Implicit_Part(example.face_velocities,dt,time,true);
    example.projection->p*=(1/dt); // unscale pressure
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after projection",0,1);
}
//#####################################################################
// Advance_To_Target_Time
//#####################################################################
template<class TV> void INCOMPRESSIBLE_INTERACTIVE_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;for(int substep=1;!done;substep++){
        // choose time step
        T dt=example.cfl>0?(example.cfl*example.incompressible.CFL(example.face_velocities)):FLT_MAX;
        if(example.mpi_grid) example.mpi_grid->Synchronize_Dt(dt);
        if(time+dt>=target_time){dt=target_time-time;done=true;}
        else if(time+2*dt>=target_time){dt=(T).5*(target_time-time);}
        std::stringstream ss;ss<<"dt is "<<dt<<std::endl;LOG::filecout(ss.str());
        ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>* advection_conservative=dynamic_cast<ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>*>(example.incompressible.advection);
        example.number_of_ghost_cells=max(3,1+example.incompressible.Real_CFL(example.face_velocities,false,false,dt));
        if(example.mpi_grid) example.mpi_grid->Synchronize_Ghost_Cells(example.number_of_ghost_cells);
        if(advection_conservative){
            advection_conservative->number_of_ghost_cells=example.number_of_ghost_cells;}

        // kinematic_update
        kinematic_evolution.Set_External_Positions(example.rigid_geometry_collection.particles.X,example.rigid_geometry_collection.particles.rotation,time);
        kinematic_evolution.Set_External_Velocities(example.rigid_geometry_collection.particles.V,example.rigid_geometry_collection.particles.angular_velocity,time,time);
        for(int i=1;i<=example.rigid_geometry_collection.kinematic_rigid_geometry.m;i++){
            RIGID_GEOMETRY<TV>& rigid_geometry=example.rigid_geometry_collection.Rigid_Geometry(i);            
            rigid_geometry.X()+=dt*rigid_geometry.V();
            rigid_geometry.Rotation()=ROTATION<TV>::From_Rotation_Vector(dt*rigid_geometry.Angular_Velocity())*rigid_geometry.Rotation();rigid_geometry.Rotation().Normalize();}
        kinematic_evolution.Set_External_Positions(example.rigid_geometry_collection.particles.X,example.rigid_geometry_collection.particles.rotation,time+dt);

        Scalar_Advance(dt,time);
        Convect(dt,time);
        Add_Forces(dt,time);
        Project(dt,time);
        time+=dt;}
}
//#####################################################################
// Function Update_Positions_From_Interactive
//#####################################################################
template<class TV> void INCOMPRESSIBLE_INTERACTIVE_DRIVER<TV>::Update_Positions_From_Interactive()
{ 
    if(visualizer->current_selection && visualizer->current_selection->type==OPENGL_SELECTION::COMPONENT_RIGID_BODIES_3D){
        OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T> *real_selection=dynamic_cast<OPENGL_SELECTION_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T>*>(visualizer->current_selection);
        example.rigid_geometry_collection.particles.X(real_selection->body_id)[1]=visualizer->opengl_world.object_location(1);
        example.rigid_geometry_collection.particles.X(real_selection->body_id)[2]=visualizer->opengl_world.object_location(2);
        example.rigid_geometry_collection.particles.X(real_selection->body_id)[3]=visualizer->opengl_world.object_location(3);
        if(rigid_component) rigid_component->Reinitialize(true,false);}
}
//#####################################################################
// Simulate_To_Frame
//#####################################################################
template<class TV> void INCOMPRESSIBLE_INTERACTIVE_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    while(frame<0 || current_frame<frame){
        if(visualizer->opengl_world.drag_current_selection) Update_Positions_From_Interactive();        
        kinematic_evolution.Get_Current_Kinematic_Keyframes(example.Time_At_Frame(current_frame+1)-time,time);        
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        if(example.mpi_grid && example.mpi_grid->threaded_grid && visualizer){
            if(scalar_component->opengl_scalar_field.upsample_scale>1){
                LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T> interpolation;
                for(typename GRID<TV>::CELL_ITERATOR iterator(scalar_component->opengl_scalar_field.grid);iterator.Valid();iterator.Next())
                    scalar_component->opengl_scalar_field.values(iterator.Cell_Index()+example.mpi_grid->threaded_grid->local_to_global_offset)=interpolation.Clamped_To_Array(example.mac_grid,example.density,iterator.Location());}
            else example.mpi_grid->threaded_grid->Sync_Scalar(example.density,scalar_component->opengl_scalar_field.values);}
            //example.mpi_grid->threaded_grid->Sync_Face_Scalar(example.face_velocities,mac_component->opengl_mac_velocity_field.face_velocities);} TODO: Fix inconsistency with opengl_mac_velocity
        //if(example.mpi_grid && example.mpi_grid->threaded_grid) pthread_barrier_wait(example.mpi_grid->threaded_grid->barr);
        if(((example.mpi_grid && example.mpi_grid->threaded_grid && example.mpi_grid->threaded_grid->tid==1) || !example.mpi_grid || !example.mpi_grid->threaded_grid) && visualizer){
            now=std::clock();
            if((now-since_last_frame)<(1./example.frame_rate)*CLOCKS_PER_SEC){
#ifdef WIN32
                Sleep((DWORD)(1000*(1./example.frame_rate)-(now-since_last_frame)/CLOCKS_PER_SEC));
#else
                sleep((unsigned int)((1./example.frame_rate)-(now-since_last_frame)/CLOCKS_PER_SEC));
#endif 
                now=std::clock();}
            since_last_frame=now;
            visualizer->Set_Frame(current_frame);}
        if(example.write_output_files) Write_Output_Files(++output_number);
        current_frame++;}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void INCOMPRESSIBLE_INTERACTIVE_DRIVER<TV>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    if(level<=example.write_substeps_level){
        example.frame_title=title;
        std::stringstream ss;ss<<"Writing substep ["<<example.frame_title<<"]: output_number="<<output_number+1<<", time="<<time<<", frame="<<current_frame<<", substep="<<substep<<std::endl;LOG::filecout(ss.str());
        Write_Output_Files(++output_number);example.frame_title="";}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void INCOMPRESSIBLE_INTERACTIVE_DRIVER<TV>::
Write_Output_Files(const int frame)
{
    FILE_UTILITIES::Create_Directory(example.output_directory);
    FILE_UTILITIES::Create_Directory(example.output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
    FILE_UTILITIES::Create_Directory(example.output_directory+"/common");
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+STRING_UTILITIES::string_sprintf("/%d/frame_title",frame),example.frame_title);
    if(frame==example.first_frame) 
        FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/first_frame",frame,"\n");
    example.Write_Output_Files(frame);
    FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/last_frame",frame,"\n");
}
//#####################################################################
template class INCOMPRESSIBLE_INTERACTIVE_DRIVER<VECTOR<float,1> >;
template class INCOMPRESSIBLE_INTERACTIVE_DRIVER<VECTOR<float,2> >;
template class INCOMPRESSIBLE_INTERACTIVE_DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INCOMPRESSIBLE_INTERACTIVE_DRIVER<VECTOR<double,1> >;
template class INCOMPRESSIBLE_INTERACTIVE_DRIVER<VECTOR<double,2> >;
template class INCOMPRESSIBLE_INTERACTIVE_DRIVER<VECTOR<double,3> >;
#endif
