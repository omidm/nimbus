//#####################################################################
// Copyright 2009-2010, Mridul Aanjaneya, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_DRIVER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_EXAMPLE.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM_FORWARD.h>

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
using namespace PhysBAM;
namespace{
template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
{
    ((INCOMPRESSIBLE_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
}
};
//#####################################################################
// Initialize
//#####################################################################
template<class TV> INCOMPRESSIBLE_DRIVER<TV>::
INCOMPRESSIBLE_DRIVER(INCOMPRESSIBLE_EXAMPLE<TV>& example)
    :example(example),kinematic_evolution(example.rigid_geometry_collection,true)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> INCOMPRESSIBLE_DRIVER<TV>::
~INCOMPRESSIBLE_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void INCOMPRESSIBLE_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void INCOMPRESSIBLE_DRIVER<TV>::
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
    if(example.mpi_grid){
        example.boundary=new BOUNDARY_MPI<GRID<TV> >(example.mpi_grid,example.boundary_scalar);}
    else example.boundary=&example.boundary_scalar;
    example.incompressible.Set_Custom_Boundary(*example.boundary);

    // setup grids and velocities
    example.incompressible.Initialize_Grids(example.mac_grid);
    example.projection->Initialize_Grid(example.mac_grid);
    example.face_velocities.Resize(example.mac_grid);
    example.face_velocities_save.Resize(example.mac_grid);
    sum_jc.Resize(example.mac_grid);
    for(typename GRID<TV>::FACE_ITERATOR iterator(example.mac_grid);iterator.Valid();iterator.Next()) sum_jc(iterator.Full_Index())=1;
    example.density.Resize(example.mac_grid.Domain_Indices(example.number_of_ghost_cells));
    example.incompressible.kinetic_energy.Resize(example.mac_grid);
    if(example.incompressible.conserve_kinetic_energy) example.incompressible.Conserve_Kinetic_Energy();
    example.Initialize_Fields();
    example.Initialize_Confinement();

    // setup laplace
    example.projection->elliptic_solver->Set_Relative_Tolerance((T)1e-8);
    example.projection->elliptic_solver->pcg.Set_Maximum_Iterations(40);
    example.projection->elliptic_solver->pcg.evolution_solver_type=krylov_solver_cg;
    example.projection->elliptic_solver->pcg.cg_restart_iterations=40;

    if(example.restart) example.Read_Output_Files(example.restart);
    if(example.restart){
        ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>* advection_conservative=dynamic_cast<ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>*>(example.incompressible.advection);
        if(advection_conservative) sum_jc=advection_conservative->sum_jc;}
    
    // setup domain boundaries
    VECTOR<VECTOR<bool,2>,TV::dimension> constant_extrapolation;constant_extrapolation.Fill(VECTOR<bool,2>::Constant_Vector(true));
    example.incompressible.boundary->Set_Constant_Extrapolation(constant_extrapolation);
    example.boundary->Set_Constant_Extrapolation(constant_extrapolation);
    example.incompressible.Set_Custom_Boundary(*example.boundary);
    example.Set_Boundary_Conditions(time); // get so CFL is correct

    if(!example.restart) example.incompressible.projection.Make_Divergence_Free(example.face_velocities,0,time);
    if(example.use_viscosity){
        PHYSBAM_FATAL_ERROR();
        example.incompressible.Use_Explicit_Part_Of_Implicit_Viscosity(false);    
        example.incompressible.Set_Viscosity(1);
        example.incompressible.Set_Variable_Viscosity(false);}

    example.Get_Scalar_Field_Sources(time);
    example.face_velocities_save=example.face_velocities;
    example.incompressible.Calculate_Kinetic_Energy(example.face_velocities,0,0);
    if(!example.restart) Write_Output_Files(example.first_frame);
}
//#####################################################################
// Scalar_Advance
//#####################################################################
template<class TV> void INCOMPRESSIBLE_DRIVER<TV>::
Scalar_Advance(const T dt,const T time)
{
    LOG::Time("Scalar Advance");
    //example.incompressible.Correct_Kinetic_Energy(example.face_velocities,dt,time);
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
template<class TV> void INCOMPRESSIBLE_DRIVER<TV>::
Convect(const T dt,const T time)
{
    //static bool first=true;
    LOG::Time("Velocity Advection");
    example.boundary->Set_Fixed_Boundary(true,0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("before convection",0,1);
    if(example.incompressible.conserve_kinetic_energy){
        //ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>* advection_conservative=dynamic_cast<ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>*>(example.incompressible.advection);
        //ARRAY<T,FACE_INDEX<TV::dimension> > sum_jc_back(example.mac_grid);
        //if(!first) sum_jc_back=advection_conservative->sum_jc;
        //advection_conservative->sum_jc=sum_jc;
        ARRAY<T,FACE_INDEX<TV::dimension> > kinetic_energy_ghost(example.mac_grid,2*example.number_of_ghost_cells+1);
        ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities_ghost(example.mac_grid,2*example.number_of_ghost_cells+1);
        example.boundary->Fill_Ghost_Cells_Face(example.mac_grid,example.face_velocities,face_velocities_ghost,dt,2*example.number_of_ghost_cells+1);
        example.boundary->Fill_Ghost_Cells_Face(example.mac_grid,example.incompressible.kinetic_energy,kinetic_energy_ghost,dt,2*example.number_of_ghost_cells+1);
        example.incompressible.advection->Update_Advection_Equation_Face(example.mac_grid,example.incompressible.kinetic_energy,kinetic_energy_ghost,face_velocities_ghost,*example.boundary,dt,time);
        example.incompressible.Advance_One_Time_Step_Convection(dt,time,example.face_velocities,example.face_velocities,example.number_of_ghost_cells);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("before correction",0,1);
        //sum_jc=advection_conservative->sum_jc;
        //advection_conservative->sum_jc=sum_jc_back;
        example.incompressible.Correct_Kinetic_Energy(example.face_velocities,dt,time);} 
    else example.incompressible.Advance_One_Time_Step_Convection(dt,time,example.face_velocities,example.face_velocities,example.number_of_ghost_cells);
    example.boundary->Set_Fixed_Boundary(false);
    //first=false;
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after convection",0,1);
}
//#####################################################################
// Add_Forces
//#####################################################################
template<class TV> void INCOMPRESSIBLE_DRIVER<TV>::
Add_Forces(const T dt,const T time)
{
    LOG::Time("Add Forces");
    //example.incompressible.Add_Energy_With_Vorticity(example.face_velocities,example.domain_boundary,dt,time,example.number_of_ghost_cells,0,&example.density);
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
template<class TV> void INCOMPRESSIBLE_DRIVER<TV>::
Project(const T dt,const T time)
{
    LOG::Time("Project");
    example.Set_Boundary_Conditions(time+dt);
    example.incompressible.Calculate_Kinetic_Energy(example.face_velocities,dt,time);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after boundary",0,1);
    example.projection->p*=dt; // rescale pressure for guess
    example.incompressible.Advance_One_Time_Step_Implicit_Part(example.face_velocities,dt,time,true);
    example.projection->p*=(1/dt); // unscale pressure
    example.incompressible.Apply_Pressure_Kinetic_Energy(example.face_velocities_save,example.face_velocities,dt,time);
    //example.incompressible.Correct_Kinetic_Energy(example.face_velocities,dt,time);
    example.face_velocities_save=example.face_velocities;
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after projection",0,1);
}
//#####################################################################
// Advance_To_Target_Time
//#####################################################################
template<class TV> void INCOMPRESSIBLE_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP","substep %d",substep);
        // choose time step
        T dt=example.cfl*example.incompressible.CFL(example.face_velocities);
        if(example.mpi_grid) example.mpi_grid->Synchronize_Dt(dt);
        if(time+dt>=target_time){dt=target_time-time;done=true;}
        else if(time+2*dt>=target_time){dt=(T).5*(target_time-time);}
        std::stringstream ss;ss<<"dt is "<<dt<<std::endl;LOG::filecout(ss.str());

        // kinematic_update
        kinematic_evolution.Set_External_Positions(example.rigid_geometry_collection.particles.X,example.rigid_geometry_collection.particles.rotation,time);
        kinematic_evolution.Set_External_Velocities(example.rigid_geometry_collection.particles.V,example.rigid_geometry_collection.particles.angular_velocity,time,time);
        for(int i=1;i<=example.rigid_geometry_collection.kinematic_rigid_geometry.m;i++){
            RIGID_GEOMETRY<TV>& rigid_geometry=example.rigid_geometry_collection.Rigid_Geometry(i);            
            rigid_geometry.X()+=dt*rigid_geometry.V();
            rigid_geometry.Rotation()=ROTATION<TV>::From_Rotation_Vector(dt*rigid_geometry.Angular_Velocity())*rigid_geometry.Rotation();rigid_geometry.Rotation().Normalize();}
        kinematic_evolution.Set_External_Positions(example.rigid_geometry_collection.particles.X,example.rigid_geometry_collection.particles.rotation,time+dt);

        if(example.order>1){
            RUNGEKUTTA<ARRAY<T,FACE_INDEX<TV::dimension> > > rungekutta_u(example.face_velocities);
            RUNGEKUTTA<ARRAY<T,TV_INT> > rungekutta_scalar(example.density);
            rungekutta_u.Set_Grid_And_Boundary_Condition(example.mac_grid,*example.boundary);
            rungekutta_u.Set_Order(example.order);rungekutta_u.Set_Time(time);rungekutta_u.Start(dt);
            rungekutta_scalar.Set_Grid_And_Boundary_Condition(example.mac_grid,*example.boundary);
            rungekutta_scalar.Set_Order(example.order);rungekutta_scalar.Set_Time(time);rungekutta_scalar.Start(dt);T rk_time=time;
            for(int rk_substep=1;rk_substep<=rungekutta_u.order;rk_substep++){
                Scalar_Advance(dt,rk_time);
                // velocity update
                if(!example.analytic_test){
                    Convect(dt,rk_time);
                    Add_Forces(dt,rk_time);
                    Project(dt,rk_time);}
                rungekutta_scalar.Main();rk_time=rungekutta_u.Main();}}
        else{
            Scalar_Advance(dt,time);
            // velocity update
            if(!example.analytic_test){
                Convect(dt,time);
                Add_Forces(dt,time);
                Project(dt,time);}}
        time+=dt;
    }
}
//#####################################################################
// Simulate_To_Frame
//#####################################################################
template<class TV> void INCOMPRESSIBLE_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    while(current_frame<frame){
        LOG::SCOPE scope("FRAME","Frame %d",current_frame+1);
        kinematic_evolution.Get_Current_Kinematic_Keyframes(example.Time_At_Frame(current_frame+1)-time,time);        
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        Write_Output_Files(++output_number);
        current_frame++;}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void INCOMPRESSIBLE_DRIVER<TV>::
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
template<class TV> void INCOMPRESSIBLE_DRIVER<TV>::
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
template class INCOMPRESSIBLE_DRIVER<VECTOR<float,1> >;
template class INCOMPRESSIBLE_DRIVER<VECTOR<float,2> >;
template class INCOMPRESSIBLE_DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INCOMPRESSIBLE_DRIVER<VECTOR<double,1> >;
template class INCOMPRESSIBLE_DRIVER<VECTOR<double,2> >;
template class INCOMPRESSIBLE_DRIVER<VECTOR<double,3> >;
#endif
