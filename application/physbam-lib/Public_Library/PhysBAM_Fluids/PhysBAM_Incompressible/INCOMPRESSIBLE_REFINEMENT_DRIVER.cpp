//#####################################################################
// Copyright 2009, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_REFINEMENT_DRIVER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_REFINEMENT_EXAMPLE.h>
using namespace PhysBAM;
namespace{
template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
{
    ((INCOMPRESSIBLE_REFINEMENT_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
}
};
//#####################################################################
// Initialize
//#####################################################################
template<class TV> INCOMPRESSIBLE_REFINEMENT_DRIVER<TV>::
INCOMPRESSIBLE_REFINEMENT_DRIVER(INCOMPRESSIBLE_REFINEMENT_EXAMPLE<TV>& example)
    :example(example),kinematic_evolution(example.rigid_geometry_collection,true)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> INCOMPRESSIBLE_REFINEMENT_DRIVER<TV>::
~INCOMPRESSIBLE_REFINEMENT_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void INCOMPRESSIBLE_REFINEMENT_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void INCOMPRESSIBLE_REFINEMENT_DRIVER<TV>::
Initialize()
{
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);

    // setup time
    if(example.restart) current_frame=example.restart;else current_frame=example.first_frame;    
    time=example.Time_At_Frame(current_frame);

    // initialize collision objects
    example.Initialize_Bodies();
    kinematic_evolution.Get_Current_Kinematic_Keyframes(1/example.frame_rate,time);
    kinematic_evolution.Set_External_Positions(example.rigid_geometry_collection.particles.X,example.rigid_geometry_collection.particles.rotation,time);
    kinematic_evolution.Set_External_Velocities(example.rigid_geometry_collection.particles.V,example.rigid_geometry_collection.particles.angular_velocity,time,time);

    // setup domain boundaries and mpi
    if(example.coarse_mpi_grid) example.coarse_mpi_grid->Initialize(example.domain_boundary);
    example.incompressible.mpi_grid=example.fine_mpi_grid;
    example.incompressible.projection.elliptic_solver->mpi_grid=example.coarse_mpi_grid;
    if(example.fine_mpi_grid){
        example.boundary=new BOUNDARY_MPI<GRID<TV> >(example.fine_mpi_grid,example.boundary_scalar);
        example.boundary_coarse=new BOUNDARY_MPI<GRID<TV> >(example.coarse_mpi_grid,example.boundary_scalar);}
    else example.boundary=example.boundary_coarse=&example.boundary_scalar;
    VECTOR<VECTOR<bool,2>,TV::dimension> constant_extrapolation;constant_extrapolation.Fill(VECTOR<bool,2>::Constant_Vector(true));
    example.incompressible.boundary->Set_Constant_Extrapolation(constant_extrapolation);
    example.boundary->Set_Constant_Extrapolation(constant_extrapolation);
    example.incompressible.Set_Custom_Boundary(*example.boundary);
    example.Initialize_MPI();

    //threads
    if(example.thread_queue){
        example.threaded_advection_scalar=new THREADED_ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T>();
        example.threaded_advection_scalar->thread_queue=example.thread_queue;
        example.threaded_advection_scalar->row_jump=1;
        example.incompressible.Set_Custom_Advection(*example.threaded_advection_scalar);}
    else{
        example.advection_scalar=new ADVECTION_SEMI_LAGRANGIAN_UNIFORM<GRID<TV>,T>();
        example.incompressible.Set_Custom_Advection(*example.advection_scalar);}

    // setup grids and velocities
    example.incompressible.Initialize_Grids(example.fine_mac_grid);
    example.projection.Initialize_Grid(example.coarse_mac_grid);
    example.coarse_face_velocities.Resize(example.coarse_mac_grid);
    example.fine_face_velocities.Resize(example.fine_mac_grid);
    example.density.Resize(example.fine_mac_grid.Domain_Indices(3));
    example.Initialize_Fields();
    example.Initialize_Confinement();

    // setup laplace
    example.projection.elliptic_solver->Set_Relative_Tolerance((T)1e-9);
    example.projection.elliptic_solver->pcg.Set_Maximum_Iterations(1000);
    example.projection.elliptic_solver->pcg.evolution_solver_type=krylov_solver_cg;
    example.projection.elliptic_solver->pcg.cg_restart_iterations=40;

    if(example.kolmogorov){
        example.turbulence.Set_Lowest_Angular_Frequency(16);
        example.turbulence.Initialize_Grid(example.fine_mac_grid);
        example.turbulence.Generate_Initial_Turbulence(time,time+(T)20/example.frame_rate);}

    example.Set_Boundary_Conditions(time); // get so CFL is correct

    output_number=current_frame;
    if(example.restart) example.Read_Output_Files(example.restart);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after read",0,1);    

    if(!example.restart || example.split_dir!="") Write_Output_Files(current_frame);
    if(example.split_dir!="") exit(0);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("after init",0,1);    
}
//#####################################################################
// Advance_To_Target_Time
//#####################################################################
void Find_Neighbors(ARRAY<VECTOR<int,3> >& neighbor_cells,int number_of_cells,VECTOR<int,3>& local_index,int sub_scale)
{
    assert(number_of_cells==8);
    for(int i=2;i<=3;i++){
        if(local_index(i-1)<=sub_scale/2){neighbor_cells(i)(i-1)--;neighbor_cells(4)(i-1)--;neighbor_cells(i+4)(i-1)--;neighbor_cells(8)(i-1)--;}
        else{neighbor_cells(i)(i-1)++;neighbor_cells(4)(i-1)++;neighbor_cells(i+4)(i-1)++;neighbor_cells(8)(i-1)++;}}
    if(local_index(3)<=sub_scale/2){neighbor_cells(5)(3)--;neighbor_cells(6)(3)--;neighbor_cells(7)(3)--;neighbor_cells(8)(3)--;}
    else{neighbor_cells(5)(3)++;neighbor_cells(6)(3)++;neighbor_cells(7)(3)++;neighbor_cells(8)(3)++;}
}
void Find_Neighbors(ARRAY<VECTOR<int,2> >& neighbor_cells,int number_of_cells,VECTOR<int,2>& local_index,int sub_scale)
{
    assert(number_of_cells==4);
    for(int i=2;i<=3;i++){
        if(local_index(i-1)<=sub_scale/2){neighbor_cells(i)(i-1)--;neighbor_cells(4)(i-1)--;}
        else{neighbor_cells(i)(i-1)++;neighbor_cells(4)(i-1)++;}}
}
template<class TV> void INCOMPRESSIBLE_REFINEMENT_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP","substep %d",substep);
        // choose time step
        T dt=example.cfl*example.incompressible.CFL(example.fine_face_velocities);
        if(example.fine_mpi_grid) example.fine_mpi_grid->Synchronize_Dt(dt);
        if(time+dt>=target_time){dt=target_time-time;done=true;}
        else if(time+2*dt>=target_time){dt=(T).5*(target_time-time);}

        // kinematic_update
        kinematic_evolution.Set_External_Positions(example.rigid_geometry_collection.particles.X,example.rigid_geometry_collection.particles.rotation,time);
        kinematic_evolution.Set_External_Velocities(example.rigid_geometry_collection.particles.V,example.rigid_geometry_collection.particles.angular_velocity,time,time);
        for(int i=1;i<=example.rigid_geometry_collection.kinematic_rigid_geometry.m;i++){
            RIGID_GEOMETRY<TV>& rigid_geometry=example.rigid_geometry_collection.Rigid_Geometry(i);            
            rigid_geometry.X()+=dt*rigid_geometry.V();
            rigid_geometry.Rotation()=ROTATION<TV>::From_Rotation_Vector(dt*rigid_geometry.Angular_Velocity())*rigid_geometry.Rotation();rigid_geometry.Rotation().Normalize();}
        kinematic_evolution.Set_External_Positions(example.rigid_geometry_collection.particles.X,example.rigid_geometry_collection.particles.rotation,time+dt);

        LOG::Time("Scalar Advance");
        // scalar update
        example.Get_Scalar_Field_Sources(time);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("before advection",0,1);
        ARRAY<T,TV_INT> density_ghost(example.fine_mac_grid.Domain_Indices(example.number_of_ghost_cells));
        example.boundary->Fill_Ghost_Cells(example.fine_mac_grid,example.density,density_ghost,dt,time,example.number_of_ghost_cells);
        example.incompressible.advection->Update_Advection_Equation_Cell(example.fine_mac_grid,example.density,density_ghost,example.fine_face_velocities,*example.boundary,dt,time);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after advection",0,1);

        // velocity update
        LOG::Time("Velocity Advection");
        example.boundary->Set_Fixed_Boundary(true,0);
        example.incompressible.Advance_One_Time_Step_Convection(dt,time,example.fine_face_velocities,example.fine_face_velocities,example.number_of_ghost_cells);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after convection",0,1);
        example.boundary->Set_Fixed_Boundary(false);
        LOG::Time("Add Forces");
        if(example.incompressible.use_force){
            example.boundary->Fill_Ghost_Cells(example.fine_mac_grid,example.density,density_ghost,dt,time,example.number_of_ghost_cells);
            example.Get_Body_Force(example.incompressible.force,density_ghost,dt,time);}
        if(example.use_coarse_forces){
            example.Preprocess_Projection(dt,time);
            ARRAY<T,FACE_INDEX<TV::dimension> > coarse_face_velocities_ghost;coarse_face_velocities_ghost.Resize(example.coarse_mac_grid,example.number_of_ghost_cells,false);
            example.boundary_coarse->Fill_Ghost_Cells_Face(example.coarse_mac_grid,example.coarse_face_velocities,coarse_face_velocities_ghost,time,example.number_of_ghost_cells);
            if(example.incompressible.vorticity_confinement || example.incompressible.use_variable_vorticity_confinement){
                ARRAY<TV,TV_INT> F(example.coarse_mac_grid.Cell_Indices(1),false);F.Fill(TV());
                ARRAY<TV,TV_INT> F_fine(example.fine_mac_grid.Cell_Indices(1),false);F_fine.Fill(TV());
                example.incompressible.Compute_Vorticity_Confinement_Force(example.coarse_mac_grid,coarse_face_velocities_ghost,F);
                if(example.incompressible.use_variable_vorticity_confinement){F*=dt*(T).5;
                    for(typename GRID<TV>::CELL_ITERATOR iterator(example.coarse_mac_grid,1);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
                        F(cell)=example.incompressible.variable_vorticity_confinement(cell)!=0?(F(cell)*example.incompressible.variable_vorticity_confinement(cell)):TV();}}
                else F*=dt*example.incompressible.vorticity_confinement*(T).5;
                if(!example.use_interpolated_vorticity) for(typename GRID<TV>::CELL_ITERATOR iterator(example.fine_mac_grid);iterator.Valid();iterator.Next()){
                    TV_INT coarse_cell=(iterator.Cell_Index()-TV_INT::All_Ones_Vector())/example.sub_scale+TV_INT::All_Ones_Vector();
                    F_fine(iterator.Cell_Index())=F(coarse_cell);}
                else for(typename GRID<TV>::CELL_ITERATOR iterator(example.fine_mac_grid);iterator.Valid();iterator.Next()){
                    TV_INT coarse_cell=(iterator.Cell_Index()-TV_INT::All_Ones_Vector())/example.sub_scale+TV_INT::All_Ones_Vector();
                    TV_INT local_index=iterator.Cell_Index()-(coarse_cell-TV_INT::All_Ones_Vector())*example.sub_scale;
                    int number_of_cells=(int)pow((T)2,TV::dimension);
                    ARRAY<TV_INT> neighbor_coarse_cells(number_of_cells);
                    ARRAY<T> areas(number_of_cells);
                    T total_area=0;
                    for(int i=1;i<=number_of_cells;i++) neighbor_coarse_cells(i)=coarse_cell;Find_Neighbors(neighbor_coarse_cells,number_of_cells,local_index,example.sub_scale);
                    for(int i=1;i<=number_of_cells;i++){
                        TV diff_vector=iterator.Location()-example.coarse_mac_grid.Center(neighbor_coarse_cells(i));
                        areas(i)=T(1);
                        for(int axis=1;axis<=TV::dimension;axis++) areas(i)*=example.coarse_mac_grid.dX(axis)-abs(diff_vector(axis));
                        total_area+=areas(i);}
                    for(int i=1;i<=number_of_cells;i++) F_fine(iterator.Cell_Index())+=(areas(i)/total_area)*F(neighbor_coarse_cells(i));}
                example.incompressible.Apply_Vorticity_Confinement_Force(example.fine_face_velocities,F_fine);}}
        else example.incompressible.Advance_One_Time_Step_Forces(example.fine_face_velocities,dt,time,false,0,example.number_of_ghost_cells);
        if(example.kolmogorov){
            while(time>example.turbulence.time_end) example.turbulence.Advance_Turbulence();
            T b=(T)1-pow(max((T)1-example.kolmogorov,(T)0),dt),fraction=(time-example.turbulence.time_start)/(example.turbulence.time_end-example.turbulence.time_start);
            for(typename GRID<TV>::FACE_ITERATOR iterator(example.fine_mac_grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();T& face_velocity=example.fine_face_velocities.Component(axis)(iterator.Face_Index());
                face_velocity=(1-b)*face_velocity+b*example.turbulence.Turbulent_Face_Velocity(axis,iterator.Location(),fraction);}}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after forces",0,1);
        kinematic_evolution.Set_External_Velocities(example.rigid_geometry_collection.particles.V,example.rigid_geometry_collection.particles.angular_velocity,time+dt,time+dt);
        example.Set_Boundary_Conditions(time+dt);
        //example.incompressible.Set_Dirichlet_Boundary_Conditions(0,0);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after boundary",0,1);
        example.Preprocess_Projection(dt,time);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after preprocess",0,1);
        LOG::Time("Project");
        example.incompressible.projection.p*=dt; // rescale pressure for guess
        example.incompressible.Advance_One_Time_Step_Implicit_Part(example.coarse_face_velocities,dt,time,false,example.boundary_coarse);
        example.incompressible.projection.p*=(1/dt); // unscale pressure
        if(example.boundary_coarse) example.boundary_coarse->Apply_Boundary_Condition_Face(example.coarse_mac_grid,example.coarse_face_velocities,time+dt);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after projection",0,1);
        LOG::Time("Postprocess Project");
        example.Postprocess_Projection(dt,time);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after postprocess",0,1);
        example.boundary->Apply_Boundary_Condition_Face(example.fine_mac_grid,example.fine_face_velocities,time+dt);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after apply boundary condition face",0,1);
        
        LOG::Stop_Time();
        
        FILE_UTILITIES::Write_To_Text_File(example.output_directory+"/common/last_substep",substep," ",dt," ",time,"\n");

        time+=dt;
    }
}
//#####################################################################
// Simulate_To_Frame
//#####################################################################
template<class TV> void INCOMPRESSIBLE_REFINEMENT_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    while(current_frame<frame){
        LOG::SCOPE scope("FRAME","Frame %d",current_frame+1);
        example.Preprocess_Frame(current_frame+1);
        kinematic_evolution.Get_Current_Kinematic_Keyframes(example.Time_At_Frame(current_frame+1)-time,time);        
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        example.Postprocess_Frame(current_frame+1);
        Write_Output_Files(++output_number);
        current_frame++;}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void INCOMPRESSIBLE_REFINEMENT_DRIVER<TV>::
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
template<class TV> void INCOMPRESSIBLE_REFINEMENT_DRIVER<TV>::
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
template class INCOMPRESSIBLE_REFINEMENT_DRIVER<VECTOR<float,2> >;
template class INCOMPRESSIBLE_REFINEMENT_DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INCOMPRESSIBLE_REFINEMENT_DRIVER<VECTOR<double,2> >;
template class INCOMPRESSIBLE_REFINEMENT_DRIVER<VECTOR<double,3> >;
#endif
