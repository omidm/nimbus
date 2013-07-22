//#####################################################################
// Copyright 2009, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_ADAPTIVE_DRIVER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_ADAPTIVE_EXAMPLE.h>

#include <PhysBAM_Tools/Read_Write/Math_Tools/READ_WRITE_RANGE.h>
using namespace PhysBAM;
namespace{
template<class TV> void Write_Substep_Helper(void* writer,const std::string& title,int substep,int level)
{
    ((INCOMPRESSIBLE_ADAPTIVE_DRIVER<TV>*)writer)->Write_Substep(title,substep,level);
}
};
//#####################################################################
// Initialize
//#####################################################################
template<class TV> INCOMPRESSIBLE_ADAPTIVE_DRIVER<TV>::
INCOMPRESSIBLE_ADAPTIVE_DRIVER(INCOMPRESSIBLE_ADAPTIVE_EXAMPLE<TV>& example)
    :example(example)
{
    DEBUG_SUBSTEPS::Set_Substep_Writer((void*)this,&Write_Substep_Helper<TV>);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> INCOMPRESSIBLE_ADAPTIVE_DRIVER<TV>::
~INCOMPRESSIBLE_ADAPTIVE_DRIVER()
{
    DEBUG_SUBSTEPS::Clear_Substep_Writer((void*)this);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void INCOMPRESSIBLE_ADAPTIVE_DRIVER<TV>::
Execute_Main_Program()
{
    Initialize();
    Simulate_To_Frame(example.last_frame);
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void INCOMPRESSIBLE_ADAPTIVE_DRIVER<TV>::
Initialize()
{
    DEBUG_SUBSTEPS::Set_Write_Substeps_Level(example.write_substeps_level);

    // setup time
    current_frame=example.first_frame;
    time=example.Time_At_Frame(current_frame);

    // setup grids and velocities
    example.incompressible.Initialize_Grids(example.mac_grid);
    example.projection.Initialize_Grid(example.mac_grid);
    example.face_velocities.Resize(example.mac_grid);
    example.density.Resize(example.mac_grid.Domain_Indices(),false);
    example.temperature.Resize(example.mac_grid.Domain_Indices(),false);
    for(typename GRID<TV>::CELL_ITERATOR iterator(example.mac_grid);iterator.Valid();iterator.Next()){
        example.face_velocities.sub_arrays(iterator.Cell_Index())=new FACE_ARRAY_ADAPTIVE<T,TV::dimension>(*example.mac_grid.sub_mac_grids(iterator.Cell_Index()));
        example.density(iterator.Cell_Index()).Resize(example.mac_grid.sub_mac_grids(iterator.Cell_Index())->Domain_Indices(3));
        example.temperature(iterator.Cell_Index()).Resize(example.mac_grid.sub_mac_grids(iterator.Cell_Index())->Domain_Indices(3));}
    example.Initialize_Fields();

    // setup laplace
    example.projection.elliptic_solver->pcg.Set_Maximum_Iterations(40);
    example.projection.elliptic_solver->pcg.Use_Modified_Incomplete_Cholesky();

    // setup domain boundaries
    VECTOR<VECTOR<bool,2>,TV::dimension> constant_extrapolation;constant_extrapolation.Fill(VECTOR<bool,2>::Constant_Vector(true));    
    example.incompressible.boundary->Set_Constant_Extrapolation(constant_extrapolation);
    example.boundary_scalar.Set_Constant_Extrapolation(constant_extrapolation);

    Write_Output_Files(example.first_frame);
    output_number=example.first_frame;
}
//#####################################################################
// Advance_To_Target_Time
//#####################################################################
template<class TV> void INCOMPRESSIBLE_ADAPTIVE_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP","substep %d",substep);
        // choose time step
        T dt=target_time-time;
        for(typename GRID<TV>::CELL_ITERATOR iterator(example.mac_grid);iterator.Valid();iterator.Next()){
            GRID<TV>& local_mac_grid=*example.mac_grid.sub_mac_grids(iterator.Cell_Index());ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities=*example.face_velocities.sub_arrays(iterator.Cell_Index());
            example.Set_Boundary_Conditions(time,local_mac_grid,local_face_velocities,0);
            example.incompressible.Initialize_Grids(local_mac_grid);
            dt=min(dt,example.cfl*example.incompressible.CFL(local_face_velocities));}
        if(time+dt>=target_time){dt=target_time-time;done=true;}
        else if(time+2*dt>=target_time){dt=(T).5*(target_time-time);}

        LOG::Time("Scalar Advance");
        // scalar update
        PHYSBAM_DEBUG_WRITE_SUBSTEP("before advection",0,1);
        //TODO: Make this a new type of advection?
        LOG::Time("Velocity Advection and Forces");
        for(typename GRID<TV>::CELL_ITERATOR iterator(example.mac_grid);iterator.Valid();iterator.Next()){
            GRID<TV>& local_mac_grid=*example.mac_grid.sub_mac_grids(iterator.Cell_Index());ARRAY<T,FACE_INDEX<TV::dimension> >& local_face_velocities=*example.face_velocities.sub_arrays(iterator.Cell_Index());
            ARRAY<T,TV_INT>& local_density=example.density(iterator.Cell_Index());ARRAY<T,TV_INT>& local_temperature=example.temperature(iterator.Cell_Index());
            ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities_ghost(local_mac_grid,example.number_of_ghost_cells);
            ARRAY<T,TV_INT> density_ghost(local_mac_grid.Domain_Indices(example.number_of_ghost_cells)),temperature_ghost(local_mac_grid.Domain_Indices(example.number_of_ghost_cells));
            example.incompressible.Initialize_Grids(local_mac_grid);
            example.Get_Scalar_Field_Sources(time,local_mac_grid,local_density);
            example.boundary_scalar.Fill_Ghost_Cells(local_mac_grid,local_density,density_ghost,dt,time,example.number_of_ghost_cells);
            example.boundary_scalar.Fill_Ghost_Cells(local_mac_grid,local_temperature,temperature_ghost,dt,time,example.number_of_ghost_cells);
            example.boundary_scalar.Fill_Ghost_Cells_Face(local_mac_grid,local_face_velocities,face_velocities_ghost,time,example.number_of_ghost_cells);
            for(typename GRID<TV>::CELL_ITERATOR local_iterator(local_mac_grid,example.number_of_ghost_cells,GRID<TV>::GHOST_REGION);local_iterator.Valid();local_iterator.Next()){
                TV_INT global_index=example.mac_grid.Clamp_To_Cell(local_iterator.Location());
                GRID<TV>& ghost_mac_grid=*example.mac_grid.sub_mac_grids(global_index);
                ARRAY<T,TV_INT>& ghost_local_density=example.density(global_index);ARRAY<T,TV_INT>& ghost_local_temperature=example.temperature(global_index);
                TV_INT ghost_local_index=ghost_mac_grid.Clamp_To_Cell(local_iterator.Location());
                density_ghost(local_iterator.Cell_Index())=ghost_local_density(ghost_local_index);temperature_ghost(local_iterator.Cell_Index())=ghost_local_temperature(ghost_local_index);}
            for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid,example.number_of_ghost_cells);local_iterator.Valid();local_iterator.Next()){
                if(!local_mac_grid.Outside(local_iterator.Location())) continue;
                TV_INT axis;axis(local_iterator.Axis())=1;
                TV_INT global_index=example.mac_grid.Clamp_To_Cell(local_iterator.Location());
                GRID<TV>& ghost_mac_grid=*example.mac_grid.sub_mac_grids(global_index);ARRAY<T,FACE_INDEX<TV::dimension> >& ghost_face_velocities=*example.face_velocities.sub_arrays(global_index);
                FACE_INDEX<TV::dimension> ghost_local_index(local_iterator.Axis(),ghost_mac_grid.Clamped_Index(local_iterator.Location())+TV_INT::All_Ones_Vector());
                for(int i=1;i<=TV::dimension;i++) if(i!=local_iterator.Axis() && ghost_local_index.index(i)>local_mac_grid.Counts()(i)) ghost_local_index.index(i)=local_mac_grid.Counts()(i);
                face_velocities_ghost(local_iterator.Full_Index())=ghost_face_velocities(ghost_local_index);}
            example.advection_scalar.Update_Advection_Equation_Cell(local_mac_grid,local_density,density_ghost,local_face_velocities,example.boundary_scalar,dt,time);
            example.advection_scalar.Update_Advection_Equation_Cell(local_mac_grid,local_temperature,temperature_ghost,local_face_velocities,example.boundary_scalar,dt,time);
            example.incompressible.Advance_One_Time_Step_Convection(dt,time,face_velocities_ghost,local_face_velocities,example.number_of_ghost_cells);
            example.incompressible.Advance_One_Time_Step_Forces(local_face_velocities,dt,time,false,0,example.number_of_ghost_cells);}
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after advection",0,1);

        LOG::Time("Project");
        example.incompressible.Initialize_Grids(example.mac_grid);
        example.Preprocess_Projection(dt,time);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after preprocess",0,1);
        example.Set_Boundary_Conditions(time+dt,example.mac_grid,example.face_velocities,&example.incompressible.projection);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after boundary",0,1);
        example.incompressible.projection.p*=dt; // rescale pressure for guess
        example.incompressible.Advance_One_Time_Step_Implicit_Part(example.face_velocities,dt,time,false);
        example.incompressible.projection.p*=(1/dt); // unscale pressure
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after projection",0,1);
        example.Postprocess_Projection(dt,time);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("after postprocess",0,1);

        time+=dt;
    }
}
//#####################################################################
// Simulate_To_Frame
//#####################################################################
template<class TV> void INCOMPRESSIBLE_ADAPTIVE_DRIVER<TV>::
Simulate_To_Frame(const int frame)
{
    while(current_frame<frame){
        LOG::SCOPE scope("FRAME","Frame %d",current_frame+1);
        Advance_To_Target_Time(example.Time_At_Frame(current_frame+1));
        Write_Output_Files(++output_number);
        current_frame++;}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void INCOMPRESSIBLE_ADAPTIVE_DRIVER<TV>::
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
template<class TV> void INCOMPRESSIBLE_ADAPTIVE_DRIVER<TV>::
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
template class INCOMPRESSIBLE_ADAPTIVE_DRIVER<VECTOR<float,2> >;
template class INCOMPRESSIBLE_ADAPTIVE_DRIVER<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INCOMPRESSIBLE_ADAPTIVE_DRIVER<VECTOR<double,2> >;
template class INCOMPRESSIBLE_ADAPTIVE_DRIVER<VECTOR<double,3> >;
#endif
