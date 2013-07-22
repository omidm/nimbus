//#####################################################################
// Copyright 2009, Michael Lentine, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_FACE_INDEX.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/FLUID_GRAVITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/INCOMPRESSIBILITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_REFINEMENT_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_EXAMPLE.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM_FORWARD.h>
using namespace PhysBAM;
//#####################################################################
// INCOMPRESSIBLE_EXAMPLE
//#####################################################################
template<class TV> INCOMPRESSIBLE_EXAMPLE<TV>::
INCOMPRESSIBLE_EXAMPLE(const STREAM_TYPE stream_type_input,const int number_of_threads,const int scaling_factor)
    :stream_type(stream_type_input),initial_time(0),first_frame(0),last_frame(100),frame_rate(24),
    restart(0),write_debug_data(false),write_output_files(true),analytic_test(false),order(1),output_directory("output"),
    number_of_ghost_cells(3),cfl((T).9),use_viscosity(false),mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),mpi_grid(0),thread_queue(number_of_threads>1?new THREAD_QUEUE(number_of_threads):0),
    projection(scaling_factor>1?(new PROJECTION_REFINEMENT_UNIFORM<GRID<TV> >(mac_grid,scaling_factor,(T)1,false,false,true,true,thread_queue)):(new PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >(mac_grid,false,false,false,false,thread_queue))),
    incompressible(mac_grid,*projection),advection_scalar(thread_queue),boundary(0),rigid_geometry_collection(this)
{
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    PROCESS_UTILITIES::Set_Backtrace(true);
    incompressible.Set_Custom_Advection(advection_scalar);
    for(int i=1;i<=TV::dimension;i++){domain_boundary(i)(1)=true;domain_boundary(i)(2)=true;}
    Initialize_Geometry_Particle();
    Initialize_Read_Write_Structures();
}
//#####################################################################
// ~INCOMPRESSIBLE_EXAMPLE
//#####################################################################
template<class TV> INCOMPRESSIBLE_EXAMPLE<TV>::
~INCOMPRESSIBLE_EXAMPLE()
{
    if(mpi_grid) delete boundary;
    if(thread_queue) delete thread_queue;
    delete projection;
}
//#####################################################################
// 
//#####################################################################
template<class TV> void INCOMPRESSIBLE_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    if(!write_output_files) return;
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    if(mpi_grid) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/global_grid",mpi_grid->global_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/grid",mac_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",mac_grid);
    if(incompressible.use_analytic_energy) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/analytic_energy",incompressible.analytic_energy);
    ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>* advection_conservative=dynamic_cast<ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>*>(incompressible.advection);
    if(advection_conservative){
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/weights_barjc",advection_conservative->sum_jc);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/weights_barjc_cell",advection_conservative->sum_jc_cell);}
    PROJECTION_REFINEMENT_UNIFORM<GRID<TV> >* projection_refinement=dynamic_cast<PROJECTION_REFINEMENT_UNIFORM<GRID<TV> >*>(projection);
    if(projection_refinement){
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/coarse_grid",projection_refinement->coarse_grid);
        if(!projection->elliptic_solver->psi_D.Valid_Index(mac_grid.Domain_Indices().max_corner)){
            FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/coarse_psi_N",projection->elliptic_solver->psi_N);
            FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/coarse_psi_D",projection->elliptic_solver->psi_D);}
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/coarse_mac_velocities",projection_refinement->coarse_face_velocities);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/coarse_beta_face",projection_refinement->beta_face);}
    if(write_debug_data){
        if(advection_conservative){
            FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/weights_i",advection_conservative->weights_to);
            FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/weights_j",advection_conservative->weights_from);}
        ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities_ghost(mac_grid,number_of_ghost_cells,false),kinetic_energy_ghost(mac_grid,number_of_ghost_cells,false);
        ARRAY<T,TV_INT> density_ghost(mac_grid.Domain_Indices(number_of_ghost_cells),false);
        incompressible.boundary->Fill_Ghost_Cells_Face(mac_grid,face_velocities,face_velocities_ghost,0,number_of_ghost_cells);
        if(incompressible.conserve_kinetic_energy) incompressible.boundary->Fill_Ghost_Cells_Face(mac_grid,incompressible.kinetic_energy,kinetic_energy_ghost,0,number_of_ghost_cells);
        incompressible.boundary->Fill_Ghost_Cells(mac_grid,density,density_ghost,(T)0,0,number_of_ghost_cells);
        if(incompressible.conserve_kinetic_energy) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/kinetic_energy",kinetic_energy_ghost);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities_ghost);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density",density_ghost);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure",incompressible.projection.p);
        if(projection->elliptic_solver->psi_D.Valid_Index(mac_grid.Domain_Indices().max_corner)){
            FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_N",projection->elliptic_solver->psi_N);
            FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",projection->elliptic_solver->psi_D);}}
    else{
        if(incompressible.conserve_kinetic_energy) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/kinetic_energy",incompressible.kinetic_energy);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density",density);}
    if(!stream_type.use_doubles)
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,float>::Write(stream_type,output_directory,frame,rigid_geometry_collection);
    else
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,double>::Write(stream_type,output_directory,frame,rigid_geometry_collection);
#else
        PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif
}
template<class TV> void INCOMPRESSIBLE_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/"+f+"/density",density);
    if(incompressible.use_analytic_energy) FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/"+f+"/analytic_energy",incompressible.analytic_energy);
    std::string filename;
    filename=output_directory+"/"+f+"/mac_velocities";
    if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading mac_velocities "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,face_velocities);LOG::filecout(ss.str());}
    filename=output_directory+"/"+f+"/pressure";
    if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading pressure "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible.projection.p);LOG::filecout(ss.str());}
    ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>* advection_conservative=dynamic_cast<ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>*>(incompressible.advection);
    if(advection_conservative){
        filename=output_directory+"/"+f+"/weights_barjc";
        if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading wjc "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,advection_conservative->sum_jc);LOG::filecout(ss.str());}
        filename=output_directory+"/"+f+"/weights_barjc_cell";
        if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading wjc cell "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,advection_conservative->sum_jc_cell);LOG::filecout(ss.str());}}
    if(!stream_type.use_doubles)
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,float>::Read(stream_type,output_directory,frame,rigid_geometry_collection);
    else
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,double>::Read(stream_type,output_directory,frame,rigid_geometry_collection);
#else
        PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif
}
//#####################################################################
template class INCOMPRESSIBLE_EXAMPLE<VECTOR<float,1> >;
template class INCOMPRESSIBLE_EXAMPLE<VECTOR<float,2> >;
template class INCOMPRESSIBLE_EXAMPLE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INCOMPRESSIBLE_EXAMPLE<VECTOR<double,1> >;
template class INCOMPRESSIBLE_EXAMPLE<VECTOR<double,2> >;
template class INCOMPRESSIBLE_EXAMPLE<VECTOR<double,3> >;
#endif
