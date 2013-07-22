//#####################################################################
// Copyright 2009, Michael Lentine, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_FACE_INDEX.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_1D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/FLUID_GRAVITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/INCOMPRESSIBILITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM_FORWARD.h>
#include <PhysBAM_Dynamics/Geometry/GENERAL_GEOMETRY_FORWARD.h>
#include <PhysBAM_Dynamics/PLS_EXAMPLE.h>
using namespace PhysBAM;
//#####################################################################
// PLS_EXAMPLE
//#####################################################################
template<class TV> PLS_EXAMPLE<TV>::
PLS_EXAMPLE(const STREAM_TYPE stream_type_input)
    :stream_type(stream_type_input),initial_time(0),first_frame(0),last_frame(100),frame_rate(24),
    write_substeps_level(-1),write_output_files(true),output_directory("output"),restart(0),
    number_of_ghost_cells(3),cfl((T).9),mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),//incompressible_fluid_collection(mac_grid),
    projection(mac_grid),particle_levelset_evolution(mac_grid,number_of_ghost_cells),incompressible(mac_grid,projection),boundary(0),
    collision_bodies_affecting_fluid(mac_grid)
{
    Initialize_Particles();Initialize_Read_Write_General_Structures();
    incompressible.Set_Custom_Advection(advection_scalar);
    for(int i=1;i<=TV::dimension;i++){domain_boundary(i)(1)=true;domain_boundary(i)(2)=true;}
    domain_boundary(2)(2)=false;
}
//#####################################################################
// ~PLS_EXAMPLE
//#####################################################################
template<class TV> PLS_EXAMPLE<TV>::
~PLS_EXAMPLE()
{
    if(mpi_grid){
        delete boundary;
        delete phi_boundary;}
}
//#####################################################################
// 
//#####################################################################
template<class TV> void PLS_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    if(!write_output_files) return;
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",mac_grid);
    if(mpi_grid) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/global_grid",mpi_grid->global_grid);
    if(PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<TV> >* refinement=dynamic_cast<PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<TV> >*>(&projection)){
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure",refinement->levelset_projection.p);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_N",refinement->levelset_projection.elliptic_solver->psi_N);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",refinement->levelset_projection.elliptic_solver->psi_D);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/coarse_mac_velocities",refinement->coarse_face_velocities);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/coarse_grid",refinement->coarse_grid);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/coarse_psi_N",refinement->elliptic_solver->psi_N);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/coarse_psi_D",refinement->elliptic_solver->psi_D);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/coarse_levelset",refinement->coarse_levelset);}
    else{
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure",incompressible.projection.p);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_N",projection.elliptic_solver->psi_N);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",projection.elliptic_solver->psi_D);}
    T_PARTICLE_LEVELSET& particle_levelset=particle_levelset_evolution.particle_levelset;
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/levelset",particle_levelset.levelset);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"positive_particles"),particle_levelset.positive_particles);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"negative_particles"),particle_levelset.negative_particles);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_positive_particles"),particle_levelset.removed_positive_particles);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_negative_particles"),particle_levelset.removed_negative_particles);
    ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>* advection_conservative=dynamic_cast<ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>*>(incompressible.advection);
    if(advection_conservative){
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/weights_barjc",advection_conservative->sum_jc);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/weights_barjc_cell",advection_conservative->sum_jc_cell);
        if(TV::dimension==2){
            FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/weights_i",advection_conservative->weights_to);
            FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/weights_j",advection_conservative->weights_from);}}
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/"+f+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);
}
template<class TV> void PLS_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    T_PARTICLE_LEVELSET& particle_levelset=particle_levelset_evolution.particle_levelset;
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/"+f+"/levelset",particle_levelset.levelset);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"positive_particles"),particle_levelset.positive_particles);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"negative_particles"),particle_levelset.negative_particles);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_positive_particles"),particle_levelset.removed_positive_particles);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_negative_particles"),particle_levelset.removed_negative_particles);
    FILE_UTILITIES::Read_From_Text_File(output_directory+"/"+f+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);
    ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>* advection_conservative=dynamic_cast<ADVECTION_CONSERVATIVE_UNIFORM<GRID<TV>,T>*>(incompressible.advection);
    if(advection_conservative){
        std::string filename;
        filename=output_directory+"/"+f+"/weights_barjc";
        if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading wjc "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,advection_conservative->sum_jc);LOG::filecout(ss.str());}
        filename=output_directory+"/"+f+"/weights_barjc_cell";
        if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading wjc cell "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,advection_conservative->sum_jc_cell);LOG::filecout(ss.str());}}
    std::string filename;
    filename=output_directory+"/"+f+"/pressure";
    if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading pressure "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible.projection.p);LOG::filecout(ss.str());}
    filename=output_directory+"/"+f+"/mac_velocities";
    if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading mac_velocities "<<filename<<std::endl;FILE_UTILITIES::Read_From_File(stream_type,filename,face_velocities);LOG::filecout(ss.str());}
}
//#####################################################################
template class PLS_EXAMPLE<VECTOR<float,2> >;
template class PLS_EXAMPLE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PLS_EXAMPLE<VECTOR<double,2> >;
template class PLS_EXAMPLE<VECTOR<double,3> >;
#endif
