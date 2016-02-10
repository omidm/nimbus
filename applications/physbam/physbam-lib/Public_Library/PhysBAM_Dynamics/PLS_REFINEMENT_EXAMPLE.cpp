//#####################################################################
// Copyright 2009, Michael Lentine, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_1D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/FLUID_GRAVITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/INCOMPRESSIBILITY.h>
#include <PhysBAM_Dynamics/Geometry/GENERAL_GEOMETRY_FORWARD.h>
#include <PhysBAM_Dynamics/PLS_REFINEMENT_EXAMPLE.h>
using namespace PhysBAM;
//#####################################################################
// PLS_REFINEMENT_EXAMPLE
//#####################################################################
template<class TV> PLS_REFINEMENT_EXAMPLE<TV>::
PLS_REFINEMENT_EXAMPLE(const STREAM_TYPE stream_type_input)
    :stream_type(stream_type_input),initial_time(0),first_frame(0),last_frame(100),frame_rate(24),write_debug_data(true),
    output_directory("output"),restart(0),number_of_ghost_cells(3),cfl((T).9),use_collidable_advection(false),gravity((T)9.8),
    fine_mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),coarse_mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),fine_mpi_grid(0),coarse_mpi_grid(0),
    projection(coarse_mac_grid),particle_levelset_evolution(fine_mac_grid,number_of_ghost_cells),incompressible(fine_mac_grid,projection),boundary(0),boundary_coarse(0),phi_boundary(0),
    //projection(coarse_mac_grid,false,false,true,true),particle_levelset_evolution(fine_mac_grid),incompressible(fine_mac_grid,projection),boundary(0),
    rigid_geometry_collection(this),collision_bodies_affecting_fluid(fine_mac_grid)
{
    Initialize_Particles();Initialize_Read_Write_General_Structures();
    incompressible.Set_Custom_Advection(advection_scalar);
    for(int i=1;i<=TV::dimension;i++){non_mpi_boundary(i)(1)=true;non_mpi_boundary(i)(2)=true;}
    for(int i=1;i<=TV::dimension;i++){domain_boundary(i)(1)=true;domain_boundary(i)(2)=true;}
    domain_boundary(2)(2)=false;
}
//#####################################################################
// ~PLS_REFINEMENT_EXAMPLE
//#####################################################################
template<class TV> PLS_REFINEMENT_EXAMPLE<TV>::
~PLS_REFINEMENT_EXAMPLE()
{
    if(fine_mpi_grid) delete boundary;
}
//#####################################################################
// 
//#####################################################################
template<class TV> void PLS_REFINEMENT_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    ARRAY<T,FACE_INDEX<TV::dimension> > fine_face_velocities_ghost(fine_mac_grid,3,false);
    boundary->Fill_Ghost_Cells_Face(fine_mac_grid,fine_face_velocities,fine_face_velocities_ghost,0,3);
    ARRAY<T,FACE_INDEX<TV::dimension> > coarse_face_velocities_ghost(coarse_mac_grid,3,false);
    if(boundary_coarse) boundary_coarse->Fill_Ghost_Cells_Face(coarse_mac_grid,coarse_face_velocities,coarse_face_velocities_ghost,0,3);
    else boundary->Fill_Ghost_Cells_Face(coarse_mac_grid,coarse_face_velocities,coarse_face_velocities_ghost,0,3);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",fine_face_velocities_ghost);
    bool split=split_dir!="",first_frame=(frame==0 || (split && frame==restart));
    if(first_frame) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",fine_mac_grid);
    if(fine_mpi_grid && first_frame) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/global_grid",fine_mpi_grid->global_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/coarse_mac_velocities",coarse_face_velocities_ghost);
    if(first_frame) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/coarse_grid",coarse_mac_grid);
    if(write_debug_data){
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/coarse_pressure",incompressible.projection.p);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/coarse_psi_N",projection.elliptic_solver->psi_N);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/coarse_psi_D",projection.elliptic_solver->psi_D);}
    T_PARTICLE_LEVELSET& particle_levelset=particle_levelset_evolution.particle_levelset;
    T_LEVELSET coarse_levelset(coarse_mac_grid,coarse_phi);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/levelset",particle_levelset.levelset);
    if(write_debug_data){
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/coarse_levelset",coarse_levelset);
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"positive_particles"),particle_levelset.positive_particles);
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"negative_particles"),particle_levelset.negative_particles);
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_positive_particles"),particle_levelset.removed_positive_particles);
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_negative_particles"),particle_levelset.removed_negative_particles);
        FILE_UTILITIES::Write_To_Text_File(output_directory+"/"+f+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);}
}
template<class TV> void PLS_REFINEMENT_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    bool split=split_dir!="";
    std::string my_output_directory=split?split_dir:output_directory;
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    T_PARTICLE_LEVELSET& particle_levelset=particle_levelset_evolution.particle_levelset;
    if(split){
        ARRAY<T,TV_INT> phi_global(fine_mpi_grid->global_grid.Domain_Indices());
        T_LEVELSET levelset_global(fine_mpi_grid->global_grid,phi_global);        
        FILE_UTILITIES::Read_From_File(stream_type,my_output_directory+"/"+f+"/levelset",levelset_global);
        for(typename GRID<TV>::CELL_ITERATOR iterator(fine_mac_grid);iterator.Valid();iterator.Next()){
            particle_levelset.levelset.phi(iterator.Cell_Index())=phi_global(fine_mpi_grid->global_grid.Clamped_Index(iterator.Location()));}}
    else FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/"+f+"/levelset",particle_levelset.levelset);
    if(write_debug_data){
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"positive_particles"),particle_levelset.positive_particles);
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"negative_particles"),particle_levelset.negative_particles);
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_positive_particles"),particle_levelset.removed_positive_particles);
        FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_negative_particles"),particle_levelset.removed_negative_particles);
        FILE_UTILITIES::Read_From_Text_File(output_directory+"/"+f+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);}
    std::string filename;
    filename=my_output_directory+"/"+f+"/coarse_pressure";
    if(FILE_UTILITIES::File_Exists(filename)){
        std::stringstream ss;ss<<"Reading pressure "<<filename<<std::endl;LOG::filecout(ss.str());
        if(split){
            ARRAY<T,TV_INT> p(coarse_mpi_grid->global_grid.Domain_Indices());
            FILE_UTILITIES::Read_From_File(stream_type,filename,p);
            for(typename GRID<TV>::CELL_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
                incompressible.projection.p(iterator.Cell_Index())=p(coarse_mpi_grid->global_grid.Clamped_Index(iterator.Location()));}}
        else FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible.projection.p);}
    filename=my_output_directory+"/"+f+"/mac_velocities";
    if(FILE_UTILITIES::File_Exists(filename)){
        std::stringstream ss;ss<<"Reading mac_velocities "<<filename<<std::endl;LOG::filecout(ss.str());
        if(split){
            ARRAY<T,FACE_INDEX<TV::dimension> > face_vel_global(fine_mpi_grid->global_grid,0,false);
            FILE_UTILITIES::Read_From_File(stream_type,filename,face_vel_global);
            for(typename GRID<TV>::FACE_ITERATOR iterator(fine_mac_grid);iterator.Valid();iterator.Next()){
                TV cell_location=iterator.Location()+fine_mac_grid.DX()/2.*TV::Axis_Vector(iterator.Axis());
                fine_face_velocities(iterator.Full_Index())=face_vel_global(FACE_INDEX<TV::dimension>(iterator.Axis(),fine_mpi_grid->global_grid.Clamped_Face_Index(cell_location)));}}
        else FILE_UTILITIES::Read_From_File(stream_type,filename,fine_face_velocities);}
    filename=my_output_directory+"/"+f+"/coarse_mac_velocities";
    if(FILE_UTILITIES::File_Exists(filename)){
        std::stringstream ss;ss<<"Reading coarse_mac_velocities "<<filename<<std::endl;LOG::filecout(ss.str());
        if(split){
            ARRAY<T,FACE_INDEX<TV::dimension> > face_vel_global(coarse_mpi_grid->global_grid,0,false);
            FILE_UTILITIES::Read_From_File(stream_type,filename,face_vel_global);
            for(typename GRID<TV>::FACE_ITERATOR iterator(coarse_mac_grid);iterator.Valid();iterator.Next()){
                TV cell_location=iterator.Location()+coarse_mac_grid.DX()/2.*TV::Axis_Vector(iterator.Axis());
                coarse_face_velocities(iterator.Full_Index())=face_vel_global(FACE_INDEX<TV::dimension>(iterator.Axis(),coarse_mpi_grid->global_grid.Clamped_Face_Index(cell_location)));}}
        else FILE_UTILITIES::Read_From_File(stream_type,filename,coarse_face_velocities);}
}
//#####################################################################
template class PLS_REFINEMENT_EXAMPLE<VECTOR<float,2> >;
template class PLS_REFINEMENT_EXAMPLE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PLS_REFINEMENT_EXAMPLE<VECTOR<double,2> >;
template class PLS_REFINEMENT_EXAMPLE<VECTOR<double,3> >;
#endif
