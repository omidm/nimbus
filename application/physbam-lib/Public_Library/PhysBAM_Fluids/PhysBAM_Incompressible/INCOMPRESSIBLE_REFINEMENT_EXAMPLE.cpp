//#####################################################################
// Copyright 2009, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_REFINEMENT_EXAMPLE.h>

using namespace PhysBAM;

//#####################################################################
// INCOMPRESSIBLE_REFINEMENT_EXAMPLE
//#####################################################################
template<class TV> INCOMPRESSIBLE_REFINEMENT_EXAMPLE<TV>::
INCOMPRESSIBLE_REFINEMENT_EXAMPLE(const STREAM_TYPE stream_type_input)
    :stream_type(stream_type_input),initial_time(0),first_frame(0),last_frame(100),frame_rate(24),write_debug_data(true),
    output_directory("output"),use_coarse_forces(false),use_interpolated_vorticity(false),restart(0),kolmogorov(0),number_of_ghost_cells(3),cfl((T).9),
    coarse_mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),fine_mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    projection(coarse_mac_grid,false,false,true,true),incompressible(fine_mac_grid,projection),advection_scalar(0),threaded_advection_scalar(0),
    boundary(0),rigid_geometry_collection(this),thread_queue(0)
{
    for(int i=1;i<=TV::dimension;i++){domain_boundary(i)(1)=true;domain_boundary(i)(2)=true;}
    Initialize_Geometry_Particle();
    Initialize_Read_Write_Structures();
}
//#####################################################################
// ~INCOMPRESSIBLE_REFINEMENT_EXAMPLE
//#####################################################################
template<class TV> INCOMPRESSIBLE_REFINEMENT_EXAMPLE<TV>::
~INCOMPRESSIBLE_REFINEMENT_EXAMPLE()
{
    if(advection_scalar) delete advection_scalar;
    if(threaded_advection_scalar) delete threaded_advection_scalar;
    if(thread_queue) delete thread_queue;    
    if(fine_mpi_grid) delete boundary;
}
//#####################################################################
// 
//#####################################################################
template<class TV> void INCOMPRESSIBLE_REFINEMENT_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    ARRAY<T,FACE_INDEX<TV::dimension> > fine_face_velocities_ghost(fine_mac_grid,3,false);
    boundary->Fill_Ghost_Cells_Face(fine_mac_grid,fine_face_velocities,fine_face_velocities_ghost,0,3);
    ARRAY<T,FACE_INDEX<TV::dimension> > coarse_face_velocities_ghost(coarse_mac_grid,3,false);
    if(boundary_coarse) boundary_coarse->Fill_Ghost_Cells_Face(coarse_mac_grid,coarse_face_velocities,coarse_face_velocities_ghost,0,3);
    else boundary->Fill_Ghost_Cells_Face(coarse_mac_grid,coarse_face_velocities,coarse_face_velocities_ghost,0,3);
    ARRAY<T,TV_INT> density_ghost(fine_mac_grid.Domain_Indices(3));
    boundary->Fill_Ghost_Cells(fine_mac_grid,density,density_ghost,0,0,3);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",fine_face_velocities_ghost);
    bool split=split_dir!="",first_frame=(frame==0 || (split && frame==restart+1));
    if(first_frame) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",fine_mac_grid);
    if(first_frame && fine_mpi_grid) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/global_grid",fine_mpi_grid->global_grid);
    if(first_frame && coarse_mpi_grid) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/global_coarse_grid",coarse_mpi_grid->global_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/coarse_mac_velocities",coarse_face_velocities_ghost);
    if(first_frame) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/coarse_grid",coarse_mac_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density",density_ghost);
    if(write_debug_data){
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/coarse_pressure",incompressible.projection.p);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_N",projection.elliptic_solver->psi_N);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",projection.elliptic_solver->psi_D);}
    if(!stream_type.use_doubles)
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,float>::Write(stream_type,output_directory,frame,rigid_geometry_collection);
    else
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,double>::Write(stream_type,output_directory,frame,rigid_geometry_collection);
#else
        PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif
}
template<class TV> void INCOMPRESSIBLE_REFINEMENT_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    bool split=split_dir!="";
    std::string my_output_directory=split?split_dir:output_directory;
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    if(split){
        ARRAY<T,TV_INT> density_global(fine_mpi_grid->global_grid.Domain_Indices());
        FILE_UTILITIES::Read_From_File(stream_type,my_output_directory+"/"+f+"/density",density_global);
        for(typename GRID<TV>::CELL_ITERATOR iterator(fine_mac_grid);iterator.Valid();iterator.Next()){
            density(iterator.Cell_Index())=density_global(fine_mpi_grid->global_grid.Clamped_Index(iterator.Location()));}}
    else FILE_UTILITIES::Read_From_File(stream_type,my_output_directory+"/"+f+"/density",density);
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
    if(!stream_type.use_doubles)
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,float>::Read(stream_type,my_output_directory,frame,rigid_geometry_collection);
    else
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,double>::Read(stream_type,my_output_directory,frame,rigid_geometry_collection);
#else
        PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif
}
//#####################################################################
template class INCOMPRESSIBLE_REFINEMENT_EXAMPLE<VECTOR<float,2> >;
template class INCOMPRESSIBLE_REFINEMENT_EXAMPLE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INCOMPRESSIBLE_REFINEMENT_EXAMPLE<VECTOR<double,2> >;
template class INCOMPRESSIBLE_REFINEMENT_EXAMPLE<VECTOR<double,3> >;
#endif
