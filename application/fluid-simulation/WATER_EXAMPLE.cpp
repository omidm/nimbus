//#####################################################################
// Copyright 2009, Michael Lentine, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_1D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/FLUID_GRAVITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Forces/INCOMPRESSIBILITY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM.h>
#include <PhysBAM_Dynamics/Geometry/GENERAL_GEOMETRY_FORWARD.h>
#include "WATER_EXAMPLE.h"
using namespace PhysBAM;
//#####################################################################
// WATER_EXAMPLE
//#####################################################################
template<class TV> WATER_EXAMPLE<TV>::
WATER_EXAMPLE(const STREAM_TYPE stream_type_input,int number_of_threads,int refine)
    :stream_type(stream_type_input),initial_time(0),first_frame(0),last_frame(100),frame_rate(24),
    write_substeps_level(-1),write_output_files(true),output_directory("output"),restart(0),number_of_ghost_cells(3),test_number(1),
    cfl(.9),mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),mpi_grid(0),//incompressible_fluid_collection(mac_grid),
    thread_queue(number_of_threads>1?new THREAD_QUEUE(number_of_threads):0),
    projection(refine>1?*new PROJECTION_FREE_SURFACE_REFINEMENT_UNIFORM<GRID<TV> >(mac_grid,particle_levelset_evolution.particle_levelset.levelset,refine):*new PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >(mac_grid,false,false,false,false,thread_queue)),
    particle_levelset_evolution(mac_grid,number_of_ghost_cells),incompressible(mac_grid,projection),boundary(0),rigid_geometry_collection(this),collision_bodies_affecting_fluid(mac_grid)
{
    Initialize_Particles();Initialize_Read_Write_General_Structures();
    incompressible.Set_Custom_Advection(advection_scalar);
    for(int i=1;i<=TV::dimension;i++){domain_boundary(i)(1)=true;domain_boundary(i)(2)=true;}
    domain_boundary(2)(2)=false;
}
//#####################################################################
// ~WATER_EXAMPLE
//#####################################################################
template<class TV> WATER_EXAMPLE<TV>::
~WATER_EXAMPLE()
{
    delete &projection;
    if(mpi_grid){
        delete boundary;
        delete phi_boundary;}
}
//#####################################################################
// Initialize_Phi
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Initialize_Phi()
{
    ARRAY<T,TV_INT>& phi=particle_levelset_evolution.phi;
    
    if(test_number==1)
    {
    	for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next())
	{
            const TV &X=iterator.Location();
            phi(iterator.Cell_Index())=X.y-(T).25;
	}
    }
    else if(test_number==2)
    {
      for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next())
      {
        const TV &X=iterator.Location();
        phi(iterator.Cell_Index())=min((T)(sqrt((X.x-0.5)*(X.x-0.5) + 0.25*(X.y-0.4)*(X.y-0.4)))-(T)0.16, (T)X.y-(T).35);
      }
    }
    else if(test_number==3)
    {
      for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next())
      {
        const TV &X=iterator.Location();
        phi(iterator.Cell_Index())=min(min(
		(T)(sqrt((X.x-0.2)*(X.x-0.2) + 0.30*(X.y-0.4)*(X.y-0.4)))-(T)0.12,
		(T)X.y-(T).10),
		(T)(sqrt((X.x-0.2)*(X.x-0.2) + 0.30*(X.y-0.85)*(X.y-0.85)))-(T)0.12);
      }
    }
    else if(test_number==4)
    {
      for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next())
      {
        const TV &X=iterator.Location();
        phi(iterator.Cell_Index())=min(min(
		(T)(sqrt((X.x-0.5)*(X.x-0.5) + 0.25*(X.y-0.4)*(X.y-0.4)))-(T)0.16,
		(T)X.y-(T).65),
		(T)(-0.5*X.y+(X.x-0.5)*(X.x-0.5)+(T).30));
      }
    }
}
//#####################################################################
// Initialize_Phi
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Initialize_Grid(TV_INT counts,RANGE<TV> domain)
{
    mac_grid.Initialize(counts,domain,true);
}
//#####################################################################
// Set_Boundary_Conditions
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Set_Boundary_Conditions(const T time)
{
    projection.elliptic_solver->psi_D.Fill(false);projection.elliptic_solver->psi_N.Fill(false);
    for(int axis=1;axis<=TV::dimension;axis++) for(int axis_side=1;axis_side<=2;axis_side++){int side=2*(axis-1)+axis_side;
        TV_INT interior_cell_offset=axis_side==1?TV_INT():-TV_INT::Axis_Vector(axis);
        TV_INT exterior_cell_offset=axis_side==1?-TV_INT::Axis_Vector(axis):TV_INT();
        TV_INT boundary_face_offset=axis_side==1?TV_INT::Axis_Vector(axis):-TV_INT::Axis_Vector(axis);
        if(domain_boundary(axis)(axis_side)){
            for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){
                TV_INT face=iterator.Face_Index()+boundary_face_offset;
                if(particle_levelset_evolution.phi(face+interior_cell_offset)<=0){
                    if(face_velocities.Component(axis).Valid_Index(face)){projection.elliptic_solver->psi_N.Component(axis)(face)=true;face_velocities.Component(axis)(face)=0;}}
                else{TV_INT cell=face+exterior_cell_offset;projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;}}}
        else for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Face_Index()+interior_cell_offset;
            projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;}}
    for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
        for(int i=1;i<=sources.m;i++){
            if(time<=3 && sources(i)->Lazy_Inside(iterator.Location())){
                projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
                if((TV::dimension==2 && iterator.Axis()==1)|| (TV::dimension==3 && iterator.Axis()==3)) face_velocities(iterator.Full_Index())=-1;
                else face_velocities(iterator.Full_Index())=0;}}
        for(int i=1;i<=rigid_geometry_collection.particles.array_collection->Size();i++){
            if(rigid_geometry_collection.particles.rigid_geometry(i)->Implicit_Geometry_Lazy_Inside(iterator.Location())){
                projection.elliptic_solver->psi_N(iterator.Full_Index())=true;
                face_velocities(iterator.Full_Index())=rigid_geometry_collection.particles.V(i)(iterator.Axis());}}}
}
//#####################################################################
// Adjust_Phi_With_Sources
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Adjust_Phi_With_Sources(const T time)
{
    if(time>3) return;
    for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){TV_INT index=iterator.Cell_Index();
        for(int i=1;i<=sources.m;i++) particle_levelset_evolution.phi(index)=min(particle_levelset_evolution.phi(index),sources(i)->Extended_Phi(iterator.Location()));}
}
//#####################################################################
// Adjust_Phi_With_Objects
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Adjust_Phi_With_Objects(const T time)
{
    T tolerance=(T)9.8/24; // dt*gravity where dt=1/24 is based on the length of a frame
    for(int id=1;id<=rigid_geometry_collection.particles.array_collection->Size();id++){
        for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
            TV_INT index=iterator.Cell_Index();TV location=mac_grid.X(index);
            if(particle_levelset_evolution.phi(index)<0 && rigid_geometry_collection.Rigid_Geometry(id).Implicit_Geometry_Extended_Value(location)<0){
                TV V_fluid;
                for(int i=1;i<=TV::dimension;i++) V_fluid(i)=(face_velocities(FACE_INDEX<TV::dimension>(i,iterator.First_Face_Index(i)))+face_velocities(FACE_INDEX<TV::dimension>(i,iterator.Second_Face_Index(i))))/2.;
                TV V_object=rigid_geometry_collection.Rigid_Geometry(id).Pointwise_Object_Velocity(location); // velocity object should be spatially varying
                TV V_relative=V_fluid-V_object;
                TV normal=rigid_geometry_collection.Rigid_Geometry(id).Implicit_Geometry_Normal(location);
                T VN=TV::Dot_Product(V_relative,normal),magnitude=V_relative.Magnitude();
                if(VN > max(tolerance,(T).1*magnitude)) particle_levelset_evolution.phi(index)=-rigid_geometry_collection.Rigid_Geometry(id).Implicit_Geometry_Extended_Value(location);}}}
}
//#####################################################################
// Extrapolate_Phi_Into_Objects
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Extrapolate_Phi_Into_Objects(const T time)
{
    for(int id=1;id<=rigid_geometry_collection.particles.array_collection->Size();id++){
        ARRAY<T,TV_INT> phi_object(mac_grid.Domain_Indices(3));
        for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next())
            phi_object(iterator.Cell_Index())=-rigid_geometry_collection.Rigid_Geometry(id).Implicit_Geometry_Extended_Value(iterator.Location());
        EXTRAPOLATION_UNIFORM<GRID<TV>,T> extrapolate(mac_grid,phi_object,particle_levelset_evolution.particle_levelset.levelset.phi,3);extrapolate.Set_Band_Width(3);extrapolate.Extrapolate();}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    if(!write_output_files) return;
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",mac_grid);
    if(mpi_grid) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/global_grid",mpi_grid->global_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure",incompressible.projection.p);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_N",projection.elliptic_solver->psi_N);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",projection.elliptic_solver->psi_D);
    T_PARTICLE_LEVELSET& particle_levelset=particle_levelset_evolution.particle_levelset;
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/levelset",particle_levelset.levelset);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"positive_particles"),particle_levelset.positive_particles);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"negative_particles"),particle_levelset.negative_particles);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_positive_particles"),particle_levelset.removed_positive_particles);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_negative_particles"),particle_levelset.removed_negative_particles);
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/"+f+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);
    if(!stream_type.use_doubles)
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,float>::Write(stream_type,output_directory,frame,rigid_geometry_collection);
    else
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,double>::Write(stream_type,output_directory,frame,rigid_geometry_collection);
#else
        PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif
}
//#####################################################################
// Read_Output_Files
//#####################################################################
template<class TV> void WATER_EXAMPLE<TV>::
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
    std::string filename;
    filename=output_directory+"/"+f+"/pressure";
    if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading pressure "<<filename<<std::endl;LOG::filecout(ss.str());FILE_UTILITIES::Read_From_File(stream_type,filename,incompressible.projection.p);}
    filename=output_directory+"/"+f+"/mac_velocities";
    if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading mac_velocities "<<filename<<std::endl;LOG::filecout(ss.str());FILE_UTILITIES::Read_From_File(stream_type,filename,face_velocities);}
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
template class WATER_EXAMPLE<VECTOR<float,2> >;
template class WATER_EXAMPLE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class WATER_EXAMPLE<VECTOR<double,2> >;
template class WATER_EXAMPLE<VECTOR<double,3> >;
#endif
