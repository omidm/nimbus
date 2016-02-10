/*
 * Minor change from original PhysBAM code "PROJECTION_EXAMPLE.cpp".
 * Author: Hang Qu <quhang@stanford.edu>
 */

#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include "projection_example.h"

using namespace PhysBAM;

template<class TV> PROJECTION_EXAMPLE<TV>::
PROJECTION_EXAMPLE(const STREAM_TYPE stream_type_input,int number_of_threads)
    :stream_type(stream_type_input),initial_time(0),first_frame(0),frame_rate(24),
    restart(0),write_debug_data(false),output_directory("output"),mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),mpi_grid(0),
    thread_queue(number_of_threads>1?new THREAD_QUEUE(number_of_threads):0),projection(mac_grid,false,false,thread_queue),boundary(0)
{
    for(int i=1;i<=TV::dimension;i++){domain_boundary(i)(1)=true;domain_boundary(i)(2)=true;}
}

template<class TV> PROJECTION_EXAMPLE<TV>::
~PROJECTION_EXAMPLE()
{
    if(mpi_grid) delete boundary;
}

template<class TV> void PROJECTION_EXAMPLE<TV>::
Set_Boundary_Conditions(const T time)
{
    projection.elliptic_solver->psi_D.Fill(false);projection.elliptic_solver->psi_N.Fill(false);
    for(int axis=1;axis<=TV::dimension;axis++) for(int axis_side=1;axis_side<=2;axis_side++){
        int side=2*(axis-1)+axis_side;
        if(domain_boundary(axis)(axis_side)){
            TV_INT interior_cell_offset=axis_side==1?TV_INT():-TV_INT::Axis_Vector(axis);    
            for(typename GRID<TV>::FACE_ITERATOR iterator(mac_grid,1,GRID<TV>::BOUNDARY_REGION,side);
		iterator.Valid();
		iterator.Next()) {
	      TV_INT cell=iterator.Face_Index()+interior_cell_offset;
              TV_INT boundary_face=axis_side==1?iterator.Face_Index()+TV_INT::Axis_Vector(axis):iterator.Face_Index()-TV_INT::Axis_Vector(axis);
              projection.elliptic_solver->psi_D(cell)=true;projection.p(cell)=0;
	    }
	}
    }
}

template<class TV> void PROJECTION_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
    if(mpi_grid) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/global_grid",mpi_grid->global_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/grid",mac_grid);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",mac_grid);
    if(write_debug_data){
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure",projection.p);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_N",projection.elliptic_solver->psi_N);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",projection.elliptic_solver->psi_D);}
}

template<class TV> void PROJECTION_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    std::string filename;
    filename=output_directory+"/"+f+"/mac_velocities";
    if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading mac_velocities "<<filename<<std::endl;LOG::filecout(ss.str());FILE_UTILITIES::Read_From_File(stream_type,filename,face_velocities);}
    filename=output_directory+"/"+f+"/pressure";
    if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading pressure "<<filename<<std::endl;LOG::filecout(ss.str());FILE_UTILITIES::Read_From_File(stream_type,filename,projection.p);}
}

template class PROJECTION_EXAMPLE<VECTOR<float,2> >;
template class PROJECTION_EXAMPLE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PROJECTION_EXAMPLE<VECTOR<double,2> >;
template class PROJECTION_EXAMPLE<VECTOR<double,3> >;
#endif
