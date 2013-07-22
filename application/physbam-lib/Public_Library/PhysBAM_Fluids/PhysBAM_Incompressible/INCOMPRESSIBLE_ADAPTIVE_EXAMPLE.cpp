//#####################################################################
// Copyright 2009, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/INCOMPRESSIBLE_ADAPTIVE_EXAMPLE.h>

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>

using namespace PhysBAM;

//#####################################################################
// INCOMPRESSIBLE_ADAPTIVE_EXAMPLE
//#####################################################################
template<class TV> INCOMPRESSIBLE_ADAPTIVE_EXAMPLE<TV>::
INCOMPRESSIBLE_ADAPTIVE_EXAMPLE(const STREAM_TYPE stream_type_input)
    :stream_type(stream_type_input),initial_time(0),first_frame(0),last_frame(100),frame_rate(24),
    output_directory("output"),number_of_ghost_cells(3),cfl((T).9),mac_grid(TV_INT(),RANGE<TV>::Unit_Box(),true),
    projection(mac_grid),incompressible(mac_grid,projection)
{
    incompressible.Set_Custom_Advection(advection_scalar);
}
//#####################################################################
// ~INCOMPRESSIBLE_ADAPTIVE_EXAMPLE
//#####################################################################
template<class TV> INCOMPRESSIBLE_ADAPTIVE_EXAMPLE<TV>::
~INCOMPRESSIBLE_ADAPTIVE_EXAMPLE()
{}
//#####################################################################
// 
//#####################################################################
template<class TV> void INCOMPRESSIBLE_ADAPTIVE_EXAMPLE<TV>::
Write_Output_Files(const int frame)
{
    //TODO: Fix this
    ARRAY<GRID<TV>,TV_INT> sub_mac_grids;
    ARRAY<ARRAY<T,FACE_INDEX<d> >,TV_INT> sub_arrays;    
    sub_mac_grids.Resize(mac_grid.Domain_Indices());
    sub_arrays.Resize(mac_grid.Domain_Indices(),false);
    for(typename GRID<TV>::CELL_ITERATOR iterator(mac_grid);iterator.Valid();iterator.Next()){
        GRID<TV>& local_mac_grid=*mac_grid.sub_mac_grids(iterator.Cell_Index());
        sub_arrays(iterator.Cell_Index()).Resize(local_mac_grid);
        sub_mac_grids(iterator.Cell_Index()).Initialize(local_mac_grid.counts,RANGE<TV>(local_mac_grid.domain),local_mac_grid.Is_MAC_Grid());
        for(typename GRID<TV>::FACE_ITERATOR local_iterator(local_mac_grid);local_iterator.Valid();local_iterator.Next())
            sub_arrays(iterator.Cell_Index())(local_iterator.Full_Index())=(*face_velocities.sub_arrays(iterator.Cell_Index()))(local_iterator.Full_Index());}

    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",static_cast<GRID<TV> >(mac_grid));
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",static_cast<ARRAY<T,FACE_INDEX<d> > >(face_velocities));
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/sub_mac_velocities",sub_arrays);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/sub_grids",sub_mac_grids);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/sub_density",density);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/sub_temperature",temperature);
}
//#####################################################################
template class INCOMPRESSIBLE_ADAPTIVE_EXAMPLE<VECTOR<float,2> >;
template class INCOMPRESSIBLE_ADAPTIVE_EXAMPLE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INCOMPRESSIBLE_ADAPTIVE_EXAMPLE<VECTOR<double,2> >;
template class INCOMPRESSIBLE_ADAPTIVE_EXAMPLE<VECTOR<double,3> >;
#endif
