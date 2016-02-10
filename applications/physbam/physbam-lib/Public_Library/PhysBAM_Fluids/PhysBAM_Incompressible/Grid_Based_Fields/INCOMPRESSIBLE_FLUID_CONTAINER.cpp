//#####################################################################
// Copyright 2011, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Grid_Based_Fields/INCOMPRESSIBLE_FLUID_CONTAINER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>::
INCOMPRESSIBLE_FLUID_CONTAINER(RIGID_GRID_COLLECTION<T_GRID>& rigid_grid_collection,int number_of_ghost_cells_input):
    rigid_grid(rigid_grid_collection),density_container(rigid_grid.grid),temperature_container(rigid_grid.grid),particle_levelset_evolution(rigid_grid.grid,number_of_ghost_cells_input),collision_bodies_affecting_fluid(rigid_grid.grid),number_of_ghost_cells(number_of_ghost_cells_input)
{}
template<class T_GRID> INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>::
INCOMPRESSIBLE_FLUID_CONTAINER(RIGID_GRID<T_GRID> rigid_grid_input)
    :rigid_grid(rigid_grid_input),density_container(rigid_grid_input.grid),temperature_container(rigid_grid_input.grid),particle_levelset_evolution(rigid_grid.grid,3),collision_bodies_affecting_fluid(rigid_grid.grid),number_of_ghost_cells(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>::
~INCOMPRESSIBLE_FLUID_CONTAINER()
{
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>::
Write_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
}
//#####################################################################
// Function Read_Output_Files
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>::
Read_Output_Files(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    std::string filename=output_directory+"/"+f+"/mac_velocities";
    std::string centered_velocity_filename=output_directory+"/"+f+"/centered_velocities";
    if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading mac_velocities "<<filename<<std::endl;LOG::filecout(ss.str());
        FILE_UTILITIES::Read_From_File(stream_type,filename,face_velocities);}
    else if(FILE_UTILITIES::File_Exists(centered_velocity_filename)){
        std::stringstream ss;ss<<"Converting from centered velocities "<<centered_velocity_filename<<" to mac_velocities"<<std::endl;LOG::filecout(ss.str());
        T_ARRAYS_VECTOR centered_velocities;
        FILE_UTILITIES::Read_From_File(stream_type,centered_velocity_filename,centered_velocities);

        TV_INT face_index,first_cell_index,second_cell_index;int axis;
        for(FACE_ITERATOR iterator(rigid_grid.grid);iterator.Valid();iterator.Next()){
            face_index=iterator.Face_Index();axis=iterator.Axis();
            first_cell_index=iterator.First_Cell_Index();second_cell_index=iterator.Second_Cell_Index();
            face_velocities.Component(axis)(face_index)=
                (centered_velocities(first_cell_index)(axis)+centered_velocities(second_cell_index)(axis))*(T).5;}}
}
//#####################################################################
// Function Initialize_Grids
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>::
Initialize_Grids(const bool initialize_phi,const bool initialize_viscosity,const bool initialize_density,const bool initialize_temperature,bool use_ghost_cells)
{
    int allocated_ghost_cells=0;
    if(use_ghost_cells)
        allocated_ghost_cells=number_of_ghost_cells;
    LOG::cout << "allocating container, ghost cells=" << allocated_ghost_cells << std::endl;
    face_velocities.Resize(rigid_grid.grid.Domain_Indices(allocated_ghost_cells));
    psi_N.Resize(rigid_grid.grid.Domain_Indices(allocated_ghost_cells));
    pressure.Resize(rigid_grid.grid.Domain_Indices(allocated_ghost_cells));
    psi_D.Resize(rigid_grid.grid.Domain_Indices(allocated_ghost_cells));
    if(initialize_phi){
        phi.Resize(rigid_grid.grid.Domain_Indices(allocated_ghost_cells));
        particle_levelset_evolution.phi.Resize(rigid_grid.grid.Domain_Indices(allocated_ghost_cells));}
    if(initialize_viscosity) viscosity.Resize(rigid_grid.grid.Domain_Indices(allocated_ghost_cells));
    if(initialize_density) density_container.Initialize_Array(allocated_ghost_cells);
    if(initialize_temperature) temperature_container.Initialize_Array(allocated_ghost_cells);
}
//#####################################################################
// Function Sync_Data
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>::
Sync_Data(INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& fluid_container,THREADED_UNIFORM_GRID<T_GRID>& threaded_grid)
{
    threaded_grid.Sync_Face_Scalar(face_velocities,fluid_container.face_velocities);    
    threaded_grid.Sync_Scalar(viscosity,fluid_container.viscosity);    
}
//#####################################################################
// Function Distribute_Data
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>::
Distribute_Data(INCOMPRESSIBLE_FLUID_CONTAINER<T_GRID>& fluid_container,THREADED_UNIFORM_GRID<T_GRID>& threaded_grid)
{
    threaded_grid.Distribute_Face_Scalar(face_velocities,fluid_container.face_velocities);    
    threaded_grid.Distribute_Scalar(viscosity,fluid_container.viscosity);    
}
//#####################################################################
template class INCOMPRESSIBLE_FLUID_CONTAINER<GRID<VECTOR<float,1> > >;
template class INCOMPRESSIBLE_FLUID_CONTAINER<GRID<VECTOR<float,2> > >;
template class INCOMPRESSIBLE_FLUID_CONTAINER<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INCOMPRESSIBLE_FLUID_CONTAINER<GRID<VECTOR<double,1> > >;
template class INCOMPRESSIBLE_FLUID_CONTAINER<GRID<VECTOR<double,2> > >;
template class INCOMPRESSIBLE_FLUID_CONTAINER<GRID<VECTOR<double,3> > >;
#endif
