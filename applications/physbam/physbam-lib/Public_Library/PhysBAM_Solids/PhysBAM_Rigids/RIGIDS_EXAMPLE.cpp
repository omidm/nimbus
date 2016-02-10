//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_EXAMPLE
//#####################################################################
#include <PhysBAM_Tools/Arrays/EXTERNAL_ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_ELEMENT_ID.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Particles/RIGIDS_PARTICLES_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Read_Write/Particles/READ_WRITE_RIGIDS_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigids_Evolution/RIGIDS_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/RIGIDS_EXAMPLE.h>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGIDS_EXAMPLE<TV>::
RIGIDS_EXAMPLE(const STREAM_TYPE stream_type,const int array_collection_type)
    :BASE((Initialize_Geometry_Particle(),Initialize_Rigids_Particles(),stream_type)),rigids_parameters(*new RIGIDS_PARAMETERS<TV>),collision_body_list(*new COLLISION_GEOMETRY_COLLECTION<TV>),
    rigid_body_collection(*new RIGID_BODY_COLLECTION<TV>(this,&collision_body_list,array_collection_type?new EXTERNAL_ARRAY_COLLECTION():new ARRAY_COLLECTION())),
    rigids_evolution(new RIGIDS_EVOLUTION<TV>(rigids_parameters,rigid_body_collection)),mpi_rigids(0)
{
    Initialize_Read_Write_Structures();
    Set_Write_Substeps_Level(-1);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGIDS_EXAMPLE<TV>::
~RIGIDS_EXAMPLE()
{
    delete rigids_evolution;
    delete &rigid_body_collection;
    delete &collision_body_list;
    delete &rigids_parameters;
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
template<class TV> void RIGIDS_EXAMPLE<TV>::
Read_Output_Files_Solids(const int frame)
{
    rigid_body_collection.Read(stream_type,output_directory,frame);
}
//#####################################################################
// Function Log_Parameters
//#####################################################################
template<class TV> void RIGIDS_EXAMPLE<TV>::
Log_Parameters() const
{
    LOG::SCOPE scope("RIGIDS_EXAMPLE parameters");
    BASE::Log_Parameters();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void RIGIDS_EXAMPLE<TV>::
Write_Output_Files(const int frame) const
{
    if(!write_output_files) return;
    FILE_UTILITIES::Create_Directory(output_directory);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    Write_Frame_Title(frame);
    rigid_body_collection.Write(stream_type,output_directory,frame);
    if(mpi_rigids)
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/partition",mpi_rigids->particles_of_partition(PARTITION_ID(mpi_rigids->rank+1)));
}
//#####################################################################
// Function Adjust_Output_Directory_For_MPI
//#####################################################################
template<class TV> template<class T_MPI> void RIGIDS_EXAMPLE<TV>::
Adjust_Output_Directory_For_MPI(const T_MPI mpi)
{
    if(mpi_rigids && mpi_rigids->number_of_processors>1){
        output_directory+=STRING_UTILITIES::string_sprintf("_NP%d/%d",mpi_rigids->number_of_processors,(mpi_rigids->rank+1));
        FILE_UTILITIES::Create_Directory(output_directory);
        FILE_UTILITIES::Create_Directory(output_directory+"/common");
        LOG::Instance()->Copy_Log_To_File(output_directory+"/common/log.txt",false);}
}
//#####################################################################
template class RIGIDS_EXAMPLE<VECTOR<float,1> >;
template class RIGIDS_EXAMPLE<VECTOR<float,2> >;
template class RIGIDS_EXAMPLE<VECTOR<float,3> >;
template void RIGIDS_EXAMPLE<VECTOR<float,1> >::Adjust_Output_Directory_For_MPI<MPI_RIGIDS<VECTOR<float,1> >*>(MPI_RIGIDS<VECTOR<float,1> >*);
template void RIGIDS_EXAMPLE<VECTOR<float,2> >::Adjust_Output_Directory_For_MPI<MPI_RIGIDS<VECTOR<float,2> >*>(MPI_RIGIDS<VECTOR<float,2> >*);
template void RIGIDS_EXAMPLE<VECTOR<float,3> >::Adjust_Output_Directory_For_MPI<MPI_RIGIDS<VECTOR<float,3> >*>(MPI_RIGIDS<VECTOR<float,3> >*);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGIDS_EXAMPLE<VECTOR<double,1> >;
template class RIGIDS_EXAMPLE<VECTOR<double,2> >;
template class RIGIDS_EXAMPLE<VECTOR<double,3> >;
template void RIGIDS_EXAMPLE<VECTOR<double,1> >::Adjust_Output_Directory_For_MPI<MPI_RIGIDS<VECTOR<double,1> >*>(MPI_RIGIDS<VECTOR<double,1> >*);
template void RIGIDS_EXAMPLE<VECTOR<double,2> >::Adjust_Output_Directory_For_MPI<MPI_RIGIDS<VECTOR<double,2> >*>(MPI_RIGIDS<VECTOR<double,2> >*);
template void RIGIDS_EXAMPLE<VECTOR<double,3> >::Adjust_Output_Directory_For_MPI<MPI_RIGIDS<VECTOR<double,3> >*>(MPI_RIGIDS<VECTOR<double,3> >*);
#endif
