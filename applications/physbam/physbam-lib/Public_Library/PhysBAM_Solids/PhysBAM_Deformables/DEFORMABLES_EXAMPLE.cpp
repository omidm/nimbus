//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLES_EXAMPLE
//#####################################################################
#include <PhysBAM_Tools/Arrays/EXTERNAL_ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLES_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformables_Evolution/DEFORMABLES_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/DEFORMABLES_EXAMPLE.h>
#include <stdexcept>
using namespace PhysBAM;
namespace PhysBAM{
void Initialize_Deformables_Particles();
void Initialize_Geometry_Particle();
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLES_EXAMPLE<TV>::
DEFORMABLES_EXAMPLE(const STREAM_TYPE stream_type,const int array_collection_type)
    :BASE((Initialize_Geometry_Particle(),Initialize_Deformables_Particles(),stream_type)),deformables_parameters(*new DEFORMABLES_PARAMETERS<TV>),collision_body_list(*new COLLISION_GEOMETRY_COLLECTION<TV>),
    deformable_body_collection(*new DEFORMABLE_BODY_COLLECTION<TV>(this,collision_body_list,array_collection_type?new EXTERNAL_ARRAY_COLLECTION():new ARRAY_COLLECTION())),
    rigid_geometry_collection(*new RIGID_GEOMETRY_COLLECTION<TV>(*new RIGID_GEOMETRY_PARTICLES<TV>(array_collection_type?new EXTERNAL_ARRAY_COLLECTION():new ARRAY_COLLECTION()),this,&collision_body_list)),
    deformables_evolution(new DEFORMABLES_EVOLUTION<TV>(deformables_parameters,deformable_body_collection,rigid_geometry_collection))
{
    Initialize_Read_Write_Structures();
    Set_Write_Substeps_Level(-1);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DEFORMABLES_EXAMPLE<TV>::
~DEFORMABLES_EXAMPLE()
{
    delete deformables_evolution;
    delete &deformable_body_collection;
    delete &rigid_geometry_collection.particles;
    delete &rigid_geometry_collection;
    delete &collision_body_list;
    delete &deformables_parameters;
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
template<class TV> void DEFORMABLES_EXAMPLE<TV>::
Read_Output_Files_Solids(const int frame)
{
    ARRAY<int> *needs_init=0,*needs_destroy=0;
    deformable_body_collection.Read(stream_type,output_directory,frame,frame,false,true);
    if(!stream_type.use_doubles)
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,float>::Read(stream_type,output_directory,frame,rigid_geometry_collection,needs_init,needs_destroy);
    else
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,double>::Read(stream_type,output_directory,frame,rigid_geometry_collection,needs_init,needs_destroy);
#else
        PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif
}
//#####################################################################
// Function Log_Parameters
//#####################################################################
template<class TV> void DEFORMABLES_EXAMPLE<TV>::
Log_Parameters() const
{
    LOG::SCOPE scope("DEFORMABLES_EXAMPLE parameters");
    BASE::Log_Parameters();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void DEFORMABLES_EXAMPLE<TV>::
Write_Output_Files(const int frame) const
{
    FILE_UTILITIES::Create_Directory(output_directory);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    Write_Frame_Title(frame);
    deformable_body_collection.Write(stream_type,output_directory,frame,-1,frame==first_frame,true);
    if(!stream_type.use_doubles)
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,float>::Write(stream_type,output_directory,frame,rigid_geometry_collection);
    else
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,double>::Write(stream_type,output_directory,frame,rigid_geometry_collection);
#else
        PHYSBAM_FATAL_ERROR("Cannot write doubles");
#endif
}
//#####################################################################
template class DEFORMABLES_EXAMPLE<VECTOR<float,1> >;
template class DEFORMABLES_EXAMPLE<VECTOR<float,2> >;
template class DEFORMABLES_EXAMPLE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DEFORMABLES_EXAMPLE<VECTOR<double,1> >;
template class DEFORMABLES_EXAMPLE<VECTOR<double,2> >;
template class DEFORMABLES_EXAMPLE<VECTOR<double,3> >;
#endif
