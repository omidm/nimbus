//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_EXAMPLE
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Read_Write/Particles/READ_WRITE_DEFORMABLES_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Read_Write/Particles/READ_WRITE_RIGIDS_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/SOLIDS_EXAMPLE.h>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLIDS_EXAMPLE<TV>::
SOLIDS_EXAMPLE(const STREAM_TYPE stream_type)
    :BASE((Initialize_Geometry_Particle(),Initialize_Rigids_Particles(),Initialize_Deformables_Particles(),stream_type)),
    solids_parameters(*new SOLIDS_PARAMETERS<TV>),solid_body_collection(*new SOLID_BODY_COLLECTION<TV>(this)),solids_evolution(new NEWMARK_EVOLUTION<TV>(solids_parameters,solid_body_collection))
{
    Initialize_Read_Write_Solids_Structures();
    Set_Write_Substeps_Level(-1);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLIDS_EXAMPLE<TV>::
~SOLIDS_EXAMPLE()
{
    delete solids_evolution;
    delete &solid_body_collection;
    delete &solids_parameters;
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
Read_Output_Files_Solids(const int frame)
{
    solid_body_collection.Read(stream_type,output_directory,frame,frame,solids_parameters.write_static_variables_every_frame,solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies,
        solids_parameters.write_deformable_body,solids_parameters.write_from_every_process);
}
//#####################################################################
// Function Log_Parameters
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
Log_Parameters() const
{
    LOG::SCOPE scope("SOLIDS_EXAMPLE parameters");
    BASE::Log_Parameters();
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
Write_Output_Files(const int frame) const
{
    FILE_UTILITIES::Create_Directory(output_directory);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    FILE_UTILITIES::Create_Directory(output_directory+"/"+f);
    FILE_UTILITIES::Create_Directory(output_directory+"/common");
    Write_Frame_Title(frame);
    solid_body_collection.Write(stream_type,output_directory,frame,first_frame,solids_parameters.write_static_variables_every_frame,
        solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies,solids_parameters.write_deformable_body,solids_parameters.write_from_every_process,
        solids_parameters.triangle_collision_parameters.output_interaction_pairs);
}
//#####################################################################
// Function Register_Options
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
Register_Options()
{
    BASE::Register_Options();
    parse_args->Add_String_Argument("-params","","parameter file");
    parse_args->Add_Double_Argument("-solidscfl",.9,"solids CFL");
    parse_args->Add_Option_Argument("-solidscg","Use CG for time integration");
    parse_args->Add_Option_Argument("-solidscr","Use CONJUGATE_RESIDUAL for time integration");
    parse_args->Add_Option_Argument("-solidssymmqmr","Use SYMMQMR for time integration");
    parse_args->Add_Double_Argument("-rigidcfl",.5,"rigid CFL");
}
//#####################################################################
// Function Parse_Late_Options
//#####################################################################
template<class TV> void SOLIDS_EXAMPLE<TV>::
Parse_Late_Options()
{
    BASE::Parse_Late_Options();
    if(parse_args->Is_Value_Set("-solidscfl")) solids_parameters.cfl=(T)parse_args->Get_Double_Value("-solidscfl");
    if(parse_args->Is_Value_Set("-rigidcfl")) solids_parameters.rigid_body_evolution_parameters.rigid_cfl=(T)parse_args->Get_Double_Value("-rigidcfl");
    if(parse_args->Is_Value_Set("-solidscg")) solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cg;
    if(parse_args->Is_Value_Set("-solidscr")) solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_cr;
    if(parse_args->Is_Value_Set("-solidssymmqmr")) solids_parameters.implicit_solve_parameters.evolution_solver_type=krylov_solver_symmqmr;
}
//#####################################################################
template class SOLIDS_EXAMPLE<VECTOR<float,1> >;
template class SOLIDS_EXAMPLE<VECTOR<float,2> >;
template class SOLIDS_EXAMPLE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SOLIDS_EXAMPLE<VECTOR<double,1> >;
template class SOLIDS_EXAMPLE<VECTOR<double,2> >;
template class SOLIDS_EXAMPLE<VECTOR<double,3> >;
#endif
