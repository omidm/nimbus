//#####################################################################
// Copyright 2009. 
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_EXAMPLE
//#####################################################################
#ifndef __SOLIDS_EXAMPLE__
#define __SOLIDS_EXAMPLE__

#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION_CALLBACKS.h>
namespace PhysBAM{

template<class TV> class SOLIDS_PARAMETERS;
template<class TV> class SOLIDS_EVOLUTION;
template<class TV> class SOLID_BODY_COLLECTION;

template<class TV>
class SOLIDS_EXAMPLE:public EXAMPLE<TV>,public EXAMPLE_FORCES_AND_VELOCITIES<TV>,public SOLIDS_EVOLUTION_CALLBACKS<TV>,public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef EXAMPLE<TV> BASE;
public:
    using BASE::output_directory;using BASE::frame_title;using BASE::stream_type;using BASE::initial_time;using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;
    using BASE::write_first_frame;using BASE::write_last_frame;using BASE::write_time;using BASE::write_substeps_level;using BASE::Set_Write_Substeps_Level;using BASE::data_directory;
    using BASE::Write_Frame_Title;using EXAMPLE_FORCES_AND_VELOCITIES<TV>::Set_External_Positions;using BASE::parse_args;

    SOLIDS_PARAMETERS<TV>& solids_parameters;
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    SOLIDS_EVOLUTION<TV>* solids_evolution; // defaults to newmark

    SOLIDS_EXAMPLE(const STREAM_TYPE stream_type);
    virtual ~SOLIDS_EXAMPLE();

//#####################################################################
    virtual void Post_Initialization(){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Preprocess_Frame(const int frame){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Postprocess_Frame(const int frame){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Preprocess_Substep(const T dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Postprocess_Substep(const T dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();} // time at start of substep
    virtual void Log_Parameters() const;

    virtual void Initialize_Bodies() {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Read_Output_Files_Solids(const int frame);
    void Write_Output_Files(const int frame) const;
    void Register_Options() PHYSBAM_OVERRIDE=0;
    void Parse_Late_Options() PHYSBAM_OVERRIDE=0;
//#####################################################################
};
}
#endif
