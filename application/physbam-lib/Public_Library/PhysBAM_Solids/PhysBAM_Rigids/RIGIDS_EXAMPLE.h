//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_EXAMPLE
//#####################################################################
#ifndef __RIGIDS_EXAMPLE__
#define __RIGIDS_EXAMPLE__

#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Parallel_Computation/MPI_RIGIDS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigids_Evolution/RIGIDS_EVOLUTION_CALLBACKS.h>
namespace PhysBAM{

template<class TV> class RIGIDS_PARAMETERS;
template<class TV> class RIGIDS_EVOLUTION;
template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> class COLLISION_GEOMETRY_COLLECTION;

template<class TV>
class RIGIDS_EXAMPLE:public EXAMPLE<TV>,public RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>,public RIGIDS_EVOLUTION_CALLBACKS<TV>,public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename TV::SPIN T_SPIN;
    typedef EXAMPLE<TV> BASE;
public:
    using BASE::output_directory;using BASE::frame_title;using BASE::stream_type;using BASE::initial_time;using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;
    using BASE::write_first_frame;using BASE::write_last_frame;using BASE::write_time;using BASE::write_substeps_level;using BASE::Set_Write_Substeps_Level;using BASE::data_directory;
    using BASE::write_output_files;using BASE::Write_Frame_Title;

    RIGIDS_PARAMETERS<TV>& rigids_parameters;
    COLLISION_GEOMETRY_COLLECTION<TV>& collision_body_list;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    RIGIDS_EVOLUTION<TV>* rigids_evolution;
    MPI_RIGIDS<TV>* mpi_rigids;
    TV_INT processes_per_dimension;

    RIGIDS_EXAMPLE(const STREAM_TYPE stream_type,const int array_collection_type=0);
    virtual ~RIGIDS_EXAMPLE();

//#####################################################################
    virtual void Post_Initialization(){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Preprocess_Frame(const int frame){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Postprocess_Frame(const int frame){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Preprocess_Substep(const T dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Postprocess_Substep(const T dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();} // time at start of substep
    // solids
    virtual void Initialize_Bodies() {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Read_Output_Files_Solids(const int frame);
    virtual void Log_Parameters() const;
    void Write_Output_Files(const int frame) const;
    template<class T_MPI> void Adjust_Output_Directory_For_MPI(const T_MPI mpi);
//#####################################################################
};
}
#endif
