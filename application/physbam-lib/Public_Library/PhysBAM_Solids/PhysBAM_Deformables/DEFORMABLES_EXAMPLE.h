//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLES_EXAMPLE
//#####################################################################
#ifndef __DEFORMABLES_EXAMPLE__
#define __DEFORMABLES_EXAMPLE__

#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Solids_Geometry_Evolution/RIGID_GEOMETRY_EXAMPLE_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformables_Evolution/DEFORMABLES_EVOLUTION_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES.h>
namespace PhysBAM{

template<class TV> class DEFORMABLES_EVOLUTION;
template<class TV> class DEOFMRABLE_BODY_COLLECTION;
template<class TV> class RIGID_GEOMETRY_COLLECTION;
template<class TV> class COLLISION_GEOMETRY_COLLECTION;
template<class TV> class DEFORMABLES_PARAMETERS;
    
template<class TV>
class DEFORMABLES_EXAMPLE:public EXAMPLE<TV>,public DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES<TV>,public RIGID_GEOMETRY_EXAMPLE_VELOCITIES<TV>,public DEFORMABLES_EVOLUTION_CALLBACKS<TV>,public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef EXAMPLE<TV> BASE;
public:
    using BASE::output_directory;using BASE::frame_title;using BASE::stream_type;using BASE::initial_time;using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;
    using BASE::write_first_frame;using BASE::write_last_frame;using BASE::write_time;using BASE::write_substeps_level;using BASE::Set_Write_Substeps_Level;using BASE::data_directory;
    using BASE::Write_Frame_Title;
    
    DEFORMABLES_PARAMETERS<TV>& deformables_parameters;
    COLLISION_GEOMETRY_COLLECTION<TV>& collision_body_list;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection;
    RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection; //kinematic,static bodies
    DEFORMABLES_EVOLUTION<TV>* deformables_evolution;

    DEFORMABLES_EXAMPLE(const STREAM_TYPE stream_type,const int array_collection_type=0);
    virtual ~DEFORMABLES_EXAMPLE();

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
//#####################################################################
};
}
#endif
