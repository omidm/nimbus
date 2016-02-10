//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_BODY_MOTION_SEQUENCE
//#####################################################################
#ifndef __READ_WRITE_BODY_MOTION_SEQUENCE__
#define __READ_WRITE_BODY_MOTION_SEQUENCE__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Dynamics/Motion/BODY_MOTION_SEQUENCE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<BODY_MOTION_SEQUENCE<T>,RW>
{
    typedef VECTOR<T,3> TV;
public:
    static void Read(std::istream& input,BODY_MOTION_SEQUENCE<T>& object)
    {Read_Binary<RW>(input,object.bone_hierarchy,object.base_position,object.time_grid,object.trajectories,object.valid,object.names,object.saved_frame);object.Update_Name_Lookup();object.ui_position.Resize(object.base_position.m);
    for(int i=1;i<=object.trajectories.m;i++) for(int j=1;j<=object.trajectories(i).counts.x;j++){
        object.trajectories(i)(j).translation=object.trajectories(i)(j).targeted_translation=FRAME<TV>(object.trajectories(i)(j).transform.t,ROTATION<TV>());
        object.trajectories(i)(j).rotation=object.trajectories(i)(j).targeted_rotation=FRAME<TV>(TV(),object.trajectories(i)(j).transform.r);}}

    static void Write(std::ostream& output,const BODY_MOTION_SEQUENCE<T>& object)
    {Write_Binary<RW>(output,object.bone_hierarchy,object.base_position,object.time_grid,object.trajectories,object.valid,object.names,object.saved_frame);}
};
}

#endif
