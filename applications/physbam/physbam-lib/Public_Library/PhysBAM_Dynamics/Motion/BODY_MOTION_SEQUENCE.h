//#####################################################################
// Copyright 2008, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __BODY_MOTION_SEQUENCE__
#define __BODY_MOTION_SEQUENCE__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#include <PhysBAM_Dynamics/Motion/BONE.h>
#include <PhysBAM_Dynamics/Motion/MOTION_SEQUENCE.h>

namespace PhysBAM{

template<class T>
class BODY_MOTION_SEQUENCE:public MOTION_SEQUENCE<T,BONE<T> >
{
    typedef VECTOR<T,3> TV;
public:
    ARRAY<ARRAY<std::string> > bone_hierarchy;
    ARRAY<BONE<T> > base_position;
    ARRAY<BONE<T> > ui_position; //internal use only - this should not be read or writen
    TV bone_axis;
    int saved_frame;

    BODY_MOTION_SEQUENCE()
        :MOTION_SEQUENCE<T,BONE<T> >()
    {}

    void Resize(int start_frame,int end_frame)
    {if(end_frame==0) end_frame=this->trajectories(1).m;
    for(int i=1;i<=this->trajectories.m;i++) for(int j=start_frame;j<=end_frame;j++) this->trajectories(i)(j-start_frame+1)=this->trajectories(i)(j);}

    void Rescale(T scaling_factor)
    {for(int i=1;i<=this->trajectories.m;i++){base_position(i).Rescale(scaling_factor);for(int j=1;j<=this->trajectories(i).m;j++) this->trajectories(i)(j).Rescale(scaling_factor);}}

    FRAME<TV> Inverse_Trans(FRAME<TV>& trans_input)
    {return FRAME<TV>(TV(-1*trans_input.t.x,-1*trans_input.t.y,-1*trans_input.t.z),ROTATION<TV>());}

    void Update_Bone_Length(int bone,int frame)
    {this->trajectories(bone)(frame+1).length=base_position(bone).length;}

    void Update_Trajectories(int bone)
    {for(int i=1;i<=this->trajectories(bone).counts.x;i++) Update_Trajectories(bone,i-1);
    ui_position(bone).targeted_transform=FRAME<TV>();ui_position(bone).targeted_translation=FRAME<TV>();ui_position(bone).targeted_rotation=FRAME<TV>();}

    void Update_Trajectories(int bone,int frame)
    {this->trajectories(bone)(frame+1).targeted_rotation=ui_position(bone).targeted_rotation*this->trajectories(bone)(frame+1).targeted_rotation;
    this->trajectories(bone)(frame+1).targeted_translation=ui_position(bone).targeted_translation*this->trajectories(bone)(frame+1).targeted_translation;}

    void Update_Base_Position()
    {for(int i=1;i<=base_position.m;i++) Update_Base_Position(i);}

    void Update_Base_Position(int bone)
    {base_position(bone).transform=base_position(bone).targeted_transform;
    base_position(bone).targeted_translation=FRAME<TV>();base_position(bone).targeted_rotation=FRAME<TV>();}

    void Update_Transforms_From_Scale(int bone,int frame,T length_change)
    {base_position(bone).targeted_transform.t+=TV((T)0,length_change,(T)0);
    for(int i=1;i<=bone_hierarchy(bone).m;i++) Update_Transforms_From_Scale(this->name_to_track_index.Get(bone_hierarchy(bone)(i)),frame,length_change);}
    
    void Update_Base_Transforms(int bone,int frame,T length_change)
    {for(int i=1;i<=bone_hierarchy(bone).m;i++) if(this->names(bone)!="Root"||bone_hierarchy(bone)(i)=="Spine1"){ //this is a hack so legs are not affected when root is scaled
        Update_Transforms_From_Scale(this->name_to_track_index.Get(bone_hierarchy(bone)(i)),frame,length_change);Update_Targeted_Transforms(bone,frame);}};

    void Update_Base_Transforms(int bone,int frame)
    {Update_Base_Transforms(bone,frame,base_position(bone).length-base_position(bone).length/ui_position(bone).length);}

    void Update_Base_Transforms(int bone)
    {for(int i=1;i<=this->trajectories(bone).m;i++) Update_Base_Transforms(bone,i-1);}

    void Update_Children_Rotation(int bone,int frame,FRAME<TV> &root_trans,FRAME<TV> &transform)
    {base_position(bone).targeted_rotation=FRAME<TV>(root_trans.t,ROTATION<TV>())*transform*Inverse_Trans(root_trans)*base_position(bone).rotation;
    Update_Targeted_Transforms(bone,frame);
    for(int i=1;i<=bone_hierarchy(bone).m;i++) Update_Children_Rotation(this->name_to_track_index.Get(bone_hierarchy(bone)(i)),frame,root_trans,transform);}

    void Update_Children_Rotation(int bone,int frame,FRAME<TV> &transform)
    {for(int i=1;i<=bone_hierarchy(bone).m;i++) 
        Update_Children_Rotation(this->name_to_track_index.Get(bone_hierarchy(bone)(i)),frame,this->trajectories(bone)(frame+1).transform,transform);}

    void Update_Children_Translation_Helper(int bone,int frame,FRAME<TV> &transform)
    {base_position(bone).targeted_translation=transform*base_position(bone).translation;
    Update_Targeted_Transforms(bone,frame);
    for(int i=1;i<=bone_hierarchy(bone).m;i++) Update_Children_Translation_Helper(this->name_to_track_index.Get(bone_hierarchy(bone)(i)),frame,transform);}

    void Update_Children_Translation(int bone,int frame,FRAME<TV> &transform)
    {for(int i=1;i<=bone_hierarchy(bone).m;i++) 
        Update_Children_Translation_Helper(this->name_to_track_index.Get(bone_hierarchy(bone)(i)),frame,transform);}

    void Update_Children_Scale(int bone,int frame,FRAME<TV> &transform)
    {for(int i=1;i<=bone_hierarchy(bone).m;i++) if(this->names(bone)!="Root"||bone_hierarchy(bone)(i)=="Spine1")
        Update_Children_Translation_Helper(this->name_to_track_index.Get(bone_hierarchy(bone)(i)),frame,transform);}

    void Update_Transforms(int bone,int frame)
    {this->trajectories(bone)(frame+1).transform=this->trajectories(bone)(frame+1).targeted_transform;}

    void Update_Targeted_Transforms(int bone,int frame)
    {base_position(bone).targeted_transform=base_position(bone).targeted_translation*base_position(bone).targeted_rotation*base_position(bone).transform;
    ui_position(bone).targeted_transform=ui_position(bone).targeted_translation*ui_position(bone).targeted_rotation;
    this->trajectories(bone)(frame+1).targeted_transform=base_position(bone).targeted_transform*this->trajectories(bone)(frame+1).targeted_translation*
        ui_position(bone).targeted_transform*this->trajectories(bone)(frame+1).targeted_rotation;}

    void Update_Lengths(int bone,int frame)
    {this->trajectories(bone)(frame+1).length=base_position(bone).length*ui_position(bone).length;}
    
    void Update_Targeted_Transforms_Bone(int bone)
    {for(int i=1;i<=this->trajectories(bone).m;i++) Update_Targeted_Transforms(bone,i-1);}
    
    void Update_Targeted_Transforms_Frame(int frame)
    {for(int i=1;i<=this->trajectories.m;i++) Update_Targeted_Transforms(i,frame);}
    
    void Update_Targeted_Transforms()
    {for(int i=1;i<=this->trajectories.m;i++) for(int j=1;j<=this->trajectories(i).m;j++) Update_Targeted_Transforms(i,j-1);}

    void Propogate_Lengths(int frame)
    {for(int i=1;i<=this->trajectories.m;i++){
        for(int j=1;j<=frame;j++) this->trajectories(i)(j).length=this->trajectories(i)(frame+1).length;
        for(int j=frame+2;j<=this->trajectories(i).counts.x;j++) this->trajectories(i)(j).length=this->trajectories(i)(frame+1).length;}}

    void Propogate_Transforms(int frame)
    {for(int i=1;i<=this->trajectories.m;i++){
        for(int j=1;j<=frame;j++)
            this->trajectories(i)(j).targeted_transform=base_position(i).targeted_transform*this->trajectories(i)(j).targeted_translation*this->trajectories(i)(j).targeted_rotation;
        for(int j=frame+2;j<=this->trajectories(i).counts.x;j++)
            this->trajectories(i)(j).targeted_transform=base_position(i).targeted_transform*this->trajectories(i)(j).targeted_translation*this->trajectories(i)(j).targeted_rotation;}}

    void Initialize(const int number_of_tracks,const GRID<TV>& time_grid_input)
    {base_position.Resize(number_of_tracks);ui_position.Resize(number_of_tracks);bone_hierarchy.Resize(number_of_tracks);MOTION_SEQUENCE<T,BONE<T> >::Initialize(number_of_tracks,time_grid_input);}

    void Set_Transform(int bone,int frame,FRAME<TV>& f)
    {this->trajectories(bone)(frame).transform=f;}

    void Set_Targeted(int bone,int frame,FRAME<TV>& f)
    {this->trajectories(bone)(frame).targeted_transform=f;}

//#####################################################################
};
}
#include <PhysBAM_Dynamics/Read_Write/Motion/READ_WRITE_BODY_MOTION_SEQUENCE.h>
#endif
