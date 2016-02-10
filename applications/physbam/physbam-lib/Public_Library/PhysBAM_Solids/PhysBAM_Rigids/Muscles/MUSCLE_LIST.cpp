//#####################################################################
// Copyright 2005-2007, Kevin Der, Eran Guendelman, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/ATTACHMENT_POINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize_Muscle_Attachments_On_Rigid_Body
//#####################################################################
template<class TV> void MUSCLE_LIST<TV>::
Initialize_Muscle_Attachments_On_Rigid_Body()
{
    muscle_attachments_on_rigid_body.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size());
    for(int i(1);i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++) muscle_attachments_on_rigid_body(i).Remove_All();
    for(int i=1;i<=muscles.m;i++){
        ATTACHMENT_POINT<TV>* attachment_1=muscles(i)->attachment_point_1;
        for(int j=1;j<=muscles(i)->via_points.m+1;j++){
            ATTACHMENT_POINT<TV>* attachment_2=j>muscles(i)->via_points.m?muscles(i)->attachment_point_2:muscles(i)->via_points(j);
            if(&attachment_1->rigid_body!=&attachment_2->rigid_body){ // store neighboring attachment points that are on different bodies
                muscle_attachments_on_rigid_body(attachment_1->rigid_body.particle_index).Append(
                   TRIPLE<int,ATTACHMENT_POINT<TV>*,ATTACHMENT_POINT<TV>*>(i,attachment_1,attachment_2));
                muscle_attachments_on_rigid_body(attachment_2->rigid_body.particle_index).Append(
                   TRIPLE<int,ATTACHMENT_POINT<TV>*,ATTACHMENT_POINT<TV>*>(i,attachment_2,attachment_1));}
            attachment_1=attachment_2;}}

#if 0
    for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++) for(int j=1;j<=muscle_attachments_on_rigid_body(i).m;j++){
        TRIPLE<int,T_ATTACHMENT_POINT*,T_ATTACHMENT_POINT*>& value=muscle_attachments_on_rigid_body(i)(j);
        {std::stringstream ss;ss<<"Body "<<i<<", muscle "<<value.x<<", pt1 ("<<value.y->Rigid_Body().id_number<<","<<value.y->object_space_position 
           <<"), pt2 ("<<value.z->Rigid_Body().id_number<<","<< value.z->object_space_position<<")"<<std::endl;LOG::filecout(ss.str());}
    }
#endif
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void MUSCLE_LIST<TV>::
Read(const STREAM_TYPE stream_type,const std::string& directory,const int frame)
{
    Clean_Memory();
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(STRING_UTILITIES::string_sprintf("%s/%d/muscle_states",directory.c_str(),frame));
    TYPED_ISTREAM typed_input=TYPED_ISTREAM(*input,stream_type);
    int num_muscles;Read_Binary(typed_input,num_muscles);//muscles.Resize(num_muscles);
    for(int i=1;i<=num_muscles;i++){
        MUSCLE<TV>* muscle=new MUSCLE<TV>(muscle_force_curve);
        muscle->Read(typed_input,rigid_body_collection);
        muscles.Append(muscle);}
    delete input;
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void MUSCLE_LIST<TV>::
Write(const STREAM_TYPE stream_type,const std::string& directory,const int frame) const
{
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(STRING_UTILITIES::string_sprintf("%s/%d/muscle_states",directory.c_str(),frame));
    TYPED_OSTREAM typed_output=TYPED_OSTREAM(*output,stream_type);
    Write_Binary(typed_output,muscles.m);
    for(int i=1;i<=muscles.m;i++) Write_Binary(typed_output,*muscles(i));
    delete output;
}
//#####################################################################
template class MUSCLE_LIST<VECTOR<float,1> >;
template class MUSCLE_LIST<VECTOR<float,2> >;
template class MUSCLE_LIST<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MUSCLE_LIST<VECTOR<double,1> >;
template class MUSCLE_LIST<VECTOR<double,2> >;
template class MUSCLE_LIST<VECTOR<double,3> >;
#endif
