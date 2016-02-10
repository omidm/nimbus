//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Eran Guendelman, Craig Schroeder, Tamar Shinar, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_DYNAMIC_LIST.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_UNDIRECTED_GRAPH.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_FRAME.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/NORMAL_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/PRISMATIC_TWIST_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/RIGID_JOINT.h>
using namespace PhysBAM;
//#####################################################################
namespace{
enum JOINT_TYPE {POINT_JOINT_TYPE=1,RIGID_JOINT_TYPE,TYPE_ANGLE_JOINT,PRISMATIC_TWIST_JOINT_TYPE,NORMAL_JOINT_TYPE};
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> JOINT_MESH<TV>::
JOINT_MESH()
    :joints(dynamic_list.Active_List())
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> JOINT_MESH<TV>::
~JOINT_MESH()
{}
//#####################################################################
// Function Add_Joint
//#####################################################################
template<class TV> JOINT_ID JOINT_MESH<TV>::
Add_Joint(JOINT<TV>* new_joint_input)
{
    JOINT_ID id=dynamic_list.Add_Element(new_joint_input);new_joint_input->Set_Id_Number(id);return id;
}
//#####################################################################
// Function Read
//#####################################################################
// assumes static variables are already read in
template<class TV> void JOINT_MESH<TV>::
Read(TYPED_ISTREAM& input,const std::string& directory,const int frame) // assumes static variables are already read in
{
    std::string prefix=STRING_UTILITIES::string_sprintf("%s/%d/joint_mesh_",directory.c_str(),frame);
    ARRAY<JOINT_ID> needs_init;
    if(!input.type.use_doubles)
        Read_Write<DYNAMIC_LIST<JOINT<TV>,JOINT_ID>,float>::Read(dynamic_list,prefix,needs_init);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else
        Read_Write<DYNAMIC_LIST<JOINT<TV>,JOINT_ID>,double>::Read(dynamic_list,prefix,needs_init);
#endif
    for(int i=1;i<=joints.m;i++){JOINT_TYPE joint_type;Read_Binary(input,joint_type);
        delete joints(i);
        ARRAY<JOINT<TV>*>& nonconst_joints=const_cast<ARRAY<JOINT<TV>*>&>(joints);        
        switch(joint_type){
            case POINT_JOINT_TYPE: nonconst_joints(i)=new POINT_JOINT<TV>();break;
            case RIGID_JOINT_TYPE: nonconst_joints(i)=new RIGID_JOINT<TV>();break;
            case TYPE_ANGLE_JOINT: nonconst_joints(i)=new ANGLE_JOINT<TV>();break;
            case PRISMATIC_TWIST_JOINT_TYPE: nonconst_joints(i)=new PRISMATIC_TWIST_JOINT<TV>();break;
            case NORMAL_JOINT_TYPE: nonconst_joints(i)=new NORMAL_JOINT<TV>();break;
            default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Invalid joint type %d",joint_type));}
        Read_Binary(input,*joints(i));}
    Read_Binary(input,undirected_graph);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void JOINT_MESH<TV>::
Write(TYPED_OSTREAM& output,const std::string& directory,const int frame) const
{
    std::string prefix=STRING_UTILITIES::string_sprintf("%s/%d/joint_mesh_",directory.c_str(),frame);
    if(!output.type.use_doubles)
        Read_Write<DYNAMIC_LIST<JOINT<TV>,JOINT_ID>,float>::Write(dynamic_list,prefix);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        Read_Write<DYNAMIC_LIST<JOINT<TV>,JOINT_ID>,double>::Write(dynamic_list,prefix);
#endif
    for(int i=1;i<=joints.m;i++){
        JOINT_TYPE joint_type;
        const std::type_info& type=typeid(*joints(i));
        if(type==typeid(POINT_JOINT<TV>)) joint_type=POINT_JOINT_TYPE;
        else if(type==typeid(RIGID_JOINT<TV>)) joint_type=RIGID_JOINT_TYPE;
        else if(type==typeid(ANGLE_JOINT<TV>)) joint_type=TYPE_ANGLE_JOINT;
        else if(type==typeid(PRISMATIC_TWIST_JOINT<TV>)) joint_type=PRISMATIC_TWIST_JOINT_TYPE;
        else if(type==typeid(NORMAL_JOINT<TV>)) joint_type=NORMAL_JOINT_TYPE;
        else PHYSBAM_NOT_IMPLEMENTED(STRING_UTILITIES::string_sprintf("Don't know how to write out a joint of type %s",type.name()));
        Write_Binary(output,joint_type,*joints(i));}
    Write_Binary(output,undirected_graph);
}
//#####################################################################
template class JOINT_MESH<VECTOR<float,1> >;
template class JOINT_MESH<VECTOR<float,2> >;
template class JOINT_MESH<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class JOINT_MESH<VECTOR<double,1> >;
template class JOINT_MESH<VECTOR<double,2> >;
template class JOINT_MESH<VECTOR<double,3> >;
#endif
