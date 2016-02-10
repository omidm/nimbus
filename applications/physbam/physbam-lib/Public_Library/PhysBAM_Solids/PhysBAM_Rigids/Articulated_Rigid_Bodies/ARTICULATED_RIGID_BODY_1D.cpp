//#####################################################################
// Copyright 2007-2008, Nipun Kwatra, Craig Schroeder, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARTICULATED_RIGID_BODY
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_1D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
ARTICULATED_RIGID_BODY(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
    :rigid_body_collection(rigid_body_collection_input),joint_mesh(*new JOINT_MESH<TV>),muscle_list(new MUSCLE_LIST<TV>(rigid_body_collection)),use_epsilon_scale(false),
    contact_level_iterations(0),shock_propagation_level_iterations(0),actuation_iterations(0),use_shock_propagation(false),do_final_pass(false),use_pd_actuators(false),
    use_muscle_actuators(false),constrain_pd_directions(false),use_krylov_prestab(false),max_iterations(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
~ARTICULATED_RIGID_BODY()
{
    delete muscle_list;
    delete &joint_mesh;
}
//#####################################################################
// Function Parent_Id
//#####################################################################
template<class T> int ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Parent_Id(JOINT_ID joint_id) const
{
    return joint_mesh.Parent_Id(joint_id);
}
//#####################################################################
// Function Child_Id
//#####################################################################
template<class T> int ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Child_Id(JOINT_ID joint_id) const
{
    return joint_mesh.Child_Id(joint_id);
}
//#####################################################################
// Function Parent
//#####################################################################
template<class T> RIGID_BODY<VECTOR<T,1> >* ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Parent(JOINT_ID joint_id)
{
    return &rigid_body_collection.Rigid_Body(Parent_Id(joint_id));
}
//#####################################################################
// Function Parent
//#####################################################################
template<class T> const RIGID_BODY<VECTOR<T,1> >* ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Parent(JOINT_ID joint_id) const
{
    return &rigid_body_collection.Rigid_Body(Parent_Id(joint_id));
}
//#####################################################################
// Function Child
//#####################################################################
template<class T> RIGID_BODY<VECTOR<T,1> >* ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Child(JOINT_ID joint_id)
{
    return &rigid_body_collection.Rigid_Body(Child_Id(joint_id));
}
//#####################################################################
// Function Child
//#####################################################################
template<class T> const RIGID_BODY<VECTOR<T,1> >* ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Child(JOINT_ID joint_id) const
{
    return &rigid_body_collection.Rigid_Body(Child_Id(joint_id));
}
//#####################################################################
// Function Read
//#####################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Read(const STREAM_TYPE stream_type,const std::string& directory,const int frame)
{
    PHYSBAM_WARNING("ARTICULATED_RIGID_BODY<VECTOR<T,1> > not implemented");
}
//#####################################################################
// Function Apply_Poststabilization
//#####################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Apply_Poststabilization(bool test_system,bool print_matrix,const bool target_pd,const bool no_global_post_stabilization_only,const bool angular_damping_only)
{
    PHYSBAM_WARNING("ARTICULATED_RIGID_BODY<VECTOR<T,1> > not implemented");
}
//#####################################################################
// Function Store_Velocities_And_Momenta
//#####################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Store_Velocities_And_Momenta()
{
    PHYSBAM_WARNING("ARTICULATED_RIGID_BODY<VECTOR<T,1> > not implemented");
}
//#####################################################################
// Function Restore_Velocities_And_Momenta
//#####################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Restore_Velocities_And_Momenta()
{
    PHYSBAM_WARNING("ARTICULATED_RIGID_BODY<VECTOR<T,1> > not implemented");
}
//#####################################################################
// Function Compute_Position_Based_State
//#####################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Compute_Position_Based_State(const T dt,const T time)
{
    PHYSBAM_WARNING("ARTICULATED_RIGID_BODY<VECTOR<T,1> > not implemented");
}
//#####################################################################
// Function Apply_Prestabilization_To_Joint
//#####################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Apply_Prestabilization_To_Joint(const JOINT_ID joint_id,const T dt,const T epsilon_scale)
{
    PHYSBAM_WARNING("ARTICULATED_RIGID_BODY<VECTOR<T,1> > not implemented");
}
//#####################################################################
// Function Compute_Desired_PD_Velocity
//#####################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Compute_Desired_PD_Velocity(const T dt,const T time)
{
    PHYSBAM_WARNING("ARTICULATED_RIGID_BODY<VECTOR<T,1> > not implemented");
}
//#####################################################################
// Function Solve_Velocities_for_PD
//#####################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Solve_Velocities_for_PD(const T time,const T dt,bool test_system,bool print_matrix)
{
    PHYSBAM_WARNING("ARTICULATED_RIGID_BODY<VECTOR<T,1> > not implemented");
}
//#####################################################################
// Function Initialize_Poststabilization_Projection
//#####################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Initialize_Poststabilization_Projection()
{
    PHYSBAM_WARNING("ARTICULATED_RIGID_BODY<VECTOR<T,1> > not implemented");
}
//#####################################################################
// Function Poststabilization_Projection
//#####################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Poststabilization_Projection(ARRAY_VIEW<TWIST<TV> > twist,const bool symmetric)
{
    PHYSBAM_WARNING("ARTICULATED_RIGID_BODY<VECTOR<T,1> > not implemented");
}
//#####################################################################
// Function Generate_Process_List_Using_Contact_Graph
//#####################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Generate_Process_List_Using_Contact_Graph(const RIGID_BODY_CONTACT_GRAPH<TV>& contact_graph)
{
    PHYSBAM_WARNING("ARTICULATED_RIGID_BODY<VECTOR<T,1> > not implemented");
}
//#####################################################################
// Function Write
//#####################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Write(const STREAM_TYPE stream_type,const std::string& directory,const int frame) const
{
    PHYSBAM_WARNING("ARTICULATED_RIGID_BODY<VECTOR<T,1> > not implemented");
}
//#####################################################################
// Function Substitute_Joint_Parent_Body
//#####################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Substitute_Joint_Parent_Body(JOINT_ID joint_id,int new_parent,const FRAME<TV>& frame)
{
    PHYSBAM_WARNING("ARTICULATED_RIGID_BODY<VECTOR<T,1> > not implemented");
}
//#####################################################################
// Function Substitute_Joint_Child_Body
//#####################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Substitute_Joint_Child_Body(JOINT_ID joint_id,int new_child,const FRAME<TV>& frame)
{
    PHYSBAM_WARNING("ARTICULATED_RIGID_BODY<VECTOR<T,1> > not implemented");
}
//#####################################################################
// Function Substitute_Joint_Parent_Body
//#####################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Substitute_Joint_Parent_Body(JOINT_ID joint_id,int new_parent)
{
    PHYSBAM_WARNING("ARTICULATED_RIGID_BODY<VECTOR<T,1> > not implemented");
}
//#####################################################################
// Function Substitute_Joint_Child_Body
//#####################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,1> >::
Substitute_Joint_Child_Body(JOINT_ID joint_id,int new_child)
{
    PHYSBAM_WARNING("ARTICULATED_RIGID_BODY<VECTOR<T,1> > not implemented");
}
//####################################################################################
template class ARTICULATED_RIGID_BODY<VECTOR<float,1> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ARTICULATED_RIGID_BODY<VECTOR<double,1> >;
#endif
