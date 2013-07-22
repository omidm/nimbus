//#####################################################################
// Copyright 2004-2007, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/ATTACHMENT_POINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h> // TODO: remove once MUSCLE.cpp exists (windows workaround)
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> ARTICULATED_RIGID_BODY<VECTOR<T,2> >::
ARTICULATED_RIGID_BODY(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
    :BASE(rigid_body_collection_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> ARTICULATED_RIGID_BODY<VECTOR<T,2> >::
~ARTICULATED_RIGID_BODY()
{}
//####################################################################################
// Compute_Position_Based_State
//####################################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,2> >::
Compute_Position_Based_State(const T dt,const T time)
{
    for(int i=1;i<=joint_mesh.joints.m;i++) if(joint_mesh.joints(i)->joint_function) joint_mesh.joints(i)->joint_function->Compute_Desired_PD_Velocity(dt,time);
    if(use_muscle_actuators){
        muscle_list->Initialize_Muscle_Attachments_On_Rigid_Body();
        // Only looks for muscles connecting parent-child across joint
        muscles_crossing_joints.Resize(joint_mesh.joints.m);
        muscle_activations.Resize(muscle_list->muscles.m);ARRAYS_COMPUTATIONS::Fill(muscle_activations,(T)0);
        for(int i=1;i<=joint_mesh.joints.m;i++){
            JOINT<TV>& joint=*joint_mesh.joints(i);
            int parent_id=Parent_Id(joint.id_number),child_id=Child_Id(joint.id_number);
            TV location=joint.Location(*Parent(joint.id_number),*Child(joint.id_number));

            muscles_crossing_joints(i).Remove_All();
            for(int j=1;j<=muscle_list->muscle_attachments_on_rigid_body(parent_id).m;j++)
                if(muscle_list->muscle_attachments_on_rigid_body(parent_id)(j).z->rigid_body.particle_index==child_id){
                    TRIPLE<int,ATTACHMENT_POINT<TV>*,ATTACHMENT_POINT<TV>*>& muscle_attachments=muscle_list->muscle_attachments_on_rigid_body(parent_id)(j);
                    TV direction=(muscle_attachments.z->Embedded_Position()-muscle_attachments.y->Embedded_Position()).Normalized();
                    T moment_arm=TV::Cross_Product(muscle_attachments.y->Embedded_Position()-location,direction).x;
                    muscles_crossing_joints(i).Append(PAIR<int,T>(muscle_attachments.x,moment_arm));
                    {std::stringstream ss;ss<<"joint "<<i<<" muscle "<<j<<" moment arm "<<moment_arm<<std::endl;LOG::filecout(ss.str());}}}
        for(int i=1;i<=muscle_list->muscles.m;i++){muscle_list->muscles(i)->Update_Segments();}}
}
//####################################################################################
// Compute_Target_PD_Angular_Impulse
//####################################################################################
template<class T> VECTOR<T,1> ARTICULATED_RIGID_BODY<VECTOR<T,2> >::
Compute_Target_PD_Angular_Impulse(const JOINT_ID joint_id)
{
    JOINT<TV>& joint=*joint_mesh(joint_id);PHYSBAM_ASSERT(joint.joint_function);
    RIGID_BODY<TV> &parent=*Parent(joint_id),&child=*Child(joint_id);
    if(parent.Has_Infinite_Inertia() && child.Has_Infinite_Inertia()) return VECTOR<T,1>((T)0);
    parent.Update_Angular_Velocity();child.Update_Angular_Velocity();
    // IMPORTANT: joint function's relative velocity is child w.r.t. parent while rigid body uses parent w.r.t. child
    VECTOR<T,1> relative_angular_velocity=RIGID_BODY<TV>::Relative_Angular_Velocity(parent,child),
        relative_angular_velocity_new=-joint.joint_function->desired_angular_velocity;
    MATRIX<T,1> I1_inv_plus_I2_inv;if(!parent.Has_Infinite_Inertia()) I1_inv_plus_I2_inv+=parent.Inertia_Tensor().Inverse();
    if(!child.Has_Infinite_Inertia()) I1_inv_plus_I2_inv+=child.Inertia_Tensor().Inverse();
    return (I1_inv_plus_I2_inv.Inverse())*(relative_angular_velocity_new-relative_angular_velocity);
}
//####################################################################################
// Post_Stabilization_With_Actuation
//####################################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,2> >::
Post_Stabilization_With_Actuation(const JOINT_ID joint_id)
{
    JOINT<TV>& joint=*joint_mesh(joint_id);
    if(joint.joint_function && joint.joint_function->active){
        RIGID_BODY<TV> &parent=*Parent(joint_id),&child=*Child(joint_id);
        if(parent.Has_Infinite_Inertia() && child.Has_Infinite_Inertia()) return;

        VECTOR<T,1> angular_impulse=Compute_Target_PD_Angular_Impulse(joint_id);
        if(use_pd_actuators){
            if(joint.impulse_accumulator) joint.impulse_accumulator->Add_Impulse(TV(),TWIST<TV>(TV(),angular_impulse));
            RIGID_BODY<TV>::Apply_Impulse(parent,child,TV(),TV(),angular_impulse);}
        else if(use_muscle_actuators && angular_impulse.Magnitude()>(T)1e-10){
            T min_activation_penalty=(T).1;
            int active_muscles=muscles_crossing_joints(joint_mesh.Joint_Index_From_Id(joint_id)).m;
            // Create A and b
            MATRIX_MXN<T> A(1,active_muscles);
            for(int i=1;i<=active_muscles;i++){const T& moment_arm=muscles_crossing_joints(joint_mesh.Joint_Index_From_Id(joint_id))(i).y;A(1,i)=moment_arm;}
            VECTOR_ND<T> b(1);b.Set_Subvector(1,angular_impulse);
            // Solve least squares assuming all muscles are active
            MATRIX_MXN<T> A_transpose_A=A.Normal_Equations_Matrix();
            // assume all relative weights are 1
            for(int i=1;i<=active_muscles;i++) A_transpose_A(i,i)+=min_activation_penalty;
            VECTOR_ND<T> A_transpose_b=A.Transpose_Times(b);
            VECTOR_ND<T> impulse_magnitudes=A_transpose_A.Cholesky_Solve(A_transpose_b);
            for(int i=1;i<=active_muscles;i++){int muscle_index=muscles_crossing_joints(joint_mesh.Joint_Index_From_Id(joint_id))(i).x;
                {std::stringstream ss;ss<<"muscle "<<i<<" impulse "<<impulse_magnitudes(i)<<std::endl;LOG::filecout(ss.str());}
                muscle_list->muscles(muscle_index)->Apply_Fixed_Impulse_At_All_Points(impulse_magnitudes(i));
                muscle_activations(muscle_index)+=impulse_magnitudes(i);}}}

    Apply_Poststabilization_To_Joint(joint_id);
}
//####################################################################################
// Solve_Velocities_for_PD
//####################################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,2> >::
Solve_Velocities_for_PD(const T time,const T dt,bool test_system,bool print_matrix)
{
    // don't do extra actuation iterations if not using actuators
    if(!use_pd_actuators && !use_muscle_actuators){Apply_Poststabilization(test_system,print_matrix);return;}
    for(int iteration=1;iteration<=actuation_iterations;iteration++) for(int k=1;k<=poststabilization_iterations;k++)
        for(int i=1;i<=joint_mesh.joints.m;i++) Post_Stabilization_With_Actuation(joint_mesh.joints(i)->id_number);
}
//####################################################################################
// Function Create_Joint_Function
//####################################################################################
template<class T> JOINT_FUNCTION<VECTOR<T,2> >* ARTICULATED_RIGID_BODY<VECTOR<T,2> >::
Create_Joint_Function(const JOINT_ID joint_id)
{
    joint_mesh(joint_id)->Set_Joint_Function(new JOINT_FUNCTION<TV>(*joint_mesh(joint_id),*Parent(joint_id),*Child(joint_id)));
    return joint_mesh(joint_id)->joint_function;
}
//####################################################################################
// Output_Articulation_Points
//####################################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,2> >::
Output_Articulation_Points(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const
{
    if(joint_mesh.joints.m==0) return;
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(STRING_UTILITIES::string_sprintf("%s/%d/arb_info",output_directory.c_str(),frame));
    TYPED_OSTREAM typed_output(*output,stream_type);
    Write_Binary(typed_output,joint_mesh.joints.m*2);
    for(int i=1;i<=joint_mesh.joints.m;i++){
        JOINT<TV>& joint=*joint_mesh.joints(i);
        const RIGID_BODY<TV> &parent=*Parent(joint.id_number),&child=*Child(joint.id_number);
        TV ap1=parent.World_Space_Point(joint.F_pj().t),ap2=child.World_Space_Point(joint.F_cj().t);
        Write_Binary(typed_output,ap1,ap2);}
    delete output;
}
//#####################################################################
template class ARTICULATED_RIGID_BODY<VECTOR<float,2> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ARTICULATED_RIGID_BODY<VECTOR<double,2> >;
#endif
