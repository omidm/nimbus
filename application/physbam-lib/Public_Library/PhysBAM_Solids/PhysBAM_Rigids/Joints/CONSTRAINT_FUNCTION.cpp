//#####################################################################
// Copyright 2007, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSTRAINT_FUNCTION
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/Robust_Functions.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/CONSTRAINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> LINEAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >::
LINEAR_CONSTRAINT_FUNCTION(const ARTICULATED_RIGID_BODY<TV>& arb,const JOINT_ID joint_id,const T dt_input,const T epsilon_scale_input,TV& location_input)
    :joint(*arb.joint_mesh(joint_id)),dt(dt_input),epsilon_scale(epsilon_scale_input),location(location_input)
{
    rigid_body[0]=arb.Parent(joint_id);rigid_body[1]=arb.Child(joint_id);
    Initialize();
}
//#####################################################################
// Function F
//#####################################################################
template<class T> VECTOR<T,2> LINEAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >::
F(const T_IMPULSE& j) const
{
    return F_Helper(j,0)-F_Helper(j,1)+one_over_m*j+c;
}
//#####################################################################
// Function Jacobian
//#####################################################################
template<class T> MATRIX<T,2> LINEAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >::
Jacobian(const T_IMPULSE& j) const
{
    return Jacobian_Helper(j,0)+Jacobian_Helper(j,1)+one_over_m; // NOTE: the plus sign here is correct
}
//#####################################################################
// Function F_Helper
//#####################################################################
template<class T> VECTOR<T,2> LINEAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >::
F_Helper(const T_IMPULSE& j,const int i) const
{
    T sign=(T)(1-2*i);T_SPIN V_w=dt_angular_velocity[i];if(!rigid_body[i]->Has_Infinite_Inertia()) V_w+=sign*(inverse_inertia_rhat_star[i]*j);
    return ROTATION<TV>::From_Rotation_Vector(V_w).Rotate(r[i]);
}
//#####################################################################
// Function Jacobian_Helper
//#####################################################################
template<class T> MATRIX<T,2> LINEAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >::
Jacobian_Helper(const T_IMPULSE& j,const int i) const
{
    if(rigid_body[i]->Has_Infinite_Inertia()) return MATRIX<T,2>();
    T sign=(T)(1-2*i);
    // rotate by d/dt e^{it} = i e^{it} = e^{i(t+pi/2)}
    return MATRIX<T,2,1>(ROTATION<TV>::From_Rotation_Vector(dt_angular_velocity[i]+(sign*inverse_inertia_rhat_star[i]*j)(1)+(T)half_pi).Rotate(r[i]))*inverse_inertia_rhat_star[i];
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void LINEAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >::
Initialize()
{
    FRAME<TV> projected_frame[2];
    for(int i=0;i<=1;i++){
        p[i]=rigid_body[i]->X();
        dt_angular_velocity[i]=dt*rigid_body[i]->Twist().angular;
        projected_frame[i]=FRAME<TV>(p[i]+dt*rigid_body[i]->Twist().linear,ROTATION<TV>::From_Rotation_Vector(dt_angular_velocity[i])*rigid_body[i]->Rotation());}
    c=projected_frame[0].t-projected_frame[1].t;
    TV desired_translation=(joint.F_jp()*projected_frame[0].Inverse()*projected_frame[1]*joint.F_cj()).t;joint.Constrain_Prismatically(desired_translation);
    TV ap[2]={rigid_body[0]->World_Space_Point(joint.F_pj()*desired_translation),rigid_body[1]->World_Space_Point(joint.F_cj().t)};
    location=(T).5*(ap[0]+ap[1]);
    one_over_m=0;
    for(int i=0;i<=1;i++){
        r[i]=ap[i]-p[i];
        if(!rigid_body[i]->Has_Infinite_Inertia()){
            MATRIX<T,1> inverse_inertia=rigid_body[i]->World_Space_Inertia_Tensor_Inverse();TV rhat=location-p[i];
            one_over_m+=(T)1/rigid_body[i]->Mass();inverse_inertia_rhat_star[i]=inverse_inertia.Times_Cross_Product_Matrix(rhat);
            metric_tensor+=inverse_inertia.Conjugate_With_Cross_Product_Matrix(rhat);}}
    metric_tensor+=one_over_m;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >::
ANGULAR_CONSTRAINT_FUNCTION(const ARTICULATED_RIGID_BODY<TV>& arb,const JOINT_ID joint_id,const T dt_input,const T epsilon_scale_input)
    :joint(*arb.joint_mesh(joint_id)),dt(dt_input),epsilon_scale(epsilon_scale_input)
{
    rigid_body[0]=arb.Parent(joint_id);
    rigid_body[1]=arb.Child(joint_id);
    Initialize();
}
//#####################################################################
// Function F
//#####################################################################
template<class T> VECTOR<T,1> ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >::
F(const T_IMPULSE& j_tau) const
{
    return wrap(F_Helper(j_tau,0)-F_Helper(j_tau,1),T_IMPULSE((T)-pi),T_IMPULSE((T)pi));
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >::
Initialize()
{
    T_SPIN q_projected_body[2];
    for(int i=0;i<=1;i++){
        if(!rigid_body[i]->Has_Infinite_Inertia()) inverse_inertia[i]=rigid_body[i]->World_Space_Inertia_Tensor_Inverse();
        dt_angular_velocity[i]=dt*rigid_body[i]->Twist().angular;
        q_projected_body[i]=dt_angular_velocity[i]+rigid_body[i]->Rotation().Angle();}
    T_SPIN desired_angles=(joint.F_jp().r*joint.F_cj().r).Rotation_Vector()-q_projected_body[0]+q_projected_body[1];joint.Constrain_Angles(desired_angles);
    q_w_old[0]=(rigid_body[0]->Rotation()*joint.F_pj().r).Rotation_Vector()+desired_angles;q_w_old[1]=(rigid_body[1]->Rotation()*joint.F_cj().r).Rotation_Vector();
    metric_tensor=inverse_inertia[0]+inverse_inertia[1];
    length_scale_squared=max(rigid_body[0]->Length_Scale_Squared(),rigid_body[1]->Length_Scale_Squared());
}
//#####################################################################
// Function F_Helper
//#####################################################################
template<class T> typename ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >::T_SPIN ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >::
F_Helper(const T_IMPULSE& j_tau,const int i) const
{
    const T sign=(T)(1-2*i);
    T_SPIN V_w=dt_angular_velocity[i];
    if(!rigid_body[i]->Has_Infinite_Inertia()) V_w+=sign*inverse_inertia[i]*j_tau;
    return V_w+q_w_old[i];
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >::
LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION(const ARTICULATED_RIGID_BODY<TV>& arb,const JOINT_ID joint_id,const T dt_input,const T epsilon_scale_input,TV& location)
    :linear(arb,joint_id,dt_input,epsilon_scale_input,location),angular(arb,joint_id,dt_input,epsilon_scale_input),dt(dt_input),epsilon_scale(epsilon_scale_input)
{
    Initialize();
}
//#####################################################################
// Function F
//#####################################################################
template<class T> VECTOR<T,3> LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >::
F(const T_IMPULSE& j) const
{
    TV jn;j.Get_Subvector(1,jn);T j_tau=j(3);
    TV f_linear[2];T_SPIN f_angular[2];T_CONSTRAINT_ERROR f_of_j;
    for(int i=0;i<=1;i++){
        T_SPIN j_total=TV::Cross_Product(rhat[i],jn)+j_tau;
        f_linear[i]=F_Linear_Helper(j_total,i);
        f_angular[i]=angular.F_Helper(j_total,i);}
    f_of_j.Set_Subvector(1,f_linear[0]-f_linear[1]+linear.one_over_m*jn+linear.c);
    f_of_j.Set_Subvector(3,f_angular[0]-f_angular[1]);
    return f_of_j;
}
//#####################################################################
// Function Jacobian
//#####################################################################
template<class T> MATRIX<T,3> LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >::
Jacobian(const T_IMPULSE& j) const
{
    TV jn;j.Get_Subvector(1,jn);T j_tau=j(3);
    MATRIX<T,3> A;
    A.Add_To_Submatrix(1,1,one_over_m_matrix);
    T_SPIN j_total[2]={TV::Cross_Product(rhat[0],jn)+j_tau,TV::Cross_Product(rhat[1],jn)+j_tau};
    for(int i=0;i<=1;i++){
        MATRIX<T,3,1> A_angular;
        A_angular.Add_To_Submatrix(1,1,Jacobian_Linear_Helper(j_total[i],i));
        A_angular(3,1)=angular.Jacobian_Helper(i).x11;
        A.Add_To_Submatrix(1,3,A_angular);
        A.Add_To_Submatrix(1,1,A_angular*rhat_star[i]);}
    return A;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >::
Initialize()
{
    for(int i=0;i<=1;i++){
        rhat[i]=linear.location-linear.p[i];
        rhat_star[i]=MATRIX<T,1,2>::Cross_Product_Matrix(rhat[i]);}
    one_over_m_matrix=linear.one_over_m*DIAGONAL_MATRIX<T,2>::Identity_Matrix();
    // compute the metric tensor
    MATRIX<T,1,2> cross_term=linear.inverse_inertia_rhat_star[0]+linear.inverse_inertia_rhat_star[1];
    metric_tensor.Add_To_Submatrix(1,1,linear.metric_tensor);
    metric_tensor.Add_To_Submatrix(1,3,cross_term.Transposed());
    metric_tensor.Add_To_Submatrix(3,1,cross_term);
    metric_tensor.Add_To_Submatrix(3,3,angular.metric_tensor);
}
//#####################################################################
// Function F_Linear_Helper
//#####################################################################
template<class T> VECTOR<T,2> LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >::
F_Linear_Helper(const T_SPIN& j,const int i) const
{
    T sign=(T)(1-2*i);T_SPIN V_w=linear.dt_angular_velocity[i];if(!linear.rigid_body[i]->Has_Infinite_Inertia()) V_w+=sign*angular.inverse_inertia[i]*j;
    return ROTATION<TV>::From_Rotation_Vector(V_w).Rotate(linear.r[i]);
}
//#####################################################################
// Function Jacobian_Linear_Helper
//#####################################################################
template<class T> MATRIX<T,2,1> LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >::
Jacobian_Linear_Helper(const T_SPIN& j,const int i) const
{
    if(linear.rigid_body[i]->Has_Infinite_Inertia()) return MATRIX<T,2,1>();
    T sign=(T)(1-2*i);
    return MATRIX<T,2,1>(ROTATION<TV>::From_Rotation_Vector(linear.dt_angular_velocity[i]+sign*angular.inverse_inertia[i]*j+(T)half_pi).Rotate(linear.r[i]))*
        angular.inverse_inertia[i];
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> LINEAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >::
LINEAR_CONSTRAINT_FUNCTION(const ARTICULATED_RIGID_BODY<TV>& arb,const JOINT_ID joint_id,const T dt_input,const T epsilon_scale_input,TV& location_input)
    :joint(*arb.joint_mesh(joint_id)),dt(dt_input),epsilon_scale(epsilon_scale_input),location(location_input)
{
    rigid_body[0]=arb.Parent(joint_id);rigid_body[1]=arb.Child(joint_id);
    Initialize();
}
//#####################################################################
// Function F
//#####################################################################
template<class T> VECTOR<T,3> LINEAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >::
F(const T_IMPULSE& j) const
{
    return F_Helper(j,0)-F_Helper(j,1)+one_over_m*j+c;
}
//#####################################################################
// Function Jacobian
//#####################################################################
template<class T> MATRIX<T,3> LINEAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >::
Jacobian(const T_IMPULSE& j) const
{
    return Jacobian_Helper(j,0)+Jacobian_Helper(j,1)+one_over_m; // NOTE: the plus sign here is correct
}
//#####################################################################
// Function F_Helper
//#####################################################################
template<class T> VECTOR<T,3> LINEAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >::
F_Helper(const T_IMPULSE& j,const int i) const
{
    T sign=(T)(1-2*i);
    T_SPIN V_w=dt_angular_velocity[i];
    if(!rigid_body[i]->Has_Infinite_Inertia()) V_w+=sign*(inverse_inertia_rhat_star[i]*j);
    return ROTATION<TV>::From_Rotation_Vector(V_w).Rotate(r[i]);
}
//#####################################################################
// Function Jacobian_Helper
//#####################################################################
template<class T> MATRIX<T,3> LINEAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >::
Jacobian_Helper(const T_IMPULSE& j,const int i) const
{
    if(rigid_body[i]->Has_Infinite_Inertia()) return MATRIX<T,3>();
    T sign=(T)(1-2*i);T_SPIN v=dt_angular_velocity[i];v+=sign*inverse_inertia_rhat_star[i]*j;
    T theta=v.Normalize(),c_theta=cos(theta),sinc_theta=sinc(theta),s_theta=sinc_theta*theta,r_dot_v=TV::Dot_Product(r[i],v);
    return (MATRIX<T,3>::Outer_Product(((r_dot_v*v-r[i])*s_theta-TV::Cross_Product(r[i],v)*c_theta),v)
        +(MATRIX<T,3>::Cross_Product_Matrix(r[i]*sinc_theta)-(MATRIX<T,3>::Outer_Product(v,r[i])+r_dot_v)*one_minus_cos_x_over_x(theta))
        *(SYMMETRIC_MATRIX<T,3>::Outer_Product(v)-1))*inverse_inertia_rhat_star[i];
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void LINEAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >::
Initialize()
{
    FRAME<TV> projected_frame[2];
    for(int i=0;i<=1;i++){
        p[i]=rigid_body[i]->X();q[i]=rigid_body[i]->Rotation();
        dt_angular_velocity[i]=dt*rigid_body[i]->Twist().angular;
        projected_frame[i]=FRAME<TV>(p[i]+dt*rigid_body[i]->Twist().linear,ROTATION<TV>::From_Rotation_Vector(dt_angular_velocity[i])*q[i]);}
    c=projected_frame[0].t-projected_frame[1].t;
    TV desired_translation=(joint.F_jp()*projected_frame[0].Inverse()*projected_frame[1]*joint.F_cj()).t;joint.Constrain_Prismatically(desired_translation);
    TV ap[2]={rigid_body[0]->World_Space_Point(joint.F_pj()*desired_translation),rigid_body[1]->World_Space_Point(joint.F_cj().t)};
    location=(T).5*(ap[0]+ap[1]);
    one_over_m=0;
    for(int i=0;i<=1;i++){
        r[i]=ap[i]-p[i];
        if(!rigid_body[i]->Has_Infinite_Inertia()){
            one_over_m+=(T)1/rigid_body[i]->Mass();
            inverse_inertia_rhat_star[i]=rigid_body[i]->World_Space_Inertia_Tensor_Inverse().Times_Cross_Product_Matrix(location-p[i]);
            metric_tensor-=inverse_inertia_rhat_star[i].Cross_Product_Matrix_Times(location-p[i]);}}
    metric_tensor+=one_over_m;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >::
ANGULAR_CONSTRAINT_FUNCTION(const ARTICULATED_RIGID_BODY<TV>& arb,const JOINT_ID joint_id,const T dt_input,const T epsilon_scale_input)
    :joint(*arb.joint_mesh(joint_id)),dt(dt_input),epsilon_scale(epsilon_scale_input)
{
    rigid_body[0]=arb.Parent(joint_id);rigid_body[1]=arb.Child(joint_id);
    Initialize();
}
//#####################################################################
// Function F
//#####################################################################
template<class T> VECTOR<T,3> ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >::
F(const T_IMPULSE& j_tau) const
{
    return (F_Helper(j_tau,0)*F_Helper(j_tau,1).Inverse()).Quaternion().v;
}
//#####################################################################
// Function F_Helper
//#####################################################################
template<class T> ROTATION<VECTOR<T,3> > ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >::
F_Helper(const T_IMPULSE& j_tau,const int i) const
{
    const T sign=(T)(1-2*i);T_SPIN V_w=dt_angular_velocity[i];if(!rigid_body[i]->Has_Infinite_Inertia()) V_w+=sign*inverse_inertia[i]*j_tau;
    return ROTATION<TV>::From_Rotation_Vector(V_w)*q_w_old[i];
}
//#####################################################################
// Function Jacobian_Helper
//#####################################################################
template<class T> MATRIX<T,3> ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >::
Jacobian_Helper(const T_IMPULSE& j_tau_i,const T_IMPULSE& j_tau_1_m_i,const int i) const
{
    ROTATION<TV> f=F_Helper(j_tau_1_m_i,1-i);MATRIX<T,3,4> m;m.Set_Column(1,-f.Quaternion().v);m.Set_Submatrix(1,2,MATRIX<T,3>::Cross_Product_Matrix(f.Quaternion().v)+f.Quaternion().s);
    return m*Jacobian_Old_Helper(j_tau_i,i);
}
//#####################################################################
// Function Jacobian_Old_Helper
//#####################################################################
template<class T> MATRIX<T,4,3> ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >::
Jacobian_Old_Helper(const T_IMPULSE& j_tau,const int i) const
{
    const T sign=(T)(1-2*i);T_SPIN V_w=dt_angular_velocity[i];if(!rigid_body[i]->Has_Infinite_Inertia()) V_w+=sign*inverse_inertia[i]*j_tau;
    T_SPIN v_w=V_w;T two_theta_w=v_w.Normalize();MATRIX<T,4,3> A;
    if(!rigid_body[i]->Has_Infinite_Inertia()){
        TV dtheta_w_dj=sign*(T).5*inverse_inertia[i]*v_w;
        T half_sinc_theta=(T).5*sinc((T).5*two_theta_w),s_theta=half_sinc_theta*two_theta_w,c_theta=cos((T).5*two_theta_w);
        MATRIX<T,3> dv_w_dj_times_s_theta=sign*half_sinc_theta*(inverse_inertia[i]-MATRIX<T,3>::Outer_Product(v_w,inverse_inertia[i]*v_w));
        TV scalar=(-q_w_old[i].Quaternion().s*s_theta-c_theta*TV::Dot_Product(q_w_old[i].Quaternion().v,v_w))*dtheta_w_dj-dv_w_dj_times_s_theta.Transpose_Times(q_w_old[i].Quaternion().v);
        MATRIX<T,3> S=MATRIX<T,3>::Outer_Product(-s_theta*q_w_old[i].Quaternion().v+c_theta*modified_b_star[i]*v_w,dtheta_w_dj)+modified_b_star[i]*dv_w_dj_times_s_theta;
        A.Set_Row(1,scalar);A.Set_Submatrix(2,1,S);}
    return A;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >::
Initialize()
{
    ROTATION<TV> q_projected_body[2];
    for(int i=0;i<=1;i++){
        if(!rigid_body[i]->Has_Infinite_Inertia()) inverse_inertia[i]=rigid_body[i]->World_Space_Inertia_Tensor_Inverse();
        dt_angular_velocity[i]=dt*rigid_body[i]->Twist().angular;
        q_projected_body[i]=ROTATION<TV>::From_Rotation_Vector(dt_angular_velocity[i])*rigid_body[i]->Rotation();}
    ROTATION<TV> q_projected_joint=(joint.F_jp().r*q_projected_body[0].Inverse()*q_projected_body[1]*joint.F_cj().r).Normalized();
    TV desired_angles=q_projected_joint.Euler_Angles();
    joint.Constrain_Angles(desired_angles);
    ROTATION<TV> q_target=ROTATION<TV>::From_Euler_Angles(desired_angles).Normalized();
    q_w_old[0]=(rigid_body[0]->Rotation()*joint.F_pj().r*q_target).Normalized();q_w_old[1]=(rigid_body[1]->Rotation()*joint.F_cj().r).Normalized();
    for(int i=0;i<=1;i++) modified_b_star[i]=q_w_old[i].Quaternion().s-MATRIX<T,3>::Cross_Product_Matrix(q_w_old[i].Quaternion().v);
    metric_tensor=inverse_inertia[0]+inverse_inertia[1];
    length_scale_squared=max(rigid_body[0]->Length_Scale_Squared(),rigid_body[1]->Length_Scale_Squared());
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >::
LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION(const ARTICULATED_RIGID_BODY<TV>& arb,const JOINT_ID joint_id,const T dt_input,const T epsilon_scale_input,TV& location)
    :linear(arb,joint_id,dt_input,epsilon_scale_input,location),angular(arb,joint_id,dt_input,epsilon_scale_input),dt(dt_input),epsilon_scale(epsilon_scale_input)
{
    Initialize();
}
//#####################################################################
// Function F
//#####################################################################
template<class T> VECTOR<T,6> LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >::
F(const T_IMPULSE& j) const
{
    TV jn,j_tau;j.Get_Subvector(1,jn);j.Get_Subvector(4,j_tau);
        TV f_linear[2];ROTATION<TV> f_angular[2];T_CONSTRAINT_ERROR f_of_j;
    for(int i=0;i<=1;i++){
        TV j_total=TV::Cross_Product(rhat[i],jn)+j_tau;
        f_linear[i]=F_Linear_Helper(j_total,i);
        f_angular[i]=angular.F_Helper(j_total,i);}
    f_of_j.Set_Subvector(1,f_linear[0]-f_linear[1]+linear.one_over_m*jn+linear.c);
    f_of_j.Set_Subvector(4,(f_angular[0]*f_angular[1].Inverse()).Quaternion().v);
    return f_of_j;
}
//#####################################################################
// Function Jacobian
//#####################################################################
template<class T> MATRIX<T,6> LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >::
Jacobian(const T_IMPULSE& j) const
{
    TV jn,j_tau;j.Get_Subvector(1,jn);j.Get_Subvector(4,j_tau);
    MATRIX<T,6> A;MATRIX<T,6,3> A_angular;
    A.Add_To_Submatrix(1,1,one_over_m_matrix);
    TV j_total[2]={TV::Cross_Product(rhat[0],jn)+j_tau,TV::Cross_Product(rhat[1],jn)+j_tau};
    for(int i=0;i<=1;i++){
        int sign=1-2*i;
        A_angular.Set_Zero_Matrix();
        A_angular.Add_To_Submatrix(1,1,Jacobian_Linear_Helper(j_total[i],i));
        A_angular.Add_To_Submatrix(4,1,angular.Jacobian_Helper(j_total[i],j_total[1-i],i)*(T)sign);
        A.Add_To_Submatrix(1,4,A_angular);
        A.Add_To_Submatrix(1,1,A_angular.Times_Cross_Product_Matrix(rhat[i]));}
    return A;
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >::
Initialize()
{
    for(int i=0;i<=1;i++) rhat[i]=linear.location-linear.p[i];
    one_over_m_matrix=linear.one_over_m*DIAGONAL_MATRIX<T,3>::Identity_Matrix();
    // compute the metric tensor
    MATRIX<T,3> cross_term=linear.inverse_inertia_rhat_star[0]+linear.inverse_inertia_rhat_star[1];
    metric_tensor.Add_To_Submatrix(1,1,linear.metric_tensor);
    metric_tensor.Add_To_Submatrix(1,4,cross_term.Transposed());
    metric_tensor.Add_To_Submatrix(4,1,cross_term);
    metric_tensor.Add_To_Submatrix(4,4,angular.metric_tensor);
}
//#####################################################################
// Function F_Linear_Helper
//#####################################################################
template<class T> VECTOR<T,3> LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >::
F_Linear_Helper(const TV& j,const int i) const
{
    T sign=(T)(1-2*i);T_SPIN V_w=linear.dt_angular_velocity[i];if(!linear.rigid_body[i]->Has_Infinite_Inertia()) V_w+=sign*angular.inverse_inertia[i]*j;
    return ROTATION<TV>::From_Rotation_Vector(V_w).Rotate(linear.r[i]);
}
//#####################################################################
// Function Jacobian_Linear_Helper
//#####################################################################
template<class T> MATRIX<T,3> LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >::
Jacobian_Linear_Helper(const TV& j,const int i) const
{
    if(linear.rigid_body[i]->Has_Infinite_Inertia()) return MATRIX<T,3>();
    T sign=(T)(1-2*i);T_SPIN v=linear.dt_angular_velocity[i];v+=sign*angular.inverse_inertia[i]*j;
    T theta=v.Normalize(),c_theta=cos(theta),sinc_theta=sinc(theta),s_theta=sinc_theta*theta,r_dot_v=TV::Dot_Product(linear.r[i],v);
    return (MATRIX<T,3>::Outer_Product(((r_dot_v*v-linear.r[i])*s_theta-TV::Cross_Product(linear.r[i],v)*c_theta),v)
        +(MATRIX<T,3>::Cross_Product_Matrix(linear.r[i]*sinc_theta)-(MATRIX<T,3>::Outer_Product(v,linear.r[i])+r_dot_v)*one_minus_cos_x_over_x(theta))
        *(SYMMETRIC_MATRIX<T,3>::Outer_Product(v)-1))*angular.inverse_inertia[i];
}
//#####################################################################
template class ANGULAR_CONSTRAINT_FUNCTION<VECTOR<float,2> >;
template class ANGULAR_CONSTRAINT_FUNCTION<VECTOR<float,3> >;
template class LINEAR_CONSTRAINT_FUNCTION<VECTOR<float,2> >;
template class LINEAR_CONSTRAINT_FUNCTION<VECTOR<float,3> >;
template class LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<VECTOR<float,2> >;
template class LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ANGULAR_CONSTRAINT_FUNCTION<VECTOR<double,2> >;
template class ANGULAR_CONSTRAINT_FUNCTION<VECTOR<double,3> >;
template class LINEAR_CONSTRAINT_FUNCTION<VECTOR<double,2> >;
template class LINEAR_CONSTRAINT_FUNCTION<VECTOR<double,3> >;
template class LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<VECTOR<double,2> >;
template class LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<VECTOR<double,3> >;
#endif
