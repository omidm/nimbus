//#####################################################################
// Copyright 2004-2007, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Dynamics/Motion/ATTACHMENT_FRAME_CONTROL_SET.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> ATTACHMENT_FRAME_CONTROL_SET<T>::
ATTACHMENT_FRAME_CONTROL_SET(ARRAY<TV>& X_input,ARRAY<ARRAY<int> >& attached_nodes_input,const int jaw_attachment_index_input)
    :coefficient_active(24),X(X_input),muscle_force(0),attached_nodes(attached_nodes_input),jaw_attachment_index(jaw_attachment_index_input),rigidity_penalty_coefficient((T)1000),
    jaw_constraint_penalty_coefficient((T)1000),max_opening_angle((T).19)
{
    Save_Controls();
    ARRAYS_COMPUTATIONS::Fill(coefficient_active,true);
    X_save=X;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> ATTACHMENT_FRAME_CONTROL_SET<T>::
~ATTACHMENT_FRAME_CONTROL_SET()
{
}
//#####################################################################
// Function Size
//#####################################################################
template<class T> int ATTACHMENT_FRAME_CONTROL_SET<T>::
Size() const
{
    return 24;
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> typename FACE_CONTROL_SET<T>::TYPE ATTACHMENT_FRAME_CONTROL_SET<T>::
Type() const
{
    return ATTACHMENT_FRAME;
}
//#####################################################################
// operator ()
//#####################################################################
template<class T> T ATTACHMENT_FRAME_CONTROL_SET<T>::
operator()(const int control_id) const
{
    return (control_id<=12)?cranium_transform(control_id):jaw_transform(control_id-12);
}
//#####################################################################
// operator ()
//#####################################################################
template<class T> T& ATTACHMENT_FRAME_CONTROL_SET<T>::
operator()(const int control_id)
{
    return (control_id<=12)?cranium_transform(control_id):jaw_transform(control_id-12);
}
//#####################################################################
// Function Identity
//#####################################################################
template<class T> T ATTACHMENT_FRAME_CONTROL_SET<T>::
Identity(const int control_id) const
{
    return (control_id<=12)?cranium_transform.Identity(control_id):jaw_transform.Identity(control_id-12);
}
//#####################################################################
// Function Maximal_Controls
//#####################################################################
template<class T> void ATTACHMENT_FRAME_CONTROL_SET<T>::
Maximal_Controls()
{
}
//#####################################################################
// Function Control_Active
//#####################################################################
template<class T> bool ATTACHMENT_FRAME_CONTROL_SET<T>::
Control_Active(const int control_id) const
{
    return coefficient_active(control_id);
}
//#####################################################################
// Function Control_Active
//#####################################################################
template<class T> bool& ATTACHMENT_FRAME_CONTROL_SET<T>::
Control_Active(const int control_id)
{
    return coefficient_active(control_id);
}
//#####################################################################
// Function Positions_Determined_Kinematically
//#####################################################################
template<class T> bool ATTACHMENT_FRAME_CONTROL_SET<T>::
Positions_Determined_Kinematically(const int control_id) const
{
    assert(1<=control_id && control_id<=24);
    return control_id<=12;
}
//#####################################################################
// Function Force_Derivative
//#####################################################################
template<class T> void ATTACHMENT_FRAME_CONTROL_SET<T>::
Force_Derivative(ARRAY<TV>& dFdl,ARRAY<TWIST<TV> >& dFrdl,const int control_id) const
{
    assert(13<=control_id && control_id<=24);
    assert(muscle_force);
    assert(Control_Active(control_id));
    assert(jaw_attachment_index);
    MATRIX<T,3> affine_differential;
    TV translation_differential;
    if(control_id<=21) affine_differential.x[control_id-13]=(T)1;
    else translation_differential[control_id-21]=(T)1;
    ARRAY<TV> dXdl(X.m);
    for(int i=1;i<=attached_nodes(jaw_attachment_index).m;i++)
        dXdl(attached_nodes(jaw_attachment_index)(i))=cranium_transform.affine_transform*(affine_differential*X(attached_nodes(jaw_attachment_index)(i))+translation_differential);
    ARRAYS_COMPUTATIONS::Fill(dFdl,TV());
    ARRAYS_COMPUTATIONS::Fill(dFrdl,TWIST<TV>());
    muscle_force->Add_Velocity_Independent_Forces(dFdl,dFrdl,0);
}
//#####################################################################
// Function Position_Derivative
//#####################################################################
template<class T> void ATTACHMENT_FRAME_CONTROL_SET<T>::
Position_Derivative(ARRAY<TV>& dXdl,const int control_id) const
{
    assert(1<=control_id && control_id<=12);
    MATRIX<T,3> affine_differential,affine_differential_transformed;
    TV translation_differential,translation_differential_transformed;
    if(control_id<=9) affine_differential.x[control_id-1]=(T)1;
    else translation_differential[control_id-9]=(T)1;
    affine_differential_transformed=affine_differential*cranium_transform.affine_transform.Inverse();
    translation_differential_transformed=translation_differential-affine_differential_transformed*cranium_transform.translation;
    for(int i=1;i<=dXdl.m;i++) dXdl(i)=affine_differential_transformed*X(i)+translation_differential_transformed;
    for(int i=1;i<=attached_nodes.m;i++) for(int j=1;j<=attached_nodes(i).m;j++)
        if(i==jaw_attachment_index) affine_differential*(jaw_transform.affine_transform*X_save(attached_nodes(i)(j))+jaw_transform.translation)+translation_differential;
        else dXdl(attached_nodes(i)(j))=affine_differential*X_save(attached_nodes(i)(j))+translation_differential;
}
//#####################################################################
// Function Set_Attachment_Positions
//#####################################################################
template<class T> void ATTACHMENT_FRAME_CONTROL_SET<T>::
Set_Attachment_Positions(ARRAY<TV>&X) const
{
    QUASI_RIGID_TRANSFORM_3D<T> jaw_transform_composite=QUASI_RIGID_TRANSFORM_3D<T>::Composite_Transform(cranium_transform,jaw_transform);
    for(int i=1;i<=attached_nodes.m;i++) for(int j=1;j<=attached_nodes(i).m;j++)
        if(i==jaw_attachment_index) X(attached_nodes(i)(j))=jaw_transform_composite.affine_transform*X_save(attached_nodes(i)(j))+jaw_transform_composite.translation;
        else X(attached_nodes(i)(j))=cranium_transform.affine_transform*X_save(attached_nodes(i)(j))+cranium_transform.translation;
}
//#####################################################################
// Function Save_Controls
//#####################################################################
template<class T> void ATTACHMENT_FRAME_CONTROL_SET<T>::
Save_Controls()
{
    cranium_transform_save=cranium_transform;
    jaw_transform_save=jaw_transform;
}
//#####################################################################
// Function Kinematically_Update_Positions
//#####################################################################
template<class T> void ATTACHMENT_FRAME_CONTROL_SET<T>::
Kinematically_Update_Positions(ARRAY<TV>&X) const
{
    QUASI_RIGID_TRANSFORM_3D<T> cranium_transform_incremental=QUASI_RIGID_TRANSFORM_3D<T>::Incremental_Transform(cranium_transform,cranium_transform_save);
    for(int i=1;i<=X.m;i++) X(i)=cranium_transform_incremental.affine_transform*X(i)+cranium_transform_incremental.translation;
}
//#####################################################################
// Function Kinematically_Update_Jacobian
//#####################################################################
template<class T> void ATTACHMENT_FRAME_CONTROL_SET<T>::
Kinematically_Update_Jacobian(ARRAY<TV>&dX) const
{
    QUASI_RIGID_TRANSFORM_3D<T> cranium_transform_incremental=QUASI_RIGID_TRANSFORM_3D<T>::Incremental_Transform(cranium_transform,cranium_transform_save);
    for(int i=1;i<=X.m;i++) dX(i)=cranium_transform_incremental.affine_transform*dX(i);
}
//#####################################################################
// Function Penalty
//#####################################################################
template<class T> T ATTACHMENT_FRAME_CONTROL_SET<T>::
Penalty() const
{
    T penalty=0;
    penalty+=rigidity_penalty_coefficient*cranium_transform.Rigidity_Penalty();
    penalty+=rigidity_penalty_coefficient*jaw_transform.Rigidity_Penalty();
    penalty+=jaw_constraint_penalty_coefficient*Jaw_Constraint_Penalty();
    return penalty;
}
//#####################################################################
// Function Penalty_Gradient
//#####################################################################
template<class T> T ATTACHMENT_FRAME_CONTROL_SET<T>::
Penalty_Gradient(const int control_id) const
{
    assert(1<=control_id && control_id<=24);
    T penalty_gradient=0;
    if(control_id<=12) penalty_gradient+=rigidity_penalty_coefficient*cranium_transform.Rigidity_Penalty_Gradient(control_id);
    if(control_id>=13) penalty_gradient+=rigidity_penalty_coefficient*jaw_transform.Rigidity_Penalty_Gradient(control_id-12);
    if(control_id>=13) penalty_gradient+=jaw_constraint_penalty_coefficient*Jaw_Constraint_Penalty_Gradient(control_id-12);
    return penalty_gradient;
}
//#####################################################################
// Function Penalty_Hessian
//#####################################################################
template<class T> T ATTACHMENT_FRAME_CONTROL_SET<T>::
Penalty_Hessian(const int control_id1,const int control_id2) const
{
    assert(1<=control_id1 && control_id1<=24 && 1<=control_id2 && control_id2<=24);
    T penalty_hessian=0;
    if(control_id1<=12 && control_id2<=12) penalty_hessian+=rigidity_penalty_coefficient*cranium_transform.Ridigity_Penalty_Hessian_Definite_Part(control_id1,control_id2);
    if(control_id1>=13 && control_id2>=13) penalty_hessian+=rigidity_penalty_coefficient*jaw_transform.Ridigity_Penalty_Hessian_Definite_Part(control_id1-12,control_id2-12);
    if(control_id1>=13 && control_id2>=13) penalty_hessian+=jaw_constraint_penalty_coefficient*Jaw_Constraint_Penalty_Hessian(control_id1-12,control_id2-12);
    return penalty_hessian;
}
//#####################################################################
// Function Jaw_Constraint_Penalty
//#####################################################################
template<class T> T ATTACHMENT_FRAME_CONTROL_SET<T>::
Jaw_Constraint_Penalty() const
{
    T penalty=0;
    penalty+=sqr(TV::Dot_Product(jaw_normal,jaw_transform.affine_transform*jaw_axis));
    penalty+=sqr(TV::Dot_Product(jaw_normal,jaw_transform.affine_transform*jaw_midpoint+jaw_transform.translation-jaw_midpoint));
    penalty+=sqr(TV::Dot_Product(jaw_axis,jaw_transform.affine_transform*jaw_midpoint+jaw_transform.translation-jaw_midpoint));
    T left_sliding_validity_measure=TV::Dot_Product(jaw_front,jaw_transform.affine_transform*(jaw_midpoint-(T).5*jaw_axis_length*jaw_axis)+jaw_transform.translation-jaw_midpoint);
    T right_sliding_validity_measure=TV::Dot_Product(jaw_front,jaw_transform.affine_transform*(jaw_midpoint+(T).5*jaw_axis_length*jaw_axis)+jaw_transform.translation-jaw_midpoint);
    T opening_angle_validity_measure=TV::Dot_Product(jaw_normal,jaw_transform.affine_transform*jaw_front);
    if(left_sliding_validity_measure<0) penalty+=sqr(left_sliding_validity_measure);
    if(left_sliding_validity_measure>jaw_sliding_length) penalty+=sqr(left_sliding_validity_measure-jaw_sliding_length);
    if(right_sliding_validity_measure<0) penalty+=sqr(right_sliding_validity_measure);
    if(right_sliding_validity_measure>jaw_sliding_length) penalty+=sqr(right_sliding_validity_measure-jaw_sliding_length);
    if(opening_angle_validity_measure>0) penalty+=sqr(opening_angle_validity_measure);
    T opening_angle_threshold=-asin(max_opening_angle)*(jaw_transform.affine_transform*jaw_front).Magnitude();
    if(opening_angle_validity_measure<opening_angle_threshold) penalty+=sqr(opening_angle_validity_measure-opening_angle_threshold);
    return penalty;
}
//#####################################################################
// Function Jaw_Constraint_Penalty_Gradient
//#####################################################################
template<class T> T ATTACHMENT_FRAME_CONTROL_SET<T>::
Jaw_Constraint_Penalty_Gradient(const int control_id) const
{
    T penalty_gradient=0;
    MATRIX<T,3> affine_transform_differential;
    if(control_id<=9) affine_transform_differential.x[control_id-1]=(T)1;
    TV translation_differential;
    if(control_id>=10) translation_differential[control_id-9]=(T)1;
    penalty_gradient+=(T)2*TV::Dot_Product(jaw_normal,jaw_transform.affine_transform*jaw_axis)
                          *TV::Dot_Product(jaw_normal,affine_transform_differential*jaw_axis);
    penalty_gradient+=(T)2*TV::Dot_Product(jaw_normal,jaw_transform.affine_transform*jaw_midpoint+jaw_transform.translation-jaw_midpoint)
                          *TV::Dot_Product(jaw_normal,affine_transform_differential*jaw_midpoint+translation_differential);
    penalty_gradient+=(T)2*TV::Dot_Product(jaw_axis,jaw_transform.affine_transform*jaw_midpoint+jaw_transform.translation-jaw_midpoint)
                          *TV::Dot_Product(jaw_axis,affine_transform_differential*jaw_midpoint+translation_differential);
    T left_sliding_validity_measure= TV::Dot_Product(jaw_front,jaw_transform.affine_transform*(jaw_midpoint-(T).5*jaw_axis_length*jaw_axis)+jaw_transform.translation-jaw_midpoint);
    T left_sliding_differential=     TV::Dot_Product(jaw_front,affine_transform_differential*(jaw_midpoint-(T).5*jaw_axis_length*jaw_axis)+translation_differential);
    T right_sliding_validity_measure=TV::Dot_Product(jaw_front,jaw_transform.affine_transform*(jaw_midpoint+(T).5*jaw_axis_length*jaw_axis)+jaw_transform.translation-jaw_midpoint);
    T right_sliding_differential=    TV::Dot_Product(jaw_front,affine_transform_differential*(jaw_midpoint+(T).5*jaw_axis_length*jaw_axis)+translation_differential);
    T opening_angle_validity_measure=TV::Dot_Product(jaw_normal,jaw_transform.affine_transform*jaw_front);
    T opening_angle_differential=    TV::Dot_Product(jaw_normal,affine_transform_differential*jaw_front);
    if(left_sliding_validity_measure<0) penalty_gradient+=(T)2*left_sliding_validity_measure*left_sliding_differential;
    if(left_sliding_validity_measure>jaw_sliding_length) penalty_gradient+=(T)2*(left_sliding_validity_measure-jaw_sliding_length)*left_sliding_differential;
    if(right_sliding_validity_measure<0) penalty_gradient+=(T)2*right_sliding_validity_measure*right_sliding_differential;
    if(right_sliding_validity_measure>jaw_sliding_length) penalty_gradient+=(T)2*(right_sliding_validity_measure-jaw_sliding_length)*right_sliding_differential;
    if(opening_angle_validity_measure>0) penalty_gradient+=(T)2*opening_angle_validity_measure*opening_angle_differential;
    T opening_angle_threshold=-asin(max_opening_angle)*(jaw_transform.affine_transform*jaw_front).Magnitude();
    if(opening_angle_validity_measure<opening_angle_threshold) penalty_gradient+=(T)2*(opening_angle_validity_measure-opening_angle_threshold)*opening_angle_differential;
    return penalty_gradient;
}
//#####################################################################
// Function Jaw_Constraint_Penalty_Hessian
//#####################################################################
template<class T> T ATTACHMENT_FRAME_CONTROL_SET<T>::
Jaw_Constraint_Penalty_Hessian(const int control_id1,const int control_id2) const
{
    T penalty_hessian=0;
    MATRIX<T,3> affine_transform1_differential,affine_transform2_differential;
    TV translation1_differential,translation2_differential;
    if(control_id1<=9) affine_transform1_differential.x[control_id1-1]=(T)1;
    if(control_id2<=9) affine_transform2_differential.x[control_id2-1]=(T)1;
    if(control_id1>=10) translation1_differential[control_id1-9]=(T)1;
    if(control_id2>=10) translation2_differential[control_id2-9]=(T)1;
    penalty_hessian+=(T)2*TV::Dot_Product(jaw_normal,affine_transform1_differential*jaw_axis)
                         *TV::Dot_Product(jaw_normal,affine_transform2_differential*jaw_axis);
    penalty_hessian+=(T)2*TV::Dot_Product(jaw_normal,affine_transform1_differential*jaw_midpoint+translation1_differential)
                         *TV::Dot_Product(jaw_normal,affine_transform2_differential*jaw_midpoint+translation2_differential);
    penalty_hessian+=(T)2*TV::Dot_Product(jaw_axis,affine_transform1_differential*jaw_midpoint+translation1_differential)
                         *TV::Dot_Product(jaw_axis,affine_transform2_differential*jaw_midpoint+translation2_differential);
    T left_sliding_validity_measure= TV::Dot_Product(jaw_front,jaw_transform.affine_transform*(jaw_midpoint-(T).5*jaw_axis_length*jaw_axis)+jaw_transform.translation-jaw_midpoint);
    T left_sliding_differential1=    TV::Dot_Product(jaw_front,affine_transform1_differential*(jaw_midpoint-(T).5*jaw_axis_length*jaw_axis)+translation1_differential);
    T left_sliding_differential2=    TV::Dot_Product(jaw_front,affine_transform2_differential*(jaw_midpoint-(T).5*jaw_axis_length*jaw_axis)+translation2_differential);
    T right_sliding_validity_measure=TV::Dot_Product(jaw_front,jaw_transform.affine_transform*(jaw_midpoint+(T).5*jaw_axis_length*jaw_axis)+jaw_transform.translation-jaw_midpoint);
    T right_sliding_differential1=   TV::Dot_Product(jaw_front,affine_transform1_differential*(jaw_midpoint+(T).5*jaw_axis_length*jaw_axis)+translation1_differential);
    T right_sliding_differential2=   TV::Dot_Product(jaw_front,affine_transform2_differential*(jaw_midpoint+(T).5*jaw_axis_length*jaw_axis)+translation2_differential);
    T opening_angle_validity_measure=TV::Dot_Product(jaw_normal,jaw_transform.affine_transform*jaw_front);
    T opening_angle_differential1=   TV::Dot_Product(jaw_normal,affine_transform1_differential*jaw_front);
    T opening_angle_differential2=   TV::Dot_Product(jaw_normal,affine_transform2_differential*jaw_front);
    if(left_sliding_validity_measure<0 || left_sliding_validity_measure>jaw_sliding_length) penalty_hessian+=(T)2*left_sliding_differential1*left_sliding_differential2;
    if(right_sliding_validity_measure<0 || right_sliding_validity_measure>jaw_sliding_length) penalty_hessian+=(T)2*right_sliding_differential1*right_sliding_differential2;
    T opening_angle_threshold=-asin(max_opening_angle)*(jaw_transform.affine_transform*jaw_front).Magnitude();
    if(opening_angle_validity_measure>0 || opening_angle_validity_measure<opening_angle_threshold) penalty_hessian+=(T)2*opening_angle_differential1*opening_angle_differential2;
    return penalty_hessian;
}
//#####################################################################
// Function Read_Jaw_Joint_From_File
//#####################################################################
template<class T> void ATTACHMENT_FRAME_CONTROL_SET<T>::
Read_Jaw_Joint_From_File(const STREAM_TYPE& stream_type,const std::string& filename)
{
    FILE_UTILITIES::Read_From_File(stream_type,filename,jaw_midpoint,jaw_normal,jaw_axis,jaw_front,jaw_axis_length,jaw_sliding_length);
}
//#####################################################################
// Function Write_Jaw_Joint_To_File
//#####################################################################
template<class T> void ATTACHMENT_FRAME_CONTROL_SET<T>::
Write_Jaw_Joint_To_File(const STREAM_TYPE& stream_type,const std::string& filename) const
{
    FILE_UTILITIES::Write_To_File(stream_type,filename,jaw_midpoint,jaw_normal,jaw_axis,jaw_front,jaw_axis_length,jaw_sliding_length);
}
//#####################################################################
// Function Read_Configuration
//#####################################################################
template<class T> void ATTACHMENT_FRAME_CONTROL_SET<T>::
Read_Configuration(const STREAM_TYPE& stream_type,std::istream& input_stream)
{
    TYPED_ISTREAM typed_input(input_stream,stream_type);
    Read_Binary(typed_input,jaw_midpoint,jaw_normal,jaw_axis,jaw_front,jaw_axis_length,jaw_sliding_length);
}
//#####################################################################
// Function Write_Configuration
//#####################################################################
template<class T> void ATTACHMENT_FRAME_CONTROL_SET<T>::
Write_Configuration(const STREAM_TYPE& stream_type,std::ostream& output_stream) const
{
    TYPED_OSTREAM typed_output(output_stream,stream_type);
    Write_Binary(typed_output,jaw_midpoint,jaw_normal,jaw_axis,jaw_front,jaw_axis_length,jaw_sliding_length);
}
//#####################################################################
// Function Set_Original_Attachment_Configuration
//#####################################################################
template<class T> void ATTACHMENT_FRAME_CONTROL_SET<T>::
Set_Original_Attachment_Configuration(const FRAME<TV>& cranium_frame_input,const FRAME<TV>& jaw_frame_input)
{
    cranium_frame_save=cranium_frame_input;
    jaw_frame_save=jaw_frame_input;
}
//#####################################################################
// Function Project_Parameters_To_Allowable_Range
//#####################################################################
template<class T> void ATTACHMENT_FRAME_CONTROL_SET<T>::
Project_Parameters_To_Allowable_Range(const bool active_controls_only)
{
    if(!active_controls_only || coefficient_active(1)) cranium_transform.Make_Rigid();
    if(!active_controls_only || coefficient_active(13)){
        jaw_transform.Make_Rigid();
        T left_sliding_parameter=TV::Dot_Product(jaw_front,jaw_transform.affine_transform*(jaw_midpoint-(T).5*jaw_axis_length*jaw_axis)+jaw_transform.translation-jaw_midpoint)/jaw_sliding_length;
        T right_sliding_parameter=TV::Dot_Product(jaw_front,jaw_transform.affine_transform*(jaw_midpoint+(T).5*jaw_axis_length*jaw_axis)+jaw_transform.translation-jaw_midpoint)/jaw_sliding_length;
        T opening_angle=-asin(TV::Dot_Product(jaw_normal,(jaw_transform.affine_transform*jaw_front).Normalized()));
        left_sliding_parameter=clamp<T>(left_sliding_parameter,0,(T)1);right_sliding_parameter=clamp<T>(right_sliding_parameter,0,(T)1);opening_angle=clamp<T>(opening_angle,0,max_opening_angle);
        TV left_rotation_point=jaw_midpoint-(T).5*jaw_axis_length*jaw_axis+jaw_sliding_length*left_sliding_parameter*jaw_front;
        TV right_rotation_point=jaw_midpoint+(T).5*jaw_axis_length*jaw_axis+jaw_sliding_length*right_sliding_parameter*jaw_front;
        T in_plane_rotation_angle=asin(TV::Triple_Product(jaw_axis,(left_rotation_point-right_rotation_point).Normalized(),jaw_normal));
        jaw_transform.affine_transform=MATRIX<T,3>::Rotation_Matrix(jaw_normal,-in_plane_rotation_angle)*MATRIX<T,3>::Rotation_Matrix(jaw_axis,-opening_angle);
        jaw_transform.translation=(T).5*(left_rotation_point+right_rotation_point)-jaw_transform.affine_transform*jaw_midpoint;}
}
//#####################################################################
// Function Set_From_Generalized_Coordinates
//#####################################################################
template<class T> void ATTACHMENT_FRAME_CONTROL_SET<T>::
Set_From_Generalized_Coordinates(const T left_sliding_parameter,const T right_sliding_parameter,const T opening_angle)
{
    TV left_rotation_point=jaw_midpoint-(T).5*jaw_axis_length*jaw_axis+jaw_sliding_length*left_sliding_parameter*jaw_front;
    TV right_rotation_point=jaw_midpoint+(T).5*jaw_axis_length*jaw_axis+jaw_sliding_length*right_sliding_parameter*jaw_front;
    T in_plane_rotation_angle=asin(TV::Triple_Product(jaw_axis,(left_rotation_point-right_rotation_point).Normalized(),jaw_normal));
    jaw_transform.affine_transform=MATRIX<T,3>::Rotation_Matrix(jaw_normal,-in_plane_rotation_angle)*MATRIX<T,3>::Rotation_Matrix(jaw_axis,-opening_angle);
    jaw_transform.translation=(T).5*(left_rotation_point+right_rotation_point)-jaw_transform.affine_transform*jaw_midpoint;
}
//#####################################################################
// Function Opening_Angle
//#####################################################################
template<class T> T ATTACHMENT_FRAME_CONTROL_SET<T>::
Opening_Angle() const
{
    return -asin(TV::Dot_Product(jaw_normal,(jaw_transform.affine_transform*jaw_front).Normalized()));
}
//#####################################################################
// Function Increment_Opening_Angle
//#####################################################################
template<class T> void ATTACHMENT_FRAME_CONTROL_SET<T>::
Increment_Opening_Angle(const T angle)
{
    TV rotated_axis=jaw_transform.affine_transform*jaw_axis;
    TV old_midpoint_rotation=jaw_transform.affine_transform*jaw_midpoint;
    MATRIX<T,3> rotation=MATRIX<T,3>::Rotation_Matrix(rotated_axis,angle);
    jaw_transform.affine_transform=rotation*jaw_transform.affine_transform;
    TV new_midpoint_rotation=jaw_transform.affine_transform*jaw_midpoint;
    jaw_transform.translation+=(old_midpoint_rotation-new_midpoint_rotation);
}
//#####################################################################
// Function Left_Condyle_Sliding_Parameter
//#####################################################################
template<class T> T ATTACHMENT_FRAME_CONTROL_SET<T>::
Left_Condyle_Sliding_Parameter() const
{
    return TV::Dot_Product(jaw_front,jaw_transform.affine_transform*(jaw_midpoint-(T).5*jaw_axis_length*jaw_axis)+jaw_transform.translation-jaw_midpoint)/jaw_sliding_length;
}
//#####################################################################
// Function Right_Condyle_Sliding_Parameter
//#####################################################################
template<class T> T ATTACHMENT_FRAME_CONTROL_SET<T>::
Right_Condyle_Sliding_Parameter() const
{
    return TV::Dot_Product(jaw_front,jaw_transform.affine_transform*(jaw_midpoint+(T).5*jaw_axis_length*jaw_axis)+jaw_transform.translation-jaw_midpoint)/jaw_sliding_length;
}
//#####################################################################
// Function Interpolate
//#####################################################################
template<class T> void ATTACHMENT_FRAME_CONTROL_SET<T>::
Interpolate(const T interpolation_fraction)
{
    cranium_transform=QUASI_RIGID_TRANSFORM_3D<T>::Interpolate(cranium_transform_save,cranium_transform,interpolation_fraction);
    jaw_transform=QUASI_RIGID_TRANSFORM_3D<T>::Interpolate(jaw_transform_save,jaw_transform,interpolation_fraction);
}
//#####################################################################
// Function Scale
//#####################################################################
template<class T> void ATTACHMENT_FRAME_CONTROL_SET<T>::
Scale(const T scale)
{
    cranium_transform=QUASI_RIGID_TRANSFORM_3D<T>::Interpolate(QUASI_RIGID_TRANSFORM_3D<T>(),cranium_transform,scale);
    jaw_transform=QUASI_RIGID_TRANSFORM_3D<T>::Interpolate(QUASI_RIGID_TRANSFORM_3D<T>(),jaw_transform,scale);
}
//#####################################################################
// Function Distance
//#####################################################################
template<class T> T ATTACHMENT_FRAME_CONTROL_SET<T>::
Distance(const ARRAY<T>& weights,int base)
{
    return weights(base+1)*QUASI_RIGID_TRANSFORM_3D<T>::Distance(cranium_transform_save,cranium_transform)+weights(base+2)*QUASI_RIGID_TRANSFORM_3D<T>::Distance(jaw_transform_save,jaw_transform);
}
//#####################################################################
// Function Print_Diagnostics
//#####################################################################
template<class T> void ATTACHMENT_FRAME_CONTROL_SET<T>::
Print_Diagnostics(std::ostream& output) const
{
    output<<"Cranium frame control diagnostics"<<std::endl;cranium_transform.Print_Diagnostics(output);
    output<<"Jaw frame control diagnostics"<<std::endl;jaw_transform.Print_Diagnostics(output);
    output<<"Jaw plane deviation angle : "<<asin(TV::Dot_Product(jaw_normal,(jaw_transform.affine_transform*jaw_axis).Normalized()))<<std::endl;
    output<<"Jaw midpoint off-plane deviation : "<<abs(TV::Dot_Product(jaw_normal,jaw_transform.affine_transform*jaw_midpoint+jaw_transform.translation-jaw_midpoint))<<std::endl;
    output<<"Jaw midpoint off-axis deviation : "<<abs(TV::Dot_Product(jaw_axis,jaw_transform.affine_transform*jaw_midpoint+jaw_transform.translation-jaw_midpoint))<<std::endl;
    output<<"Left condyle sliding parameter : "<<Left_Condyle_Sliding_Parameter()<<std::endl;
    output<<"Right condyle sliding parameter : "<<Right_Condyle_Sliding_Parameter()<<std::endl;
    output<<"Opening angle : "<<Opening_Angle()<<std::endl;
    output<<"Jaw constraint penalty at current configuration : "<<Jaw_Constraint_Penalty()<<std::endl;
}
template class ATTACHMENT_FRAME_CONTROL_SET<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ATTACHMENT_FRAME_CONTROL_SET<double>;
#endif
