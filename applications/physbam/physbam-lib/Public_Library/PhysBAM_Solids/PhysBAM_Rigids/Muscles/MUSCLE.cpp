//#####################################################################
// Copyright 2005-2007, Kevin Der, Eran Guendelman, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/ATTACHMENT_POINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_FORCE_CURVE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MUSCLE<TV>::
MUSCLE(const MUSCLE_FORCE_CURVE<T>& force_curve_input)
    :force_curve(force_curve_input),attachment_point_1(0),attachment_point_2(0),
    optimal_length(0),peak_force(0),pennation_angle(0),tendon_slack_length(0),max_shortening_velocity(0),peak_force_times_cos_pennation_angle(0),
    length_normalization(0),velocity_normalization(0),tendon_normalization(0),one_over_peak_force(0),cos_pennation_angle((T)1),one_over_cos_pennation_angle(1)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MUSCLE<TV>::
~MUSCLE() // assume we own the constrained points!
{
    via_points.Delete_Pointers_And_Clean_Memory();
    muscle_segments.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV> void MUSCLE<TV>::
Clean_Memory()
{
    delete attachment_point_1;
    attachment_point_1=0;
    delete attachment_point_2;
    attachment_point_2=0;
    via_points.Delete_Pointers_And_Clean_Memory();
    muscle_segments.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void MUSCLE<TV>::
Initialize()
{
    ARRAY<T_MUSCLE_SEGMENT_DATA> segment_data;
    for(int i=1;i<=via_points.m+1;i++) {segment_data.Append(T_MUSCLE_SEGMENT_DATA(MUSCLE_SEGMENT<TV>::LINEAR_SEGMENT,
        ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::CURVE_NONE,ARRAY<T>()));}
    Initialize(segment_data);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void MUSCLE<TV>::
Initialize(ARRAY<T_MUSCLE_SEGMENT_DATA>& segment_data)
{
    assert(segment_data.m==via_points.m+1);
    if(!via_points.m) muscle_segments.Append(Create_Muscle_Segment(attachment_point_1,attachment_point_2,segment_data(1).x,segment_data(1).y,segment_data(1).z));
    else{
        muscle_segments.Append(Create_Muscle_Segment(attachment_point_1,via_points(1),segment_data(1).x,segment_data(1).y,segment_data(1).z));
        for(int i=1;i<via_points.m;i++){muscle_segments.Append(Create_Muscle_Segment(via_points(i),via_points(i+1),
            segment_data(i+1).x,segment_data(i+1).y,segment_data(i+1).z));}
        muscle_segments.Append(Create_Muscle_Segment(via_points(via_points.m),attachment_point_2,segment_data(via_points.m+1).x,segment_data(via_points.m+1).y,segment_data(via_points.m+1).z));}
    for(int i=1;i<=muscle_segments.m;i++){muscle_segments(i)->Initialize();}
}
//#####################################################################
// Function Create_Muscle_Segment
//#####################################################################
template<class TV> MUSCLE_SEGMENT<TV>* MUSCLE<TV>::
Create_Muscle_Segment(ATTACHMENT_POINT<TV>* point_1,ATTACHMENT_POINT<TV>* point_2,const T_MUSCLE_SEGMENT_TYPE& segment_type,
    const T_MUSCLE_SEGMENT_CURVE_TYPE& curve_type,const ARRAY<T>& parameters)
{
    MUSCLE_SEGMENT<TV>* muscle_segment_to_return;
    switch(segment_type){
      case MUSCLE_SEGMENT<TV>::NUM_MUSCLE_SEGMENT_TYPE: 
        muscle_segment_to_return=new ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>(curve_type,point_1,point_2,parameters(1),parameters(2),parameters(3),parameters(4),
            parameters(5));break;
      case MUSCLE_SEGMENT<TV>::LINEAR_SEGMENT: 
        muscle_segment_to_return=new MUSCLE_SEGMENT<TV>(point_1,point_2);break;
      case MUSCLE_SEGMENT<TV>::ANALYTIC_SURFACE_SEGMENT: default: 
        muscle_segment_to_return=new ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>(curve_type,point_1,point_2,parameters(1),parameters(2),parameters(3),parameters(4),
            parameters(5));break;}
    return muscle_segment_to_return;
}
//#####################################################################
// Function Total_Length
//#####################################################################
template<class TV> typename TV::SCALAR MUSCLE<TV>::
Total_Length() const
{
    if(!via_points.m) return (attachment_point_2->Embedded_Position()-attachment_point_1->Embedded_Position()).Magnitude();
    T length=(via_points(1)->Embedded_Position()-attachment_point_1->Embedded_Position()).Magnitude();
    for(int i=2;i<=via_points.m;i++) length+=(via_points(i)->Embedded_Position()-via_points(i-1)->Embedded_Position()).Magnitude();
    length+=(attachment_point_2->Embedded_Position()-via_points(via_points.m)->Embedded_Position()).Magnitude(); 
    return length;
}
//#####################################################################
// Function Total_Velocity
//#####################################################################
template<class TV> typename TV::SCALAR MUSCLE<TV>::
Total_Velocity() const
{
    T total_velocity=0;TV position2,velocity2,path;
    TV position1=attachment_point_1->Embedded_Position(),velocity1=attachment_point_1->Embedded_Velocity();
    for(int i=1;i<=via_points.m;i++){
        position2=via_points(i)->Embedded_Position();velocity2=via_points(i)->Embedded_Velocity();
        total_velocity+=TV::Dot_Product(velocity1-velocity2,(position1-position2).Normalized());
        position1=position2;velocity1=velocity2;}
    return total_velocity+TV::Dot_Product(velocity1-attachment_point_2->Embedded_Velocity(),(position1-attachment_point_2->Embedded_Position()).Normalized());
}
//#####################################################################
// Function Calculate_Activation
//#####################################################################
template<class TV> typename TV::SCALAR MUSCLE<TV>::
Calculate_Activation(const T desired_force) const
{
    T total_length=Total_Length(),total_velocity=Total_Velocity();
    T tendon_length=tendon_slack_length+force_curve.Tendon_Length(desired_force*one_over_peak_force)*tendon_slack_length;
    tendon_length=min(tendon_length,total_length); // clamp tendon length... careful for optimization
    T muscle_length=total_length-tendon_length;
    T desired_active_force=desired_force-Passive_Force(muscle_length);
    if(desired_active_force<=0){{std::stringstream ss;ss<<"\n****------******** ACTIVE_FORCE=0 ***********---------***********\n"<<std::endl;LOG::filecout(ss.str());}return 0;}
    T peak_active_force_quasistatic=Peak_Active_Force_Quasistatic(muscle_length);
    // since velocity scaling could no more than double the quasistatic force, we clamp to 1 if desired active force is too high (this also catches peak_active_force_quasistatic==0 case)
    if(desired_active_force >= 2*peak_active_force_quasistatic){{std::stringstream ss;ss<<"\n****************** CLAMPING ACTIVATION to 1 *********************************\n"<<std::endl;LOG::filecout(ss.str());}return 1;}
    T activation_times_velocity_scale=desired_active_force/peak_active_force_quasistatic;
    // stiffnesses computed using derivatives of force curves w.r.t. length (assume 0 acceleration)
    T tendon_stiffness=Tendon_Force_Slope(tendon_length);
    T muscle_stiffness=Passive_Force_Slope(muscle_length)+Peak_Active_Force_Quasistatic_Slope(muscle_length)*activation_times_velocity_scale;
    T muscle_velocity;
    assert(tendon_stiffness+muscle_stiffness>0); // TODO: hopefully this never even gets very close to 0!
    muscle_velocity=total_velocity*tendon_stiffness/(muscle_stiffness+tendon_stiffness);
    T peak_active_force=peak_active_force_quasistatic*Velocity_Scale(muscle_velocity);
    if(desired_active_force >= peak_active_force){{std::stringstream ss;ss<<"\n****************** CLAMPING ACTIVATION to 1 *********************************\n"<<std::endl;LOG::filecout(ss.str());}return 1;}
    return desired_active_force/peak_active_force;
}
//#####################################################################
// Function Update_Segments
//#####################################################################
template<class TV> void MUSCLE<TV>::
Update_Segments()
{
    for(int i=1;i<=muscle_segments.m;i++){
        muscle_segments(i)->Update_Parameters();
        muscle_segments(i)->Update_Frame();}
}
//#####################################################################
// Function Set_Segment_Activations
//#####################################################################
template<class TV> void MUSCLE<TV>::
Set_Segment_Activations(const T activation) 
{
    // all segments have the same activation currently
    for(int i=1;i<=muscle_segments.m;i++){muscle_segments(i)->Set_Current_Activation(activation);}
}
//#####################################################################
// Function Apply_Fixed_Impulse_At_All_Points
//#####################################################################
template<class TV> void MUSCLE<TV>::
Apply_Fixed_Impulse_At_All_Points(const T impulse)
{
    if(!via_points.m){
        TV F=impulse*(attachment_point_2->Embedded_Position()-attachment_point_1->Embedded_Position()).Normalized();
        attachment_point_1->Apply_Impulse(F);attachment_point_2->Apply_Impulse(-F);}
    else{
        TV F=impulse*(via_points(1)->Embedded_Position()-attachment_point_1->Embedded_Position()).Normalized();
        attachment_point_1->Apply_Impulse(F);via_points(1)->Apply_Impulse(-F);
        for(int i=1;i<=via_points.m-1;i++){
            F=impulse*(via_points(i+1)->Embedded_Position()-via_points(i)->Embedded_Position()).Normalized();
            via_points(i)->Apply_Impulse(F);via_points(i+1)->Apply_Impulse(-F);}
        F=impulse*(attachment_point_2->Embedded_Position()-via_points(via_points.m)->Embedded_Position()).Normalized();
        via_points(via_points.m)->Apply_Impulse(F);attachment_point_2->Apply_Impulse(-F);}
}
//#####################################################################
// Function Read_Constrained_Point
//#####################################################################
template<class TV> ATTACHMENT_POINT<TV>* MUSCLE<TV>::
Read_Constrained_Point(TYPED_ISTREAM& input_stream,RIGID_BODY_COLLECTION<TV>& rigid_body_collection)
{
    TV object_space_position;int id;Read_Binary(input_stream,object_space_position,id);
    // set-up a rigid body binding with a null deformable particle (usable only through the binding callbacks)
    return new ATTACHMENT_POINT<TV>(rigid_body_collection.Rigid_Body(id),object_space_position);
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void MUSCLE<TV>::
Read(TYPED_ISTREAM& input_stream,RIGID_BODY_COLLECTION<TV>& rigid_body_collection)
{
    Read_Binary(input_stream,optimal_length,peak_force,pennation_angle,tendon_slack_length,max_shortening_velocity,name);
    Set_Optimal_Length(optimal_length);Set_Peak_Force(peak_force);Set_Pennation_Angle(pennation_angle);Set_Tendon_Slack_Length(tendon_slack_length);Set_Max_Shortening_Velocity(max_shortening_velocity);
    // Delete_Pointers_And_Clean_Memory(); // TODO: figure out who should own the attachment points (who should delete them)
    attachment_point_1=Read_Constrained_Point(input_stream,rigid_body_collection);attachment_point_2=Read_Constrained_Point(input_stream,rigid_body_collection);
    int num_via_points;Read_Binary(input_stream,num_via_points);via_points.Resize(num_via_points);
    for(int i=1;i<=via_points.m;i++) via_points(i)=Read_Constrained_Point(input_stream,rigid_body_collection);
    int num_muscle_segments;Read_Binary(input_stream,num_muscle_segments);muscle_segments.Resize(num_muscle_segments);
    for(int i=1;i<=muscle_segments.m;i++) muscle_segments(i)=MUSCLE_SEGMENT<TV>::Create_From_Input(input_stream);
    if(!via_points.m){muscle_segments(1)->point_1=attachment_point_1;
    muscle_segments(1)->point_2=attachment_point_2;}
    else{muscle_segments(1)->point_1=attachment_point_1;muscle_segments(1)->point_2=via_points(1);
        for(int i=1;i<via_points.m;i++){muscle_segments(i+1)->point_1=via_points(i);muscle_segments(i+1)->point_2=via_points(i+1);}
        muscle_segments(muscle_segments.m)->point_1=via_points(via_points.m);muscle_segments(muscle_segments.m)->point_2=attachment_point_2;}
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void MUSCLE<TV>::
Write(TYPED_OSTREAM& output_stream) const
{
    Write_Binary(output_stream,optimal_length,peak_force,pennation_angle,tendon_slack_length,max_shortening_velocity,name);
    Write_Constrained_Point(output_stream,attachment_point_1);
    Write_Constrained_Point(output_stream,attachment_point_2);
    Write_Binary(output_stream,via_points.m);for(int i=1;i<=via_points.m;i++) Write_Constrained_Point(output_stream,via_points(i));
    Write_Binary(output_stream,muscle_segments.m);for(int i=1;i<=muscle_segments.m;i++){muscle_segments(i)->Write(output_stream);}
}
//#####################################################################
template class MUSCLE<VECTOR<float,3> >;
template MUSCLE<VECTOR<float,1> >::MUSCLE(MUSCLE_FORCE_CURVE<float> const&);
template void MUSCLE<VECTOR<float,1> >::Set_Segment_Activations(float);
template void MUSCLE<VECTOR<float,1> >::Update_Segments();
template void MUSCLE<VECTOR<float,2> >::Apply_Fixed_Impulse_At_All_Points(float);
template MUSCLE<VECTOR<float,2> >::MUSCLE(MUSCLE_FORCE_CURVE<float> const&);
template void MUSCLE<VECTOR<float,2> >::Set_Segment_Activations(float);
template void MUSCLE<VECTOR<float,2> >::Update_Segments();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MUSCLE<VECTOR<double,3> >;
template MUSCLE<VECTOR<double,1> >::MUSCLE(MUSCLE_FORCE_CURVE<double> const&);
template void MUSCLE<VECTOR<double,1> >::Set_Segment_Activations(double);
template void MUSCLE<VECTOR<double,1> >::Update_Segments();
template void MUSCLE<VECTOR<double,2> >::Apply_Fixed_Impulse_At_All_Points(double);
template MUSCLE<VECTOR<double,2> >::MUSCLE(MUSCLE_FORCE_CURVE<double> const&);
template void MUSCLE<VECTOR<double,2> >::Set_Segment_Activations(double);
template void MUSCLE<VECTOR<double,2> >::Update_Segments();
#endif
