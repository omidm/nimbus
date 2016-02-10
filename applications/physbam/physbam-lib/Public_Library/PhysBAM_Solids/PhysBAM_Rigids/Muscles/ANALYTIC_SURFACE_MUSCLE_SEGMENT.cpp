//#####################################################################
// Copyright 2006-2007, Kevin Der, Ranjitha Kumar, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/QUEUE.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/ANALYTIC_SURFACE_MUSCLE_SEGMENT.h>
#include <cassert>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
ANALYTIC_SURFACE_MUSCLE_SEGMENT()
{}
//#####################################################################
// Constructor
//#####################################################################
template<class T> ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
ANALYTIC_SURFACE_MUSCLE_SEGMENT(CURVE_TYPE curve_type_input,ATTACHMENT_POINT<TV>* point_1_input,ATTACHMENT_POINT<TV>* point_2_input,
    const T curve_thickness_input,const T curve_offset_thickness_input,const T tendon_slack_length_input,const T tendon_fraction_1_input,const T tendon_fraction_2_input)
    :MUSCLE_SEGMENT<TV>(point_1_input,point_2_input),curve_type(curve_type_input),initial_curve_thickness(curve_thickness_input),
    initial_curve_offset_thickness(curve_offset_thickness_input),curve_thickness(curve_thickness_input),curve_offset_thickness(curve_offset_thickness_input),
    tendon_fraction_1(tendon_fraction_1_input),tendon_fraction_2(tendon_fraction_2_input),tendon_slack_length(tendon_slack_length_input),activation_factor(0)
{
    segment_type=MUSCLE_SEGMENT<TV>::ANALYTIC_SURFACE_SEGMENT;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
~ANALYTIC_SURFACE_MUSCLE_SEGMENT()
{}
//#####################################################################
// Function Curve_Length
//#####################################################################
template<class T> T ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
Curve_Length() const
{
    T curve_length=Length()*(1-tendon_fraction_1-tendon_fraction_2)*activation_factor;
    return curve_length;
}
//#####################################################################
// Function Tendon_Length
//#####################################################################
template<class T> T ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
Tendon_Length() const
{
    return Length()-Curve_Length();
}
//#####################################################################
// Function Set_Current_Activation
//#####################################################################
template<class T> void ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
Set_Current_Activation(const T activation)
{
    MUSCLE_SEGMENT<TV>::Set_Current_Activation(activation);
    activation_factor=0;
    for(int i=1;i<=min(activation_memory,activations->Size());i++){T activation=min(max((T)0,(*activations)(i)),(T)1);activation_factor+=log(1+activation);}
    activation_factor/=min(activation_memory,activations->Size());activation_factor=1-activation_factor;
    assert(activation_factor<=1);assert(activation_factor>=0);
}
//#####################################################################
// Function Compute_Volume -- only volume of the bulge
//#####################################################################
template<class T> T ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
Compute_Volume() const
{
    if(curve_type==CURVE_ELLIPSE)
        return (T)pi*Tendon_Length()*sqr(curve_offset_thickness)+Curve_Length()*(T)pi/6*(4*sqr(curve_thickness)+3*curve_thickness*curve_offset_thickness*(T)pi+6*sqr(curve_offset_thickness));
    else{assert(curve_type==CURVE_COSINE);
        return (T)pi*Curve_Length()/2*(3*sqr(curve_thickness)+4*curve_thickness*curve_offset_thickness+2*sqr(curve_offset_thickness));}
}
//#####################################################################
// Function Scaled_Volume
//#####################################################################
template<class T> T ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
Scaled_Volume(const T segment_endpoint) const
{
    assert(curve_type==CURVE_COSINE);
    T subst=curve_thickness/curve_offset_thickness+1; 
    return segment_endpoint*(1+2*sqr(subst))+sin(segment_endpoint)*(cos(segment_endpoint)-4*subst);
}
//#####################################################################
// Function Scaled_Volume_Derivative
//#####################################################################
template<class T> T ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
Scaled_Volume_Derivative(const T segment_endpoint) const
{
    assert(curve_type==CURVE_COSINE);
    T subst=curve_thickness/curve_offset_thickness+1; 
    return 2*sqr(subst-cos(segment_endpoint));
}
//#####################################################################
// Function Update_Parameters
//#####################################################################
template<class T> void ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
Update_Parameters()
{   
    T curve_length=Curve_Length(),tendon_length=Tendon_Length(),length_ratio=tendon_length/curve_length;
    
    if(curve_type==CURVE_ELLIPSE){
        T discriminant=9*sqr(initial_curve_offset_thickness)*(T)pi_squared-16*(3*sqr(initial_curve_offset_thickness)*(2+2*length_ratio)-6*initial_volume/(curve_length*(T)pi));
        if(discriminant<0){curve_thickness=0;curve_offset_thickness=sqrt(2*initial_volume/((T)pi*curve_length)/(2+2*length_ratio));}
        else{T new_curve_thickness=2*(3*sqr(initial_curve_offset_thickness)*(2+2*length_ratio)-6*initial_volume/(curve_length*(T)pi))/(-3*initial_curve_offset_thickness*(T)pi-sqrt(discriminant));
            if(new_curve_thickness>=0) {curve_thickness=new_curve_thickness;curve_offset_thickness=initial_curve_offset_thickness;}
            else{curve_thickness=0;curve_offset_thickness=sqrt(2*initial_volume/((T)pi*curve_length)/(2+2*length_ratio));}}}
    else{assert(curve_type==CURVE_COSINE);
        T quad_term=2*(initial_volume/((T)pi*sqr(initial_curve_offset_thickness)*curve_length)-1),discriminant=4+3*quad_term;
        if(discriminant<0){curve_thickness=0;curve_offset_thickness=sqrt(initial_volume/((T)pi*curve_length));}
        else{T new_curve_thickness=(quad_term*curve_offset_thickness)/(2+sqrt(discriminant));
            if(new_curve_thickness>=0){curve_thickness=new_curve_thickness;curve_offset_thickness=initial_curve_offset_thickness;}
            else{curve_thickness=0;curve_offset_thickness=sqrt(initial_volume/((T)pi*curve_length));}}
        for(int i=0;i<num_segments_over_2;i++){
            T prev_endpoint,curr_endpoint=fractional_segment_endpoints(i);
            do{prev_endpoint=curr_endpoint;curr_endpoint=curr_endpoint-(Scaled_Volume(curr_endpoint)-initial_segment_volumes(i))/Scaled_Volume_Derivative(curr_endpoint);}
            while(abs(curr_endpoint-prev_endpoint)>=.01);
            fractional_segment_endpoints(i)=curr_endpoint;fractional_segment_endpoints(2*num_segments_over_2-i)=2*(T)pi-curr_endpoint;}}
}
//#####################################################################
// Function Get_Elongation_At_Local_Position
//#####################################################################
template<class T> T ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
Get_Elongation_At_Local_Position(const TV& initial_local_position,const int segment_number) const
{
    if(segment_number<1||segment_number>(2*num_segments_over_2)) return 0;
    T initial_segment_width=initial_curve_length/(2*num_segments_over_2),segment_width;
    if(segment_number==1) segment_width=fractional_segment_endpoints(segment_number)*Curve_Length()/(2*(T)pi);
    else segment_width=(fractional_segment_endpoints(segment_number)-fractional_segment_endpoints(segment_number-1))*Curve_Length()/(2*(T)pi);
    return segment_width/initial_segment_width;
}
//#####################################################################
// Function Get_Fractional_Curve_Value
//#####################################################################
template<class T> T ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
Get_Fractional_Curve_Value(const T fraction,const bool initial) const
{
    assert(fraction>=(T)0&&fraction<=(T)1);
    T curve_thickness_var=initial?initial_curve_thickness:curve_thickness,curve_offset_thickness_var=initial?initial_curve_offset_thickness:curve_offset_thickness;
    T curve_length_fraction=activation_factor*(1-tendon_fraction_1-tendon_fraction_2);
    T tendon_fraction_1_scaled=(tendon_fraction_1+tendon_fraction_2<small_number?(T)0.5:tendon_fraction_1/(tendon_fraction_1+tendon_fraction_2))*(1-curve_length_fraction);
    T tendon_fraction_2_scaled=(tendon_fraction_1+tendon_fraction_2<small_number?(T)0.5:tendon_fraction_2/(tendon_fraction_1+tendon_fraction_2))*(1-curve_length_fraction);
    if(fraction<=tendon_fraction_1_scaled||fraction>=(1-tendon_fraction_2_scaled)) return curve_offset_thickness_var;
    T curve_fraction=(fraction-tendon_fraction_1_scaled)/curve_length_fraction;
    if(curve_type==CURVE_ELLIPSE) return curve_offset_thickness_var+sqrt(sqr(curve_thickness_var)*(1-sqr(2*curve_fraction-1)));
    else{assert(curve_type==CURVE_COSINE);return -curve_thickness_var*cos((T)two_pi*curve_fraction)+curve_thickness_var+curve_offset_thickness_var;}
}
//#####################################################################
// Function Analytic_Inside_Test
//#####################################################################
template<class T> bool ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
Analytic_Inside_Test(const TV& local_position) const
{
    if(local_position(1)<=0 || local_position(1)>=Length()) return false;
    T distance_from_axis=sqrt(sqr(local_position(2))+sqr(local_position(3)));
    T curve_height=Get_Fractional_Curve_Value(local_position(1)/Length(),false);
    return distance_from_axis<=curve_height;
}
//#####################################################################
// Function Initialize_Inside_Particles
//#####################################################################
template<class T> void ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
Initialize_Inside_Particles(const TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{
    GEOMETRY_PARTICLES<TV>& particles=tetrahedralized_volume.particles;
    inside_particle_rest_positions.Resize(particles.array_collection->Size());inside_particle_segments.Resize(particles.array_collection->Size());
    for(int t=1;t<=tetrahedralized_volume.mesh.elements.m;t++){
        for(int v=1;v<=4;v++){int node=tetrahedralized_volume.mesh.elements(t)(v);
            inside_particle_rest_positions(node)=frame.Inverse()*particles.X(node);
            inside_particle_segments(node)=(int)((2*(inside_particle_rest_positions(node)(1)-tendon_fraction_1*Length())*num_segments_over_2)/Curve_Length()+1);
            if(!Analytic_Inside_Test(inside_particle_rest_positions(node))) break;}}
}
//#####################################################################
// Function Get_Local_Positions_For_Particles
//#####################################################################
template<class T> void ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
Get_Local_Positions_For_Particles(const int m,const int n,GEOMETRY_PARTICLES<TV>& particles)
{
    assert(particles.array_collection->Size()==m*n+2);
    T length=Length(),dtheta=(T)two_pi/n;
    for(int i=1;i<=m;i++){
        T x_fraction=(i-1)/(T)(m-1);T distance_from_axis=Get_Fractional_Curve_Value(x_fraction,false);
        for(int j=1;j<=n;j++) particles.X(j+(i-1)*n)=TV(x_fraction*length,distance_from_axis*sin((j-1)*dtheta),distance_from_axis*cos((j-1)*dtheta));}
    particles.X(m*n+1)=TV(0,0,0);particles.X(m*n+2)=TV(length,0,0);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
Initialize()
{
    BASE::Initialize();for(int i=1;i<=activation_memory;i++) Set_Current_Activation(0);
    initial_volume=Compute_Volume();initial_length=Length();initial_curve_length=Curve_Length();
    T subst=0;num_segments_over_2=100;
    initial_segment_volumes.Resize(2*num_segments_over_2);fractional_segment_endpoints.Resize(2*num_segments_over_2);
    for(int i=1;i<=2*num_segments_over_2;i++){subst+=(T)pi/num_segments_over_2;fractional_segment_endpoints(i)=subst;initial_segment_volumes(i)=Scaled_Volume(subst);}
}
//#####################################################################
// Function Read_And_Set_Parameters
//#####################################################################
template<class T> void ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
Read_And_Set_Parameters(TYPED_ISTREAM& input)
{
    BASE::Read_And_Set_Parameters(input);Read_Binary(input,curve_type);
    Read_Binary(input,initial_curve_thickness,initial_curve_offset_thickness,curve_thickness,curve_offset_thickness,initial_volume,initial_length);
    Read_Binary(input,tendon_fraction_1,tendon_fraction_2,tendon_slack_length,activation_factor);
}
//#####################################################################
// Function Write_Parameters
//#####################################################################
template<class T> void ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
Write_Parameters(TYPED_OSTREAM& output) const
{
    BASE::Write_Parameters(output);Write_Binary(output,curve_type);
    Write_Binary(output,initial_curve_thickness,initial_curve_offset_thickness,curve_thickness,curve_offset_thickness,initial_volume,initial_length);
    Write_Binary(output,tendon_fraction_1,tendon_fraction_2,tendon_slack_length,activation_factor);
}
//#####################################################################
// Function Maximum_Radius
//#####################################################################
template<class T> T ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
Maximum_Radius() const
{
    switch(curve_type){
        case CURVE_ELLIPSE:
            return curve_offset_thickness+curve_thickness;
            break;
        default:
            assert(curve_type==CURVE_COSINE);
            return curve_offset_thickness+2*curve_thickness;}
    return 0;
}
//#####################################################################
// Function Create
//#####################################################################
template<class T> ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>* ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
Create()
{
    return new ANALYTIC_SURFACE_MUSCLE_SEGMENT();
}
//#####################################################################
// Function Name
//#####################################################################
template<class T> std::string ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::
Static_Name()
{
    return "analytic_surface_muscle_segment";
}
//#####################################################################
template class ANALYTIC_SURFACE_MUSCLE_SEGMENT<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ANALYTIC_SURFACE_MUSCLE_SEGMENT<double>;
#endif
