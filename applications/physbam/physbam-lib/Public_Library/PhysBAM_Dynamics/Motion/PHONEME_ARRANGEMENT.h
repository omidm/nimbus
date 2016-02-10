//#####################################################################
// Copyright 2005-2006, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PHONEME_ARRANGEMENT__
#define __PHONEME_ARRANGEMENT__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_STATE.h>
#include <PhysBAM_Dynamics/Motion/FACE_CONTROL_PARAMETERS.h>
#include <PhysBAM_Dynamics/Motion/PHONEME_SEGMENT.h>
namespace PhysBAM{

template<class T>
class PHONEME_ARRANGEMENT
{
    typedef VECTOR<T,3> TV;
public:
    ARRAY<PHONEME_SEGMENT<T> > list;
    const INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,VECTOR_ND<T> >* interpolation;
    typedef T (*BLENDING_CALLBACK_TYPE)(const PHONEME_SEGMENT<T>&,const T);
    BLENDING_CALLBACK_TYPE blending_callback;
    RANGE<VECTOR<T,1> >* valid_segment;
    FACE_CONTROL_PARAMETERS<T>* face_control_parameters;

    PHONEME_ARRANGEMENT()
        :interpolation(0),valid_segment(0),face_control_parameters(0)
    {Set_Blending_Callback();}

    ~PHONEME_ARRANGEMENT()
    {}

    static T Linear_Falloff_Blending(const PHONEME_SEGMENT<T>& phoneme_segment,const T time)
    {if(time<phoneme_segment.time_segment.min_corner.x-phoneme_segment.leading_context_duration) return 0;
    else if(time<phoneme_segment.time_segment.min_corner.x) return (phoneme_segment.leading_context_duration-phoneme_segment.time_segment.min_corner.x+time)/phoneme_segment.leading_context_duration;
    else if(time<phoneme_segment.time_segment.max_corner.x) return (T)1;
    else if(time<phoneme_segment.time_segment.max_corner.x+phoneme_segment.trailing_context_duration) 
        return (phoneme_segment.trailing_context_duration+phoneme_segment.time_segment.max_corner.x-time)/phoneme_segment.trailing_context_duration;
    else return 0;}

    static T Cubic_Falloff_Blending(const PHONEME_SEGMENT<T>& phoneme_segment,const T time)
    {if(time<phoneme_segment.time_segment.min_corner.x-phoneme_segment.leading_context_duration) return 0;
    else if(time<phoneme_segment.time_segment.min_corner.x) return cube((phoneme_segment.leading_context_duration-phoneme_segment.time_segment.min_corner.x+time)/phoneme_segment.leading_context_duration);
    else if(time<phoneme_segment.time_segment.max_corner.x) return (T)1;
    else if(time<phoneme_segment.time_segment.max_corner.x+phoneme_segment.trailing_context_duration) 
        return cube((phoneme_segment.trailing_context_duration+phoneme_segment.time_segment.max_corner.x-time)/phoneme_segment.trailing_context_duration);
    else return 0;}

    void Set_Blending_Callback(BLENDING_CALLBACK_TYPE blending_callback_input=PHONEME_ARRANGEMENT<T>::Linear_Falloff_Blending)
    {blending_callback=blending_callback_input;}

    void Set_Custom_Interpolation(const INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,VECTOR_ND<T> >* interpolation_input)
    {interpolation=interpolation_input;}

    void Set_Face_Control_Parameters(FACE_CONTROL_PARAMETERS<T>* face_control_parameters_input)
    {face_control_parameters=face_control_parameters_input;}

    template<class RW>
    void Read_Samples(const std::string& base_directory)
    {for(int i=1;i<=list.m;i++) list(i).template Read_Sample<RW>(base_directory);}

    void Update_Valid_Segment()
    {if(!valid_segment)valid_segment=new RANGE<VECTOR<T,1> >;(*valid_segment)=RANGE<VECTOR<T,1> >(FLT_MAX,-FLT_MAX);for(int i=1;i<=list.m;i++)valid_segment->Enlarge_To_Include_Box(list(i).Context_Segment());}

    void Segments_Intersecting_Time_Instance(const T time,ARRAY<int>& segment_list) const
    {segment_list.Remove_All();for(int i=1;i<=list.m;i++)if(list(i).Context_Segment().Inside(VECTOR<T,1>(time),0))segment_list.Append(i);}

    template <class RW>
    void Interpolate_Rigid_Body_State(RIGID_BODY_STATE<VECTOR<T,3> >& rigid_body_state,const std::string phoneme_data_directory,const std::string rigid_transform_filename_prefix,const T time)
    {ARRAY<int> segment_list;Segments_Intersecting_Time_Instance(time,segment_list);
    switch(segment_list.m){
      case 0: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("No value to interpolate from at specified time instance %d",time));
      case 1: list(segment_list(1)).template Interpolate_Rigid_Body_State<RW>(rigid_body_state,phoneme_data_directory,rigid_transform_filename_prefix,time);return;
      case 2:
          {ARRAY<RIGID_BODY_STATE<VECTOR<T,3> > > rigid_body_states(2);for(int i=1;i<=2;i++)list(segment_list(i)).template Interpolate_Rigid_Body_State<RW>(rigid_body_states(i),phoneme_data_directory,rigid_transform_filename_prefix,time);
          ARRAY<T> blending_weights(2);for(int i=1;i<=2;i++)blending_weights(i)=(*blending_callback)(list(segment_list(i)),time);
          T blending_fraction=blending_weights(2)/(blending_weights(1)+blending_weights(2));
          rigid_body_state.time=time;
          rigid_body_state.frame=FRAME<TV>::Interpolation(rigid_body_states(1).frame,rigid_body_states(2).frame,blending_fraction);
          rigid_body_state.velocity=((T)1-blending_fraction)*rigid_body_states(1).velocity+blending_fraction*rigid_body_states(2).velocity;
          rigid_body_state.angular_velocity=((T)1-blending_fraction)*rigid_body_states(1).angular_velocity+blending_fraction*rigid_body_states(2).angular_velocity;return;}
      default: PHYSBAM_FATAL_ERROR("Interpolation between 3 or more phoneme segments is not supported");}}

    template <class RW>
    void Interpolate_Positions(ARRAY<VECTOR<T,3> >& position,const std::string phoneme_data_directory,const std::string position_filename_prefix,const T time)
    {ARRAY<int> segment_list;Segments_Intersecting_Time_Instance(time,segment_list);
    switch(segment_list.m){
      case 0: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("No value to interpolate from at specified time instance %d",time));
      case 1: list(segment_list(1)).template Interpolate_Positions<RW>(position,phoneme_data_directory,position_filename_prefix,time);return;
      case 2:
          {ARRAY<ARRAY<VECTOR<T,3> > > positions(2);for(int i=1;i<=2;i++)list(segment_list(i)).template Interpolate_Positions<RW>(positions(i),phoneme_data_directory,position_filename_prefix,time);
          ARRAY<T> blending_weights(2);for(int i=1;i<=2;i++)blending_weights(i)=(*blending_callback)(list(segment_list(i)),time);
          T blending_fraction=blending_weights(2)/(blending_weights(1)+blending_weights(2));
          if(!position.m) position.Resize(positions(1).m);
          ARRAY<VECTOR<T,3> >::copy((T)1-blending_fraction,positions(1),blending_fraction,positions(2),position);return;}
      default: PHYSBAM_FATAL_ERROR("Interpolation between 3 or more phoneme segments is not supported");}}

    VECTOR_ND<T> Controls(const T time) const
    {assert(interpolation);ARRAY<int> segment_list;Segments_Intersecting_Time_Instance(time,segment_list);
    switch(segment_list.m){
      case 0: 
          {PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("No value to interpolate from at specified time instance %d",time));}
      case 1: 
          {face_control_parameters->Set(list(segment_list(1)).Controls(time));face_control_parameters->Scale(list(segment_list(1)).scaling);
          VECTOR_ND<T> controls;face_control_parameters->Get(controls);return controls;}
      case 2:
          {ARRAY<VECTOR_ND<T> ,VECTOR<int,1> > blending_controls(1,2);
          for(int i=1;i<=2;i++){
              blending_controls(i)=list(segment_list(i)).Controls(time);
              face_control_parameters->Set(blending_controls(i));
              face_control_parameters->Scale(list(segment_list(i)).scaling);
              face_control_parameters->Get(blending_controls(i));}
          ARRAY<T> blending_weights(2);for(int i=1;i<=2;i++)blending_weights(i)=(*blending_callback)(list(segment_list(i)),time);
          T blending_fraction=blending_weights(2)/(blending_weights(1)+blending_weights(2));
          GRID<VECTOR<T,1> > blending_grid(2,0,(T)1);
          return interpolation->Clamped_To_Array(blending_grid,blending_controls,blending_fraction);}
      default: PHYSBAM_FATAL_ERROR("Interpolation between 3 or more phoneme segments is not supported");}}

    T Discrepancy_Between_Inputs(const T time,const ARRAY<T>& weights) const
    {assert(interpolation);assert(face_control_parameters);ARRAY<int> segment_list;Segments_Intersecting_Time_Instance(time,segment_list);
    switch(segment_list.m){
      case 0: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("No value to interpolate from at specified time instance %d",time));
      case 1: return 0;
      case 2:
        // Fill the face_control_parameters twice, with appropriate scaling, and get good distance
          {const PHONEME_SEGMENT<T>& seg1=list(segment_list(1)),&seg2=list(segment_list(2));
          VECTOR_ND<T> controls1=seg1.Controls(time),controls2=seg2.Controls(time);
          face_control_parameters->Set(controls1);
          if(seg1.scaling!=1)face_control_parameters->Scale(seg1.scaling);
          face_control_parameters->Save_Controls();face_control_parameters->Set(controls2);
          if(seg2.scaling!=1)face_control_parameters->Scale(seg2.scaling);
          return (seg2.scaling+seg1.scaling+1/seg1.scaling+1/seg2.scaling)*face_control_parameters->Distance(weights);}
      default: PHYSBAM_FATAL_ERROR("Interpolation between 3 or more phoneme segments is not supported");}}
    
    VECTOR_ND<T> Optimal_Scaling(const GRID<VECTOR<T,1> > sampling_grid,const T stiffness,int* samples_used=0)
    {if(samples_used)*samples_used=0;
    MATRIX_MXN<T> normal_equations_matrix(list.m,list.m);VECTOR_ND<T> normal_equations_rhs(list.m);
    for(int sample=1;sample<=sampling_grid.counts.x;sample++){
        T time=sampling_grid.Axis_X(sample,1);ARRAY<int> segment_list;Segments_Intersecting_Time_Instance(time,segment_list);
        if(segment_list.m==2){
            if(samples_used)(*samples_used)++;
            int i=segment_list(1),j=segment_list(2);
            VECTOR_ND<T> ci=list(i).Controls(time),cj=list(j).Controls(time);
            normal_equations_matrix(i,i)+=VECTOR_ND<T>::Dot_Product(ci,ci);
            normal_equations_matrix(i,j)-=VECTOR_ND<T>::Dot_Product(ci,cj);
            normal_equations_matrix(j,i)-=VECTOR_ND<T>::Dot_Product(cj,ci);
            normal_equations_matrix(j,j)+=VECTOR_ND<T>::Dot_Product(cj,cj);}
        for(int segment=1;segment<=segment_list.m;segment++){
            int i=segment_list(segment);VECTOR_ND<T> ci=list(i).Controls(time);
            normal_equations_matrix(i,i)+=sqr(stiffness)*VECTOR_ND<T>::Dot_Product(ci,ci);
            normal_equations_rhs(i)+=sqr(stiffness)*VECTOR_ND<T>::Dot_Product(ci,ci);}}
    return normal_equations_matrix.PLU_Solve(normal_equations_rhs);}

    template<class RW>
    void Read(std::istream& input_stream)
    {list.template Read<RW>(input_stream);}

    template<class RW> 
    void Write(std::ostream& output_stream) const
    {list.template Write<RW>(output_stream);}

//#####################################################################
};
}
#endif
