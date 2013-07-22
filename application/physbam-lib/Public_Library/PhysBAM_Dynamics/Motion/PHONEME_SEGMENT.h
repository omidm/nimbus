//#####################################################################
// Copyright 2005-2006, Avi Robinson-Mosher, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PHONEME_SEGMENT__
#define __PHONEME_SEGMENT__

#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_STATE.h>
#include <PhysBAM_Dynamics/Motion/PHONEME.h>
#include <cstdio>
#include <cstring>
namespace PhysBAM{

#define ALTERNATE_PHONEMES 1

template<class T>
class PHONEME_SEGMENT
{
public:
    PHONEME<T>* phoneme_sample;
    std::string phoneme_sample_name;
    RANGE<VECTOR<T,1> > time_segment;
    T leading_context_duration,trailing_context_duration;
    T scaling;
#if ALTERNATE_PHONEMES
    T peak; // middle of PHONEME
    T time_scaling;
#endif
private:
    bool owns_data;
public:
    PHONEME_SEGMENT()
        :phoneme_sample(0),scaling(1),owns_data(false)
    {}

    PHONEME_SEGMENT(const PHONEME_SEGMENT<T>& phoneme_segment_input)
    {*this=phoneme_segment_input;scaling=phoneme_segment_input.scaling;owns_data=false;}

    ~PHONEME_SEGMENT()
    {if(owns_data)delete phoneme_sample;}

    void Set_Sample(PHONEME<T>& phoneme_sample_input,const bool owns_data_input=false)
    {phoneme_sample=&phoneme_sample_input;owns_data=owns_data_input;time_scaling=time_segment.Size()/phoneme_sample->time_length;if(time_scaling<0){*((int*)0)=1; }}

    T World_Time_To_Sample_Time(T time)
    {T offset=time-peak;return offset*time_scaling+phoneme_sample->time_length/2;}

    void Set_Custom_Interpolation(const INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,VECTOR_ND<T> >* interpolation_input)
    {assert(phoneme_sample);phoneme_sample->Set_Custom_Interpolation(interpolation_input);}

    template <class RW>
    void Interpolate_Rigid_Body_State(RIGID_BODY_STATE<VECTOR<T,3> >& rigid_body_state,const std::string phoneme_data_directory,const std::string rigid_transform_filename_prefix,const T time)
    {assert(phoneme_sample);char parse_buffer[1000];char name[100],previous_name[100],next_name[100],sequence_name[100];
    int first_frame_sequence,last_frame_sequence,first_frame_context,last_frame_context,first_frame_phoneme,last_frame_phoneme;
    strcpy(parse_buffer,phoneme_sample_name.c_str());for(int i=0;parse_buffer[i]!='\0';i++)if(parse_buffer[i]=='.'||parse_buffer[i]=='-')parse_buffer[i]=' ';
    sscanf(parse_buffer,"%s %s %s %s %d %d %d %d %d %d",
           name,previous_name,next_name,sequence_name,&first_frame_sequence,&last_frame_sequence,&first_frame_context,&last_frame_context,&first_frame_phoneme,&last_frame_phoneme);
    phoneme_sample->template Interpolate_Rigid_Body_State<RW>(rigid_body_state,STRING_UTILITIES::string_sprintf("%s/%s/%d-%d",phoneme_data_directory.c_str(),sequence_name,first_frame_sequence,last_frame_sequence),
                                                              rigid_transform_filename_prefix,first_frame_phoneme,World_Time_To_Sample_Time(time));
    rigid_body_state.time=time;rigid_body_state.velocity*=1/time_scaling;rigid_body_state.angular_velocity*=1/time_scaling;}

    template <class RW>
    void Interpolate_Positions(ARRAY<VECTOR<T,3> >& positions,const std::string phoneme_data_directory,const std::string position_filename_prefix,const T time)
    {assert(phoneme_sample);char parse_buffer[1000];char name[100],previous_name[100],next_name[100],sequence_name[100];
    int first_frame_sequence,last_frame_sequence,first_frame_context,last_frame_context,first_frame_phoneme,last_frame_phoneme;
    strcpy(parse_buffer,phoneme_sample_name.c_str());for(int i=0;parse_buffer[i]!='\0';i++)if(parse_buffer[i]=='.'||parse_buffer[i]=='-')parse_buffer[i]=' ';
    sscanf(parse_buffer,"%s %s %s %s %d %d %d %d %d %d",
           name,previous_name,next_name,sequence_name,&first_frame_sequence,&last_frame_sequence,&first_frame_context,&last_frame_context,&first_frame_phoneme,&last_frame_phoneme);
    phoneme_sample->template Interpolate_Positions<RW>(positions,STRING_UTILITIES::string_sprintf("%s/%s/%d-%d",phoneme_data_directory.c_str(),sequence_name,first_frame_sequence,last_frame_sequence),
                                                       position_filename_prefix,first_frame_phoneme,World_Time_To_Sample_Time(time));}

    VECTOR_ND<T> Controls(const T time) const
    {assert(phoneme_sample);GRID<VECTOR<T,1> > time_grid(phoneme_sample->frame_length,peak-(phoneme_sample->frame_length/(T)2)*phoneme_sample->time_length/phoneme_sample->frame_length*time_scaling,
                                                 peak+(phoneme_sample->frame_length/(T)2)*phoneme_sample->time_length/phoneme_sample->frame_length*time_scaling);
    return phoneme_sample->Controls(clamp(time,time_segment.min_corner.x-leading_context_duration,time_segment.max_corner.x+trailing_context_duration),time_grid);}

    RANGE<VECTOR<T,1> > Context_Segment() const 
    {return RANGE<VECTOR<T,1> >(time_segment.min_corner.x-leading_context_duration,time_segment.max_corner.x+trailing_context_duration);}

    template<class RW>
    void Read_Sample(const std::string& base_directory)
    {if(phoneme_sample)delete phoneme_sample;phoneme_sample=new PHONEME<T>();FILE_UTILITIES::Read_From_File<RW>(base_directory+"/"+phoneme_sample_name,*phoneme_sample);owns_data=true;}

    template<class RW>
    void Read(std::istream& input_stream)
    {Read_Binary<RW>(input_stream,phoneme_sample_name,time_segment,leading_context_duration,trailing_context_duration,scaling,peak,time_scaling);}

    template<class RW>
    void Read_Legacy(std::istream& input_stream)
    {Read_Binary<RW>(input_stream,phoneme_sample_name,time_segment,leading_context_duration,trailing_context_duration,scaling);}

    template<class RW> 
    void Write(std::ostream& output_stream) const
    {Write_Binary<RW>(output_stream,phoneme_sample_name,time_segment,leading_context_duration,trailing_context_duration,scaling,peak,time_scaling);}
//#####################################################################
};
}
#endif
