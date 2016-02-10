//#####################################################################
// Copyright 2005-2006, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PHONEME_CORPUS__
#define __PHONEME_CORPUS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Dynamics/Motion/FACE_CONTROL_PARAMETERS.h>
#include <PhysBAM_Dynamics/Motion/PHONEME.h>
#include <PhysBAM_Dynamics/Motion/PHONEME_ARRANGEMENT.h>

// Store phonemes; given list of phonemes, produce PHONEME_ARRANGEMENT
namespace PhysBAM{

template<class T>
class PHONEME_CORPUS
{
public:
    HASHTABLE<std::string,ARRAY<PHONEME_SEGMENT<T>*> *> string_to_phoneme_index;
    ARRAY<std::string> phoneme_labels;
    HASHTABLE<std::string,VECTOR_ND<T>*> phoneme_stiffness;
    const INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,VECTOR_ND<T> >* interpolation;
    FACE_CONTROL_PARAMETERS<T>* face_control_parameters;
    ARRAY<T> weights;

    PHONEME_CORPUS()
        :interpolation(0),face_control_parameters(0)
    {}
    
    void Set_Custom_Interpolation(const INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,VECTOR_ND<T> >* interpolation_input)
    {interpolation=interpolation_input;}

    void Set_Face_Control_Parameters(FACE_CONTROL_PARAMETERS<T>* face_control_parameters_input)
    {face_control_parameters=face_control_parameters_input;}

    ARRAY<PHONEME_SEGMENT<T>*>* Get_Phonemes(const std::string& phoneme)
    {ARRAY<PHONEME_SEGMENT<T>*>* phonemes=0;
    if(!string_to_phoneme_index.Get(phoneme,phonemes)) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Phoneme %s not found",phoneme.c_str()));
    return phonemes;}

    void Add_Phoneme_Segment(const std::string& phoneme_file)
    {PHONEME_SEGMENT<T>* phoneme_segment=new PHONEME_SEGMENT<T>();phoneme_segment->phoneme_sample_name=FILE_UTILITIES::Get_Short_Name(phoneme_file);
    phoneme_segment->template Read_Sample<T>(FILE_UTILITIES::Get_Base_Directory_Name(phoneme_file));Add_Phoneme_Segment(phoneme_segment);}

    void Add_Phoneme_Segment(PHONEME_SEGMENT<T>* phoneme_segment)
    {ARRAY<PHONEME_SEGMENT<T>*>* phonemes=0;const std::string& phoneme_name=phoneme_segment->phoneme_sample->name;
    if(!string_to_phoneme_index.Get(phoneme_name,phonemes)){phonemes=new ARRAY<PHONEME_SEGMENT<T>*>();string_to_phoneme_index.Insert(phoneme_name,phonemes);phoneme_labels.Append(phoneme_name);}
    phonemes->Append_Unique(phoneme_segment);phoneme_segment->Set_Custom_Interpolation(interpolation);}

//#####################################################################
    void Get_Next_Kind_Of_Phoneme(PHONEME_SEGMENT<T>& phoneme_segment);
    void Get_Previous_Kind_Of_Phoneme(PHONEME_SEGMENT<T>& phoneme_segment);
    void Get_Next_Phoneme(PHONEME_SEGMENT<T>& phoneme_segment);
    void Get_Previous_Phoneme(PHONEME_SEGMENT<T>& phoneme_segment);
    void Gather_Phoneme_Statistics();
    void Build_Corpus(const std::string& corpus_directory);
    PHONEME_ARRANGEMENT<T>* Produce_Phoneme_Arrangement(const ARRAY<std::string>& phonemes,const ARRAY<T>& starting_times,const ARRAY<T>& lengths,const ARRAY<T>& leading_times,const ARRAY<T>& trailing_times,const T decay=.9995,const int step=20000,const T Tzero=200);
    PHONEME_ARRANGEMENT<T>* Brute_Force_Phoneme_Arrangement(const ARRAY<std::string>& phonemes,const ARRAY<T>& starting_times,const ARRAY<T>& lengths,const ARRAY<T>& leading_times,const ARRAY<T>& trailing_times);
private:
    static bool Brute_Force_Phoneme_Arrangement_Helper(const ARRAY<ARRAY<PHONEME_SEGMENT<T>*> >& phoneme_candidates,const int position,PHONEME_ARRANGEMENT<T>& phoneme_arrangement,T &best_result,ARRAY<int>& indices,const ARRAY<T>& weights);
    static T Phoneme_Metric_Helper_Scaled(PHONEME_ARRANGEMENT<T>& phoneme_arrangement,const ARRAY<T>& scalings,const ARRAY<T>& weights);
    static T Get_Optimal_Scalings(PHONEME_ARRANGEMENT<T>& phoneme_arrangement,const ARRAY<T>& weights);
    static void Gather_Phoneme_Statistics_Helper(PHONEME_CORPUS<T>* phoneme_corpus,ARRAY<PHONEME_SEGMENT<T>*>* node);
    T Phoneme_Metric_Helper(PHONEME_ARRANGEMENT<T>& phoneme_arrangement,const ARRAY<T>& weights);
//#####################################################################
};
}
#endif
