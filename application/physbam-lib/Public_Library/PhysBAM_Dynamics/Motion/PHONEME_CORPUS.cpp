//#####################################################################
// Copyright 2005-2006, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Dynamics/Motion/PHONEME_CORPUS.h>
#if defined(__linux__) || defined(__CYGWIN__) || defined(__APPLE__)
#include <dirent.h>
#include <sys/types.h>
#endif
using namespace PhysBAM;

//#####################################################################
// Function Gather_Phoneme_Statistics_Helper
//#####################################################################
template<class T> void PHONEME_CORPUS<T>::
Gather_Phoneme_Statistics_Helper(PHONEME_CORPUS<T>* phoneme_corpus,ARRAY<PHONEME_SEGMENT<T>*>* node)
{
    ARRAY<PHONEME_SEGMENT<T>*>& phonemes=*node;
    VECTOR_ND<T> *errors=new VECTOR_ND<T>;
    T time_length=0;
    for(int i=1;i<=phonemes.m;i++)time_length+=phonemes(i)->phoneme_sample->time_length;
    time_length/=phonemes.m;
    int frame_length=(int)(time_length*120);
    ARRAY<VECTOR_ND<T> > means(frame_length);
    for(int i=1;i<=phonemes.m;i++){
        phonemes(i)->phoneme_sample->time_length=time_length;
        for(int t=0;t<frame_length;t++){
            VECTOR_ND<T> controls=phonemes(i)->Controls(t/(T)120);
            if(!means(t+1).n)means(t+1).Resize(controls.n);
            means(t+1)+=controls;}}
    VECTOR_ND<T> per_attribute_mean(means(1).n);
    for(int t=1;t<=frame_length;t++){
        means(t)/=(T)phonemes.m;
        for(int i=1;i<=per_attribute_mean.n;i++) per_attribute_mean(i)+=abs(means(t)(i));}
    per_attribute_mean/=(T)frame_length;
    errors->Resize(means(1).n);
    // get deviations
    for(int i=1;i<=phonemes.m;i++)for(int t=0;t<frame_length;t++){VECTOR_ND<T> controls=phonemes(i)->Controls(t/(T)120);for(int j=1;j<=(*errors).n;j++)(*errors)(j)+=pow(controls(j)-means(t+1)(j),2);}
    std::stringstream ss;
    ss<<"Stiffnesses for phoneme "<<phonemes(1)->phoneme_sample->name<<std::endl;
    for(int i=1;i<=(*errors).n;i++){(*errors)(i)=sqrt((*errors)(i)/phonemes.m)/abs(per_attribute_mean(i));ss<<i<<": "<<(*errors)(i)<<std::endl;}
    phoneme_corpus->phoneme_stiffness.Set(phonemes(1)->phoneme_sample->name,errors);
    LOG::filecout(ss.str());
}
//#####################################################################
// Function Gather_Phoneme_Statistics
//#####################################################################
template<class T> void PHONEME_CORPUS<T>::
Gather_Phoneme_Statistics()
{
    for(HASHTABLE_ITERATOR<std::string,ARRAY<PHONEME_SEGMENT<T>*>*> iterator(string_to_phoneme_index);iterator.Valid();iterator.Next())
        Gather_Phoneme_Statistics_Helper(this,iterator.Data());
}
//#####################################################################
// Function Get_Next_Kind_Of_Phoneme
//#####################################################################
template<class T> void PHONEME_CORPUS<T>::
Get_Next_Kind_Of_Phoneme(PHONEME_SEGMENT<T>& phoneme_segment)
{
    int index=0;
    for(int i=1;i<=phoneme_labels.m;i++) if(!phoneme_segment.phoneme_sample->name.compare(phoneme_labels(i))){index=i;break;}
    if(!index) return;
    index=(index%phoneme_labels.m)+1;
    ARRAY<PHONEME_SEGMENT<T>*>& phoneme_segments=*Get_Phonemes(phoneme_labels(index));
    phoneme_segment.Set_Sample(*phoneme_segments(1)->phoneme_sample);
    phoneme_segment.phoneme_sample_name=phoneme_segments(1)->phoneme_sample_name;
}
//#####################################################################
// Function Get_Previous_Kind_Of_Phoneme
//#####################################################################
template<class T> void PHONEME_CORPUS<T>::
Get_Previous_Kind_Of_Phoneme(PHONEME_SEGMENT<T>& phoneme_segment)
{
    int index=0;
    for(int i=1;i<=phoneme_labels.m;i++) if(!phoneme_segment.phoneme_sample->name.compare(phoneme_labels(i))){index=i;break;}
    if(!index) return;
    index=(index==1?phoneme_labels.m:index-1);
    ARRAY<PHONEME_SEGMENT<T>*>& phoneme_segments=*Get_Phonemes(phoneme_labels(index));
    phoneme_segment.Set_Sample(*phoneme_segments(1)->phoneme_sample);
    phoneme_segment.phoneme_sample_name=phoneme_segments(1)->phoneme_sample_name;
}
//#####################################################################
// Function Get_Next_Phoneme
//#####################################################################
template<class T> void PHONEME_CORPUS<T>::
Get_Next_Phoneme(PHONEME_SEGMENT<T>& phoneme_segment)
{
    ARRAY<PHONEME_SEGMENT<T>*>& phoneme_segments=*Get_Phonemes(phoneme_segment.phoneme_sample->name);
    int index=0;
    for(int i=1;i<=phoneme_segments.m;i++) if(!phoneme_segment.phoneme_sample_name.compare(phoneme_segments(i)->phoneme_sample_name)){index=i;break;}
    if(!index) return;
    index=(index%phoneme_segments.m)+1;
    phoneme_segment.Set_Sample(*phoneme_segments(index)->phoneme_sample);
    phoneme_segment.phoneme_sample_name=phoneme_segments(index)->phoneme_sample_name;
}
//#####################################################################
// Function Get_Next_Phoneme
//#####################################################################
template<class T> void PHONEME_CORPUS<T>::
Get_Previous_Phoneme(PHONEME_SEGMENT<T>& phoneme_segment)
{
    ARRAY<PHONEME_SEGMENT<T>*>& phoneme_segments=*Get_Phonemes(phoneme_segment.phoneme_sample->name);
    int index=0;
    for(int i=1;i<=phoneme_segments.m;i++) if(!phoneme_segment.phoneme_sample_name.compare(phoneme_segments(i)->phoneme_sample_name)){index=i;break;}
    if(!index) return;
    index=(index==1?phoneme_segments.m:index-1);
    phoneme_segment.Set_Sample(*phoneme_segments(index)->phoneme_sample);
    phoneme_segment.phoneme_sample_name=phoneme_segments(index)->phoneme_sample_name;
}
//#####################################################################
// Function Build_Corpus
//#####################################################################
template<class T> void PHONEME_CORPUS<T>::
Build_Corpus(const std::string& corpus_directory)
{
#if defined(__linux__) || defined(__CYGWIN__)
    DIR *dir;
    dir=opendir(corpus_directory.c_str());
    if(!dir) PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Can't open %s for read",corpus_directory.c_str()));
    dirent *entry;
    while((entry=readdir(dir))){
        std::string filename(entry->d_name);
        if(FILE_UTILITIES::Is_Phoneme_File(filename)) Add_Phoneme_Segment(corpus_directory+"/"+filename);}
    closedir(dir);
    
    for(HASHTABLE_ITERATOR<std::string,ARRAY<PHONEME_SEGMENT<T>*>*> iterator(string_to_phoneme_index);iterator.Valid();iterator.Next()){
        ARRAY<PHONEME_SEGMENT<T>*>& phoneme_segment_list=*iterator.Data();
        for(int i=1;i<=phoneme_segment_list.m;i++) phoneme_segment_list(i)->Set_Custom_Interpolation(interpolation);}

    Sort(phoneme_labels);
#else
    LOG::filecout("THIS PHONEME_CORPUS FUNCTION IS NOT DEFINED EXCEPT UNDER LINUX!\n");
#endif
}

//#####################################################################
// Function Phoneme_Metric_Helper
//#####################################################################
template<class T> T PHONEME_CORPUS<T>::
Phoneme_Metric_Helper(PHONEME_ARRANGEMENT<T>& phoneme_arrangement,const ARRAY<T>& weights)
{
    T total=0;
    // sample at ~1 frame
    phoneme_arrangement.Update_Valid_Segment();
    for(T t=phoneme_arrangement.valid_segment->min_corner.x;t<=phoneme_arrangement.valid_segment->max_corner.x;t+=(T).00833333){
        total+=phoneme_arrangement.Discrepancy_Between_Inputs(t,weights);}
        //if(isnan(total)){//panic
        //total+=phoneme_arrangement.Discrepancy_Between_Inputs(t,weights);
        //return FLT_MAX;}}
    T time_scaling_penalty=0;
    int count=0;
    for(int i=1;i<=phoneme_arrangement.list.m;i++){
        PHONEME_SEGMENT<T>& phoneme_segment=phoneme_arrangement.list(i);
        if(phoneme_segment.phoneme_sample->frame_length>2){
#if ALTERNATE_PHONEMES
            T ratio=phoneme_segment.time_scaling;
#else
            T ratio=(phoneme_segment.time_segment.Length()/phoneme_segment.phoneme_sample->time_length);
#endif
            count++;
            time_scaling_penalty+=(ratio+1/ratio);}}
    //time_scaling_penalty=(count?time_scaling_penalty/count:1);
    return total+100*time_scaling_penalty;
}
//#####################################################################
// Function Phoneme_Metric_Helper
//#####################################################################
template<class T> T PHONEME_CORPUS<T>::
Phoneme_Metric_Helper_Scaled(PHONEME_ARRANGEMENT<T>& phoneme_arrangement,const ARRAY<T>& scalings,const ARRAY<T>& weights)
{
    // go through and get squared error sum
    for(int i=1;i<=scalings.m;i++) phoneme_arrangement.list(i).scaling=scalings(i);
    T total=0;
    // sample at ~1 frame
    phoneme_arrangement.Update_Valid_Segment();
    for(T t=phoneme_arrangement.valid_segment->min_corner.x;t<=phoneme_arrangement.valid_segment->max_corner.x;t+=(T).00833333){
        total+=phoneme_arrangement.Discrepancy_Between_Inputs(t,weights);}
        //if(isnan(total)){//panic
        //    return FLT_MAX;
            //total+=phoneme_arrangement.Discrepancy_Between_Inputs(t); *((int*)0)=1;
        //}}
    return total;
}
//#####################################################################
// Function Get_Optimal_Scalings
//#####################################################################
template<class T> T PHONEME_CORPUS<T>::
Get_Optimal_Scalings(PHONEME_ARRANGEMENT<T>& phoneme_arrangement,const ARRAY<T>& weights)
{
    /* Set scaling; try result.  Use FACE_CONTROL_PARAMETERS to compute distances in a scaling-aware way.
       Use FACE_CONTROL_PARAMETERS to do SCALING; have it (return scaled values?  scale itself?)
       Point-by-point: normal parameters: simply scale from zero.
                       rigid body params:

    */
    // Scaling search
    // go through and get squared error sum
    ARRAY<T> optimal_scalings(phoneme_arrangement.list.m);
    ARRAYS_COMPUTATIONS::Fill(optimal_scalings,(T)1);
    T result=Phoneme_Metric_Helper_Scaled(phoneme_arrangement,optimal_scalings,weights);
    T tolerance=(T)1e-4;
    T change=tolerance+1;
    const T step_size_threshold=(T)1e-5;
    while(change>tolerance){
        change=result;
        for(int i=1;i<=optimal_scalings.m;i++){
            T location=optimal_scalings(i);
            T step_size=(T).2;
            T left,right;
            left=right=result+1;
            while(step_size>step_size_threshold){
                if(left>=result&&right>=result){
                    step_size/=2;
                    optimal_scalings(i)=location-step_size;
                    left=Phoneme_Metric_Helper_Scaled(phoneme_arrangement,optimal_scalings,weights);
                    optimal_scalings(i)=location+step_size;
                    right=Phoneme_Metric_Helper_Scaled(phoneme_arrangement,optimal_scalings,weights);}
                else if(left<=result){
                    right=result;
                    result=left;
                    location-=step_size;
                    optimal_scalings(i)=location-step_size;
                    left=Phoneme_Metric_Helper_Scaled(phoneme_arrangement,optimal_scalings,weights);}
                else{
                    left=result;
                    result=right;
                    location+=step_size;
                    optimal_scalings(i)=location+step_size;
                    right=Phoneme_Metric_Helper_Scaled(phoneme_arrangement,optimal_scalings,weights);}}
            optimal_scalings(i)=location;
            result=Phoneme_Metric_Helper_Scaled(phoneme_arrangement,optimal_scalings,weights);}
        change=abs(change-result);}
    for(int i=1;i<=optimal_scalings.m;i++) {std::stringstream ss;ss<<"Scaling "<<i<<": "<<optimal_scalings(i)<<std::endl;LOG::filecout(ss.str());}
    return result;
}
//#####################################################################
// Function Brute_Force_Phoneme_Arrangement_Helper
//#####################################################################
template<class T> bool PHONEME_CORPUS<T>::
Brute_Force_Phoneme_Arrangement_Helper(const ARRAY<ARRAY<PHONEME_SEGMENT<T>*> >& phoneme_candidates,const int position,PHONEME_ARRANGEMENT<T>& phoneme_arrangement,T &best_result,ARRAY<int>& indices,const ARRAY<T>& weights)
{
    const T frame_rate=(T)120.0;
    const T stiffness=(T)1.0;
    GRID<VECTOR<T,1> > sampling_grid(1+(int)(phoneme_arrangement.valid_segment->Size()*frame_rate),*phoneme_arrangement.valid_segment);
    bool better=false;
    for(int i=1;i<=phoneme_candidates(position).m;i++){
        phoneme_arrangement.list(position).Set_Sample(*phoneme_candidates(position)(i)->phoneme_sample);
        phoneme_arrangement.list(position).phoneme_sample_name=phoneme_candidates(position)(i)->phoneme_sample_name;
        T result;
        if(position==indices.m) {//LOG::cout<<"Trying "<<std::endl;for(int j=1;j<=indices.m;j++)LOG::cout<<phoneme_arrangement.list(j).phoneme_sample_name<<"
                                 //";LOG::cout<<std::endl;
            VECTOR_ND<T> scalings=phoneme_arrangement.Optimal_Scaling(sampling_grid,stiffness);
            ARRAY<T> scalings_list(scalings.n);
            for(int j=1;j<=scalings.n;j++) scalings_list(j)=scalings(j);
            result=Phoneme_Metric_Helper_Scaled(phoneme_arrangement,scalings_list,weights);
            //result=Get_Optimal_Scalings(phoneme_arrangement);
            if(result<best_result){better=true;best_result=result;indices(position)=i;}}
        else if(Brute_Force_Phoneme_Arrangement_Helper(phoneme_candidates,position+1,phoneme_arrangement,best_result,indices,weights)){better=true;indices(position)=i;}}
    return better;
}
//#####################################################################
// Function Brute_Force_Phoneme_Arrangement
//#####################################################################
template<class T> PHONEME_ARRANGEMENT<T>* PHONEME_CORPUS<T>::
Brute_Force_Phoneme_Arrangement(const ARRAY<std::string>& phonemes,const ARRAY<T>& starting_times,const ARRAY<T>& lengths,const ARRAY<T>& leading_times,const ARRAY<T>& trailing_times)
{
    assert(face_control_parameters);assert(string_to_phoneme_index.Size());
    ARRAY<ARRAY<PHONEME_SEGMENT<T>*> > phoneme_candidates;
    ARRAY<int> indices(phonemes.m),best_indices(phonemes.m);
    T best_result=FLT_MAX;
    PHONEME_ARRANGEMENT<T>* phoneme_arrangement=new PHONEME_ARRANGEMENT<T>;
    phoneme_arrangement->Set_Custom_Interpolation(interpolation);
    phoneme_arrangement->Set_Face_Control_Parameters(face_control_parameters);
    phoneme_arrangement->list.Resize(phonemes.m);
    for(int i=1;i<=phoneme_arrangement->list.m;i++){
        PHONEME_SEGMENT<T>& phoneme_segment=phoneme_arrangement->list(i);
#if ALTERNATE_PHONEMES
        phoneme_segment.peak=starting_times(i)+lengths(i)/2;
#endif
        phoneme_segment.time_segment.min_corner.x=starting_times(i);
        phoneme_segment.time_segment.max_corner.x=phoneme_segment.time_segment.min_corner.x+lengths(i);
        phoneme_segment.leading_context_duration=leading_times(i);
        phoneme_segment.trailing_context_duration=trailing_times(i);}
    for(int i=1;i<=phonemes.m;i++){phoneme_candidates.Append(*Get_Phonemes(phonemes(i)));if(!phoneme_candidates(i).m) return 0;indices(i)=1;}
    phoneme_arrangement->Update_Valid_Segment();
    // Yankee-search!
    Brute_Force_Phoneme_Arrangement_Helper(phoneme_candidates,1,*phoneme_arrangement,best_result,indices,weights);
    for(int i=1;i<=phoneme_arrangement->list.m;i++){
        std::stringstream ss;ss<<"Phoneme "<<phoneme_candidates(i)(indices(i))->phoneme_sample_name<<" at index "<<indices(i)<<std::endl;LOG::filecout(ss.str());
        PHONEME_SEGMENT<T>& phoneme_segment=phoneme_arrangement->list(i);
        phoneme_segment.Set_Sample(*phoneme_candidates(i)(indices(i))->phoneme_sample);
        phoneme_segment.phoneme_sample_name=phoneme_candidates(i)(indices(i))->phoneme_sample_name;}
    Get_Optimal_Scalings(*phoneme_arrangement,weights);
    return phoneme_arrangement;
}
//#####################################################################
// Function Produce_Phoneme_Arrangement
//#####################################################################
template<class T> PHONEME_ARRANGEMENT<T>* PHONEME_CORPUS<T>::
Produce_Phoneme_Arrangement(const ARRAY<std::string>& phonemes,const ARRAY<T>& starting_times,const ARRAY<T>& lengths,const ARRAY<T>& leading_times,const ARRAY<T>& trailing_times,const T decay,const int steps,const T Tzero)
{
    // for the moment, no rescaling, just match
    ARRAY<ARRAY<PHONEME_SEGMENT<T>*> > phoneme_candidates;
    ARRAY<int> indices(phonemes.m);

    PHONEME_ARRANGEMENT<T>* phoneme_arrangement=new PHONEME_ARRANGEMENT<T>();
    phoneme_arrangement->Set_Custom_Interpolation(interpolation);
    phoneme_arrangement->Set_Face_Control_Parameters(face_control_parameters);
    phoneme_arrangement->list.Resize(phonemes.m);
    for(int i=1;i<=phonemes.m;i++){phoneme_candidates.Append(*Get_Phonemes(phonemes(i)));if(!phoneme_candidates(i).m) return 0;indices(i)=1;}
    for(int i=1;i<=phoneme_arrangement->list.m;i++){
        PHONEME_SEGMENT<T>& phoneme_segment=phoneme_arrangement->list(i);
        phoneme_segment.Set_Sample(*phoneme_candidates(i)(1)->phoneme_sample);
        phoneme_segment.phoneme_sample_name=phoneme_candidates(i)(1)->phoneme_sample_name;
#if ALTERNATE_PHONEMES
        phoneme_segment.peak=starting_times(i)+lengths(i)/2;
#endif
        phoneme_segment.time_segment.min_corner.x=starting_times(i);
        phoneme_segment.time_segment.max_corner.x=phoneme_segment.time_segment.min_corner.x+lengths(i);
        phoneme_segment.leading_context_duration=leading_times(i);
        phoneme_segment.trailing_context_duration=trailing_times(i);}
    

    RANDOM_NUMBERS<T> rn;
    //T frame_rate=(T)120.0;
    //T stiffness=(T)1.0;
    phoneme_arrangement->Update_Valid_Segment();
    //GRID_1D<T> sampling_grid(1+(int)(phoneme_arrangement->valid_segment->Length()*frame_rate),*phoneme_arrangement->valid_segment);
    //VECTOR_ND<T> scalings=phoneme_arrangement->Optimal_Scaling(sampling_grid,stiffness);
    //ARRAY<T> scalings_list(scalings.n);
    //for(int i=1;i<=scalings.n;i++) scalings_list(i)=scalings(i);
    //T last_result=Phoneme_Metric_Helper_Scaled(*phoneme_arrangement,scalings_list);
    //T last_result=Get_Optimal_Scalings(*phoneme_arrangement);
    T last_result=Phoneme_Metric_Helper(*phoneme_arrangement,weights);
    T Temperature=Tzero;
    int time=0;
    while(time<=steps){
        int index=rn.Get_Uniform_Integer(1,indices.m);
        int phoneme_index=rn.Get_Uniform_Integer(1,phoneme_candidates(index).m);
        int old_phoneme_index=indices(index);
        PHONEME_SEGMENT<T>& phoneme_segment=phoneme_arrangement->list(index);
        phoneme_segment.Set_Sample(*phoneme_candidates(index)(phoneme_index)->phoneme_sample);
        phoneme_segment.phoneme_sample_name=phoneme_candidates(index)(phoneme_index)->phoneme_sample_name;
        //LOG::cout<<"Trying "<<std::endl;for(int j=1;j<=indices.m;j++)LOG::cout<<phoneme_arrangement->list(j).phoneme_sample_name<<" ";LOG::cout<<std::endl;
    
        // Testing this
        /*scalings=phoneme_arrangement->Optimal_Scaling(sampling_grid,stiffness);
        for(int i=1;i<=scalings.n;i++) scalings_list(i)=scalings(i);
        T result=Phoneme_Metric_Helper_Scaled(*phoneme_arrangement,scalings_list);*/
        T result=Phoneme_Metric_Helper(*phoneme_arrangement,weights);
        //T result=Get_Optimal_Scalings(*phoneme_arrangement);
        //LOG::cout<<"Temperature is "<<Temperature<<", last was "<<last_result<<"... "<<result<<std::endl;
        if(time%100==0){std::stringstream ss;ss<<"At time step "<<time<<" temperature is "<<Temperature<<", last was "<<last_result<<"... "<<result<<std::endl;LOG::filecout(ss.str());}
        if(rn.Get_Uniform_Number((T)0,(T)1)<exp((last_result-result)/Temperature)){indices(index)=phoneme_index;last_result=result;}
        else{phoneme_segment.Set_Sample(*phoneme_candidates(index)(old_phoneme_index)->phoneme_sample);
        phoneme_segment.phoneme_sample_name=phoneme_candidates(index)(old_phoneme_index)->phoneme_sample_name;}
        Temperature*=decay;
        time++;}
        
    /*scalings=phoneme_arrangement->Optimal_Scaling(sampling_grid,stiffness);
    LOG::cout<<scalings<<std::endl;
    T rescale=scalings.Sum()/scalings.Magnitude_Squared();
    for(int i=1;i<=phoneme_arrangement->list.m;i++){
        PHONEME_SEGMENT<T>& phoneme_segment=phoneme_arrangement->list(i);
        phoneme_segment.scaling*=rescale;}*/
            
    return phoneme_arrangement;
        
    // coordinate descent isn't actually a good approach, but try it for a first cut
/*    T tolerance=1e-5;
    T change=tolerance+1;
    while(change>tolerance){
        change=Phoneme_Metric_Helper(phoneme_candidates,indices,starting_times,lengths,leading_times,trailing_times);
        for(int i=1;i<=phonemes.m;i++){
            int best_index=indices(i);
            T result=Phoneme_Metric_Helper(phoneme_candidates,indices,starting_times,lengths,leading_times,trailing_times);
            for(int j=1;j<=phoneme_candidates(i).m;j++){
                if(j==indices(i)) continue;
                ARRAY<int> local=indices;
                local(i)=j;
                T local_result=Phoneme_Metric_Helper(phoneme_candidates,local,starting_times,lengths,leading_times,trailing_times);
                if(local_result<result){result=local_result;best_index=j;}}
            indices(i)=best_index;}
        change=abs(change-Phoneme_Metric_Helper(phoneme_candidates,indices,starting_times,lengths,leading_times,trailing_times));}
    PHONEME_ARRANGEMENT<T> *phoneme_arrangement=new PHONEME_ARRANGEMENT<T>;
    phoneme_arrangement->list.Resize(phonemes.m);
    for(int i=1;i<=phoneme_arrangement->list.m;i++){
        LOG::cout<<"Phoneme "<<phoneme_candidates(i)(indices(i))->phoneme_sample_name<<" at index "<<indices(i)<<std::endl;
        PHONEME_SEGMENT<T>& phoneme_segment=phoneme_arrangement->list(i);
        phoneme_segment=PHONEME_SEGMENT<T>(*phoneme_candidates(i)(indices(i)));
        phoneme_segment.time_segment.min_corner.x=starting_times(i);
        phoneme_segment.time_segment.max_corner.x=phoneme_segment.time_segment.min_corner.x+lengths(i);
        phoneme_segment.leading_context_duration=leading_times(i);
        phoneme_segment.trailing_context_duration=trailing_times(i);}
        return phoneme_arrangement;*/
}

template class PHONEME_CORPUS<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PHONEME_CORPUS<double>;
#endif
