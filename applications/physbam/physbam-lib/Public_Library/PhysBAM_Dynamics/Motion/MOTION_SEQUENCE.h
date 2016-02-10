//#####################################################################
// Copyright 2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MOTION_SEQUENCE__
#define __MOTION_SEQUENCE__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#include <fstream>

namespace PhysBAM{

template<class T,class T2>
class MOTION_SEQUENCE
{
    typedef VECTOR<T,1> TV;
public:
    int max_size;
    GRID<TV> time_grid;
    ARRAY<ARRAY<T2,VECTOR<int,1> > > trajectories;
    ARRAY<ARRAY<bool,VECTOR<int,1> > > valid;
    ARRAY<std::string> names;

    //INTERPOLATION_UNIFORM<GRID<TV>,T2>* interpolation;
    //LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T2> default_interpolation;
    HASHTABLE<std::string,int> name_to_track_index;

    MOTION_SEQUENCE()
        :max_size(100)
        //:interpolation(&default_interpolation)
    {}

    MOTION_SEQUENCE(const MOTION_SEQUENCE<T,T2>& motion_input)
    {
        //if(motion_input.interpolation==&motion_input.default_interpolation) interpolation=&default_interpolation;
        //else interpolation=motion_input.interpolation;
        //interpolation=&default_interpolation;
        *this=motion_input;
    }

    MOTION_SEQUENCE<T,T2>& operator=(const MOTION_SEQUENCE<T,T2>& motion_input)
    {time_grid=motion_input.time_grid;
    trajectories=motion_input.trajectories;
    valid=motion_input.valid;
    names=motion_input.names;
    name_to_track_index=motion_input.name_to_track_index;
    //if(motion_input.interpolation==&motion_input.default_interpolation) interpolation=&default_interpolation;
    //else interpolation=motion_input.interpolation;
    return *this;}

    void Copy(const VECTOR<int,1>& time_index,const MOTION_SEQUENCE<T,T2>& motion_input,const VECTOR<int,1>& time_index_input)
    {for(int i=1;i<=trajectories.m;i++){trajectories(i)(time_index)=motion_input.trajectories(i)(time_index_input);valid(i)(time_index)=motion_input.valid(i)(time_index_input);}}

    void Initialize(const int number_of_tracks,const GRID<TV>& time_grid_input)
    {time_grid=time_grid_input;trajectories.Resize(number_of_tracks);valid.Resize(number_of_tracks);names.Resize(number_of_tracks);
    for(int i=1;i<=number_of_tracks;i++){trajectories(i).Resize(1,time_grid.counts(1));valid(i).Resize(1,time_grid.counts(1));}}
    
    void Resize(const GRID<TV>& time_grid_input)
    {time_grid=time_grid_input;for(int i=1;i<=trajectories.m;i++){trajectories(i).Resize(1,time_grid.counts(1));valid(i).Resize(1,time_grid.counts(1));}}

    void Rescale(const GRID<TV>& time_grid_input)
    {MOTION_SEQUENCE<T,T2> motion_sequence=*this;motion_sequence.Resize(time_grid_input);
    for(int i=1;i<=names.m;i++) for(typename GRID<TV>::CELL_ITERATOR iterator(time_grid_input);iterator.Valid();iterator.Next()){VECTOR<int,1> cell=iterator.Cell_Index();
        int lower=time_grid.Cell(iterator.Location()-.5*time_grid.dX(1),1)(1);
        T alpha=((iterator.Location()-time_grid.X(VECTOR<int,1>(lower)))/time_grid.dX)(1);
        //int lower=(int)((T)cell(1)/(T)time_grid_input.counts(1)*(time_grid.counts(1)-1)+1);
        //T alpha=(T)cell(1)/(T)time_grid_input.counts(1)*(time_grid.counts(1)-1)+1-lower;
        if(lower+1>time_grid.counts(1)){alpha=1;lower=time_grid.counts(1)-1;}
        if(lower<=0){alpha=0;lower=1;}
        VECTOR<T,3> alpha_vector=alpha?(trajectories(i)(lower+1).t):(VECTOR<T,3>());
        motion_sequence.trajectories(i)(cell).t=trajectories(i)(lower).t*(1-alpha)+alpha*alpha_vector;
        ROTATION<VECTOR<T,3> > alpha_r=alpha?(trajectories(i)(lower+1).r):(trajectories(i)(lower).r);
        motion_sequence.trajectories(i)(cell).r=ROTATION<VECTOR<T,3> >::Spherical_Linear_Interpolation(trajectories(i)(lower).r,alpha_r,alpha);}
    *this=motion_sequence;}

    void Normalize(const int normalize=3)
    {for(int i=1;i<=time_grid.counts(1);i++){VECTOR<T,3> root=trajectories(normalize)(i).t;for(int j=1;j<=names.m;j++) trajectories(Track_Index(names(j)))(i).t-=root;}}

    void Append_Trajectories(const ARRAY<T2>& data,const int normalize=0)
    {int size=time_grid.counts(1)+1;
    if(size>max_size){size=max_size/2;
        for(int i=1;i<size;i++) for(int j=1;j<trajectories.m;j++){
            trajectories(j)(i)=trajectories(j)(i+size+1);
            valid(j)(i)=valid(j)(i+size+1);}}
    VECTOR<int,1> counts(size);
    GRID<TV> time_grid_new(counts,time_grid.domain);
    Resize(time_grid_new);
    for(int i=1;i<=names.m;i++) trajectories(Track_Index(names(i)))(size)=data(i);
    if(normalize){VECTOR<T,3> root=trajectories(normalize)(size).t;
        for(int i=1;i<=names.m;i++) trajectories(Track_Index(names(i)))(size).t-=root;}}

    T Norm(const MOTION_SEQUENCE<T,FRAME<VECTOR<T,3> > >& motion_input)
    {T distance=0;for(int n=1;n<=names.m;n++) for(int m=1;m<=time_grid.counts(1);m++){
        T diff=(trajectories(n)(m).t-motion_input.trajectories(n)(m).t).Magnitude();
        distance+=diff*diff;}
    return distance;}

    T Norm(const MOTION_SEQUENCE<T,ARRAY<VECTOR<T,3> > >& motion_input)
    {T distance=0;for(int n=1;n<=names.m;n++) for(int m=1;m<=time_grid.counts(1);m++){
        for(int i=1;i<=trajectories(n)(m).m;i++){
            T diff=(trajectories(n)(m)(i)-motion_input.trajectories(n)(m)(i)).Magnitude();
            distance+=diff*diff;}}
    return distance;}

    void Write_Ply(ARRAY<VECTOR<T,3>,int>& points,const char* filename,VECTOR<int,3> color)
    {
        std::fstream outFile(filename,std::ios::out);
        outFile<<"ply\nformat ascii 1.0\nelement vertex "<<names.m<<std::endl;
        outFile<<"property float x\n"<<"property float y\n"<<"property float z\n";
        outFile<<"property uchar red\n"<<"property uchar green\n"<<"property uchar blue\n";
        outFile<<"end_header\n";
        for(int i=1;i<=names.m;i++){
            outFile<<points(i)(1)<<" "<<points(i)(2)<<" "<<points(i)(3)<<" ";
            outFile<<color(1)<<" "<<color(2)<<" "<<color(3)<<std::endl;}
        outFile.close();
    }

    T Frame_Matching_Cost(const int frame, const MOTION_SEQUENCE<T,ARRAY<VECTOR<T,3> > >& motion_input, const int input_frame, bool write_debug_data=false)
    {
        T cost=0,x_bar1=0,x_bar2=0,y_bar1=0,y_bar2=0;

        ARRAY<VECTOR<T,3>,int> points1(names.m),points2(names.m),points3(names.m);
        for(int i=1;i<=names.m;i++) for(int j=1;j<=trajectories(i)(frame).m;j++){
            points1(i)=trajectories(i)(frame)(j);points2(i)=motion_input.trajectories(i)(frame)(j);
            x_bar1+=trajectories(i)(frame)(j).x;x_bar2+=motion_input.trajectories(i)(input_frame)(j).x;
            y_bar1+=trajectories(i)(frame)(j).y;y_bar2+=motion_input.trajectories(i)(input_frame)(j).y;}

        T numerator=x_bar2*y_bar1-x_bar1*y_bar2,denominator=-x_bar1*x_bar2-y_bar1*y_bar2;
        numerator/=names.m;denominator/=names.m;

        for(int i=1;i<=names.m;i++) for(int j=1;j<=trajectories(i)(frame).m;j++){
            T x1=trajectories(i)(frame)(j).x,x2=motion_input.trajectories(i)(input_frame)(j).x;
            T y1=trajectories(i)(frame)(j).y,y2=motion_input.trajectories(i)(input_frame)(j).y;
            numerator+=(x1*y2-x2*y1);
            denominator+=x1*x2+y1*y2;}

        T theta=atan2(numerator,denominator);
        T x_0=(x_bar1-x_bar2*cos(theta)-y_bar2*sin(theta))/names.m;
        T y_0=(y_bar1+x_bar2*sin(theta)-y_bar2*cos(theta))/names.m;
        VECTOR<T,3> translation(x_0,y_0,0);
        ROTATION<VECTOR<T,3> > rotation(-theta,VECTOR<T,3>(0,0,1));

        for(int i=1;i<=names.m;i++)for(int j=1;j<=trajectories(i)(frame).m;j++){ 
            points3(i)=rotation.Rotate(motion_input.trajectories(i)(input_frame)(j))+translation;
            cost+=(trajectories(i)(frame)(j)-rotation.Rotate(motion_input.trajectories(i)(input_frame)(j))-translation).Magnitude_Squared();}

        if(write_debug_data){
            Write_Ply(points1,"points1.ply",VECTOR<int,3>(255,0,0));
            Write_Ply(points2,"points2.ply",VECTOR<int,3>(0,255,0));
            Write_Ply(points3,"points3.ply",VECTOR<int,3>(0,0,255));}

        return cost;
    
    }

    T Frame_Matching_Cost(const int frame, const MOTION_SEQUENCE<T,FRAME<VECTOR<T,3> > >& motion_input, const int input_frame, bool write_debug_data=false)
    {
        T cost=0,x_bar1=0,x_bar2=0,y_bar1=0,y_bar2=0;

        ARRAY<VECTOR<T,3>,int> points1(names.m),points2(names.m),points3(names.m);
        for(int i=1;i<=names.m;i++){
            points1(i)=trajectories(i)(frame).t;points2(i)=motion_input.trajectories(i)(frame).t;
            x_bar1+=trajectories(i)(frame).t.x;x_bar2+=motion_input.trajectories(i)(input_frame).t.x;
            y_bar1+=trajectories(i)(frame).t.y;y_bar2+=motion_input.trajectories(i)(input_frame).t.y;}

        T numerator=x_bar2*y_bar1-x_bar1*y_bar2,denominator=-x_bar1*x_bar2-y_bar1*y_bar2;
        numerator/=names.m;denominator/=names.m;

        for(int i=1;i<=names.m;i++){
            T x1=trajectories(i)(frame).t.x,x2=motion_input.trajectories(i)(input_frame).t.x;
            T y1=trajectories(i)(frame).t.y,y2=motion_input.trajectories(i)(input_frame).t.y;
            numerator+=(x1*y2-x2*y1);
            denominator+=x1*x2+y1*y2;}

        T theta=atan2(numerator,denominator);
        T x_0=(x_bar1-x_bar2*cos(theta)-y_bar2*sin(theta))/names.m;
        T y_0=(y_bar1+x_bar2*sin(theta)-y_bar2*cos(theta))/names.m;
        VECTOR<T,3> translation(x_0,y_0,0);
        ROTATION<VECTOR<T,3> > rotation(-theta,VECTOR<T,3>(0,0,1));

        for(int i=1;i<=names.m;i++){ 
            points3(i)=rotation.Rotate(motion_input.trajectories(i)(input_frame).t)+translation;
            cost+=(trajectories(i)(frame).t-rotation.Rotate(motion_input.trajectories(i)(input_frame).t)-translation).Magnitude_Squared();}

        if(write_debug_data){
            Write_Ply(points1,"points1.ply",VECTOR<int,3>(255,0,0));
            Write_Ply(points2,"points2.ply",VECTOR<int,3>(0,255,0));
            Write_Ply(points3,"points3.ply",VECTOR<int,3>(0,0,255));}

        return cost;
    }

    int Display_Alignment(const MOTION_SEQUENCE<T,T2>& motion_input,ARRAY<T,VECTOR<int,2> >& costs,const T skip_cost)
    {
        int m=time_grid.counts(1),n=motion_input.time_grid.counts(1);
        ARRAY<int,int> source_alignment,target_alignment;int i=m+1,j=n+1;

        while(i>=2 && j>=2){
            T match_cost_current=Frame_Matching_Cost(i-1,motion_input,j-1);
            if(costs(i,j)==costs(i-1,j-1)+match_cost_current){
                source_alignment.Append(i-1);target_alignment.Append(j-1);i--;j--;}
            else if(costs(i,j)==costs(i-1,j)+skip_cost){
                source_alignment.Append(i-1);target_alignment.Append(0);i--;}
            else if(costs(i,j)==costs(i,j-1)+skip_cost/2){
                source_alignment.Append(0);target_alignment.Append(j-1);j--;}}
        while(i>=2){
            source_alignment.Append(i-1);target_alignment.Append(0);i--;}
        while(j>=2){
            source_alignment.Append(0);target_alignment.Append(j-1);j--;}

        std::cout<<"Sequence"<<"\t"<<"Gesture\n";
        for(int k=source_alignment.m;k>=1;k--) std::cout<<source_alignment(k)<<"\t\t"<<target_alignment(k)<<std::endl;
        return source_alignment.m;
    }

    T Match_Sequence(const MOTION_SEQUENCE<T,T2>& motion_input,bool display_alignment=false)
    {
        int m=time_grid.counts(1),n=motion_input.time_grid.counts(1);
        T skip_cost=.1;     // manually adjusted parameter
        ARRAY<T,VECTOR<int,2> > costs(1,m+1,1,n+1);

        // initialization
        costs(1,1)=0;
        for(int i=2;i<=m+1;i++) costs(i,1)=costs(i-1,1)+skip_cost;
        for(int j=2;j<=n+1;j++) costs(1,j)=costs(1,j-1)+skip_cost;

        for(int i=2;i<=m+1;i++) for(int j=2;j<=n+1;j++){
            T match_cost_current=Frame_Matching_Cost(i-1,motion_input,j-1);
            costs(i,j)=min(costs(i-1,j-1)+match_cost_current,min(costs(i-1,j)+skip_cost,costs(i,j-1)+skip_cost/2));}

        if(display_alignment) Display_Alignment(motion_input,costs,skip_cost);
        return costs(m+1,n+1)/(m+n);
    }

    MOTION_SEQUENCE<T,FRAME<VECTOR<T,3> > > Velocity()
    {MOTION_SEQUENCE<T,FRAME<VECTOR<T,3> > > vel;
    vel.Initialize(names.m,GRID<TV>(time_grid.counts-VECTOR<int,TV::dimension>::All_Ones_Vector()*2,time_grid.domain));vel.names=names;
    for(int i=1;i<=names.m;i++) for(int j=2;j<=time_grid.counts(1)-1;j++) vel.trajectories(i)(j-1).t=(trajectories(i)(j+1).t-trajectories(i)(j-1).t)/2.;
    return vel;}

    MOTION_SEQUENCE<T,FRAME<VECTOR<T,3> > > Smoothed()
    {MOTION_SEQUENCE<T,FRAME<VECTOR<T,3> > > blurred;blurred.Initialize(names.m,time_grid);blurred.names=names;T scale=(1./17.);
    for(int name=1;name<=names.m;name++){
        blurred.trajectories(name)(1).t=scale*(trajectories(name)(1).t*12+trajectories(name)(3).t+trajectories(name)(2).t*4);
        blurred.trajectories(name)(2).t=scale*(trajectories(name)(1).t*5+trajectories(name)(2).t*7+trajectories(name)(4).t+trajectories(name)(3).t*4);
        blurred.trajectories(name)(time_grid.counts(1)).t=scale*(trajectories(name)(time_grid.counts(1)).t*12+trajectories(name)(time_grid.counts(1)-2).t+trajectories(name)(time_grid.counts(1)-1).t*4);
        blurred.trajectories(name)(time_grid.counts(1)-1).t=scale*(trajectories(name)(time_grid.counts(1)).t*5+trajectories(name)(time_grid.counts(1)-1).t*7+trajectories(name)(time_grid.counts(1)-3).t+trajectories(name)(time_grid.counts(1)-2).t*4);
        for(int i=3;i<=time_grid.counts(1)-2;i++) blurred.trajectories(name)(i).t=scale*(trajectories(name)(i-2).t+trajectories(name)(i-1).t*4+trajectories(name)(i).t*7+trajectories(name)(i+2).t+trajectories(name)(i+1).t*4);}
    return blurred;}

    MOTION_SEQUENCE<T,FRAME<VECTOR<T,3> > > Smoothed_Velocity()
    {MOTION_SEQUENCE<T,FRAME<VECTOR<T,3> > > blurred=Smoothed();
    return blurred.Velocity();}

    T Get_Score(const MOTION_SEQUENCE<T,T2>& motion_input)
    {return Match_Sequence(motion_input);}

    T2 X(const int marker,const T t) const
    {return Frame_X(marker,(int)((t-time_grid.domain.min_corner.x)/(time_grid.domain.max_corner.x-time_grid.domain.min_corner.x)*(time_grid.counts.x-1)+1));}
    //{if(time_grid.m==1) return trajectories(marker)(1);}
    //else return interpolation->Clamped_To_Array(time_grid,trajectories(marker),VECTOR<T,1>(t));}

    const T2& Frame_X(const int marker,const int frame) const
    {return trajectories(marker)(frame);}

    T2& Frame_X(const int marker,const int frame)
    {return trajectories(marker)(frame);}

    //void Set_Custom_Interpolation(const INTERPOLATION_UNIFORM<T,T2,GRID<VECTOR<T,1> > >* interpolation_input)
    //{interpolation=interpolation_input;}

    int Track_Index(const std::string& name)
    {int index;
    if(name_to_track_index.Get(name,index)) return index;
    else return -1;} 

    void Update_Name_Lookup()
    {name_to_track_index.Clean_Memory();for(int i=1;i<=names.m;i++) name_to_track_index.Insert(names(i),i);}

    void Set_Frame_Rate(const T frame_rate,const T initial_time=0)
    {time_grid.Initialize(time_grid.counts,initial_time,initial_time+frame_rate*(time_grid.counts(1)-1));}

    template<class RW> 
    void Read(std::istream& input_stream)
    {Read_Binary<RW>(input_stream,time_grid,trajectories,valid,names);Update_Name_Lookup();}

    template<class RW> 
    void Write(std::ostream& output_stream) const
    {Write_Binary<RW>(output_stream,time_grid,trajectories,valid,names);}

//#####################################################################
};
}
#include <PhysBAM_Dynamics/Read_Write/Motion/READ_WRITE_MOTION_SEQUENCE.h>
#endif
