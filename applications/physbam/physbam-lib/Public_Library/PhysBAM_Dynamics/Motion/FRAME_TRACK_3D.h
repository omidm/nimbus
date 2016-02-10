//#####################################################################
// Copyright 2005-2007, Jared Go, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __FRAME_TRACK__
#define __FRAME_TRACK__
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/wrap.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>

namespace PhysBAM{

template<class T>
class FRAME_TRACK_3D
{
    typedef VECTOR<T,3> TV;
public:
    std::string name;
    GRID<TV> time_grid;
    ARRAY<FRAME<TV> ,VECTOR<int,1> > trajectory;
    bool periodic;

    FRAME_TRACK_3D(const int samples,const T first_time,const T last_time)
        :time_grid(samples,first_time,last_time),periodic(false)
    {trajectory.Resize(time_grid);}    

    FRAME<TV> Frame(T t)
    {if(periodic) t=wrap(t,time_grid.xmin,time_grid.xmax);
    else t=clamp(t,time_grid.xmin,time_grid.xmax);
    int index=time_grid.Clamped_Index_End_Minus_One(VECTOR<T,1>(t)).x;
    return FRAME<TV>::Interpolation(trajectory(index),trajectory(index+1),(t-time_grid.x(index))*time_grid.one_over_dx);}

    VECTOR<T,3> Angular_Velocity(T t)
    {if(periodic) t=wrap(t,time_grid.xmin,time_grid.xmax);
    else if(t<time_grid.xmin || t>time_grid.xmax) return VECTOR<T,3>();
    int index=time_grid.Clamped_Index_End_Minus_One(VECTOR<T,1>(t)).x;
    return time_grid.one_over_dx*(trajectory(index+1).r*trajectory(index).r.Inverse()).Rotation_Vector();}

//#####################################################################
};
}
#endif
