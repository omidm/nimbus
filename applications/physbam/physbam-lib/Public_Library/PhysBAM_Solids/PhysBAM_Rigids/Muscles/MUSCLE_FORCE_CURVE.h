//#####################################################################
// Copyright 2005, Ron Fedkiw, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MUSCLE_FORCE_CURVE
//#####################################################################
#ifndef __MUSCLE_FORCE_CURVE__
#define __MUSCLE_FORCE_CURVE__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
#include <fstream>

namespace PhysBAM{

template<class T>
class MUSCLE_FORCE_CURVE:public NONCOPYABLE
{
public:
    GRID<VECTOR<T,1> > passive_force_grid,passive_force_slope_grid;
    ARRAY<T,VECTOR<int,1> > passive_force,passive_force_slope;
    GRID<VECTOR<T,1> > tendon_force_grid,tendon_force_slope_grid,tendon_length_grid;
    ARRAY<T,VECTOR<int,1> > tendon_force,tendon_force_slope,tendon_length;
    GRID<VECTOR<T,1> > active_force_grid,active_force_slope_grid;
    ARRAY<T,VECTOR<int,1> > active_force,active_force_slope;
    GRID<VECTOR<T,1> > velocity_grid;
    ARRAY<T,VECTOR<int,1> > velocity_curve;
private:
    INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,T>* passive_interpolation; 
    INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,T>* active_interpolation;
    INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,T>* tendon_force_interpolation;
    INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,T>* velocity_interpolation;
    LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,T> default_interpolation;

public:
    MUSCLE_FORCE_CURVE();
    virtual ~MUSCLE_FORCE_CURVE();

    void Set_Custom_Passive_Interpolation(INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,T>* passive_interpolation_input)
    {passive_interpolation=passive_interpolation_input;}

    void Set_Custom_Active_Interpolation(INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,T>* active_interpolation_input)
    {active_interpolation=active_interpolation_input;}

    void Set_Custom_Tendon_Force_Interpolation(INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,T>* tendon_force_interpolation_input)
    {tendon_force_interpolation=tendon_force_interpolation_input;}

    T Passive_Force(const T length) const // assumes normalized length
    {return passive_interpolation->Clamped_To_Array(passive_force_grid,passive_force,VECTOR<T,1>(length));}

    T Passive_Force_Slope(const T length) const // assumes normalized length
    {return passive_interpolation->Clamped_To_Array(passive_force_slope_grid,passive_force_slope,VECTOR<T,1>(length));}

    T Active_Force(const T length) const // assumes normalized length
    {return active_interpolation->Clamped_To_Array(active_force_grid,active_force,VECTOR<T,1>(length));}

    T Active_Force_Slope(const T length) const // assumes normalized length
    {return active_interpolation->Clamped_To_Array(active_force_slope_grid,active_force_slope,VECTOR<T,1>(length));}

    T Velocity_Scale(const T velocity) const // assumes normalized length and velocity
    {return velocity_interpolation->Clamped_To_Array(velocity_grid,velocity_curve,VECTOR<T,1>(velocity));}

    T Tendon_Force(const T length) const // assumes normalized length
    {return tendon_force_interpolation->Clamped_To_Array(tendon_force_grid,tendon_force,VECTOR<T,1>(length));}

    T Tendon_Force_Slope(const T length) const // assumes normalized length
    {return tendon_force_interpolation->Clamped_To_Array(tendon_force_slope_grid,tendon_force_slope,VECTOR<T,1>(length));}

    T Tendon_Length(const T force) const // assumes normalized force
    {return tendon_force_interpolation->Clamped_To_Array(tendon_length_grid,tendon_length,VECTOR<T,1>(force));}

//#####################################################################
    void Initialize(const std::string& data_directory);
    void Load_Data(const std::string& prefix,GRID<VECTOR<T,1> >& grid,ARRAY<T,VECTOR<int,1> >& values);
private:
    void Compute_Slopes(const GRID<VECTOR<T,1> >& grid,const ARRAY<T,VECTOR<int,1> >& values,GRID<VECTOR<T,1> >& slope_grid,ARRAY<T,VECTOR<int,1> >& slopes);
//#####################################################################
};
}
#endif
