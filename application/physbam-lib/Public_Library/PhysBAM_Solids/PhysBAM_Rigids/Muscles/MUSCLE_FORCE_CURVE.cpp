//#####################################################################
// Copyright 2005, Ron Fedkiw, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_FORCE_CURVE.h>
using namespace PhysBAM;
template<class T> MUSCLE_FORCE_CURVE<T>::
MUSCLE_FORCE_CURVE()
    :passive_interpolation(&default_interpolation),active_interpolation(&default_interpolation),tendon_force_interpolation(&default_interpolation),velocity_interpolation(&default_interpolation)
{
}
template<class T> MUSCLE_FORCE_CURVE<T>::
~MUSCLE_FORCE_CURVE()
{
}
template<class T> void MUSCLE_FORCE_CURVE<T>::
Initialize(const std::string& data_directory)
{
    Load_Data(data_directory+"/Muscle_Curves/passive_cubic",passive_force_grid,passive_force);
    Load_Data(data_directory+"/Muscle_Curves/active_cubic",active_force_grid,active_force);
    Load_Data(data_directory+"/Muscle_Curves/tendon_cubic",tendon_force_grid,tendon_force);
    velocity_grid.Initialize(1001,-2,2);velocity_curve.Resize(velocity_grid.Domain_Indices());
    for(int v=1;v<=1001;v++){T velocity=velocity_grid.Axis_X(v,1);velocity_curve(v)=(velocity<-1)?0:(T).54*atan((T)5.69*velocity+(T).51)+(T).745;}
    tendon_length_grid.Initialize(1001,0,3.5);tendon_length.Resize(tendon_length_grid.Domain_Indices());
    LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,1> >,T>::Compute_Inverse_Map(tendon_force_grid,tendon_force,tendon_length_grid,tendon_length);
    Compute_Slopes(passive_force_grid,passive_force,passive_force_slope_grid,passive_force_slope);
    Compute_Slopes(active_force_grid,active_force,active_force_slope_grid,active_force_slope);
    Compute_Slopes(tendon_force_grid,tendon_force,tendon_force_slope_grid,tendon_force_slope);
}
/*function above from "Model the leg cycling movement with neural oscillator", Dingguo Zhang; Kuanyi Zhu; Hang Zheng;Systems, Man and Cybernetics, 2004 IEEE International Conference on Volume 1,  10-13 Oct. 2004 Page(s):740 - 744 vol.1 */
template<class T> void MUSCLE_FORCE_CURVE<T>::
Load_Data(const std::string& prefix,GRID<VECTOR<T,1> >& grid,ARRAY<T,VECTOR<int,1> >& values)
{
    {std::istream* input=FILE_UTILITIES::Safe_Open_Input(prefix+".grid",false);
    *input>>grid.counts.x>>grid.domain.min_corner.x>>grid.domain.max_corner.x;
    grid.Initialize(grid.counts.x,grid.domain.min_corner.x,grid.domain.max_corner.x);
    values.Resize(grid.Domain_Indices());
    delete input;}
    std::istream* input(FILE_UTILITIES::Safe_Open_Input(prefix+".values",false));
    for(int i=1;i<=values.counts.x;i++) *input>>values(i);
    delete input;
}
template<class T> void MUSCLE_FORCE_CURVE<T>::
Compute_Slopes(const GRID<VECTOR<T,1> >& grid,const ARRAY<T,VECTOR<int,1> >& values,GRID<VECTOR<T,1> >& slope_grid,ARRAY<T,VECTOR<int,1> >& slopes)
{
    assert(!grid.Is_MAC_Grid());
    slope_grid=grid.Get_MAC_Grid();
    slopes.Resize(1,slope_grid.counts.x);
    for(int i=1;i<=slopes.counts.x;i++) slopes(i)=(values(i+1)-values(i))*grid.one_over_dX.x;
}
template class MUSCLE_FORCE_CURVE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MUSCLE_FORCE_CURVE<double>;
#endif
