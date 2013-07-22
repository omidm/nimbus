//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COLLISION_AWARE_INDEX_MAP
//##################################################################### 
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COLLISION_AWARE_INDEX_MAP.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/GENERALIZED_FLUID_MASS.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GENERALIZED_FLUID_MASS<TV>::
GENERALIZED_FLUID_MASS(const COLLISION_AWARE_INDEX_MAP<TV>& index_map_input,const T_FACE_ARRAYS_SCALAR& beta_input,const ARRAY<T>& constrained_beta_input)
    :index_map(index_map_input),beta(beta_input),constrained_beta(constrained_beta_input)
{
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void GENERALIZED_FLUID_MASS<TV>::
Compute()
{
    one_over_fluid_mass_at_faces.Resize(index_map.Number_Faces());
    for(int i=1;i<=index_map.indexed_faces.m;i++)
        one_over_fluid_mass_at_faces(i)=beta(index_map.indexed_faces(i));
    for(int i=1;i<=index_map.indexed_constraints.m;++i)
        one_over_fluid_mass_at_faces(index_map.indexed_faces.m+i)=constrained_beta(i);
    one_over_fluid_mass_at_faces*=index_map.grid.One_Over_Cell_Size();
}
//#####################################################################
// Function Second_Order_Mass_Correction
//#####################################################################
template<class TV> void GENERALIZED_FLUID_MASS<TV>::
Second_Order_Mass_Correction(const ARRAY<T,TV_INT>& phi)
{
    for(int i=1;i<=index_map.indexed_faces.m;i++){
        FACE_INDEX<d> face=index_map.indexed_faces(i);
        T phi1=phi(face.First_Cell_Index()),phi2=phi(face.Second_Cell_Index());
        if((phi1>0)==(phi2>0)) continue;
        T theta=LEVELSET_UTILITIES<T>::Theta(phi1,phi2);
        if(phi1>0) theta=1-theta; // theta relative to phi2, since it is inside.
        if(theta<(T)1e-8) theta=(T)1e-8;
        one_over_fluid_mass_at_faces(i)/=theta;}
}
//#####################################################################
// Function Compute_Two_Phase_Pressure_Jump
//#####################################################################
template<class TV> void GENERALIZED_FLUID_MASS<TV>::
Compute_For_Two_Phase_Pressure_Jump(const ARRAY<T,TV_INT>& phi,const ARRAY<T,TV_INT>& density)
{
    one_over_fluid_mass_at_faces.Resize(index_map.Number_Faces());
    for(int i=1;i<=index_map.indexed_faces.m;i++){
        FACE_INDEX<d> face=index_map.indexed_faces(i);
        T density1=density(face.First_Cell_Index()),density2=density(face.Second_Cell_Index());
        T phi1=phi(face.First_Cell_Index()),phi2=phi(face.Second_Cell_Index());
        T theta;
        if((phi1>0)==(phi2>0)) theta=.5;
        else theta=LEVELSET_UTILITIES<T>::Theta(phi1,phi2);
        one_over_fluid_mass_at_faces(i)=Inverse(density1*theta+density2*(1-theta));}
    one_over_fluid_mass_at_faces*=index_map.grid.One_Over_Cell_Size();
}
//#####################################################################
// Function Inverse_Times
//#####################################################################
template<class TV> void GENERALIZED_FLUID_MASS<TV>::
Inverse_Times(const VECTOR_ND<T>& faces_in,VECTOR_ND<T>& faces_out) const
{
    assert(faces_in.n==faces_out.n && faces_in.n==one_over_fluid_mass_at_faces.n);
    for(int i=1;i<=faces_out.n;i++) faces_out(i)=one_over_fluid_mass_at_faces(i)*faces_in(i);
}
//#####################################################################
// Function Inverse_Times_Add
//#####################################################################
template<class TV> void GENERALIZED_FLUID_MASS<TV>::
Inverse_Times_Add(const VECTOR_ND<T>& faces_in,VECTOR_ND<T>& faces_out) const
{
    assert(faces_in.n==faces_out.n && faces_in.n==one_over_fluid_mass_at_faces.n);
    for(int i=1;i<=faces_out.n;i++) faces_out(i)+=one_over_fluid_mass_at_faces(i)*faces_in(i);
}
//#####################################################################
// Function Print_Each_Matrix
//#####################################################################
template<class TV> void GENERALIZED_FLUID_MASS<TV>::
Print_Each_Matrix(int n) const
{
    OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("Bi-%i.txt",n).c_str()).Write("Bi",one_over_fluid_mass_at_faces);
}
//#####################################################################
// Function Add_Raw_Matrix
//#####################################################################
template<class TV> void GENERALIZED_FLUID_MASS<TV>::
Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const
{
    for(int i=1;i<=one_over_fluid_mass_at_faces.n;i++)
        data.Append(TRIPLE<int,int,T>(i,i,one_over_fluid_mass_at_faces(i)));
}
//#####################################################################
template class GENERALIZED_FLUID_MASS<VECTOR<float,1> >;
template class GENERALIZED_FLUID_MASS<VECTOR<float,2> >;
template class GENERALIZED_FLUID_MASS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GENERALIZED_FLUID_MASS<VECTOR<double,1> >;
template class GENERALIZED_FLUID_MASS<VECTOR<double,2> >;
template class GENERALIZED_FLUID_MASS<VECTOR<double,3> >;
#endif
