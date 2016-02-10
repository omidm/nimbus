//#####################################################################
// Copyright 2004-2005, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Dynamics/Motion/QUASI_RIGID_TRANSFORM_3D.h>
using namespace PhysBAM;
template<class T> QUASI_RIGID_TRANSFORM_3D<T>::
//#####################################################################
// Constructor
//#####################################################################
QUASI_RIGID_TRANSFORM_3D()
    :affine_transform(MATRIX<T,3>::Identity_Matrix())
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> QUASI_RIGID_TRANSFORM_3D<T>::
QUASI_RIGID_TRANSFORM_3D(const MATRIX<T,3>& affine_transform_input,const TV& translation_input)
    :affine_transform(affine_transform_input),translation(translation_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> QUASI_RIGID_TRANSFORM_3D<T>::
~QUASI_RIGID_TRANSFORM_3D()
{
}
//#####################################################################
// Function Rigidity_Penalty
//#####################################################################
template<class T> T QUASI_RIGID_TRANSFORM_3D<T>::
Rigidity_Penalty() const
{
    SYMMETRIC_MATRIX<T,3> strain=affine_transform.Normal_Equations_Matrix()-1;
    return strain.Frobenius_Norm_Squared();
}
//#####################################################################
// Function Rigidity_Penalty_Gradient
//#####################################################################
template<class T> T QUASI_RIGID_TRANSFORM_3D<T>::
Rigidity_Penalty_Gradient(const int i) const
{
    assert(1<=i && i<=12);
    if(i>9) return 0;
    MATRIX<T,3> gradient=affine_transform*(affine_transform.Normal_Equations_Matrix()-1);
    return (T)4*gradient.x[i-1];
}
//#####################################################################
// Function Ridigity_Penalty_Hessian_Definite_Part
//#####################################################################
template<class T> T QUASI_RIGID_TRANSFORM_3D<T>::
Ridigity_Penalty_Hessian_Definite_Part(const int i,const int j) const
{
    assert(1<=i && i<=12  &&  1<=j && j<=12);
    if(i>9 || j>9) return 0;
    MATRIX<T,3> U,V;DIAGONAL_MATRIX<T,3> Sigma;
    affine_transform.Fast_Singular_Value_Decomposition(U,Sigma,V);
    DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3> rotated_hessian;
    rotated_hessian.x1111=(T)12*sqr(Sigma.x11)-(T)4;
    rotated_hessian.x2222=(T)12*sqr(Sigma.x22)-(T)4;
    rotated_hessian.x3333=(T)12*sqr(Sigma.x33)-(T)4;
    rotated_hessian.x2211=0;
    rotated_hessian.x3311=0;
    rotated_hessian.x3322=0;
    rotated_hessian.x2121=(T)4*(sqr(Sigma.x11)+sqr(Sigma.x22)-(T)1);
    rotated_hessian.x3131=(T)4*(sqr(Sigma.x11)+sqr(Sigma.x33)-(T)1);
    rotated_hessian.x3232=(T)4*(sqr(Sigma.x22)+sqr(Sigma.x33)-(T)1);
    rotated_hessian.x2112=(T)4*Sigma.x11*Sigma.x22;
    rotated_hessian.x3113=(T)4*Sigma.x11*Sigma.x33;
    rotated_hessian.x3223=(T)4*Sigma.x22*Sigma.x33;
    rotated_hessian.Enforce_Definiteness();
    MATRIX<T,3> dF1,dF2;
    dF1.x[i-1]=(T)1;
    dF2.x[j-1]=(T)1;
    dF1=U.Transpose_Times(dF1)*V;
    dF2=U.Transpose_Times(dF2)*V;
    return dF1.Transpose_Times(rotated_hessian.Differential(dF2)).Trace();
}
//#####################################################################
// Function Incremental_Transform
//#####################################################################
template<class T> QUASI_RIGID_TRANSFORM_3D<T> QUASI_RIGID_TRANSFORM_3D<T>::
Incremental_Transform(const QUASI_RIGID_TRANSFORM_3D<T>& target_transform,const QUASI_RIGID_TRANSFORM_3D<T>& initial_transform)
{
    MATRIX<T,3> affine_transform_incremental=target_transform.affine_transform*initial_transform.affine_transform.Inverse();
    TV translation_incremental=target_transform.translation-affine_transform_incremental*initial_transform.translation;
    return QUASI_RIGID_TRANSFORM_3D<T>(affine_transform_incremental,translation_incremental);
}
//#####################################################################
// Function Composite_Transform
//#####################################################################
template<class T> QUASI_RIGID_TRANSFORM_3D<T> QUASI_RIGID_TRANSFORM_3D<T>::
Composite_Transform(const QUASI_RIGID_TRANSFORM_3D<T>& master_transform,const QUASI_RIGID_TRANSFORM_3D<T>& slave_transform)
{
    MATRIX<T,3> affine_transform_composite=master_transform.affine_transform*slave_transform.affine_transform;
    TV translation_composite=master_transform.affine_transform*slave_transform.translation+master_transform.translation;
    return QUASI_RIGID_TRANSFORM_3D<T>(affine_transform_composite,translation_composite);
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> FRAME<VECTOR<T,3> > QUASI_RIGID_TRANSFORM_3D<T>::
Frame()
{
    MATRIX<T,3> U,V;DIAGONAL_MATRIX<T,3> Sigma;
    affine_transform.Fast_Singular_Value_Decomposition(U,Sigma,V);
    return FRAME<TV>(translation,ROTATION<TV>(U.Times_Transpose(V)));
}
//#####################################################################
// Function Interpolate
//#####################################################################
template<class T> QUASI_RIGID_TRANSFORM_3D<T> QUASI_RIGID_TRANSFORM_3D<T>::
Interpolate(const QUASI_RIGID_TRANSFORM_3D<T>& initial,const QUASI_RIGID_TRANSFORM_3D<T>& final,const T interpolation_fraction)
{
    QUASI_RIGID_TRANSFORM_3D<T> result;
    result.translation=LINEAR_INTERPOLATION<T,TV>::Linear(initial.translation,final.translation,interpolation_fraction);
    ROTATION<TV> initial_rotation(initial.affine_transform),final_rotation(final.affine_transform);
    result.affine_transform=ROTATION<TV>::Spherical_Linear_Interpolation(initial_rotation,final_rotation,interpolation_fraction).Rotation_Matrix();
    return result;
}
//#####################################################################
// Function Distance
//#####################################################################
template<class T> T QUASI_RIGID_TRANSFORM_3D<T>::
Distance(const QUASI_RIGID_TRANSFORM_3D<T>& first,const QUASI_RIGID_TRANSFORM_3D<T>& second)
{
    ROTATION<TV> first_transform=ROTATION<TV>(first.affine_transform),second_transform=ROTATION<TV>(second.affine_transform);
    T larger_angle=max(first_transform.Angle(),second_transform.Angle());T larger_magnitude=max(first.translation.Magnitude(),second.translation.Magnitude());
    return (first_transform.Inverse()*second_transform).Angle()/(larger_angle?larger_angle:1)+(first.translation-second.translation).Magnitude()/(larger_magnitude?larger_magnitude:1);
}
//#####################################################################
// Function Make_Rigid
//#####################################################################
template<class T> void QUASI_RIGID_TRANSFORM_3D<T>::
Make_Rigid()
{
    MATRIX<T,3> U,V;
    DIAGONAL_MATRIX<T,3> Sigma;
    affine_transform.Fast_Singular_Value_Decomposition(U,Sigma,V);
    affine_transform=U.Times_Transpose(V);
}
//#####################################################################
// Function Print_Diagnostics
//#####################################################################
template<class T> void QUASI_RIGID_TRANSFORM_3D<T>::
Print_Diagnostics(std::ostream& output) const
{
    output<<"Transformation matrix :"<<std::endl<<affine_transform;
    MATRIX<T,3> U,V;DIAGONAL_MATRIX<T,3> Sigma;affine_transform.Fast_Singular_Value_Decomposition(U,Sigma,V);
    output<<"Singular values : "<<Sigma.x11<<" , "<<Sigma.x22<<" , "<<Sigma.x33<<std::endl;
    output<<"Translation : "<<translation<<std::endl;
    output<<"Rigidity penalty at current configuration : "<<Rigidity_Penalty()<<std::endl;
}
template class QUASI_RIGID_TRANSFORM_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class QUASI_RIGID_TRANSFORM_3D<double>;
#endif
