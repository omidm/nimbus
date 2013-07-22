//#####################################################################
// Copyright 2004-2005, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUASI_RIGID_TRANSFORM_3D
//#####################################################################
#ifndef __QUASI_RIGID_TRANSFORM_3D__
#define __QUASI_RIGID_TRANSFORM_3D__

#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE.h>
namespace PhysBAM{

template<class T>
class QUASI_RIGID_TRANSFORM_3D
{
    typedef VECTOR<T,3> TV;
public:
    MATRIX<T,3> affine_transform;
    TV translation;

    QUASI_RIGID_TRANSFORM_3D();
    QUASI_RIGID_TRANSFORM_3D(const MATRIX<T,3>& affine_transform_input,const TV& translation_input);
    ~QUASI_RIGID_TRANSFORM_3D();

    T operator()(const int i) const
    {assert(1<=i && i<=12);return i<=9?affine_transform.x[i-1]:translation[i-9];}

    T& operator()(const int i)
    {assert(1<=i && i<=12);return i<=9?affine_transform.x[i-1]:translation[i-9];}

    T Identity(const int i) const
    {assert(1<=i && i<=12);return (T)((i==1 || i==5 || i==9)?1:0);}

//#####################################################################
    T Rigidity_Penalty() const;
    T Rigidity_Penalty_Gradient(const int i) const;
    T Ridigity_Penalty_Hessian_Definite_Part(const int i,const int j) const;
    static QUASI_RIGID_TRANSFORM_3D<T> Incremental_Transform(const QUASI_RIGID_TRANSFORM_3D<T>& target_transform,const QUASI_RIGID_TRANSFORM_3D<T>& initial_transform);
    static QUASI_RIGID_TRANSFORM_3D<T> Composite_Transform(const QUASI_RIGID_TRANSFORM_3D<T>& master_transform,const QUASI_RIGID_TRANSFORM_3D<T>& slave_transform);
    FRAME<TV> Frame();
    static QUASI_RIGID_TRANSFORM_3D<T> Interpolate(const QUASI_RIGID_TRANSFORM_3D<T>& initial,const QUASI_RIGID_TRANSFORM_3D<T>& final,const T interpolation_fraction);
    static T Distance(const QUASI_RIGID_TRANSFORM_3D<T>& first,const QUASI_RIGID_TRANSFORM_3D<T>& second);
    void Make_Rigid();
    void Print_Diagnostics(std::ostream& output) const;
//#####################################################################
};
}
#endif
