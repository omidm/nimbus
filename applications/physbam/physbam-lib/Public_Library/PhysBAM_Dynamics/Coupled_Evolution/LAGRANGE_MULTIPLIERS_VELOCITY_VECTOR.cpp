//#####################################################################
// Copyright 2006-2008, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Craig Schroeder, Andrew Selle, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR<TV>::
LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR(VECTOR_ND<T>& lagrange_multipliers_input,GENERALIZED_VELOCITY<TV>& solid_velocity_input)
    :lagrange_multipliers(lagrange_multipliers_input),solid_velocity(solid_velocity_input)
{
}
//#####################################################################
// Operator =
//#####################################################################
template<class TV> LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR<TV>& LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR<TV>::
operator=(const LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR& v)
{
    lagrange_multipliers=v.lagrange_multipliers;
    solid_velocity=v.solid_velocity;
    return *this;
}
//#####################################################################
// Operator +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR<TV>::
operator+=(const BASE& bv)
{
    const LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR& v=debug_cast<const LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR&>(bv);
    lagrange_multipliers+=v.lagrange_multipliers;
    solid_velocity+=v.solid_velocity;
    return *this;
}
//#####################################################################
// Operator -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR<TV>::
operator-=(const BASE& bv)
{
    const LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR& v=debug_cast<const LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR&>(bv);
    lagrange_multipliers-=v.lagrange_multipliers;
    solid_velocity-=v.solid_velocity;
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR<TV>::
operator*=(const T a)
{
    lagrange_multipliers*=a;
    solid_velocity*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR<TV>::
Copy(const T c,const BASE& bv)
{
    const LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR& v=debug_cast<const LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR&>(bv);
    assert(v.lagrange_multipliers.n==lagrange_multipliers.n);
    VECTOR_ND<T>::Copy(c,v.lagrange_multipliers,lagrange_multipliers);
    solid_velocity.Copy(c,v.solid_velocity);
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR<TV>::
Copy(const T c1,const BASE& bv1,const BASE& bv2)
{
    const LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR& v1=debug_cast<const LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR&>(bv1);
    const LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR& v2=debug_cast<const LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR&>(bv2);
    assert(v1.lagrange_multipliers.n==v2.lagrange_multipliers.n && lagrange_multipliers.n==v1.lagrange_multipliers.n);
    VECTOR_ND<T>::Copy(c1,v1.lagrange_multipliers,v2.lagrange_multipliers,lagrange_multipliers);
    solid_velocity.Copy(c1,v1.solid_velocity,v2.solid_velocity);
}
//#####################################################################
template class LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR<VECTOR<float,1> >;
template class LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR<VECTOR<float,2> >;
template class LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR<VECTOR<double,1> >;
template class LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR<VECTOR<double,2> >;
template class LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR<VECTOR<double,3> >;
#endif
