//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_VECTOR.h>
using namespace PhysBAM;
template<class TV> ARTICULATED_VECTOR<TV>::
ARTICULATED_VECTOR()
{
}
template<class TV> ARTICULATED_VECTOR<TV>::
~ARTICULATED_VECTOR()
{
}
template<class TV> const ARTICULATED_VECTOR<TV>& ARTICULATED_VECTOR<TV>::
operator= (const ARTICULATED_VECTOR& bv)
{
    v=bv.v;
    return *this;
}
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& ARTICULATED_VECTOR<TV>::
operator+=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    v+=debug_cast<const ARTICULATED_VECTOR&>(bv).v;
    return *this;
}
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& ARTICULATED_VECTOR<TV>::
operator-=(const KRYLOV_VECTOR_BASE<T>& bv)
{
    v-=debug_cast<const ARTICULATED_VECTOR&>(bv).v;
    return *this;
}
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& ARTICULATED_VECTOR<TV>::
operator*=(const T a)
{
    v*=a;
    return *this;
}
template<class TV> void ARTICULATED_VECTOR<TV>::
Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv)
{
    v=c*debug_cast<const ARTICULATED_VECTOR&>(bv).v;
}
template<class TV> void ARTICULATED_VECTOR<TV>::
Copy(const T c1,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2)
{
    v=c1*debug_cast<const ARTICULATED_VECTOR&>(bv1).v+debug_cast<const ARTICULATED_VECTOR&>(bv2).v;
}
template<class TV> int ARTICULATED_VECTOR<TV>::
Raw_Size() const
{
    return Value(v.m)*TWIST<TV>::dimension;
}
template<class TV> typename TV::SCALAR& ARTICULATED_VECTOR<TV>::
Raw_Get(int i)
{
    int o=(i-1)%TWIST<TV>::dimension+1,n=(i-1)/TWIST<TV>::dimension+1;
    if(o<=TV::dimension) return v(JOINT_ID(n)).linear(o);
    return v(JOINT_ID(n)).angular(o-TV::dimension);
}
template class ARTICULATED_VECTOR<VECTOR<float,1> >;
template class ARTICULATED_VECTOR<VECTOR<float,2> >;
template class ARTICULATED_VECTOR<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ARTICULATED_VECTOR<VECTOR<double,1> >;
template class ARTICULATED_VECTOR<VECTOR<double,2> >;
template class ARTICULATED_VECTOR<VECTOR<double,3> >;
#endif
