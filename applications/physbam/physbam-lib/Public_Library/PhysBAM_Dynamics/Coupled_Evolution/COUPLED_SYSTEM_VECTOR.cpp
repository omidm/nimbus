//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COUPLED_SYSTEM_VECTOR
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COUPLED_SYSTEM_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COUPLED_SYSTEM_VECTOR<TV>::
COUPLED_SYSTEM_VECTOR()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COUPLED_SYSTEM_VECTOR<TV>::
~COUPLED_SYSTEM_VECTOR()
{
}
//#####################################################################
// Operator =
//#####################################################################
template<class TV> COUPLED_SYSTEM_VECTOR<TV>& COUPLED_SYSTEM_VECTOR<TV>::
operator=(const COUPLED_SYSTEM_VECTOR& v)
{
    pressure=v.pressure;
    lambda=v.lambda;
    force_coefficients=v.force_coefficients;
    viscous_force_coefficients=v.viscous_force_coefficients;
    return *this;
}
//#####################################################################
// Operator +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& COUPLED_SYSTEM_VECTOR<TV>::
operator+=(const BASE& bv)
{
    const COUPLED_SYSTEM_VECTOR& v=debug_cast<const COUPLED_SYSTEM_VECTOR&>(bv);
    pressure+=v.pressure;
    lambda+=v.lambda;
    force_coefficients+=v.force_coefficients;
    viscous_force_coefficients+=v.viscous_force_coefficients;
    return *this;
}
//#####################################################################
// Operator -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& COUPLED_SYSTEM_VECTOR<TV>::
operator-=(const BASE& bv)
{
    const COUPLED_SYSTEM_VECTOR& v=debug_cast<const COUPLED_SYSTEM_VECTOR&>(bv);
    pressure-=v.pressure;
    lambda-=v.lambda;
    force_coefficients-=v.force_coefficients;
    viscous_force_coefficients-=v.viscous_force_coefficients;
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& COUPLED_SYSTEM_VECTOR<TV>::
operator*=(const T a)
{
    pressure*=a;
    lambda*=a;
    force_coefficients*=a;
    viscous_force_coefficients*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void COUPLED_SYSTEM_VECTOR<TV>::
Copy(const T c,const BASE& bv)
{
    const COUPLED_SYSTEM_VECTOR& v=debug_cast<const COUPLED_SYSTEM_VECTOR&>(bv);
    assert(v.pressure.n==pressure.n);
    VECTOR_ND<T>::Copy(c,v.pressure,pressure);
    lambda=c*v.lambda;
    force_coefficients=c*v.force_coefficients;
    viscous_force_coefficients=c*v.viscous_force_coefficients;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void COUPLED_SYSTEM_VECTOR<TV>::
Copy(const T c1,const BASE& bv1,const BASE& bv2)
{
    const COUPLED_SYSTEM_VECTOR& v1=debug_cast<const COUPLED_SYSTEM_VECTOR&>(bv1);
    const COUPLED_SYSTEM_VECTOR& v2=debug_cast<const COUPLED_SYSTEM_VECTOR&>(bv2);
    assert(v1.pressure.n==v2.pressure.n && pressure.n==v1.pressure.n);
    VECTOR_ND<T>::Copy(c1,v1.pressure,v2.pressure,pressure);
    lambda=c1*v1.lambda+v2.lambda;
    force_coefficients=c1*v1.force_coefficients+v2.force_coefficients;
    viscous_force_coefficients=c1*v1.viscous_force_coefficients+v2.viscous_force_coefficients;
}
//#####################################################################
// Function Print
//#####################################################################
template<class TV> void COUPLED_SYSTEM_VECTOR<TV>::
Print() const
{
    std::stringstream ss;
    // Flat print
    for(int i=1;i<=pressure.n;i++)
        ss<<pressure(i)<<" ";
    for(COUPLING_CONSTRAINT_ID i(1);i<=lambda.Size();i++)
        ss<<lambda(i)<<" ";
    for(FORCE_AGGREGATE_ID i(1);i<=force_coefficients.Size();i++)
        ss<<force_coefficients(i)<<" ";
    for(VISCOUS_FORCE_ID i(1);i<=viscous_force_coefficients.Size();i++)
        ss<<viscous_force_coefficients(i)<<" ";
    LOG::filecout(ss.str());
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int COUPLED_SYSTEM_VECTOR<TV>::
Raw_Size() const
{
    return pressure.n+Value(force_coefficients.m)+Value(lambda.m)+Value(viscous_force_coefficients.m);
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& COUPLED_SYSTEM_VECTOR<TV>::
Raw_Get(int i)
{
    if(i<=pressure.n) return pressure(i);
    i-=pressure.n;
    int l=Value(lambda.m);
    if(i<=l) return lambda(COUPLING_CONSTRAINT_ID(i));
    int f=Value(force_coefficients.m);
    if(i<=l+f) return force_coefficients(FORCE_AGGREGATE_ID(i-l));
    return viscous_force_coefficients(VISCOUS_FORCE_ID(i-l-f));
}
//#####################################################################
template class COUPLED_SYSTEM_VECTOR<VECTOR<float,1> >;
template class COUPLED_SYSTEM_VECTOR<VECTOR<float,2> >;
template class COUPLED_SYSTEM_VECTOR<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COUPLED_SYSTEM_VECTOR<VECTOR<double,1> >;
template class COUPLED_SYSTEM_VECTOR<VECTOR<double,2> >;
template class COUPLED_SYSTEM_VECTOR<VECTOR<double,3> >;
#endif
