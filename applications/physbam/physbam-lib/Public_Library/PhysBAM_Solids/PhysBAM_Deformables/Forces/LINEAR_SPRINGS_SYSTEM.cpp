//#####################################################################
// Copyright 2010, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_SPRING_SYSTEM
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS_SYSTEM.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LINEAR_SPRINGS_SYSTEM<TV>::
LINEAR_SPRINGS_SYSTEM(SPARSE_MATRIX_FLAT_MXN<T>& J_input)
    :KRYLOV_SYSTEM_BASE<typename TV::SCALAR>(false,true),J(J_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LINEAR_SPRINGS_SYSTEM<TV>::
~LINEAR_SPRINGS_SYSTEM()
{
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void LINEAR_SPRINGS_SYSTEM<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const
{
}
//#####################################################################
// Function Force
//#####################################################################
template<class TV> void LINEAR_SPRINGS_SYSTEM<TV>::
Force(const VECTOR_T& V,VECTOR_T& F) const
{
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void LINEAR_SPRINGS_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const
{
    const VECTOR_T& V=debug_cast<const VECTOR_T&>(BV);VECTOR_T& F=debug_cast<VECTOR_T&>(BF);
    VECTOR_ND<T> Jx(F.X.n);
    J.Times(V.X,Jx);J.Transpose_Times(Jx,F.X);
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void LINEAR_SPRINGS_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& BV) const
{
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double LINEAR_SPRINGS_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV1,const KRYLOV_VECTOR_BASE<T>& BV2) const
{
    const VECTOR_T& V1=debug_cast<const VECTOR_T&>(BV1),&V2=debug_cast<const VECTOR_T&>(BV2);
    double inner_product=VECTOR_ND<T>::Dot_Product_Double_Precision(V1.X,V2.X);
    return inner_product;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const
{
    const VECTOR_T& R=debug_cast<const VECTOR_T&>(BR);
    T convergence_norm_squared=sqr(R.X.Maximum_Magnitude());
    T convergence_norm=sqrt(convergence_norm_squared);
    return convergence_norm;
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void LINEAR_SPRINGS_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const
{
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void LINEAR_SPRINGS_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LINEAR_SPRINGS_SYSTEM_VECTOR<TV>::
LINEAR_SPRINGS_SYSTEM_VECTOR(const int size)
{
    X.Resize(size);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LINEAR_SPRINGS_SYSTEM_VECTOR<TV>::
~LINEAR_SPRINGS_SYSTEM_VECTOR()
{
}
//#####################################################################
// Operator =
//#####################################################################
template<class TV> LINEAR_SPRINGS_SYSTEM_VECTOR<TV>& LINEAR_SPRINGS_SYSTEM_VECTOR<TV>::
operator=(const LINEAR_SPRINGS_SYSTEM_VECTOR& v)
{
    X=v.X;
    return *this;
}
//#####################################################################
// Operator +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& LINEAR_SPRINGS_SYSTEM_VECTOR<TV>::
operator+=(const BASE& bv)
{
    const LINEAR_SPRINGS_SYSTEM_VECTOR& v=debug_cast<const LINEAR_SPRINGS_SYSTEM_VECTOR&>(bv);
    X+=v.X;
    return *this;
}
//#####################################################################
// Operator -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& LINEAR_SPRINGS_SYSTEM_VECTOR<TV>::
operator-=(const BASE& bv)
{
    const LINEAR_SPRINGS_SYSTEM_VECTOR& v=debug_cast<const LINEAR_SPRINGS_SYSTEM_VECTOR&>(bv);
    X-=v.X;
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& LINEAR_SPRINGS_SYSTEM_VECTOR<TV>::
operator*=(const T a)
{
    X*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void LINEAR_SPRINGS_SYSTEM_VECTOR<TV>::
Copy(const T c,const BASE& bv)
{
    const LINEAR_SPRINGS_SYSTEM_VECTOR& v=debug_cast<const LINEAR_SPRINGS_SYSTEM_VECTOR&>(bv);
    assert(v.X.n==X.n);
    VECTOR_ND<T>::Copy(c,v.X,X);
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void LINEAR_SPRINGS_SYSTEM_VECTOR<TV>::
Copy(const T c1,const BASE& bv1,const BASE& bv2)
{
    const LINEAR_SPRINGS_SYSTEM_VECTOR& v1=debug_cast<const LINEAR_SPRINGS_SYSTEM_VECTOR&>(bv1);
    const LINEAR_SPRINGS_SYSTEM_VECTOR& v2=debug_cast<const LINEAR_SPRINGS_SYSTEM_VECTOR&>(bv2);
    assert(v1.X.n==v2.X.n && X.n==v1.X.n);
    VECTOR_ND<T>::Copy(c1,v1.X,v2.X,X);
}
//#####################################################################
// Function Print
//#####################################################################
template<class TV> void LINEAR_SPRINGS_SYSTEM_VECTOR<TV>::
Print() const
{
    // Flat print
    for(int i=1;i<=X.n;i++)
    {std::stringstream ss;ss<<X(i)<<" ";LOG::filecout(ss.str());}
}
template<class TV> int LINEAR_SPRINGS_SYSTEM_VECTOR<TV>::
Raw_Size() const
{
    return X.n;
}
template<class TV> typename TV::SCALAR& LINEAR_SPRINGS_SYSTEM_VECTOR<TV>::
Raw_Get(int i)
{
    return X(i);
}
//#####################################################################
template class LINEAR_SPRINGS_SYSTEM<VECTOR<float,1> >;
template class LINEAR_SPRINGS_SYSTEM<VECTOR<float,2> >;
template class LINEAR_SPRINGS_SYSTEM<VECTOR<float,3> >;
template class LINEAR_SPRINGS_SYSTEM_VECTOR<VECTOR<float,1> >;
template class LINEAR_SPRINGS_SYSTEM_VECTOR<VECTOR<float,2> >;
template class LINEAR_SPRINGS_SYSTEM_VECTOR<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LINEAR_SPRINGS_SYSTEM<VECTOR<double,1> >;
template class LINEAR_SPRINGS_SYSTEM<VECTOR<double,2> >;
template class LINEAR_SPRINGS_SYSTEM<VECTOR<double,3> >;
template class LINEAR_SPRINGS_SYSTEM_VECTOR<VECTOR<double,1> >;
template class LINEAR_SPRINGS_SYSTEM_VECTOR<VECTOR<double,2> >;
template class LINEAR_SPRINGS_SYSTEM_VECTOR<VECTOR<double,3> >;
#endif
