//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/LEVELSET_VISCOSITY_UNIFORM_SYSTEM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LEVELSET_VISCOSITY_UNIFORM_SYSTEM<TV>::
LEVELSET_VISCOSITY_UNIFORM_SYSTEM(LEVELSET_INDEX_MAP_UNIFORM<TV>& index_map_input)
    :KRYLOV_SYSTEM_BASE<T>(false,true),poisson(index_map_input),scale(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LEVELSET_VISCOSITY_UNIFORM_SYSTEM<TV>::
~LEVELSET_VISCOSITY_UNIFORM_SYSTEM()
{
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void LEVELSET_VISCOSITY_UNIFORM_SYSTEM<TV>::
Compute(int axis,T scale_input)
{
    scale=scale_input;
    poisson.Compute(axis);
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void LEVELSET_VISCOSITY_UNIFORM_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const
{
    const VECTOR_T& vx=debug_cast<const VECTOR_T&>(x);
    VECTOR_T& vr=debug_cast<VECTOR_T&>(result);
    poisson.P.Times(vx.v,vr.v);
    vr.v=vx.v-scale*vr.v;
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double LEVELSET_VISCOSITY_UNIFORM_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    const VECTOR_T& vx=debug_cast<const VECTOR_T&>(x);
    const VECTOR_T& vy=debug_cast<const VECTOR_T&>(y);
    return vx.v.Dot_Product(vx.v,vy.v);
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR LEVELSET_VISCOSITY_UNIFORM_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const
{
    const VECTOR_T& vx=debug_cast<const VECTOR_T&>(x);
    return vx.v.Max_Abs();
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void LEVELSET_VISCOSITY_UNIFORM_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void LEVELSET_VISCOSITY_UNIFORM_SYSTEM<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void LEVELSET_VISCOSITY_UNIFORM_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Add_Constant_Part
//#####################################################################
template<class TV> void LEVELSET_VISCOSITY_UNIFORM_SYSTEM<TV>::
Add_Constant_Part(KRYLOV_VECTOR_BASE<T>& x) const
{
    VECTOR_T& vx=debug_cast<VECTOR_T&>(x);
    vx.v+=poisson.b*scale;
}
template class LEVELSET_VISCOSITY_UNIFORM_SYSTEM<VECTOR<float,1> >;
template class LEVELSET_VISCOSITY_UNIFORM_SYSTEM<VECTOR<float,2> >;
template class LEVELSET_VISCOSITY_UNIFORM_SYSTEM<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET_VISCOSITY_UNIFORM_SYSTEM<VECTOR<double,1> >;
template class LEVELSET_VISCOSITY_UNIFORM_SYSTEM<VECTOR<double,2> >;
template class LEVELSET_VISCOSITY_UNIFORM_SYSTEM<VECTOR<double,3> >;
#endif
