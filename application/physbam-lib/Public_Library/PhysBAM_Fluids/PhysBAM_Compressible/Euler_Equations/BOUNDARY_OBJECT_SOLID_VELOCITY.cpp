//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Fill solid cells with the solid velocity, rather than some reflected value.
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/BOUNDARY_OBJECT_SOLID_VELOCITY.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER.h>

using namespace PhysBAM;
//#####################################################################
// Function Apply_Neumann_Boundary_Condition
//#####################################################################
template<class T_GRID> void BOUNDARY_OBJECT_SOLID_VELOCITY<T_GRID>::
Apply_Neumann_Boundary_Condition(TV_DIMENSION& u_1d,const T neumann_face_velocity,const int axis)
{
    T rho=u_1d(1);
    if(rho){
        TV velocity=EULER<T_GRID>::Get_Velocity(u_1d);velocity(axis)=neumann_face_velocity;
        T internal_energy=EULER<T_GRID>::e(u_1d);
        u_1d=EULER<T_GRID>::Get_Euler_State_From_rho_velocity_And_internal_energy(rho,velocity,internal_energy);}
}
template<class T_GRID> void BOUNDARY_OBJECT_SOLID_VELOCITY<T_GRID>::
Apply_Neumann_Boundary_Condition(T_ARRAYS_ELEMENT& u_1d,const TV& object_velocity,const T unused)
{
    T rho=u_1d(1);
    if(rho){
        T internal_energy=EULER<T_GRID>::e(u_1d);
        u_1d=EULER<T_GRID>::Get_Euler_State_From_rho_velocity_And_internal_energy(rho,object_velocity,internal_energy);}
}
//#####################################################################
template class BOUNDARY_OBJECT_SOLID_VELOCITY<GRID<VECTOR<float,1> > >;
template class BOUNDARY_OBJECT_SOLID_VELOCITY<GRID<VECTOR<float,2> > >;
template class BOUNDARY_OBJECT_SOLID_VELOCITY<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BOUNDARY_OBJECT_SOLID_VELOCITY<GRID<VECTOR<double,1> > >;
template class BOUNDARY_OBJECT_SOLID_VELOCITY<GRID<VECTOR<double,2> > >;
template class BOUNDARY_OBJECT_SOLID_VELOCITY<GRID<VECTOR<double,3> > >;
#endif
