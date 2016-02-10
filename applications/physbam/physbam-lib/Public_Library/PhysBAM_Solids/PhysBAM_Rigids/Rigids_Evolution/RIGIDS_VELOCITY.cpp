//#####################################################################
// Copyright 2006-2009, Zhaosheng Bao, Elliot English, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Avi Robinson-Mosher, Craig Schroeder, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_VELOCITY
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigids_Evolution/RIGIDS_VELOCITY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGIDS_VELOCITY<TV>::
RIGIDS_VELOCITY(ARRAY_VIEW<TWIST<TV> > rigid_V_full,RIGID_BODY_COLLECTION<TV>& rigid_body_collection)
    :rigid_V(rigid_V_full,rigid_body_collection.dynamic_rigid_body_particles),
    kinematic_and_static_rigid_V(rigid_V_full,rigid_body_collection.static_and_kinematic_rigid_bodies)
{}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGIDS_VELOCITY<TV>::
RIGIDS_VELOCITY(ARRAY_VIEW<TWIST<TV> > rigid_V_full,ARRAY<int>& dynamic_rigid_body_particles,ARRAY<int>& static_and_kinematic_rigid_bodies)
    :rigid_V(rigid_V_full,dynamic_rigid_body_particles),kinematic_and_static_rigid_V(rigid_V_full,static_and_kinematic_rigid_bodies)
{}
//#####################################################################
// Operator +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& RIGIDS_VELOCITY<TV>::
operator+=(const BASE& bv)
{
    const RIGIDS_VELOCITY& v=debug_cast<const RIGIDS_VELOCITY&>(bv);
    rigid_V+=v.rigid_V;
    kinematic_and_static_rigid_V+=v.kinematic_and_static_rigid_V;
    return *this;
}
//#####################################################################
// Operator -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& RIGIDS_VELOCITY<TV>::
operator-=(const BASE& bv)
{
    const RIGIDS_VELOCITY& v=debug_cast<const RIGIDS_VELOCITY&>(bv);
    rigid_V-=v.rigid_V;
    kinematic_and_static_rigid_V-=v.kinematic_and_static_rigid_V;
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& RIGIDS_VELOCITY<TV>::
operator*=(const T a)
{
    rigid_V*=a;
    kinematic_and_static_rigid_V*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void RIGIDS_VELOCITY<TV>::
Copy(const T c,const BASE& bv)
{
    const RIGIDS_VELOCITY& v=debug_cast<const RIGIDS_VELOCITY&>(bv);
    rigid_V=c*v.rigid_V;kinematic_and_static_rigid_V=c*v.kinematic_and_static_rigid_V;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void RIGIDS_VELOCITY<TV>::
Copy(const T c1,const BASE& bv1,const BASE& bv2)
{
    const RIGIDS_VELOCITY& v1=debug_cast<const RIGIDS_VELOCITY&>(bv1),&v2=debug_cast<const RIGIDS_VELOCITY&>(bv2);
    rigid_V=c1*v1.rigid_V+v2.rigid_V;kinematic_and_static_rigid_V=c1*v1.kinematic_and_static_rigid_V+v2.kinematic_and_static_rigid_V;
}
//#####################################################################
// Function Pack
//#####################################################################
template<class TV> void RIGIDS_VELOCITY<TV>::
Pack(VECTOR_ND<T> &velocities) const
{
    PHYSBAM_ASSERT(velocities.Size()==Raw_Size());
    int index=0;
    for(int i=1;i<=rigid_V.Size();i++){
        for(int j=1;j<=TV::dimension;j++)
            velocities(++index)=rigid_V(i).linear(j);
        for(int j=1;j<=T_SPIN::dimension;j++)
            velocities(++index)=rigid_V(i).angular(j);}
}
//#####################################################################
// Function Unpack
//#####################################################################
template<class TV> void RIGIDS_VELOCITY<TV>::
Unpack(VECTOR_ND<T> &velocities)
{
    PHYSBAM_ASSERT(velocities.Size()==Raw_Size());
    int index=0;
    for(int i=1;i<=rigid_V.Size();i++){
        for(int j=1;j<=TV::dimension;j++)
            rigid_V(i).linear(j)=velocities(++index);
        for(int j=1;j<=T_SPIN::dimension;j++)
            rigid_V(i).angular(j)=velocities(++index);}
}
//#####################################################################
// Function Unpack_And_Add
//#####################################################################
template<class TV> void RIGIDS_VELOCITY<TV>::
Unpack_And_Add(VECTOR_ND<T> &velocities)
{
    PHYSBAM_ASSERT(velocities.Size()==Raw_Size());
    int index=0;
    for(int i=1;i<=rigid_V.Size();i++){
        for(int j=1;j<=TV::dimension;j++)
            rigid_V(i).linear(j)+=velocities(++index);
        for(int j=1;j<=T_SPIN::dimension;j++)
            rigid_V(i).angular(j)+=velocities(++index);}
}
template<class TV> int RIGIDS_VELOCITY<TV>::
Raw_Size() const
{
    return rigid_V.Size()*TWIST<TV>::dimension;
}
template<class TV> typename TV::SCALAR& RIGIDS_VELOCITY<TV>::
Raw_Get(int i)
{
    int o=(i-1)%TWIST<TV>::dimension+1,n=(i-1)/TWIST<TV>::dimension+1;
    if(o<=TV::dimension) return rigid_V(n).linear(o);
    return rigid_V(n).angular(o-TV::dimension);
}
//#####################################################################
template class RIGIDS_VELOCITY<VECTOR<float,1> >;
template class RIGIDS_VELOCITY<VECTOR<float,2> >;
template class RIGIDS_VELOCITY<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGIDS_VELOCITY<VECTOR<double,1> >;
template class RIGIDS_VELOCITY<VECTOR<double,2> >;
template class RIGIDS_VELOCITY<VECTOR<double,3> >;
#endif
