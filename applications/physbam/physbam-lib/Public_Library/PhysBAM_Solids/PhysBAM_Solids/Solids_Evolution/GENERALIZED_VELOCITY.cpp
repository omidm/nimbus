//#####################################################################
// Copyright 2006-2008, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Craig Schroeder, Andrew Selle, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_MXN
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GENERALIZED_VELOCITY<TV>::
GENERALIZED_VELOCITY(ARRAY_VIEW<TV> V_full,ARRAY_VIEW<TWIST<TV> > rigid_V_full,const SOLID_BODY_COLLECTION<TV>& solid_body_collection)
    :V(V_full,solid_body_collection.deformable_body_collection.dynamic_particles),rigid_V(rigid_V_full,solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles),
    kinematic_and_static_rigid_V(rigid_V_full,solid_body_collection.rigid_body_collection.static_and_kinematic_rigid_bodies)
{}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GENERALIZED_VELOCITY<TV>::
GENERALIZED_VELOCITY(ARRAY_VIEW<TV> V_full,ARRAY<int>& dynamic_particles,ARRAY_VIEW<TWIST<TV> > rigid_V_full,ARRAY<int>& dynamic_rigid_body_particles,
    ARRAY<int>& static_and_kinematic_rigid_bodies)
    :V(V_full,dynamic_particles),rigid_V(rigid_V_full,dynamic_rigid_body_particles),
    kinematic_and_static_rigid_V(rigid_V_full,static_and_kinematic_rigid_bodies)
{}
template<class TV> GENERALIZED_VELOCITY<TV>::
~GENERALIZED_VELOCITY()
{
}
//#####################################################################
// Operator +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& GENERALIZED_VELOCITY<TV>::
operator+=(const BASE& bv)
{
    const GENERALIZED_VELOCITY& v=debug_cast<const GENERALIZED_VELOCITY&>(bv);
    V+=v.V;
    rigid_V+=v.rigid_V;
    kinematic_and_static_rigid_V+=v.kinematic_and_static_rigid_V;
    return *this;
}
//#####################################################################
// Operator -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& GENERALIZED_VELOCITY<TV>::
operator-=(const BASE& bv)
{
    const GENERALIZED_VELOCITY& v=debug_cast<const GENERALIZED_VELOCITY&>(bv);
    V-=v.V;
    rigid_V-=v.rigid_V;
    kinematic_and_static_rigid_V-=v.kinematic_and_static_rigid_V;
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& GENERALIZED_VELOCITY<TV>::
operator*=(const T a)
{
    V*=a;
    rigid_V*=a;
    kinematic_and_static_rigid_V*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void GENERALIZED_VELOCITY<TV>::
Copy(const T c,const BASE& bv)
{
    const GENERALIZED_VELOCITY& v=debug_cast<const GENERALIZED_VELOCITY&>(bv);
    V=c*v.V;rigid_V=c*v.rigid_V;kinematic_and_static_rigid_V=c*v.kinematic_and_static_rigid_V;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void GENERALIZED_VELOCITY<TV>::
Copy(const T c1,const BASE& bv1,const BASE& bv2)
{
    const GENERALIZED_VELOCITY& v1=debug_cast<const GENERALIZED_VELOCITY&>(bv1),&v2=debug_cast<const GENERALIZED_VELOCITY&>(bv2);
    V=c1*v1.V+v2.V;rigid_V=c1*v1.rigid_V+v2.rigid_V;kinematic_and_static_rigid_V=c1*v1.kinematic_and_static_rigid_V+v2.kinematic_and_static_rigid_V;
}
//#####################################################################
// Function Pack
//#####################################################################
template<class TV> void GENERALIZED_VELOCITY<TV>::
Pack(VECTOR_ND<T> &velocities) const
{
    PHYSBAM_ASSERT(velocities.Size()==Raw_Size());
    int index=0;
    for(int i=1;i<=V.Size();i++)
        for(int j=1;j<=TV::dimension;j++)
            velocities(++index)=V.array(i)(j);
    for(int i=1;i<=rigid_V.Size();i++){
        for(int j=1;j<=TV::dimension;j++)
            velocities(++index)=rigid_V(i).linear(j);
        for(int j=1;j<=T_SPIN::dimension;j++)
            velocities(++index)=rigid_V(i).angular(j);}
}
//#####################################################################
// Function Unpack
//#####################################################################
template<class TV> void GENERALIZED_VELOCITY<TV>::
Unpack(VECTOR_ND<T> &velocities)
{
    PHYSBAM_ASSERT(velocities.Size()==Raw_Size());
    int index=0;
    for(int i=1;i<=V.Size();i++)
        for(int j=1;j<=TV::dimension;j++)
            V.array(i)(j)=velocities(++index);
    for(int i=1;i<=rigid_V.Size();i++){
        for(int j=1;j<=TV::dimension;j++)
            rigid_V(i).linear(j)=velocities(++index);
        for(int j=1;j<=T_SPIN::dimension;j++)
            rigid_V(i).angular(j)=velocities(++index);}
}
//#####################################################################
// Function Unpack_And_Add
//#####################################################################
template<class TV> void GENERALIZED_VELOCITY<TV>::
Unpack_And_Add(VECTOR_ND<T> &velocities)
{
    PHYSBAM_ASSERT(velocities.Size()==Raw_Size());
    int index=0;
    for(int i=1;i<=V.Size();i++)
        for(int j=1;j<=TV::dimension;j++)
            V.array(i)(j)+=velocities(++index);
    for(int i=1;i<=rigid_V.Size();i++){
        for(int j=1;j<=TV::dimension;j++)
            rigid_V(i).linear(j)+=velocities(++index);
        for(int j=1;j<=T_SPIN::dimension;j++)
            rigid_V(i).angular(j)+=velocities(++index);}
}
//#####################################################################
// Function Raw_Size
//#####################################################################
template<class TV> int GENERALIZED_VELOCITY<TV>::
Raw_Size() const
{
    return V.Size()*TV::dimension+rigid_V.Size()*TWIST<TV>::dimension;
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class TV> typename TV::SCALAR& GENERALIZED_VELOCITY<TV>::
Raw_Get(int i)
{
    if(i<=V.Size()*TV::dimension) return V((i-1)/TV::dimension+1)((i-1)%TV::dimension+1);
    i-=V.Size()*TV::dimension;
    int o=(i-1)%TWIST<TV>::dimension+1,n=(i-1)/TWIST<TV>::dimension+1;
    if(o<=TV::dimension) return rigid_V(n).linear(o);
    return rigid_V(n).angular(o-TV::dimension);
}
//#####################################################################
// Function Exchange
//#####################################################################
template<class TV> void GENERALIZED_VELOCITY<TV>::
Exchange(GENERALIZED_VELOCITY<TV>& gv)
{
    V.array.Exchange(gv.V.array);
    rigid_V.array.Exchange(gv.rigid_V.array);
    kinematic_and_static_rigid_V.array.Exchange(gv.kinematic_and_static_rigid_V.array);
}
//#####################################################################
template class GENERALIZED_VELOCITY<VECTOR<float,1> >;
template class GENERALIZED_VELOCITY<VECTOR<float,2> >;
template class GENERALIZED_VELOCITY<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GENERALIZED_VELOCITY<VECTOR<double,1> >;
template class GENERALIZED_VELOCITY<VECTOR<double,2> >;
template class GENERALIZED_VELOCITY<VECTOR<double,3> >;
#endif
