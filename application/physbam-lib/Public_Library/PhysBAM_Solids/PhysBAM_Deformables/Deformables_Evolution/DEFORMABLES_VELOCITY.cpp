//#####################################################################
// Copyright 2006-2009, Zhaosheng Bao, Elliot English, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Avi Robinson-Mosher, Craig Schroeder, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLES_VELOCITY
//#####################################################################
#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformables_Evolution/DEFORMABLES_VELOCITY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLES_VELOCITY<TV>::
DEFORMABLES_VELOCITY(ARRAY_VIEW<TV> V_full,ARRAY_VIEW<TWIST<TV> > rigid_V_full,DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection,RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection)
    :V(V_full,deformable_body_collection.dynamic_particles),rigid_V(rigid_V_full,*new ARRAY<int>(IDENTITY_ARRAY<int>(rigid_geometry_collection.particles.array_collection->Size())))
{}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLES_VELOCITY<TV>::
DEFORMABLES_VELOCITY(ARRAY_VIEW<TV> V_full,DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection)
    :V(V_full,deformable_body_collection.dynamic_particles),rigid_V(ARRAY_VIEW<TWIST<TV> >(0,0),*new ARRAY<int>())
{}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> DEFORMABLES_VELOCITY<TV>::
DEFORMABLES_VELOCITY(ARRAY_VIEW<TV> V_full,ARRAY<int>& dynamic_particles,ARRAY_VIEW<TWIST<TV> > rigid_V_full)
    :V(V_full,dynamic_particles),rigid_V(rigid_V_full,*new ARRAY<int>(IDENTITY_ARRAY<int>(rigid_V_full.Size())))
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> DEFORMABLES_VELOCITY<TV>::
~DEFORMABLES_VELOCITY()
{delete &rigid_V.indices;}
//#####################################################################
// Operator +=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& DEFORMABLES_VELOCITY<TV>::
operator+=(const BASE& bv)
{
    const DEFORMABLES_VELOCITY& v=debug_cast<const DEFORMABLES_VELOCITY&>(bv);
    V+=v.V;
    rigid_V+=v.rigid_V;
    return *this;
}
//#####################################################################
// Operator -=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& DEFORMABLES_VELOCITY<TV>::
operator-=(const BASE& bv)
{
    const DEFORMABLES_VELOCITY& v=debug_cast<const DEFORMABLES_VELOCITY&>(bv);
    V-=v.V;
    rigid_V-=v.rigid_V;
    return *this;
}
//#####################################################################
// Operator *=
//#####################################################################
template<class TV> KRYLOV_VECTOR_BASE<typename TV::SCALAR>& DEFORMABLES_VELOCITY<TV>::
operator*=(const T a)
{
    V*=a;
    rigid_V*=a;
    return *this;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void DEFORMABLES_VELOCITY<TV>::
Copy(const T c,const BASE& bv)
{
    const DEFORMABLES_VELOCITY& v=debug_cast<const DEFORMABLES_VELOCITY&>(bv);
    V=c*v.V;rigid_V=c*v.rigid_V;
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void DEFORMABLES_VELOCITY<TV>::
Copy(const T c1,const BASE& bv1,const BASE& bv2)
{
    const DEFORMABLES_VELOCITY& v1=debug_cast<const DEFORMABLES_VELOCITY&>(bv1),&v2=debug_cast<const DEFORMABLES_VELOCITY&>(bv2);
    V=c1*v1.V+v2.V;rigid_V=c1*v1.rigid_V+v2.rigid_V;
}
//#####################################################################
// Function Pack
//#####################################################################
template<class TV> void DEFORMABLES_VELOCITY<TV>::
Pack(VECTOR_ND<T> &velocities) const
{
    PHYSBAM_ASSERT(velocities.Size()==Raw_Size());
    int index=0;
    for(int i=1;i<=V.Size();i++)
        for(int j=1;j<=TV::dimension;j++)
            velocities(++index)=V.array(i)(j);
}
//#####################################################################
// Function Unpack
//#####################################################################
template<class TV> void DEFORMABLES_VELOCITY<TV>::
Unpack(VECTOR_ND<T> &velocities)
{
    PHYSBAM_ASSERT(velocities.Size()==Raw_Size());
    int index=0;
    for(int i=1;i<=V.Size();i++)
        for(int j=1;j<=TV::dimension;j++)
            V.array(i)(j)=velocities(++index);
}
//#####################################################################
// Function Unpack_And_Add
//#####################################################################
template<class TV> void DEFORMABLES_VELOCITY<TV>::
Unpack_And_Add(VECTOR_ND<T> &velocities)
{
    PHYSBAM_ASSERT(velocities.Size()==Raw_Size());
    int index=0;
    for(int i=1;i<=V.Size();i++)
        for(int j=1;j<=TV::dimension;j++)
            V.array(i)(j)+=velocities(++index);
}
template<class TV> int DEFORMABLES_VELOCITY<TV>::
Raw_Size() const
{
    return V.Size()*TV::dimension+rigid_V.Size()*TWIST<TV>::dimension;
}
template<class TV> typename TV::SCALAR& DEFORMABLES_VELOCITY<TV>::
Raw_Get(int i)
{
    if(i<=V.Size()*TV::dimension) return V((i-1)/TV::dimension+1)((i-1)%TV::dimension+1);
    i-=V.Size()*TV::dimension;
    int o=(i-1)%TWIST<TV>::dimension+1,n=(i-1)/TWIST<TV>::dimension+1;
    if(o<=TV::dimension) return rigid_V(n).linear(o);
    return rigid_V(n).angular(o-TV::dimension);
}
//#####################################################################
template class DEFORMABLES_VELOCITY<VECTOR<float,1> >;
template class DEFORMABLES_VELOCITY<VECTOR<float,2> >;
template class DEFORMABLES_VELOCITY<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DEFORMABLES_VELOCITY<VECTOR<double,1> >;
template class DEFORMABLES_VELOCITY<VECTOR<double,2> >;
template class DEFORMABLES_VELOCITY<VECTOR<double,3> >;
#endif
