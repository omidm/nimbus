//#####################################################################
// Copyright 2006-2009, Zhaosheng Bao, Elliot English, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Avi Robinson-Mosher, Craig Schroeder, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GENERALIZED_VELOCITY
//#####################################################################
#ifndef __GENERALIZED_VELOCITY__
#define __GENERALIZED_VELOCITY__

#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
namespace PhysBAM{
template<class TV> class SOLID_BODY_COLLECTION;
//#####################################################################
// Class GENERALIZED_VELOCITY
//#####################################################################
template<class TV>
class GENERALIZED_VELOCITY: public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef KRYLOV_VECTOR_BASE<T> BASE;
public:
    INDIRECT_ARRAY<ARRAY_VIEW<TV> > V;
    INDIRECT_ARRAY<ARRAY_VIEW<TWIST<TV> > > rigid_V;
    INDIRECT_ARRAY<ARRAY_VIEW<TWIST<TV> > > kinematic_and_static_rigid_V;  // TODO: Go away.

    GENERALIZED_VELOCITY(ARRAY_VIEW<TV> V_full,ARRAY_VIEW<TWIST<TV> > rigid_V_full,const SOLID_BODY_COLLECTION<TV>& solid_body_collection);
    GENERALIZED_VELOCITY(ARRAY_VIEW<TV> V_full,ARRAY<int>& dynamic_particles,ARRAY_VIEW<TWIST<TV> > rigid_V_full,ARRAY<int>& dynamic_rigid_body_particles,
        ARRAY<int>& static_and_kinematic_rigid_bodies);
    ~GENERALIZED_VELOCITY();

    BASE& operator+=(const BASE& bv) PHYSBAM_OVERRIDE;
    BASE& operator-=(const BASE& bv) PHYSBAM_OVERRIDE;
    BASE& operator*=(const T a) PHYSBAM_OVERRIDE;
    void Copy(const T c,const BASE& bv) PHYSBAM_OVERRIDE;
    void Copy(const T c1,const BASE& bv1,const BASE& bv2) PHYSBAM_OVERRIDE;
    void Pack(VECTOR_ND<T> &velocities) const;
    void Unpack(VECTOR_ND<T> &velocities);
    void Unpack_And_Add(VECTOR_ND<T> &velocities);
    int Raw_Size() const PHYSBAM_OVERRIDE;
    T& Raw_Get(int i) PHYSBAM_OVERRIDE;
    void Exchange(GENERALIZED_VELOCITY<TV>& gv);
};
}
#endif
