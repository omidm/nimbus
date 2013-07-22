//#####################################################################
// Copyright 2009;
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARTICULATED_SYSTEM
//#####################################################################
#ifndef __ARTICULATED_SYSTEM__
#define __ARTICULATED_SYSTEM__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_ID.h>
namespace PhysBAM{

template<class T> class KRYLOV_VECTOR_BASE;
template<class TV> class ARTICULATED_RIGID_BODY;
template<class TV> class ARTICULATED_VECTOR;
template<class TV>
class ARTICULATED_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef KRYLOV_SYSTEM_BASE<T> BASE;
public:
    ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body;
    ARRAY<MATRIX<T,TV::m>,JOINT_ID> prismatic_projection;
    ARRAY<MATRIX<T,TV::SPIN::m>,JOINT_ID> angular_projection;
    ARRAY<TV,JOINT_ID> location;
    ARRAY<MATRIX<T,TWIST<TV>::dimension>,JOINT_ID> effective_mass;
    mutable ARRAY<TWIST<TV> > intermediate_twists;
    bool break_loops;
    ARRAY<bool,JOINT_ID> keep_joint;
    ARTICULATED_VECTOR<TV>* internal_x;

    ARTICULATED_SYSTEM(ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body_input);

    virtual ~ARTICULATED_SYSTEM();

    void Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const PHYSBAM_OVERRIDE;
    void Gather(const ARRAY_VIEW<const TWIST<TV>,JOINT_ID> x,ARRAY_VIEW<TWIST<TV> > y) const;
    void Inverse_Mass(ARRAY_VIEW<TWIST<TV> > x) const;
    void Scatter(const ARRAY_VIEW<const TWIST<TV> > x,ARRAY_VIEW<TWIST<TV>,JOINT_ID> y) const;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const PHYSBAM_OVERRIDE;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE;
    void Project(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE;
    void Initialize();
    void Kinetic_Energy() const;
protected:
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
