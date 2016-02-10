//#####################################################################
// Copyright 2006-2008, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Craig Schroeder, Andrew Selle, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_BACKWARD_EULER_SYSTEM
//#####################################################################
#ifndef __RIGIDS_BACKWARD_EULER_SYSTEM__
#define __RIGIDS_BACKWARD_EULER_SYSTEM__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_POLICY.h>
namespace PhysBAM{
//#####################################################################
// Class RIGIDS_MASS
//#####################################################################
template<class TV> class RIGIDS_EVOLUTION;
template<class TV> class RIGIDS_VELOCITY;
template<class TV> class RIGID_BODY_COLLECTION;
template<class TV,bool b> class RIGID_BODY_MASS;
template<class TV> class ARTICULATED_RIGID_BODY;
template<class TV>
class RIGIDS_MASS
{
    typedef typename TV::SCALAR T;
public:
    INDIRECT_ARRAY<ARRAY_VIEW<T> > rigid_mass;
    INDIRECT_ARRAY<ARRAY_VIEW<typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR> > rigid_inertia_tensor;
    ARRAY<RIGID_BODY_MASS<TV,true> > world_space_rigid_mass;
    ARRAY<RIGID_BODY_MASS<TV,true> > world_space_rigid_mass_inverse;

    RIGIDS_MASS(RIGID_BODY_COLLECTION<TV>& rigid_body_collection);

    void Initialize_World_Space_Masses(const RIGID_BODY_COLLECTION<TV>& rigid_body_collection);
    void Inverse_Multiply(const RIGIDS_VELOCITY<TV>& V,RIGIDS_VELOCITY<TV>& F) const;
};
//#####################################################################
// Class RIGIDS_BACKWARD_EULER_SYSTEM
//#####################################################################
// assumes all solids_forces are linear in velocity, with a symmetric negative definite Jacobian.
template<class TV>
class RIGIDS_BACKWARD_EULER_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;

public:
    typedef RIGIDS_VELOCITY<TV> VECTOR_T;

    RIGIDS_EVOLUTION<TV>& rigids_evolution;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    T dt,current_velocity_time,current_position_time;
    ARTICULATED_RIGID_BODY<TV>* arb;
    bool velocity_update;
    int project_nullspace_frequency;
    struct PROJECTION_DATA
    {
        PROJECTION_DATA(RIGID_BODY_COLLECTION<TV>& rigid_body_collection)
            :mass(rigid_body_collection),project_nullspace_counter(0)
        {}
        RIGIDS_MASS<TV> mass;
        mutable int project_nullspace_counter;
    };
    PROJECTION_DATA projection_data;

    RIGIDS_BACKWARD_EULER_SYSTEM(RIGIDS_EVOLUTION<TV>& rigids_evolution_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const T dt,const T current_velocity_time,
        const T current_position_time,ARTICULATED_RIGID_BODY<TV>* arb_input,const bool velocity_update_input,const bool enforce_poststabilization_in_cg);

    virtual ~RIGIDS_BACKWARD_EULER_SYSTEM();

    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE
    {}

//#####################################################################
    void Force(const VECTOR_T& V,VECTOR_T& F) const;
    void Multiply(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& F) const PHYSBAM_OVERRIDE;
    void Set_Global_Boundary_Conditions(VECTOR_T& V,ARRAY<TV>& rigid_X_save,ARRAY<ROTATION<TV> >& rigid_rotation_save,
        ARRAY<TWIST<TV> >& rigid_velocity_save,ARRAY<typename TV::SPIN>& rigid_angular_momentum_save,bool test_system,bool print_matrix) const;// TODO: The meaning of this function has changed.
    void Project(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& V1,const KRYLOV_VECTOR_BASE<T>& V2) const PHYSBAM_OVERRIDE;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
