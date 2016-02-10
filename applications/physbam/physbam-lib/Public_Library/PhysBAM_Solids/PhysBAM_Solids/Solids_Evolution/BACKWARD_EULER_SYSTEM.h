//#####################################################################
// Copyright 2006-2008, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Craig Schroeder, Andrew Selle, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BACKWARD_EULER_SYSTEM
//#####################################################################
#ifndef __BACKWARD_EULER_SYSTEM__
#define __BACKWARD_EULER_SYSTEM__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_MASS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
namespace PhysBAM{
//#####################################################################
// Class GENERALIZED_MASS
//#####################################################################
template<class TV> class SOLIDS_EVOLUTION;
template<class TV>
class GENERALIZED_MASS
{
    typedef typename TV::SCALAR T;
public:
    INDIRECT_ARRAY<ARRAY_VIEW<T> > mass;
    INDIRECT_ARRAY<ARRAY_VIEW<T> > one_over_mass;
    INDIRECT_ARRAY<ARRAY_VIEW<T> > rigid_mass;
    INDIRECT_ARRAY<ARRAY_VIEW<typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR> > rigid_inertia_tensor;
    ARRAY<RIGID_BODY_MASS<TV,true> > world_space_rigid_mass;
    ARRAY<RIGID_BODY_MASS<TV,true> > world_space_rigid_mass_inverse_full;
    INDIRECT_ARRAY<ARRAY<RIGID_BODY_MASS<TV,true> > > world_space_rigid_mass_inverse;

    GENERALIZED_MASS(SOLID_BODY_COLLECTION<TV>& solid_body_collection);
    ~GENERALIZED_MASS();

    void Initialize_World_Space_Masses(const SOLID_BODY_COLLECTION<TV>& solid_body_collection);
    void Inverse_Multiply(const GENERALIZED_VELOCITY<TV>& V,GENERALIZED_VELOCITY<TV>& F,bool include_static) const;
};
//#####################################################################
// Class BACKWARD_EULER_SYSTEM
//#####################################################################
// assumes all solids_forces are linear in velocity, with a symmetric negative definite Jacobian.
template<class TV>
class BACKWARD_EULER_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;

public:
    typedef GENERALIZED_VELOCITY<TV> VECTOR_T;

    SOLIDS_EVOLUTION<TV>& solids_evolution;
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    T dt,current_velocity_time,current_position_time;
    ARTICULATED_RIGID_BODY<TV>* arb;
    MPI_SOLIDS<TV>* mpi_solids;
    TRIANGLE_REPULSIONS<TV>* repulsions;
    bool velocity_update;
    int project_nullspace_frequency;
    struct PROJECTION_DATA
    {
        PROJECTION_DATA(SOLID_BODY_COLLECTION<TV>& solid_body_collection)
            :mass(solid_body_collection),project_nullspace_counter(0)
        {}
        GENERALIZED_MASS<TV> mass;
        ARRAY<POINT_FACE_REPULSION_PAIR<TV> > point_face_pairs;
        ARRAY<EDGE_EDGE_REPULSION_PAIR<TV> > edge_edge_pairs;
        ARRAY<PRECOMPUTE_PROJECT_POINT_FACE<TV> > point_face_precomputed;
        ARRAY<PRECOMPUTE_PROJECT_EDGE_EDGE<TV> > edge_edge_precomputed;
        mutable int project_nullspace_counter;
    };
    PROJECTION_DATA projection_data;

    BACKWARD_EULER_SYSTEM(SOLIDS_EVOLUTION<TV>& solids_evolution_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection,const T dt,const T current_velocity_time,
        const T current_position_time,ARTICULATED_RIGID_BODY<TV>* arb_input,TRIANGLE_REPULSIONS<TV>* repulsions_input,MPI_SOLIDS<TV>* mpi_solids,const bool velocity_update_input);

    virtual ~BACKWARD_EULER_SYSTEM();

    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE
    {}

//#####################################################################
    void Force(const VECTOR_T& V,VECTOR_T& F) const;
    void Multiply(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& F) const PHYSBAM_OVERRIDE;
    void Set_Global_Boundary_Conditions(VECTOR_T& V,ARRAY<TV>& X_save,ARRAY<TV>& rigid_X_save,ARRAY<ROTATION<TV> >& rigid_rotation_save,
        ARRAY<TWIST<TV> >& rigid_velocity_save,ARRAY<typename TV::SPIN>& rigid_angular_momentum_save,ARRAY<TV>& V_save,bool test_system,
        bool print_matrix) const;// TODO: The meaning of this function has changed.
    void Project(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& V1,const KRYLOV_VECTOR_BASE<T>& V2) const PHYSBAM_OVERRIDE;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
