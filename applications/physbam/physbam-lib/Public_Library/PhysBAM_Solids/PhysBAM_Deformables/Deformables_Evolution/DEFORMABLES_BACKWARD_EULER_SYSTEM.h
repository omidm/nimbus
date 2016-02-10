//#####################################################################
// Copyright 2006-2008, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Craig Schroeder, Andrew Selle, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLES_BACKWARD_EULER_SYSTEM
//#####################################################################
#ifndef __DEFORMABLES_BACKWARD_EULER_SYSTEM__
#define __DEFORMABLES_BACKWARD_EULER_SYSTEM__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformables_Evolution/DEFORMABLES_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformables_Evolution/DEFORMABLES_VELOCITY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
namespace PhysBAM{
//#####################################################################
// Class RIGIDS_MASS
//#####################################################################
template<class TV> class RIGIDS_EVOLUTION;
template<class TV>
class DEFORMABLES_MASS
{
    typedef typename TV::SCALAR T;
public:
    INDIRECT_ARRAY<ARRAY_VIEW<T> > mass;
    INDIRECT_ARRAY<ARRAY_VIEW<T> > one_over_mass;

    DEFORMABLES_MASS(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection);

    void Inverse_Multiply(const DEFORMABLES_VELOCITY<TV>& V,DEFORMABLES_VELOCITY<TV>& F) const;
};
//#####################################################################
// Class DEFORMABLES_BACKWARD_EULER_SYSTEM
//#####################################################################
// assumes all solids_forces are linear in velocity, with a symmetric negative definite Jacobian.
template<class TV> struct POINT_FACE_REPULSION_PAIR;
template<class TV> struct EDGE_EDGE_REPULSION_PAIR;
template<class TV> struct PRECOMPUTE_PROJECT_POINT_FACE;
template<class TV> struct PRECOMPUTE_PROJECT_EDGE_EDGE;
template<class TV>
class DEFORMABLES_BACKWARD_EULER_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;

public:
    typedef DEFORMABLES_VELOCITY<TV> VECTOR_T;

    DEFORMABLES_EVOLUTION<TV>& deformables_evolution;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection;
    T dt,current_velocity_time,current_position_time;
    TRIANGLE_REPULSIONS<TV>* repulsions;
    MPI_SOLIDS<TV>* mpi_solids;
    bool velocity_update;
    int project_nullspace_frequency;
    struct PROJECTION_DATA
    {
        PROJECTION_DATA(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection)
            :mass(deformable_body_collection),project_nullspace_counter(0)
        {}
        DEFORMABLES_MASS<TV> mass;
        ARRAY<POINT_FACE_REPULSION_PAIR<TV> > point_face_pairs;
        ARRAY<EDGE_EDGE_REPULSION_PAIR<TV> > edge_edge_pairs;
        ARRAY<PRECOMPUTE_PROJECT_POINT_FACE<TV> > point_face_precomputed;
        ARRAY<PRECOMPUTE_PROJECT_EDGE_EDGE<TV> > edge_edge_precomputed;
        mutable int project_nullspace_counter;
    };
    PROJECTION_DATA projection_data;

    DEFORMABLES_BACKWARD_EULER_SYSTEM(DEFORMABLES_EVOLUTION<TV>& deformables_evolution_input,DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection,const T dt,
        const T current_velocity_time,const T current_position_time,MPI_SOLIDS<TV>* mpi_solids_input,TRIANGLE_REPULSIONS<TV>* repulsions_input,const bool velocity_update_input);

    virtual ~DEFORMABLES_BACKWARD_EULER_SYSTEM();

    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE
    {}

//#####################################################################
    void Force(const VECTOR_T& V,VECTOR_T& F) const;
    void Multiply(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& F) const PHYSBAM_OVERRIDE;
    void Set_Global_Boundary_Conditions(VECTOR_T& V) const;
    void Project(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& V1,const KRYLOV_VECTOR_BASE<T>& V2) const PHYSBAM_OVERRIDE;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
