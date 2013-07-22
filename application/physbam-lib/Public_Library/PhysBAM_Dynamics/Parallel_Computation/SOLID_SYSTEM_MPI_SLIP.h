//#####################################################################
// Copyright 2008-2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_SYSTEM_MPI_SLIP
//#####################################################################
#ifndef __SOLID_SYSTEM_MPI_SLIP__
#define __SOLID_SYSTEM_MPI_SLIP__
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE_SYSTEM.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
namespace PhysBAM{
//#####################################################################
// Class SOLID_SYSTEM_MPI_SLIP
//#####################################################################
template<class TV> class MPI_SOLID_FLUID_SLIP;
template<class TV> class NEWMARK_EVOLUTION;
template<class TV>
class SOLID_SYSTEM_MPI_SLIP:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef GENERALIZED_VELOCITY<TV> VECTOR_T;
    typedef KRYLOV_SYSTEM_BASE<T> BASE;
public:

    BACKWARD_EULER_SYSTEM<TV>& solid_system;
    MPI_SOLID_FLUID_SLIP<TV>* mpi_solid_fluid;
    ARRAY<ARRAY<int> >& coupled_deformable_particle_indices;
    mutable ARRAY<ARRAY<TV> > recv_fluid_V_boundary_arrays;
    mutable ARRAY<ARRAY<TWIST<TV> > > recv_fluid_rigid_V_boundary_arrays;
    GENERALIZED_MASS<TV>* solid_mass;
    mutable VECTOR_ND<T> packed_solid_velocities;

    NEWMARK_EVOLUTION<TV>& newmark_evolution;

    SOLID_SYSTEM_MPI_SLIP(const bool use_preconditioner_input,BACKWARD_EULER_SYSTEM<TV>& solid_system_input,MPI_SOLID_FLUID_SLIP<TV>* mpi_solid_fluid_input,
        ARRAY<ARRAY<int> >& coupled_deformable_particle_indices_input,NEWMARK_EVOLUTION<TV>& newmark_evolution_input,const int rigid_V_size);

    virtual ~SOLID_SYSTEM_MPI_SLIP();

    static void Copy(T c,const VECTOR_T& c1,const VECTOR_T& c2,VECTOR_T& v)
    {v.Copy(c,c1,c2);}

    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BR) const PHYSBAM_OVERRIDE
    {const VECTOR_T& V=debug_cast<const VECTOR_T&>(BV);VECTOR_T& R=debug_cast<VECTOR_T&>(BR);
    R.Copy(1,V);}

    void Project(KRYLOV_VECTOR_BASE<T>& BV) const PHYSBAM_OVERRIDE
    {VECTOR_T& V=debug_cast<VECTOR_T&>(BV);solid_system.Project(V);}

//#####################################################################
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& BV) const PHYSBAM_OVERRIDE;
    void Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const PHYSBAM_OVERRIDE;
    void Send_Generalized_Velocity_To_Fluid(const GENERALIZED_VELOCITY<TV>& V) const;
    void Get_Generalized_Velocity_From_Fluid(GENERALIZED_VELOCITY<TV>& V) const;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV1,const KRYLOV_VECTOR_BASE<T>& BV2) const PHYSBAM_OVERRIDE;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const PHYSBAM_OVERRIDE;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& BV) const PHYSBAM_OVERRIDE {}
//#####################################################################
};
}
#endif
