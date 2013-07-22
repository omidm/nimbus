//#####################################################################
// Copyright 2008, Avi Robinson-Mosher, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLID_SYSTEM_MPI
//#####################################################################
#ifndef __SOLID_SYSTEM_MPI__
#define __SOLID_SYSTEM_MPI__
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
namespace PhysBAM{
//#####################################################################
// Class SOLID_SYSTEM_MPI
//#####################################################################
template<class TV> class MPI_SOLID_FLUID;
template<class TV> class NEWMARK_EVOLUTION;
template<class TV>
class SOLID_SYSTEM_MPI:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;typedef typename RIGID_BODY_POLICY<TV>::WORLD_SPACE_INERTIA_TENSOR T_INERTIA_TENSOR;
    typedef typename MATRIX_POLICY<TV>::DIAGONAL_MATRIX T_DIAGONAL_MATRIX;

public:
    typedef GENERALIZED_VELOCITY<TV> VECTOR_T;

    BACKWARD_EULER_SYSTEM<TV>& solid_system;
    ARRAY<T_DIAGONAL_MATRIX>& fluid_mass;
    ARRAY<T_DIAGONAL_MATRIX>& rigid_body_fluid_mass;
    ARRAY<T_DIAGONAL_MATRIX> modified_mass,one_over_modified_mass;
    ARRAY<T_DIAGONAL_MATRIX> modified_world_space_rigid_mass,modified_world_space_rigid_mass_inverse;
    ARRAY<T_INERTIA_TENSOR> modified_world_space_rigid_inertia_tensor;
    ARRAY<T_INERTIA_TENSOR> modified_world_space_rigid_inertia_tensor_inverse;
    MPI_SOLID_FLUID<TV>* mpi_solid_fluid;
    ARRAY<ARRAY<int> >& coupled_deformable_particle_indices;
    mutable ARRAY<ARRAY<TV> > recv_fluid_V_boundary_arrays;
    mutable ARRAY<ARRAY<TWIST<TV> > > recv_fluid_rigid_V_boundary_arrays;

    NEWMARK_EVOLUTION<TV>& newmark_evolution;

    SOLID_SYSTEM_MPI(BACKWARD_EULER_SYSTEM<TV>& solid_system_input,ARRAY<T_DIAGONAL_MATRIX>& fluid_mass_input,ARRAY<T_DIAGONAL_MATRIX>& rigid_body_fluid_mass_input,
        ARRAY<T_INERTIA_TENSOR>& modified_world_space_rigid_inertia_tensor_input,MPI_SOLID_FLUID<TV>* mpi_solid_fluid_input,
        ARRAY<ARRAY<int> >& coupled_deformable_particle_indices_input,NEWMARK_EVOLUTION<TV>& newmark_evolution_input,const int rigid_V_size,bool precondition=true);

    ~SOLID_SYSTEM_MPI();

    const T Solid_Sign() const
    {return (T)-1;}

    static void Copy(T c,const VECTOR_T& c1,const VECTOR_T& c2,VECTOR_T& v)
    {v.Copy(c,c1,c2);}

    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE
    {R=V;}

    void Project(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE
    {solid_system.Project(V);}

    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE
    {return solid_system.Convergence_Norm(R);}

    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE {}

//#####################################################################
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& BV) const PHYSBAM_OVERRIDE;
    void Send_Generalized_Velocity_To_Fluid(const GENERALIZED_VELOCITY<TV>& V) const;
    void Get_Generalized_Velocity_From_Fluid(GENERALIZED_VELOCITY<TV>& V) const;
    void Multiply(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& F) const PHYSBAM_OVERRIDE;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& V1,const KRYLOV_VECTOR_BASE<T>& V2) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
