//#####################################################################
// Copyright 2011, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_PARTICLE_COUPLED_SYSTEM
//#####################################################################
#ifndef __MPI_PARTICLE_COUPLED_SYSTEM__
#define __MPI_PARTICLE_COUPLED_SYSTEM__

#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/DOT_PRODUCT.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PARTICLES.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>

namespace PhysBAM
{

template<class TV>
class MPI_PARTICLE_COUPLED_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef KRYLOV_VECTOR_WRAPPER<T,VECTOR_ND<T>&> KRYLOV_VECTOR_T;

public:
    const SPARSE_MATRIX_FLAT_MXN<T>& K;
    const VECTOR_ND<T>& Mi;
    VECTOR_ND<T>& temp_v;
    VECTOR<int,2> diagonal_range;
    const ARRAY<int>& filter;//0: interior particles; 1: ghost; 2: external
    ARRAY<int> backward_filter;// filters on other processors
    const ARRAY<int>& mapping;//from vertex index to particle index
    const ARRAY<int>& back_mapping;//from particle index to vertex index
    MPI_PARTICLES<GRID<TV> >& mpi_particles;
    ARRAY<int> interior_indices,cross_over_indices;
    ARRAY<ARRAY<int> > send_indices,recv_indices,send_indices_scalar,recv_indices_scalar;

    MPI_PARTICLE_COUPLED_SYSTEM(const SPARSE_MATRIX_FLAT_MXN<T>& K_input,const VECTOR_ND<T>& Mi_input,VECTOR_ND<T>& temp_v_input,const VECTOR<int,2>& diagonal_range_input,const ARRAY<int>& filter_input,const ARRAY<int>& mapping_input,const ARRAY<int>& back_mapping_input,MPI_PARTICLES<GRID<TV> >& mpi_particles_input);
    void Fill_Vector_Ghost_Values(VECTOR_ND<T>& v)const;
    void Fill_Scalar_Ghost_Values(VECTOR_ND<T>& v)const;
    void Multiply(const KRYLOV_VECTOR_BASE<T>& bp,KRYLOV_VECTOR_BASE<T>& bresult) const PHYSBAM_OVERRIDE;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& bx,const KRYLOV_VECTOR_BASE<T>& by) const PHYSBAM_OVERRIDE;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bx) const PHYSBAM_OVERRIDE;
//#####################################################################
    void Project(KRYLOV_VECTOR_BASE<T>& p) const PHYSBAM_OVERRIDE{}
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& p) const PHYSBAM_OVERRIDE{Project(p);}
    virtual void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& bx) const{}
};
}
#endif
