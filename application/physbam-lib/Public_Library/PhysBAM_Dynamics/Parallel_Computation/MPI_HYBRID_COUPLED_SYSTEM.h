//#####################################################################
// Copyright 2012, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_HYBRID_COUPLED_SYSTEM
//#####################################################################
#ifndef __MPI_HYBRID_COUPLED_SYSTEM__
#define __MPI_HYBRID_COUPLED_SYSTEM__

#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/DOT_PRODUCT.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>

namespace PhysBAM
{

template<class TV>
class MPI_HYBRID_COUPLED_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef KRYLOV_VECTOR_WRAPPER<T,VECTOR_ND<T>&> KRYLOV_VECTOR_T;

protected:
    ARRAY<int> neighbor_ranks;
    ARRAY<TV_INT> neighbor_directions;
    ARRAY<ARRAY<int> > send_indices,recv_indices;
    
public:
    const SPARSE_MATRIX_FLAT_NXN<T>& A;
    ARRAY<int> interior;
    MPI_GRID<GRID<TV> >& mpi_grid;

    MPI_HYBRID_COUPLED_SYSTEM(const SPARSE_MATRIX_FLAT_NXN<T>& A_input,const ARRAY<TV>& positions,const ARRAY<int>& global_ids,MPI_GRID<GRID<TV> >& mpi_grid,bool use_diagonal_preconditioner=false);
    void Fill_Ghost_Values(VECTOR_ND<T>& v)const;
    void Multiply(const KRYLOV_VECTOR_BASE<T>& bp,KRYLOV_VECTOR_BASE<T>& bresult) const PHYSBAM_OVERRIDE;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& bx,const KRYLOV_VECTOR_BASE<T>& by) const PHYSBAM_OVERRIDE;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bx) const PHYSBAM_OVERRIDE;
    void Project(KRYLOV_VECTOR_BASE<T>& p) const PHYSBAM_OVERRIDE{}
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& p) const PHYSBAM_OVERRIDE{Project(p);}
    virtual void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& bx) const{}
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& br,KRYLOV_VECTOR_BASE<T>& bz) const PHYSBAM_OVERRIDE // solve Mz=r
    {
        const KRYLOV_VECTOR_T& rr=debug_cast<const KRYLOV_VECTOR_T&>(br);KRYLOV_VECTOR_T& zz=debug_cast<KRYLOV_VECTOR_T&>(bz);
        for(int i=1;i<=interior.m;i++){int j=interior(i); zz.v(j)=rr.v(j)/A.A(A.diagonal_index(j)).a;}
        Fill_Ghost_Values(zz.v);
    }
//#####################################################################
protected:
    TV_INT Get_Direction(const TV& position)
    {
        TV_INT direction;
        for(int axis=1;axis<=TV::m;axis++) for(int axis_side=1;axis_side<=2;axis_side++)
            if(mpi_grid.Neighbor(axis,axis_side)){ //half open if it has neighbor
                if(axis_side==1 && position(axis)<mpi_grid.local_grid.domain.min_corner(axis)) direction(axis)--;
                else if(axis_side==2 && position(axis)>=mpi_grid.local_grid.domain.max_corner(axis)) direction(axis)++;}
            else{
                if(axis_side==1 && position(axis)<mpi_grid.local_grid.domain.min_corner(axis)) direction(axis)--;
                else if(axis_side==2 && position(axis)>mpi_grid.local_grid.domain.max_corner(axis)) direction(axis)++;}
        return direction;
    }
};
}
#endif
