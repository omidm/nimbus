//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_TETRAHEDRALIZED_VOLUME
//#####################################################################
//
// Solves -laplace u = f, with Neumann conditions on all boundary faces.
// f must sum to 0 for compatibility.
//
//#####################################################################  
#ifndef __LAPLACE_TETRAHEDRALIZED_VOLUME__
#define __LAPLACE_TETRAHEDRALIZED_VOLUME__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_PDE_Linear/LAPLACE.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
namespace PhysBAM{

template<class T> class SPARSE_MATRIX_FLAT_NXN;
template<class T> class TETRAHEDRALIZED_VOLUME;

template<class T>
class LAPLACE_TETRAHEDRALIZED_VOLUME:public LAPLACE<T>
{
public:
    using LAPLACE<T>::tolerance;
    
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume;
    ARRAY<T>& u;
    ARRAY<T> f;
    PCG_SPARSE<T> pcg;

    LAPLACE_TETRAHEDRALIZED_VOLUME(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume_input,ARRAY<T>& u_input)
        :tetrahedralized_volume(tetrahedralized_volume_input),u(u_input)
    {}

    virtual ~LAPLACE_TETRAHEDRALIZED_VOLUME()
    {}

//#####################################################################
    void Initialize_Tetrahedralized_Volume();
    virtual void Solve(const T time);
    void Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& b,const bool enforce_compatability=true,const bool recompute_preconditioner=true);
    virtual void Find_A(SPARSE_MATRIX_FLAT_NXN<T>& A);
    virtual void Find_b(VECTOR_ND<T>& b);
    virtual void Add_Scalar_Term(SPARSE_MATRIX_FLAT_NXN<T>& A,const T s);
//#####################################################################
};
}
#endif

