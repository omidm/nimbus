//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_FLUID_GRADIENT_BASE
//#####################################################################
#ifndef __MATRIX_FLUID_GRADIENT_BASE__
#define __MATRIX_FLUID_GRADIENT_BASE__
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>

namespace PhysBAM{

template<class TV> class COLLISION_AWARE_INDEX_MAP;
template<class TV> class GRID;

template<class TV>
class MATRIX_FLUID_GRADIENT_BASE:public NONCOPYABLE,public SYSTEM_MATRIX_BASE<typename TV::SCALAR>
{
    enum WORKAROUND {d=TV::dimension};
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,d> TV_INT;

public:
    COLLISION_AWARE_INDEX_MAP<TV>& index_map;
    SPARSE_MATRIX_FLAT_MXN<T> gradient;

    struct GHOST_GRADIENT_ENTRY
    {
        TV_INT index;
        int face;
        T weight;
    };
    ARRAY<GHOST_GRADIENT_ENTRY> ghost_gradient;
    struct INTERFACE_ENTRY
    {
        int face;
        T weight;
    };
    ARRAY<INTERFACE_ENTRY> interface_gradient;

    // TODO: the only reason index_map isn't const is because Array_View doesn't work for const quite yet
    MATRIX_FLUID_GRADIENT_BASE(COLLISION_AWARE_INDEX_MAP<TV>& index_map_input);
    virtual ~MATRIX_FLUID_GRADIENT_BASE();

//#####################################################################
    virtual void Compute(const ARRAY<bool,FACE_INDEX<d> >& psi_N_domain_boundary)=0;
    void Collect_Maxabs_Velocity(const VECTOR_ND<T>& faces,VECTOR_ND<T>& cells) const;
    void Times_Add(const VECTOR_ND<T>& cells,VECTOR_ND<T>& faces) const;
    void Times(const VECTOR_ND<T>& cells,VECTOR_ND<T>& faces) const;
    void Transpose_Times_Add(const VECTOR_ND<T>& faces,VECTOR_ND<T>& cells) const;
    void Transpose_Times(const VECTOR_ND<T>& faces,VECTOR_ND<T>& cells) const;
    void Test_Matrix() const;
    void Print_Each_Matrix(int n) const;
    void Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
