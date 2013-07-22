//#####################################################################
// Copyright 2010, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_VISCOUS_FORCES
//#####################################################################
#ifndef __MATRIX_VISCOUS_FORCES__
#define __MATRIX_VISCOUS_FORCES__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/VISCOUS_FORCE_ID.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>

namespace PhysBAM{
template<class TV> class COLLISION_AWARE_INDEX_MAP;
template<class TV> class GRID;

template<class TV>
class MATRIX_VISCOUS_FORCES:public NONCOPYABLE,public SYSTEM_MATRIX_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::dimension> TV_INT;
    enum WORKAROUND{d=TV::dimension};
private:
    const GRID<TV>& grid;
    const COLLISION_AWARE_INDEX_MAP<TV>& index_map;
    struct ENTRY
    {
        T weight;
        int face_index;
        VISCOUS_FORCE_ID viscous_id;
        
        ENTRY()
        {}
        
        ENTRY(T weight_input,int face_index_input,VISCOUS_FORCE_ID id)
            :weight(weight_input),face_index(face_index_input),viscous_id(id)
        {}
    };

    ARRAY<ENTRY> entries;
    VISCOUS_FORCE_ID last_id;

public:
    MATRIX_VISCOUS_FORCES(const COLLISION_AWARE_INDEX_MAP<TV>& index_map_input);

//#####################################################################
    void Compute(const T dt,const ARRAY<bool,FACE_INDEX<d> >& psi_N,T mu);
    void Times_Add(const VECTOR_ND<T>& velocities,ARRAY<T,VISCOUS_FORCE_ID>& viscous_force_coefficients) const;
    void Times(const VECTOR_ND<T>& velocities,ARRAY<T,VISCOUS_FORCE_ID>& viscous_force_coefficients) const;
    void Transpose_Times_Add(const ARRAY<T,VISCOUS_FORCE_ID>& viscous_force_coefficients,VECTOR_ND<T>& velocities) const;
    void Transpose_Times(const ARRAY<T,VISCOUS_FORCE_ID>& viscous_force_coefficients,VECTOR_ND<T>& velocities) const;
    VISCOUS_FORCE_ID Viscous_Forces_Size() const;
    void Test_Matrix() const;
    void Print_Each_Matrix(int n) const;
    void Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
