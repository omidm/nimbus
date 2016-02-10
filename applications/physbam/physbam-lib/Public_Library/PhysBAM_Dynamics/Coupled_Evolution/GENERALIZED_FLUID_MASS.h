//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GENERALIZED_FLUID_MASS
//#####################################################################
#ifndef __GENERALIZED_FLUID_MASS__
#define __GENERALIZED_FLUID_MASS__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>

namespace PhysBAM{
template<class TV> class COLLISION_AWARE_INDEX_MAP;
template<class TV> class GRID;

template<class TV>
class GENERALIZED_FLUID_MASS:public NONCOPYABLE,public SYSTEM_MATRIX_BASE<typename TV::SCALAR>
{
    enum WORKAROUND {d=TV::dimension};
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;

    const COLLISION_AWARE_INDEX_MAP<TV>& index_map;
    const T_FACE_ARRAYS_SCALAR& beta; // 1/rho
    const ARRAY<T>& constrained_beta;

public:
    VECTOR_ND<T> one_over_fluid_mass_at_faces;

    GENERALIZED_FLUID_MASS(const COLLISION_AWARE_INDEX_MAP<TV>& index_map_input,const T_FACE_ARRAYS_SCALAR& beta_input,const ARRAY<T>& constrained_beta_input);

    void Compute();
    void Inverse_Times(const VECTOR_ND<T>& faces_in,VECTOR_ND<T>& faces_out) const;
    void Inverse_Times_Add(const VECTOR_ND<T>& faces_in,VECTOR_ND<T>& faces_out) const;
    void Second_Order_Mass_Correction(const ARRAY<T,TV_INT>& phi);    
    void Compute_For_Two_Phase_Pressure_Jump(const ARRAY<T,TV_INT>& phi,const ARRAY<T,TV_INT>& density);

    void Print_Each_Matrix(int n) const;
    void Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const PHYSBAM_OVERRIDE;
};
}
#endif
