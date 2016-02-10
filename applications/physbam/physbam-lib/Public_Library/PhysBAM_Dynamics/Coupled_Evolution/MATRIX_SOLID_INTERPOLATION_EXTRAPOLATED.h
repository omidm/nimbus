//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED
//#####################################################################
#ifndef __MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED__
#define __MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED__
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COUPLING_CONSTRAINT_ID.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_SOLID_INTERPOLATION_BASE.h>

namespace PhysBAM{

template<class TV>
class MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED:public MATRIX_SOLID_INTERPOLATION_BASE<TV>
{
    typedef typename TV::SCALAR T;
    using MATRIX_SOLID_INTERPOLATION_BASE<TV>::iterator_info;

public:
    ARRAY<int> entries;

    MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info);
    virtual ~MATRIX_SOLID_INTERPOLATION_EXTRAPOLATED();

//#####################################################################
    COUPLING_CONSTRAINT_ID Number_Of_Constraints() const PHYSBAM_OVERRIDE;
    void Compute(const int ghost_cells) PHYSBAM_OVERRIDE;
    void Times_Add(const GENERALIZED_VELOCITY<TV>& solids,ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints) const PHYSBAM_OVERRIDE;
    void Transpose_Times_Add(const ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints,GENERALIZED_VELOCITY<TV>& solids) const PHYSBAM_OVERRIDE;
    void Print_Each_Matrix(int n,GENERALIZED_VELOCITY<TV>& G) const PHYSBAM_OVERRIDE;
    void Add_Diagonal(ARRAY<T,COUPLING_CONSTRAINT_ID>& diagonal,const GENERALIZED_MASS<TV>& solid_mass) const PHYSBAM_OVERRIDE;
    void Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
