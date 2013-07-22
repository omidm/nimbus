//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_SOLID_INTERPOLATION
//#####################################################################
#ifndef __MATRIX_SOLID_INTERPOLATION__
#define __MATRIX_SOLID_INTERPOLATION__
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COUPLING_CONSTRAINT_ID.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/MATRIX_SOLID_INTERPOLATION_BASE.h>

namespace PhysBAM{

template<class TV>
class MATRIX_SOLID_INTERPOLATION:public MATRIX_SOLID_INTERPOLATION_BASE<TV>
{
    typedef typename TV::SCALAR T;
    using MATRIX_SOLID_INTERPOLATION_BASE<TV>::iterator_info;

    struct WEIGHT_HELPER
    {
        int index;
        T weight;
        WEIGHT_HELPER(int i,T w):index(i),weight(w) {}
        bool operator<(const WEIGHT_HELPER& w) const {return index<w.index;}
        void Merge(const WEIGHT_HELPER& w){assert(index==w.index);weight+=w.weight;}
    };

    struct DEFORMABLE_WEIGHT:public WEIGHT_HELPER
    {
        DEFORMABLE_WEIGHT(int i=0,T w=0):WEIGHT_HELPER(i,w) {}
    };

    struct RIGID_WEIGHT:public WEIGHT_HELPER
    {
        TV radius;
        RIGID_WEIGHT(int i=0,T w=0,const TV& r=TV()):WEIGHT_HELPER(i,w),radius(r) {}
    };

    struct ROW
    {
        int axis;
        TV normal;
        ARRAY<DEFORMABLE_WEIGHT> deformable_weights;
        ARRAY<RIGID_WEIGHT> rigid_weights;
    };

public:
    ARRAY<ROW,COUPLING_CONSTRAINT_ID> rows;

    MATRIX_SOLID_INTERPOLATION(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info);
    virtual ~MATRIX_SOLID_INTERPOLATION();

//#####################################################################
    virtual COUPLING_CONSTRAINT_ID Number_Of_Constraints() const PHYSBAM_OVERRIDE;
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
