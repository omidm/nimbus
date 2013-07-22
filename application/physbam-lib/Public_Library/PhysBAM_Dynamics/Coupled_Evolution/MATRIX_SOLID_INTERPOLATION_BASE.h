//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_SOLID_INTERPOLATION_BASE
//#####################################################################
#ifndef __MATRIX_SOLID_INTERPOLATION_BASE__
#define __MATRIX_SOLID_INTERPOLATION_BASE__
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COUPLING_CONSTRAINT_ID.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>

namespace PhysBAM{

template<class TV> class GENERALIZED_VELOCITY;
template<class TV> class GENERALIZED_MASS;
template<class TV> class GRID;
template<class TV> class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO;

template<class TV>
class MATRIX_SOLID_INTERPOLATION_BASE:public NONCOPYABLE,public SYSTEM_MATRIX_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
protected:
    const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& iterator_info;

public:
    MATRIX_SOLID_INTERPOLATION_BASE(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info);
    virtual ~MATRIX_SOLID_INTERPOLATION_BASE();

    int V_size,rigid_V_size;
    const ARRAY<int>* V_indices,*rigid_V_indices;

//#####################################################################
    virtual COUPLING_CONSTRAINT_ID Number_Of_Constraints() const=0;
    virtual void Compute(const int ghost_cells)=0;
    virtual void Times_Add(const GENERALIZED_VELOCITY<TV>& solids,ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints) const=0;
    void Times(const GENERALIZED_VELOCITY<TV>& solids,ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints) const;
    virtual void Transpose_Times_Add(const ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints,GENERALIZED_VELOCITY<TV>& solids) const=0;
    void Transpose_Times(const ARRAY<T,COUPLING_CONSTRAINT_ID>& constraints,GENERALIZED_VELOCITY<TV>& solids) const;
    void Test_Matrix(int number_particles,int number_rigid_particles) const;
    virtual void Print_Each_Matrix(int n,GENERALIZED_VELOCITY<TV>& G) const=0;
    virtual void Add_Diagonal(ARRAY<T,COUPLING_CONSTRAINT_ID>& diagonal,const GENERALIZED_MASS<TV>& solid_mass) const=0;
    void Store_Maps(const GENERALIZED_VELOCITY<TV>& G);
//#####################################################################
};
}
#endif
