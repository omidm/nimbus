//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_TO_SOLID_INTERPOLATION_BASE
//#####################################################################
#ifndef __FLUID_TO_SOLID_INTERPOLATION_BASE__
#define __FLUID_TO_SOLID_INTERPOLATION_BASE__
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COUPLING_CONSTRAINT_ID.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>

namespace PhysBAM{

template<class TV> class GENERALIZED_VELOCITY;
template<class TV> class GENERALIZED_MASS;
template<class TV> class GRID;
template<class TV> class COLLISION_AWARE_INDEX_MAP;

template<class TV>
class FLUID_TO_SOLID_INTERPOLATION_BASE:public NONCOPYABLE,public SYSTEM_MATRIX_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
protected:
    const COLLISION_AWARE_INDEX_MAP<TV>& index_map;

public:
    FLUID_TO_SOLID_INTERPOLATION_BASE(const COLLISION_AWARE_INDEX_MAP<TV>& map);
    virtual ~FLUID_TO_SOLID_INTERPOLATION_BASE();

    int V_size;
    const ARRAY<int>* V_indices;

//#####################################################################
    virtual void Compute(const int ghost_cells)=0;
    virtual void Times_Add(const VECTOR_ND<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& solid_velocity) const=0;
    void Times(const VECTOR_ND<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& solid_velocity) const;
    virtual void Transpose_Times_Add(const GENERALIZED_VELOCITY<TV>& solid_force,VECTOR_ND<T>& fluid_force) const=0;
    void Transpose_Times(const GENERALIZED_VELOCITY<TV>& solid_force,VECTOR_ND<T>& fluid_force) const;
    void Test_Matrix(int number_fluid_faces,int number_particles,int number_rigid_particles) const;
    virtual void Print_Each_Matrix(int n,int fluid_faces,GENERALIZED_VELOCITY<TV>& G) const=0;
    void Store_Maps(const GENERALIZED_VELOCITY<TV>& G);
//#####################################################################
};
}
#endif
