//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_SOLID_FORCES
//#####################################################################
#ifndef __MATRIX_SOLID_FORCES__
#define __MATRIX_SOLID_FORCES__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/FORCE_AGGREGATE_ID.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SYSTEM_MATRIX_HELPER.h>

namespace PhysBAM{
template<class TV> class SOLID_BODY_COLLECTION;
template<class TV> class GENERALIZED_VELOCITY;

template<class TV>
class MATRIX_SOLID_FORCES:public NONCOPYABLE,public SYSTEM_MATRIX_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
public:
    const SOLID_BODY_COLLECTION<TV>& solid_body_collection;
private:
    T dt;
    T current_velocity_time;

    ARRAY<int> force_dof_counts;
    FORCE_AGGREGATE_ID total_force_dof;

public:
    MATRIX_SOLID_FORCES(const SOLID_BODY_COLLECTION<TV>& solid_body_collection);

//#####################################################################
    void Compute(const T dt_input,const T current_velocity_time_input);
    void Times_Add(const GENERALIZED_VELOCITY<TV>& solid_velocities,ARRAY<T,FORCE_AGGREGATE_ID>& force_coefficients) const;
    void Times(const GENERALIZED_VELOCITY<TV>& solid_velocities,ARRAY<T,FORCE_AGGREGATE_ID>& force_coefficients) const;
    void Transpose_Times_Add(const ARRAY<T,FORCE_AGGREGATE_ID>& force_coefficients,GENERALIZED_VELOCITY<TV>& solid_velocities) const;
    void Transpose_Times(const ARRAY<T,FORCE_AGGREGATE_ID>& force_coefficients,GENERALIZED_VELOCITY<TV>& solid_velocities) const;
    FORCE_AGGREGATE_ID Velocity_Dependent_Forces_Size() const;
    void Test_Matrix() const;
    void Print_Each_Matrix(int n,const GENERALIZED_VELOCITY<TV>& V,T sqrt_dt) const;
    void Add_Raw_Matrix(ARRAY<TRIPLE<int,int,T> >& data) const PHYSBAM_OVERRIDE;
private:
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V_full,ARRAY_VIEW<const TWIST<TV> > rigid_V_full,ARRAY_VIEW<T,FORCE_AGGREGATE_ID> aggregate,const T time) const;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T,FORCE_AGGREGATE_ID> aggregate,ARRAY_VIEW<TV> V_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const;
//#####################################################################
};
}
#endif
