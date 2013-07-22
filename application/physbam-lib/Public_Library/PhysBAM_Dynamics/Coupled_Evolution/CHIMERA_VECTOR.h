//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAGRANGE_MULTIPLIERS_VELOCITY_VECTOR
//#####################################################################
#ifndef __CHIMERA_VECTOR__
#define __CHIMERA_VECTOR__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
namespace PhysBAM{
template<class T> class VECTOR_ND;
template<class TV> class GENERALIZED_VELOCITY;

template<class TV>
class CHIMERA_VECTOR : public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef KRYLOV_VECTOR_BASE<T> BASE;
public:
    VECTOR_ND<T> pressure;
    VECTOR_ND<T> kinematic_force;

    CHIMERA_VECTOR(VECTOR_ND<T>& pressure_input,VECTOR_ND<T>& kinematic_force_input):
        pressure(pressure_input),kinematic_force(kinematic_force_input) {}
    CHIMERA_VECTOR(int n_pressure,int n_kinematic_force):
        pressure(n_pressure),kinematic_force(n_kinematic_force) {}
    CHIMERA_VECTOR() {}

    virtual ~CHIMERA_VECTOR() {}

    CHIMERA_VECTOR& operator=(const CHIMERA_VECTOR& bv)
    {
        const CHIMERA_VECTOR<TV>& v=debug_cast<const CHIMERA_VECTOR<TV>&>(bv);
        pressure=v.pressure;
        kinematic_force=v.kinematic_force;
        return *this;
    }

    BASE& operator+=(const BASE& bv) PHYSBAM_OVERRIDE
    {
        const CHIMERA_VECTOR<TV>& v=debug_cast<const CHIMERA_VECTOR<TV>&>(bv);
        pressure+=v.pressure;
        kinematic_force+=v.kinematic_force;
        return *this;
    }
    BASE& operator-=(const BASE& bv) PHYSBAM_OVERRIDE
    {
        const CHIMERA_VECTOR<TV>& v=debug_cast<const CHIMERA_VECTOR<TV>&>(bv);
        pressure-=v.pressure;
        kinematic_force-=v.kinematic_force;
        return *this;
    }
    BASE& operator*=(const T a) PHYSBAM_OVERRIDE
    {
        pressure*=a;
        kinematic_force*=a;
        return *this;
    }
    void Copy(const T c,const BASE& bv) PHYSBAM_OVERRIDE
    {
        const CHIMERA_VECTOR<TV>& v=debug_cast<const CHIMERA_VECTOR<TV>&>(bv);
        VECTOR_ND<T>::Copy(c,v.pressure,pressure);
        VECTOR_ND<T>::Copy(c,v.kinematic_force,kinematic_force);
    }
    void Copy(const T c1,const BASE& bv1,const BASE& bv2) PHYSBAM_OVERRIDE
    {
        const CHIMERA_VECTOR<TV>& v1=debug_cast<const CHIMERA_VECTOR<TV>&>(bv1);
        const CHIMERA_VECTOR<TV>& v2=debug_cast<const CHIMERA_VECTOR<TV>&>(bv2);
        VECTOR_ND<T>::Copy(c1,v1.pressure,v2.pressure,pressure);
        VECTOR_ND<T>::Copy(c1,v1.kinematic_force,v2.kinematic_force,kinematic_force);
    }
    void Print() const
    {
    }
    int Raw_Size() const PHYSBAM_OVERRIDE
    {
        return pressure.Size()+kinematic_force.Size();
    }
    T& Raw_Get(int i)
    {
        if(i<=pressure.Size())
            return pressure(i);
        else
            return kinematic_force(i-(pressure.Size()));
    }
};
}
#endif
