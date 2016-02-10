//#####################################################################
// Copyright 2006-2008, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Craig Schroeder, Andrew Selle, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PRESSURE_VELOCITY_VECTOR
//#####################################################################
#ifndef __PRESSURE_VELOCITY_VECTOR__
#define __PRESSURE_VELOCITY_VECTOR__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
namespace PhysBAM{
//#####################################################################
// Class PRESSURE_VELOCITY_VECTOR
//#####################################################################
template<class TV>
class PRESSURE_VELOCITY_VECTOR: public KRYLOV_VECTOR_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef KRYLOV_VECTOR_BASE<typename TV::SCALAR> BASE;
public:
    GENERALIZED_VELOCITY<TV>& solid_velocity;
    ARRAY<VECTOR_ND<T> >& pressure;

    PRESSURE_VELOCITY_VECTOR(GENERALIZED_VELOCITY<TV>& solid_velocity_input,ARRAY<VECTOR_ND<T> >& pressure_input)
        :solid_velocity(solid_velocity_input),pressure(pressure_input)
    {}

    PRESSURE_VELOCITY_VECTOR& operator=(const PRESSURE_VELOCITY_VECTOR& v)
    {solid_velocity=v.solid_velocity;pressure=v.pressure;return *this;}

    BASE& operator+=(const BASE& bv) PHYSBAM_OVERRIDE
    {const PRESSURE_VELOCITY_VECTOR& v=debug_cast<const PRESSURE_VELOCITY_VECTOR&>(bv);solid_velocity+=v.solid_velocity;pressure+=v.pressure;return *this;}

    BASE& operator-=(const BASE& bv) PHYSBAM_OVERRIDE
    {const PRESSURE_VELOCITY_VECTOR& v=debug_cast<const PRESSURE_VELOCITY_VECTOR&>(bv);solid_velocity-=v.solid_velocity;pressure-=v.pressure;return *this;}

    BASE& operator*=(const T a) PHYSBAM_OVERRIDE
    {solid_velocity*=a;pressure*=a;return *this;}

    void Copy(const T c,const BASE& bv) PHYSBAM_OVERRIDE
    {const PRESSURE_VELOCITY_VECTOR& v=debug_cast<const PRESSURE_VELOCITY_VECTOR&>(bv);
    assert(v.pressure.m==pressure.m);
    solid_velocity.Copy(c,v.solid_velocity);for(int i=1;i<=v.pressure.m;i++) VECTOR_ND<T>::Copy(c,v.pressure(i),pressure(i));}

    void Copy(const T c1,const BASE& bv1,const BASE& bv2) PHYSBAM_OVERRIDE
    {const PRESSURE_VELOCITY_VECTOR& v1=debug_cast<const PRESSURE_VELOCITY_VECTOR&>(bv1),&v2=debug_cast<const PRESSURE_VELOCITY_VECTOR&>(bv2);
    assert(v1.pressure.m==v2.pressure.m && pressure.m==v1.pressure.m);
    solid_velocity.Copy(c1,v1.solid_velocity,v2.solid_velocity);
    for(int i=1;i<=v1.pressure.m;i++) VECTOR_ND<T>::Copy(c1,v1.pressure(i),v2.pressure(i),pressure(i));}

    int Raw_Size() const PHYSBAM_OVERRIDE
    {int n=solid_velocity.Raw_Size();for(int i=1;i<=pressure.m;i++) n+=pressure(i).n;return n;}

    T& Raw_Get(int i) PHYSBAM_OVERRIDE
    {int n=solid_velocity.Raw_Size();if(i<=n) return solid_velocity.Raw_Get(i);i-=n;
    for(int j=1;j<=pressure.m;j++){if(i<=pressure(j).n) return pressure(j)(i);i-=pressure(j).n;}
    PHYSBAM_FATAL_ERROR();}
};
}
#endif
