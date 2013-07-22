//#####################################################################
// Copyright 2011,Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_PARTICLES
//#####################################################################
#ifndef __FLUID_PARTICLES__
#define __FLUID_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD.h>
#include <PhysBAM_Tools/Point_Clouds_Computations/EULER_STEP.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
namespace PhysBAM{

template<class TV>
class FLUID_PARTICLES:public PARTICLES<TV>
{
    typedef PARTICLES<TV> BASE;
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    using BASE::array_collection;using BASE::id;using BASE::X;using BASE::V;using BASE::Store_Id;using BASE::Store_Velocity;

    ARRAY_VIEW<TV> F,momentum_residual;
    ARRAY_VIEW<VECTOR<T,3> > color;
    ARRAY_VIEW<T> volume,volume_residual;

    ARRAY_VIEW<int> on_boundary,type;
    ARRAY_VIEW<TV> normal;
    ARRAY_VIEW<TV_INT> node;
    ARRAY_VIEW<T> age,phi;

    FLUID_PARTICLES();
    virtual ~FLUID_PARTICLES();

    int Add_New_Particle()
    {
        int p=array_collection->Add_Element();
        V(p)=TV();
        color(p)=VECTOR<T,3>((T)0,(T)0,(T)1);
        type(p)=0;
        node(p)=TV_INT();
        on_boundary(p)=0;
        age(p)=0;
        return p;
    }

    int Size()
    {return array_collection->Size();}

    void Euler_Step_Position(const T dt)
    {X+=V*dt;}

    void Euler_Step_Velocity(const T dt)
    {V+=F*dt;}
};
}
#endif
