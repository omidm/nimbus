//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_REMOVED_PARTICLES
//#####################################################################
#ifndef __PARTICLE_LEVELSET_REMOVED_PARTICLES__
#define __PARTICLE_LEVELSET_REMOVED_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Point_Clouds_Computations/EULER_STEP.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
namespace PhysBAM{

template<class TV>
class PARTICLE_LEVELSET_REMOVED_PARTICLES:public CLONEABLE<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>,PARTICLE_LEVELSET_PARTICLES<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>,PARTICLE_LEVELSET_PARTICLES<TV> > BASE;
public:
    using BASE::X;using BASE::array_collection;

    ARRAY_VIEW<TV> V;
    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* next; //This should always be 0

    //PARTICLE_LEVELSET_REMOVED_PARTICLES(ARRAY_COLLECTION& array_collection_input);
    PARTICLE_LEVELSET_REMOVED_PARTICLES();
    ~PARTICLE_LEVELSET_REMOVED_PARTICLES();

    template<class T_INDICES>
    void Euler_Step(const T_INDICES& indices,const ARRAY<TV>& F,const ARRAY<T>& mass,const T dt)
    {POINT_CLOUDS_COMPUTATIONS::Euler_Step(indices,V,F,mass,dt);}

    template<class T_INDICES>
    void Euler_Step(const T_INDICES& indices,const T dt)
    {POINT_CLOUDS_COMPUTATIONS::Euler_Step(indices,X,V,dt);}

    void Euler_Step_Position(const T dt)
    {POINT_CLOUDS_COMPUTATIONS::Euler_Step(X,V,dt);}

//#####################################################################
};
}
#endif
