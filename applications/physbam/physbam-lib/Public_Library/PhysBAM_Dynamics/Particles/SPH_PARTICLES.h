//#####################################################################
// Copyright 2004-2009, Michael Lentine, Frank Losasso, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPH_PARTICLES
//#####################################################################
#ifndef __SPH_PARTICLES__
#define __SPH_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD.h>
#include <PhysBAM_Tools/Point_Clouds_Computations/EULER_STEP.h>
namespace PhysBAM{

template<class TV>
class SPH_PARTICLES:public CLONEABLE<SPH_PARTICLES<TV>,POINT_CLOUD<TV> >
{
    typedef CLONEABLE<SPH_PARTICLES<TV>,POINT_CLOUD<TV> > BASE;
    typedef typename TV::SCALAR T;
public:
    using BASE::array_collection;

    ARRAY_VIEW<TV> V;

    SPH_PARTICLES();
    virtual ~SPH_PARTICLES();

    template<class T_INDICES>
    void Euler_Step(const T_INDICES& indices,const ARRAY<TV>& X,const T dt)
    {POINT_CLOUDS_COMPUTATIONS::Euler_Step(indices,X,V.array,dt);}

    template<class T_INDICES>
    void Euler_Step(const T_INDICES& indices,const ARRAY<TV>& F,const ARRAY<T>& mass,const T dt)
    {POINT_CLOUDS_COMPUTATIONS::Euler_Step(indices,V.array,F,mass,dt);}
};
}
#endif
