//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Nipun Kwatra, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_PARTICLES
//#####################################################################
#ifndef __PARTICLE_LEVELSET_PARTICLES__
#define __PARTICLE_LEVELSET_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD.h>
namespace PhysBAM{

template<class TV>
class PARTICLE_LEVELSET_PARTICLES:public CLONEABLE<PARTICLE_LEVELSET_PARTICLES<TV>,POINT_CLOUD<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<PARTICLE_LEVELSET_PARTICLES<TV>,POINT_CLOUD<TV> > BASE;
public:
    using BASE::array_collection;

    ARRAY_VIEW<unsigned short> quantized_collision_distance;
    ARRAY_VIEW<T> age;
    ARRAY_VIEW<T> radius;

    PARTICLE_LEVELSET_PARTICLES<TV>* next;

    //PARTICLE_LEVELSET_PARTICLES(ARRAY_COLLECTION& array_collection_input);
    PARTICLE_LEVELSET_PARTICLES();
    virtual ~PARTICLE_LEVELSET_PARTICLES();

    int Number() const
    {return next?array_collection->Size()+next->Number():array_collection->Size();}

//#####################################################################
};
}
#endif
