//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VORTICITY_PARTICLES
//#####################################################################
#ifndef __VORTICITY_PARTICLES__
#define __VORTICITY_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD.h>
namespace PhysBAM{

template<class TV>
class VORTICITY_PARTICLES:public CLONEABLE<VORTICITY_PARTICLES<TV>,POINT_CLOUD<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<VORTICITY_PARTICLES<TV>,POINT_CLOUD<TV> > BASE;
public:
    using BASE::array_collection;

    ARRAY_VIEW<typename TV::SPIN> vorticity;
    ARRAY_VIEW<T> radius;

    //VORTICITY_PARTICLES(ARRAY_COLLECTION* array_collection_input);
    VORTICITY_PARTICLES();
    virtual ~VORTICITY_PARTICLES();
};
}
#endif
