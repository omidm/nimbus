#ifndef __CODIMENSIONAL_FLUID_PARTICLES__
#define __CODIMENSIONAL_FLUID_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD.h>

namespace PhysBAM{

template<class TV>
class CODIMENSIONAL_FLUID_PARTICLES:public CLONEABLE<CODIMENSIONAL_FLUID_PARTICLES<TV>,POINT_CLOUD<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<CODIMENSIONAL_FLUID_PARTICLES<TV>,POINT_CLOUD<TV> > BASE;
public:
    using BASE::array_collection;
    CODIMENSIONAL_FLUID_PARTICLES();
    virtual ~CODIMENSIONAL_FLUID_PARTICLES();
};
}

#endif
