#ifndef NIMBUS_APPLICATION_WATER_ALTERNATE_COARSE_WATER_SOURCES_H_
#define NIMBUS_APPLICATION_WATER_ALTERNARE_COARSE_WATER_SOURCES_H_

#include "application/water_alternate_coarse/app_utils.h"
#include "application/water_alternate_coarse/water_example.h"
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM {

    typedef application::T T;

    void Add_Source(WATER_EXAMPLE<VECTOR<T,1> >* example);
    void Add_Source(WATER_EXAMPLE<VECTOR<T,2> >* example);
    void Add_Source(WATER_EXAMPLE<VECTOR<T,3> >* example);

} // namespace PhysBAM

#endif  // NIMBUS_APPLICATION_WATER_ALTERNATE_COARSE_WATER_SOURCES_H_
