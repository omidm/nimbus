//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_FLUIDS_CALLBACKS
//#####################################################################
#ifndef __SOLIDS_FLUIDS_CALLBACKS__
#define __SOLIDS_FLUIDS_CALLBACKS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>

namespace PhysBAM{
template<class TV> class TWIST;
template<class TV> class RIGIDS_NEWMARK_COLLISION_CALLBACKS;

template <class TV>
class SOLIDS_FLUIDS_CALLBACKS
{
    typedef typename TV::SCALAR T;
public:
    SOLIDS_FLUIDS_CALLBACKS()
    {}

    virtual ~SOLIDS_FLUIDS_CALLBACKS()
    {}

//#####################################################################
    virtual void Save_Fluid_Pressure_On_Solids(const ARRAY<TV>& F,const ARRAY<TWIST<TV> >& F_twist) {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Solid_Fluid_Fracture_Callback(const T dt, const T time,RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>& rigids_evolution_callbacks) {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
//#####################################################################
};
}
#endif
