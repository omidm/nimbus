//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLES_EVOLUTION_CALLBACKS
//#####################################################################
#ifndef __DEFORMABLES_EVOLUTION_CALLBACKS__
#define __DEFORMABLES_EVOLUTION_CALLBACKS__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <cfloat>
namespace PhysBAM{

template<class TV>
class DEFORMABLES_EVOLUTION_CALLBACKS
{
    typedef typename TV::SCALAR T;
public:
    DEFORMABLES_EVOLUTION_CALLBACKS()
    {}

    virtual ~DEFORMABLES_EVOLUTION_CALLBACKS()
    {}

//#####################################################################
    virtual void Update_Deformables_Parameters(const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Preprocess_Solids_Substep(const T time,const int substep){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Postprocess_Solids_Substep(const T time,const int substep){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Self_Collisions_Begin_Callback(const T time,const int substep){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Apply_Constraints(const T dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Limit_Solids_Dt(T& dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual T Constraints_CFL(){return FLT_MAX;}
//#####################################################################
};
}
#endif
