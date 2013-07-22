//#####################################################################
// Copyright 2007, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRACTURE_CALLBACKS
//##################################################################### 
#ifndef __FRACTURE_CALLBACKS__
#define __FRACTURE_CALLBACKS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
namespace PhysBAM{



template<class TV>
class FRACTURE_CALLBACKS
{
    typedef typename TV::SCALAR T;

public:

    FRACTURE_CALLBACKS()
    {}

    virtual ~FRACTURE_CALLBACKS()
    {}


//#####################################################################
    virtual TV Perturb(const TV& position) const {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return TV();}
    virtual void Perturb(ARRAY_VIEW<const TV> positions,ARRAY_VIEW<TV> perturbations,const FRAME<TV>& frame) const {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Grain_Boundary_Weakness_Multiplier(ARRAY_VIEW<const TV> positions,ARRAY_VIEW<T> seed_weakness_multipliers) const {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Get_Breakability(ARRAY_VIEW<const TV> positions,ARRAY_VIEW<bool> is_breakable) const {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Node_Regions(ARRAY_VIEW<const TV> positions,ARRAY_VIEW<int> node_regions,const FRAME<TV>& frame) const {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual bool Has_Crack_Map() const {return false;};
//#####################################################################
};  
}
#endif
