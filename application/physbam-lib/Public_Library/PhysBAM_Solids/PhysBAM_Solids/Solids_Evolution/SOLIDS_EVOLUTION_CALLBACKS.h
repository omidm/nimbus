//#####################################################################
// Copyright 2004-2008, Ronald Fedkiw, Geoffrey Irving, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_EVOLUTION_CALLBACKS
//#####################################################################
#ifndef __SOLIDS_EVOLUTION_CALLBACKS__
#define __SOLIDS_EVOLUTION_CALLBACKS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <cfloat>
namespace PhysBAM{

class RIGID_BODY_COLLISION_MANAGER;

template<class TV>
class SOLIDS_EVOLUTION_CALLBACKS
{
    typedef typename TV::SCALAR T;
public:
    SOLIDS_EVOLUTION_CALLBACKS()
    {}

    virtual ~SOLIDS_EVOLUTION_CALLBACKS()
    {}

//#####################################################################
    virtual void Update_Solids_Parameters(const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Self_Collisions_Begin_Callback(const T time,const int substep){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Preprocess_Solids_Substep(const T time,const int substep){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Postprocess_Solids_Substep(const T time,const int substep){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Align_Deformable_Bodies_With_Rigid_Bodies(){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Apply_Constraints(const T dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual T Constraints_CFL(){return FLT_MAX;}
    virtual void Limit_Solids_Dt(T& dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Pre_Advance_Cluster_Fracture(const T& dt, const T& time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Post_Advance_Cluster_Fracture(const T& dt, const T& time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual void Filter_Velocities(const T dt,const T time,const bool velocity_update){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual bool Get_Solid_Source_Velocities(ARRAY<int>& deformable_simplices,ARRAY<T>& deformable_simplex_forces,ARRAY<PAIR<int,int> >& rigid_simplices,ARRAY<T>& rigid_simplex_forces,TV& orientation,const T time)
    {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return false;}
//#####################################################################
};
}
#endif
