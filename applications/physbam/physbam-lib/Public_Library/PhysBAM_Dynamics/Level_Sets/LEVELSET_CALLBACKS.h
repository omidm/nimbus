//#####################################################################
// Copyright 2003-2005, Ronald Fedkiw, Geoffrey Irving, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_CALLBACKS
//##################################################################### 
#ifndef __LEVELSET_CALLBACKS__
#define __LEVELSET_CALLBACKS__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <PhysBAM_Dynamics/Particles/PARTICLES_FORWARD.h>
namespace PhysBAM{

template<class T_GRID> struct GRID_ARRAYS_POLICY;

template<class T_GRID>
class LEVELSET_CALLBACKS
{   
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET_MULTIPLE T_LEVELSET_MULTIPLE;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
public:
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename TV::SCALAR T;
    LEVELSET_CALLBACKS()
    {}

    virtual ~LEVELSET_CALLBACKS()
    {}

//#####################################################################
    virtual void Get_Levelset_Velocity(const T_GRID& grid,T_LEVELSET& levelset,T_FACE_ARRAYS_SCALAR& face_velocity,const T time=0) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Get_Levelset_Velocity(const T_GRID& grid,T_LEVELSET_MULTIPLE& levelset,T_FACE_ARRAYS_SCALAR& face_velocity,const T time=0) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,
        const T dt,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual bool Adjust_Particle_For_Objects(TV& X,TV& V,const T r, const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time)
        {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return true;} // return false if particle should be deleted
//#####################################################################
};
}
#endif
