//#####################################################################
// Copyright 2006, Nipun Kwatra, Frank Losasso, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPH_CALLBACKS
//##################################################################### 
#ifndef __SPH_CALLBACKS__
#define __SPH_CALLBACKS__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLES_FORWARD.h>
namespace PhysBAM{

template<class T_GRID> struct GRID_ARRAYS_POLICY;

template<class T_GRID>
class SPH_CALLBACKS
{    
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS FACE_ARRAYS;
    typedef typename REBIND<FACE_ARRAYS,bool>::TYPE FACE_ARRAYS_BOOL;
public:

    SPH_CALLBACKS()
    {}

    virtual ~SPH_CALLBACKS()
    {}

//#####################################################################
    virtual void Adjust_SPH_Particle_For_Domain_Boundaries(SPH_PARTICLES<TV>& particles,const int index,TV& V,const T dt,const T time)const{PHYSBAM_WARN_IF_NOT_OVERRIDDEN();}
    virtual bool Adjust_SPH_Particle_For_Objects(SPH_PARTICLES<TV>& particles,const int index,TV& V,const T dt,const T time)const{PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return true;} // return false if particle should be deleted
    virtual void Do_Something_With_Density(const T_GRID &grid,const typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR &cell_weight)const{}
    virtual T Target_Density_Factor(const TV& location,const T time)const{PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return 1;}
//#####################################################################
};
}
#endif
