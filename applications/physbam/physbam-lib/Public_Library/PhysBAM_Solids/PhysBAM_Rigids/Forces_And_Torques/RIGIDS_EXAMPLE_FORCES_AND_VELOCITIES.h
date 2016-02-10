//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES
//#####################################################################
#ifndef __RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES__
#define __RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Geometry/Solids_Geometry_Evolution/RIGID_GEOMETRY_EXAMPLE_VELOCITIES.h>
namespace PhysBAM{

template<class TV> class TWIST;
template<class TV> class ROTATION;
template<class TV> class FRAME;
template<class TV>
class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES:public RIGID_GEOMETRY_EXAMPLE_VELOCITIES<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
public:

    RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES()
    {}

    virtual ~RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES();

//#####################################################################
    virtual void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time);
    virtual void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time); // or zero out components of their velocities
    virtual void Set_Rigid_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated);
//#####################################################################
};
}
#endif
