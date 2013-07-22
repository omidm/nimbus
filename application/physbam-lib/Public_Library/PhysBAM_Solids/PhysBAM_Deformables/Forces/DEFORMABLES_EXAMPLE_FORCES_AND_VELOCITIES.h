//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES
//#####################################################################
#ifndef __DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES__
#define __DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
namespace PhysBAM{

template<class TV>
class DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES
{
    typedef typename TV::SCALAR T;
public:

    DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES()
    {}

    virtual ~DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES();

//#####################################################################
    virtual void Add_External_Forces(ARRAY_VIEW<TV> F,const T time);
    virtual void Set_External_Positions(ARRAY_VIEW<TV> X,const T time); // set external positions
    virtual void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time);
    virtual void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time); // zero out entries corresponding to external positions
    virtual void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time); // or zero out components of their velocities
    virtual void Set_Deformable_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated);
//#####################################################################
};
}
#endif
