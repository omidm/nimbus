//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXAMPLE_FORCES_AND_VELOCITIES
//#####################################################################
#ifndef __EXAMPLE_FORCES_AND_VELOCITIES__
#define __EXAMPLE_FORCES_AND_VELOCITIES__

#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES.h>
namespace PhysBAM{

template<class TV>
class EXAMPLE_FORCES_AND_VELOCITIES:public RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>,public DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES<TV> DEFORMABLES_BASE;
    typedef RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV> RIGIDS_BASE;
public:
    using DEFORMABLES_BASE::Add_External_Forces;using DEFORMABLES_BASE::Set_External_Positions;using DEFORMABLES_BASE::Set_External_Velocities;
    using DEFORMABLES_BASE::Zero_Out_Enslaved_Position_Nodes;using DEFORMABLES_BASE::Zero_Out_Enslaved_Velocity_Nodes;using DEFORMABLES_BASE::Set_Deformable_Particle_Is_Simulated;
    using RIGIDS_BASE::Add_External_Forces;using RIGIDS_BASE::Set_External_Positions;using RIGIDS_BASE::Set_External_Velocities;using RIGIDS_BASE::Set_Kinematic_Positions;
    using RIGIDS_BASE::Set_Kinematic_Velocities;using RIGIDS_BASE::Zero_Out_Enslaved_Velocity_Nodes;using RIGIDS_BASE::Set_Rigid_Particle_Is_Simulated;

    EXAMPLE_FORCES_AND_VELOCITIES()
    {}

    virtual ~EXAMPLE_FORCES_AND_VELOCITIES();

//#####################################################################
    // solids
    virtual void Advance_One_Time_Step_Begin_Callback(const T dt,const T time);
    virtual void Advance_One_Time_Step_End_Callback(const T dt,const T time);
    virtual void Set_PD_Targets(const T dt,const T time);
    virtual void Update_Time_Varying_Material_Properties(const T time);
    virtual void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt); // adjust velocity with external impulses
    virtual void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt); // adjust velocity with external impulses
    virtual void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt); // adjust velocity of single node with external impulse
//#####################################################################
};
}
#endif
