//#####################################################################
// Copyright 2003-2006, Ronald Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_EVOLUTION
//##################################################################### 
#ifndef __PARTICLE_LEVELSET_EVOLUTION__
#define __PARTICLE_LEVELSET_EVOLUTION__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <cassert>
namespace PhysBAM{

template<class T>
class PARTICLE_LEVELSET_EVOLUTION:public NONCOPYABLE
{
protected:
    bool use_particle_levelset;
public:
    T time;
    T cfl_number;
    bool use_reinitialization,use_fmm;
    int reseeding_frequency;
    bool use_frozen_velocity;
    int runge_kutta_order_particles,runge_kutta_order_levelset;
    bool track_mass;
    double initial_mass;

    PARTICLE_LEVELSET_EVOLUTION()
    {
        Set_Time();
        Set_CFL_Number();
        Use_Fast_Marching_Method();
        Use_Particle_Levelset(true);
        Set_Reseeding_Frequency();
        Use_Frozen_Velocity();
        Set_Runge_Kutta_Order_Levelset();Set_Runge_Kutta_Order_Particles();
    }

    virtual ~PARTICLE_LEVELSET_EVOLUTION()
    {}

    void Set_Time(const T time_input=0)
    {time=time_input;}

    virtual void Set_CFL_Number(const T cfl_number_input=.5)
    {cfl_number=cfl_number_input;}

    void Use_Reinitialization()
    {use_reinitialization=true;use_fmm=false;}

    void Use_Fast_Marching_Method()
    {use_fmm=true;use_reinitialization=false;}

    void Use_Particle_Levelset(const bool use_particle_levelset_input)
    {use_particle_levelset=use_particle_levelset_input;}

    void Set_Reseeding_Frequency(const int number_of_frames=20)
    {reseeding_frequency=number_of_frames;}

    void Use_Frozen_Velocity(const bool use_frozen_velocity_input=true)
    {use_frozen_velocity=use_frozen_velocity_input;}

    void Set_Runge_Kutta_Order_Levelset(const int order=1)
    {runge_kutta_order_levelset=order;assert(order >= 1 && order <= 3);}
   
    void Set_Runge_Kutta_Order_Particles(const int order=2)
    {runge_kutta_order_particles=order;assert(order >= 1 && order <= 3);}

//#####################################################################
};
}
#endif
