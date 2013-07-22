//#####################################################################
// Copyright 2011, Mridul Aanjaneya, Linhai Qiu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_GRID_COLLECTION
//#####################################################################
#ifndef __RIGID_GRID_COLLECTION__
#define __RIGID_GRID_COLLECTION__

#include <PhysBAM_Geometry/Geometry_Particles/GRID_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Geometry_Particles/RIGID_GRID_PARTICLES.h>

namespace PhysBAM{
class STREAM_TYPE;

template<class T_GRID> class RIGID_GRID;
template<class T_GRID> struct ALLOCATE_HELPER_RIGID_GRID{virtual RIGID_GRID<T_GRID>* Create(int index=0)=0;virtual ~ALLOCATE_HELPER_RIGID_GRID(){}};

template<class T_GRID>
class RIGID_GRID_COLLECTION:public NONCOPYABLE
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SPIN T_SPIN;
public:
    RIGID_GRID_PARTICLES<T_GRID>& particles;
    ALLOCATE_HELPER_RIGID_GRID<T_GRID>* allocate_helper;
    ARRAY<int> static_rigid_grid,kinematic_rigid_grid;
    bool owns_particles;
   
    RIGID_GRID_COLLECTION(RIGID_GRID_PARTICLES<T_GRID>& particles_input,ALLOCATE_HELPER_RIGID_GRID<T_GRID>* allocate_helper_input=0);
    RIGID_GRID_COLLECTION(ALLOCATE_HELPER_RIGID_GRID<T_GRID>* allocate_helper_input=0);
    virtual ~RIGID_GRID_COLLECTION();

    RIGID_GRID<T_GRID>& Rigid_Grid(const int index)
    {return *particles.rigid_grid(index);}

    const RIGID_GRID<T_GRID>& Rigid_Grid(const int index) const
    {return *particles.rigid_grid(index);}

    void Deactivate_Grid(const int p)
    {assert(Exists(p) && Rigid_Grid(p).index>0);Rigid_Grid(p).index=-Rigid_Grid(p).index;}

    void Reactivate_Grid(const int p)
    {assert(Exists(p) && Rigid_Grid(p).index<0);Rigid_Grid(p).index=-Rigid_Grid(p).index;}

    RIGID_GRID<T_GRID>* New_Grid(int index);

//#####################################################################
    bool Exists(const int particle) const;
    bool Is_Active(const int particle) const;
    void Update_Kinematic_Particles();
    int Add_Rigid_Grid(TV X_input=TV(),ROTATION<TV> rotation_input=ROTATION<TV>(), TV V_input=TV(), T_SPIN angular_input=T_SPIN());
    int Add_Rigid_Grid(RIGID_GRID<T_GRID>* rigid_grid,TV X_input=TV(),ROTATION<TV> rotation_input=ROTATION<TV>(), TV V_input=TV(), T_SPIN angular_input=T_SPIN());
//#####################################################################
};
}
#endif
