//#####################################################################
// Copyright 2011, Mridul Aanjaneya, Linhai Qiu. 
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_GRID_COLLECTION
//#####################################################################
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GRID_COLLECTION.h>
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#endif
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GRID.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
namespace PhysBAM{
template<class T_GRID>
struct ALLOCATE_GRID_HELPER:public ALLOCATE_HELPER_RIGID_GRID<T_GRID>
{
    RIGID_GRID_COLLECTION<T_GRID>& collection;
    ALLOCATE_GRID_HELPER(RIGID_GRID_COLLECTION<T_GRID>& collection_input): collection(collection_input) {}
    RIGID_GRID<T_GRID>* Create(int index=0) PHYSBAM_OVERRIDE {return new RIGID_GRID<T_GRID>(collection,index);}
    virtual ~ALLOCATE_GRID_HELPER(){}
};
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> RIGID_GRID_COLLECTION<T_GRID>::
RIGID_GRID_COLLECTION(RIGID_GRID_PARTICLES<T_GRID>& particles_input,ALLOCATE_HELPER_RIGID_GRID<T_GRID>* allocate_helper_input)
    :particles(particles_input),allocate_helper(allocate_helper_input),owns_particles(false)
{
    if(!allocate_helper) allocate_helper=new ALLOCATE_GRID_HELPER<T_GRID>(*this);
}
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> RIGID_GRID_COLLECTION<T_GRID>::
RIGID_GRID_COLLECTION(ALLOCATE_HELPER_RIGID_GRID<T_GRID>* allocate_helper_input)
    :particles(*new RIGID_GRID_PARTICLES<T_GRID>()),allocate_helper(allocate_helper_input),owns_particles(true)
{
    if(!allocate_helper) allocate_helper=new ALLOCATE_GRID_HELPER<T_GRID>(*this);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> RIGID_GRID_COLLECTION<T_GRID>::
~RIGID_GRID_COLLECTION()
{
    delete allocate_helper;
}
//#####################################################################
// Function Exists
//#####################################################################
template<class T_GRID> bool RIGID_GRID_COLLECTION<T_GRID>::
Exists(const int particle) const
{
    return particle>0 && particle<=particles.array_collection->Size() && particles.rigid_grid(particle);
}
//#####################################################################
// Function Exists
//#####################################################################
template<class T_GRID> bool RIGID_GRID_COLLECTION<T_GRID>::
Is_Active(const int particle) const
{
    return Exists(particle) && particles.rigid_grid(particle)->index>0;
}
//#####################################################################
// Function Update_Kinematic_Particles
//#####################################################################
template<class T_GRID> void RIGID_GRID_COLLECTION<T_GRID>::
Update_Kinematic_Particles()
{
    static_rigid_grid.Remove_All();kinematic_rigid_grid.Remove_All();
    for(int p=1;p<=particles.array_collection->Size();p++) if(Is_Active(p)){RIGID_GRID<T_GRID>& rigid_grid=Rigid_Grid(p);
        if(rigid_grid.is_static) static_rigid_grid.Append(p); else kinematic_rigid_grid.Append(p);}
}
//#####################################################################
// Function Add_Rigid_Grid
//#####################################################################
template<class T_GRID> int RIGID_GRID_COLLECTION<T_GRID>::
Add_Rigid_Grid(TV X_input,ROTATION<TV> rotation_input, TV V_input, T_SPIN angular_input)
{
    RIGID_GRID<T_GRID>* rigid_grid=new RIGID_GRID<T_GRID>(*this);
    return Add_Rigid_Grid(rigid_grid,X_input,rotation_input,V_input,angular_input);
}
//#####################################################################
// Function Add_Rigid_Grid
//#####################################################################
template<class T_GRID> int RIGID_GRID_COLLECTION<T_GRID>::
Add_Rigid_Grid(RIGID_GRID<T_GRID>* rigid_grid,TV X_input,ROTATION<TV> rotation_input, TV V_input, T_SPIN angular_input)
{
    particles.rotation(rigid_grid->index)=rotation_input;
    particles.V(rigid_grid->index)=V_input;
    particles.angular(rigid_grid->index)=angular_input;
    particles.X(rigid_grid->index)=X_input;
    return rigid_grid->index;
}
//#####################################################################
// Function New_Grid
//#####################################################################
template<class T_GRID> RIGID_GRID<T_GRID>* RIGID_GRID_COLLECTION<T_GRID>::
New_Grid(int index)
{
    return allocate_helper->Create(index);
}
//#####################################################################
template class RIGID_GRID_COLLECTION<GRID<VECTOR<float,1> > >;
template class RIGID_GRID_COLLECTION<GRID<VECTOR<float,2> > >;
template class RIGID_GRID_COLLECTION<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_GRID_COLLECTION<GRID<VECTOR<double,1> > >;
template class RIGID_GRID_COLLECTION<GRID<VECTOR<double,2> > >;
template class RIGID_GRID_COLLECTION<GRID<VECTOR<double,3> > >;
#endif
}
