//#####################################################################
// Copyright 2011, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_GRID_PARTICLES
//#####################################################################
#ifndef __RIGID_GRID_PARTICLES__
#define __RIGID_GRID_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GRID_PARTICLES_FORWARD.h>

namespace PhysBAM{

template<class T_GRID> class RIGID_GRID;

template<class T_GRID>
class RIGID_GRID_PARTICLES:public CLONEABLE<RIGID_GRID_PARTICLES<T_GRID>,GEOMETRY_PARTICLES<typename T_GRID::VECTOR_T> >
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;typedef typename TV::SPIN T_SPIN;
    typedef CLONEABLE<RIGID_GRID_PARTICLES<T_GRID>,GEOMETRY_PARTICLES<TV> > BASE;
public:
    using BASE::array_collection;using BASE::X;using BASE::V;

    ARRAY_VIEW<RIGID_GRID<T_GRID>*> rigid_grid;
    ARRAY_VIEW<ROTATION<TV> > rotation;
    ARRAY_VIEW<T_SPIN> angular;

    RIGID_GRID_PARTICLES(ARRAY_COLLECTION* array_collection_input)
        :rigid_grid(0,0),rotation(0,0),angular(0,0)
    {delete array_collection;array_collection=array_collection_input;Initialize_Array_Collection();}

    RIGID_GRID_PARTICLES()
        :rigid_grid(0,0),rotation(0,0),angular(0,0)
    {Initialize_Array_Collection();}

    void Initialize_Array_Collection()
    {POINT_CLOUD<TV>::Initialize_Array_Collection();GEOMETRY_PARTICLES<TV>::Store_Velocity(true);
    array_collection->Add_Array(ATTRIBUTE_ID_RIGID_GRID,&rigid_grid);array_collection->Add_Array(ATTRIBUTE_ID_ROTATION,&rotation);
    array_collection->Add_Array(ATTRIBUTE_ID_ANGULAR,&angular);}
    
    ~RIGID_GRID_PARTICLES()
    {}

//#####################################################################
    void Resize(const int new_size);
    void Remove_Grid(const int p);
    void Clean_Memory();
    void Delete_All_Particles();
//#####################################################################
};
}
#endif
