//#####################################################################
// Copyright 2011, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Geometry_Particles/RIGID_GRID_PARTICLES.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GRID.h>
namespace PhysBAM{
//#####################################################################
// Resize
//#####################################################################
template<class T_GRID> void RIGID_GRID_PARTICLES<T_GRID>::
Resize(const int new_size)
{
    for(int p=new_size+1;p<=array_collection->Size();p++)
        if(rigid_grid(p)) Remove_Grid(p);
    array_collection->Resize(new_size);
}
//#####################################################################
// Remove_Grid
//#####################################################################
template<class T_GRID> void RIGID_GRID_PARTICLES<T_GRID>::
Remove_Grid(const int p)
{
    delete rigid_grid(p);
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class T_GRID> void RIGID_GRID_PARTICLES<T_GRID>::
Clean_Memory()
{
    for(int p=1;p<=array_collection->Size();p++) if(rigid_grid(p)) Remove_Grid(p);
    array_collection->Clean_Memory();
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class T_GRID> void RIGID_GRID_PARTICLES<T_GRID>::
Delete_All_Particles()
{
    for(int p=1;p<=array_collection->Size();p++) if(rigid_grid(p)) Remove_Grid(p);
    array_collection->Delete_All_Elements();
}
template class RIGID_GRID_PARTICLES<GRID<VECTOR<float,1> > >;
template class RIGID_GRID_PARTICLES<GRID<VECTOR<float,2> > >;
template class RIGID_GRID_PARTICLES<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_GRID_PARTICLES<GRID<VECTOR<double,1> > >;
template class RIGID_GRID_PARTICLES<GRID<VECTOR<double,2> > >;
template class RIGID_GRID_PARTICLES<GRID<VECTOR<double,3> > >;
#endif
}
