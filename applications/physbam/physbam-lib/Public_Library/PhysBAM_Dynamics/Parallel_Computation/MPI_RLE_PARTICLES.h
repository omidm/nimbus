#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPI_RLE_PARTICLES
//#####################################################################
#ifndef __MPI_RLE_PARTICLES__
#define __MPI_RLE_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_RLE_GRID.h>
namespace PhysBAM{

template<class TV> class PARTICLE_LEVELSET_PARTICLES;
template<class TV> class PARTICLE_LEVELSET_REMOVED_PARTICLES;
template<class T_GRID> class PARTICLE_LEVELSET_RLE;

template<class T_GRID,class T_PARTICLES> void Exchange_Boundary_Particles(MPI_RLE_GRID<T_GRID>& mpi_grid,T_PARTICLES& particles,const typename T_GRID::SCALAR ghost_distance=0);
template<class T_GRID,class T_PARTICLES> void Exchange_Boundary_Particles(MPI_RLE_GRID<T_GRID>& mpi_grid,PARTICLE_LEVELSET_RLE<T_GRID>& particle_levelset,ARRAY<T_PARTICLES*>& particles,
    T_PARTICLES* particles_in_long_cells,const int bandwidth);

}
#endif
#endif
