//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __GEOMETRY_BOUNDARY_POLICY__
#define __GEOMETRY_BOUNDARY_POLICY__

#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_RLE_Boundaries/BOUNDARY_POLICY_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_POLICY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/BOUNDARY_FORWARD.h>
namespace PhysBAM{

template<class TV> class GRID;
template<class T> class QUADTREE_GRID;
template<class T> class OCTREE_GRID;
template<class T> class RLE_GRID_2D;
template<class T> class RLE_GRID_3D;

template<class T_GRID>
struct GEOMETRY_BOUNDARY_POLICY
{
    typedef PhysBAM::BOUNDARY_REFLECTION_UNIFORM<T_GRID,typename T_GRID::SCALAR> BOUNDARY_REFLECTION;
    typedef PhysBAM::BOUNDARY_PHI_WATER<T_GRID> BOUNDARY_PHI_WATER; 
    typedef PhysBAM::BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<T_GRID> BOUNDARY_MAC_GRID_SOLID_WALL_SLIP;
};

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template<class T>
struct GEOMETRY_BOUNDARY_POLICY<OCTREE_GRID<T> >
{
private:
    typedef OCTREE_GRID<T> T_GRID;
public:
    typedef BOUNDARY_DYADIC<T_GRID,T> BOUNDARY_REFLECTION;
    typedef PhysBAM::BOUNDARY_PHI_WATER<T_GRID> BOUNDARY_PHI_WATER;
    typedef BOUNDARY_DYADIC<T_GRID,T> BOUNDARY_MAC_GRID_SOLID_WALL_SLIP; // TODO: write dyadic slip
};

template<class T>
struct GEOMETRY_BOUNDARY_POLICY<QUADTREE_GRID<T> >
{
private:
    typedef QUADTREE_GRID<T> T_GRID;
public:
    typedef BOUNDARY_DYADIC<T_GRID,T> BOUNDARY_REFLECTION;
    typedef PhysBAM::BOUNDARY_PHI_WATER<T_GRID> BOUNDARY_PHI_WATER;
    typedef BOUNDARY_DYADIC<T_GRID,T> BOUNDARY_MAC_GRID_SOLID_WALL_SLIP; // TODO: write dyadic slip
};
#endif

#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template<class T>
struct GEOMETRY_BOUNDARY_POLICY<RLE_GRID_2D<T> >
{
private:
    typedef RLE_GRID_2D<T> T_GRID;
public:
    typedef BOUNDARY_RLE<T_GRID,T> BOUNDARY_REFLECTION; // TODO: write rle reflection boundaries
    typedef BOUNDARY_RLE<T_GRID,T> BOUNDARY_PHI_WATER; // TODO: write rle phi water boundaries 
    typedef BOUNDARY_RLE<T_GRID,T> BOUNDARY_MAC_GRID_SOLID_WALL_SLIP; // TODO: write rle slip boundary conditions
};

template<class T>
struct GEOMETRY_BOUNDARY_POLICY<RLE_GRID_3D<T> >
{
private:
    typedef RLE_GRID_3D<T> T_GRID;
public:
    typedef BOUNDARY_RLE<T_GRID,T> BOUNDARY_REFLECTION; // TODO: write rle reflection boundaries
    typedef BOUNDARY_RLE<T_GRID,T> BOUNDARY_PHI_WATER; // TODO: write rle phi water boundaries 
    typedef BOUNDARY_RLE<T_GRID,T> BOUNDARY_MAC_GRID_SOLID_WALL_SLIP; // TODO: write rle slip boundary conditions
};
#endif
}

#endif
