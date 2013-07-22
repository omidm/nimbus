//#####################################################################
// Copyright 2009, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INCOMPRESSIBLE_POLICY__
#define __INCOMPRESSIBLE_POLICY__

#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_FORWARD.h>

namespace PhysBAM{

template<class TV> class GRID;
template<class T> class QUADTREE_GRID;
template<class T> class OCTREE_GRID;
template<class T> class RLE_GRID_2D;
template<class T> class RLE_GRID_3D;

template<class T_GRID> struct INCOMPRESSIBLE_POLICY;

template<class TV> struct INCOMPRESSIBLE_POLICY<GRID<TV> >
{
    typedef INCOMPRESSIBLE_UNIFORM<GRID<TV> > INCOMPRESSIBLE;
    typedef PROJECTION_DYNAMICS_UNIFORM<GRID<TV> > PROJECTION;
};

template<class T> struct INCOMPRESSIBLE_POLICY<QUADTREE_GRID<T> >
{
    typedef INCOMPRESSIBLE_QUADTREE<T> INCOMPRESSIBLE;
    typedef PROJECTION_DYADIC<QUADTREE_GRID<T> > PROJECTION;
};

template<class T> struct INCOMPRESSIBLE_POLICY<OCTREE_GRID<T> >
{
    typedef INCOMPRESSIBLE_OCTREE<T> INCOMPRESSIBLE;
    typedef PROJECTION_DYADIC<OCTREE_GRID<T> > PROJECTION;
};

template<class T> struct INCOMPRESSIBLE_POLICY<RLE_GRID_2D<T> >
{
    typedef INCOMPRESSIBLE_RLE<RLE_GRID_2D<T> > INCOMPRESSIBLE;
    typedef PROJECTION_RLE<RLE_GRID_2D<T> > PROJECTION;
};

template<class T> struct INCOMPRESSIBLE_POLICY<RLE_GRID_3D<T> >
{
    typedef INCOMPRESSIBLE_RLE<RLE_GRID_3D<T> > INCOMPRESSIBLE;
    typedef PROJECTION_RLE<RLE_GRID_3D<T> > PROJECTION;
};

//#####################################################################
}
#endif
