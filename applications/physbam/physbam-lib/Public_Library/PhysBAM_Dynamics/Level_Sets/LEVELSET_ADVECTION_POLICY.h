//#####################################################################
// Copyright 2009, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_ADVECTION_POLICY 
//#####################################################################
#ifndef __LEVELSET_ADVECTION_POLICY__
#define __LEVELSET_ADVECTION_POLICY__


namespace PhysBAM{

template<class T_GRID> class LEVELSET_ADVECTION_POLICY;

template<class T_GRID> class LEVELSET_ADVECTION_1D;
template<class T_GRID> class LEVELSET_ADVECTION_2D;
template<class T_GRID> class LEVELSET_ADVECTION_3D;
template<class T_GRID> class LEVELSET_ADVECTION_RLE;
template<class T_GRID> class LEVELSET_ADVECTION_DYADIC;
template<class T_GRID> class LEVELSET_ADVECTION_MULTIPLE_UNIFORM;
template<class T_GRID> class FAST_LEVELSET_ADVECTION;
template<class TV> class GRID;
template<class T> class QUADTREE_GRID;
template<class T> class OCTREE_GRID;
template<class T> class RLE_GRID_1D;
template<class T> class RLE_GRID_2D;
template<class T,int d> class VECTOR;

template<class T> struct LEVELSET_ADVECTION_POLICY<GRID<VECTOR<T,1> > >
{
    typedef LEVELSET_ADVECTION_1D<T> LEVELSET_ADVECTION;
    typedef FAST_LEVELSET_ADVECTION<GRID<VECTOR<T,1> > > FAST_LEVELSET_ADVECTION_T;
    typedef LEVELSET_ADVECTION_MULTIPLE_UNIFORM<GRID<VECTOR<T,1> > > LEVELSET_ADVECTION_MULTIPLE;
};

template<class T> struct LEVELSET_ADVECTION_POLICY<GRID<VECTOR<T,2> > >
{
    typedef LEVELSET_ADVECTION_2D<T> LEVELSET_ADVECTION;
    typedef FAST_LEVELSET_ADVECTION<GRID<VECTOR<T,2> > > FAST_LEVELSET_ADVECTION_T;
    typedef LEVELSET_ADVECTION_MULTIPLE_UNIFORM<GRID<VECTOR<T,2> > > LEVELSET_ADVECTION_MULTIPLE;
};

template<class T> struct LEVELSET_ADVECTION_POLICY<GRID<VECTOR<T,3> > >
{
    typedef LEVELSET_ADVECTION_3D<T> LEVELSET_ADVECTION;
    typedef FAST_LEVELSET_ADVECTION<GRID<VECTOR<T,3> > > FAST_LEVELSET_ADVECTION_T;
    typedef LEVELSET_ADVECTION_MULTIPLE_UNIFORM<GRID<VECTOR<T,3> > > LEVELSET_ADVECTION_MULTIPLE;
};

template<class T> struct LEVELSET_ADVECTION_POLICY<QUADTREE_GRID<T> >
{
    typedef LEVELSET_ADVECTION_DYADIC<QUADTREE_GRID<T> > LEVELSET_ADVECTION;
    typedef FAST_LEVELSET_ADVECTION<QUADTREE_GRID<T> > FAST_LEVELSET_ADVECTION_T;
};

template<class T> struct LEVELSET_ADVECTION_POLICY<OCTREE_GRID<T> >
{
    typedef LEVELSET_ADVECTION_DYADIC<OCTREE_GRID<T> > LEVELSET_ADVECTION;
    typedef FAST_LEVELSET_ADVECTION<OCTREE_GRID<T> > FAST_LEVELSET_ADVECTION_T;
};

template<class T> struct LEVELSET_ADVECTION_POLICY<RLE_GRID_1D<T> >
{
    typedef LEVELSET_ADVECTION_RLE<RLE_GRID_1D<T> > LEVELSET_ADVECTION;
    typedef FAST_LEVELSET_ADVECTION<RLE_GRID_1D<T> > FAST_LEVELSET_ADVECTION_T;
};

template<class T> struct LEVELSET_ADVECTION_POLICY<RLE_GRID_2D<T> >
{
    typedef LEVELSET_ADVECTION_RLE<RLE_GRID_2D<T> > LEVELSET_ADVECTION;
    typedef FAST_LEVELSET_ADVECTION<RLE_GRID_2D<T> > FAST_LEVELSET_ADVECTION_T;
};

}
#endif
