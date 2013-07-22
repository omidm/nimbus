//#####################################################################
// Copyright 2005-2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __LINEAR_INTERPOLATION_CHIMERA__
#define __LINEAR_INTERPOLATION_CHIMERA__

#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>

namespace PhysBAM{

template<class T_GRID> class CELL_LOOKUP_CHIMERA;
template<class T_GRID> class LAPLACE_CHIMERA_GRID_MPI;

template<class T_GRID>
class LINEAR_INTERPOLATION_CHIMERA
{
    typedef LAPLACE_CHIMERA_GRID_MPI<T_GRID> T_LAPLACE_GRID;
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef FACE_INDEX<TV::dimension> D_FACE_INDEX;
    typedef PAIR<int,TV_INT> GRID_CELL_INDEX;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension>::SIMPLEX T_SIMPLEX;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
public:
    T_LAPLACE_GRID& laplace_grid;
    ARRAY<VECTOR<GRID_CELL_INDEX,TV::dimension+1> > simplex_indices;
    ARRAY<VECTOR<TV,TV::dimension+1> > simplex_vertices;
    BOX_HIERARCHY<TV> hierarchy;

    LINEAR_INTERPOLATION_CHIMERA(T_LAPLACE_GRID& laplace_grid_input);
    ~LINEAR_INTERPOLATION_CHIMERA() {}

    void Clean_Memory()
    {
        simplex_indices.Clean_Memory();
        simplex_vertices.Clean_Memory();
        hierarchy.Clean_Memory();
    }

    void Construct_Delaunay_Simplices();
    bool Cell_Intersect_Local_Domain(GRID_CELL_INDEX cell);
    T From_Cell_Centers(const CELL_LOOKUP_CHIMERA<T_GRID>& u,const TV& X) const;
//#####################################################################
};
}
#endif
