//#####################################################################
// Copyright 2005-2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __VORONOI_INTERPOLATION_CHIMERA__
#define __VORONOI_INTERPOLATION_CHIMERA__

#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>

namespace PhysBAM{

template<class T_GRID> class LAPLACE_CHIMERA_GRID_MPI;

template<class T_GRID>
class VORONOI_INTERPOLATION_CHIMERA
{
    typedef LAPLACE_CHIMERA_GRID_MPI<T_GRID> T_LAPLACE_GRID;
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef FACE_INDEX<TV::dimension> D_FACE_INDEX;
    typedef PAIR<int,TV_INT> GRID_CELL_INDEX;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension>::SIMPLEX T_SIMPLEX;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_FACE_ARRAYS_VECTOR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef MATRIX<T,TV::dimension,TV::dimension> D_MATRIX;
public:
    T_LAPLACE_GRID& laplace_grid;
    
    ARRAY<TV> voronoi_face_centroids;
    ARRAY<TV> voronoi_face_centroid_values;
    ARRAY<ARRAY<TV> > voronoi_face_vertex_values;
    
    ARRAY<ARRAY<TV> > cell_centroids;
    ARRAY<ARRAY<TV> > cell_centroid_values;
    ARRAY<HASHTABLE<D_FACE_INDEX,TV> > face_centroid_values;
    ARRAY<HASHTABLE<TV_INT,TV> > node_values;

    //ARRAY<HASHTABLE<D_FACE_INDEX,ARRAY<TV> >;
    const ARRAY<T_FACE_ARRAYS_SCALAR*>& face_values;
    const ARRAY<T>& voronoi_face_values;
    
    VORONOI_INTERPOLATION_CHIMERA(T_LAPLACE_GRID& laplace_grid_input,const ARRAY<T_FACE_ARRAYS_SCALAR*>& face_values_input,const ARRAY<T>& voronoi_face_values_input);
    ~VORONOI_INTERPOLATION_CHIMERA() {}
    
    TV Compute_Least_Squares_Vector(const TV& x,const T search_radius,const ARRAY<T_FACE_ARRAYS_SCALAR*>& face_values,const ARRAY<T>& voronoi_face_values) const;
    GRID_CELL_INDEX Nearest_Cell(const TV& x) const;
    TV From_Voronoi_Faces(const TV& x) const;
    TV From_Voronoi_Cell_Faces(const TV& x) const;
//#####################################################################
};
}
#endif
