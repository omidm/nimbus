#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic_Computations/VORTICITY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_OCTREE_HELPER.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/EXTRAPOLATION_DYADIC.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_OCTREE.h>
using namespace PhysBAM;
//#####################################################################
// Function Compute_Vorticity_Confinement_Force
//#####################################################################
template<class T> void INCOMPRESSIBLE_OCTREE<T>::
Compute_Vorticity_Confinement_Force(OCTREE_GRID<T>& grid,const ARRAY<VECTOR<T,3> >& V_ghost,ARRAY<VECTOR<T,3> >& F)
{
    int i;ARRAY<VECTOR<T,3> > vorticity(grid.number_of_nodes,false);VORTICITY_DYADIC<T>::Vorticity(grid,V_ghost,vorticity);
    ARRAY<T> vorticity_magnitude(grid.number_of_nodes);for(i=1;i<=grid.number_of_nodes;i++) vorticity_magnitude(i)=vorticity(i).Magnitude();
    T minimum_cell_size=grid.Minimum_Edge_Length();T one_over_two_dx=1/(2*minimum_cell_size);VECTOR<T,3> x_offset(minimum_cell_size,0,0),y_offset(0,minimum_cell_size,0),
        z_offset(0,0,minimum_cell_size);
    PHYSBAM_NOT_IMPLEMENTED(); // The following loop won't work since it touches the extreme boundary nodes, and the looks further out.  Also, it should probably be looping over cells.
    for(DYADIC_GRID_ITERATOR_NODE<OCTREE_GRID<T> > iterator(grid,grid.Map_All_Nodes());iterator.Valid();iterator.Next()){
        VECTOR<T,3> X=iterator.Location();
        VECTOR<T,3> vortex_normal_vector((LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,vorticity_magnitude,X+x_offset)
            -LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,vorticity_magnitude,X-x_offset))*one_over_two_dx,
            (LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,vorticity_magnitude,X+y_offset)
            -LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,vorticity_magnitude,X-y_offset))*one_over_two_dx,
            (LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,vorticity_magnitude,X+z_offset)
            -LINEAR_INTERPOLATION_OCTREE_HELPER<T>::Interpolate_Nodes(grid,vorticity_magnitude,X-z_offset))*one_over_two_dx);
        vortex_normal_vector.Normalize();
        F(i)=VECTOR<T,3>::Cross_Product(vortex_normal_vector,vorticity(i));}
}
//#####################################################################
// Function Compute_Aggregate_Face_Velocities
//#####################################################################
template<class T> T INCOMPRESSIBLE_OCTREE<T>::
Compute_Aggregate_Face_Velocities(void* nothing,const T v1,const T v2,const T v3,const T v4)
{
    return (T).25*(v1+v2+v3+v4);
}
//#####################################################################
// Function Interpolate_Face_Velocities_Valid_To_Direct_Children
//#####################################################################
template<class T> void INCOMPRESSIBLE_OCTREE<T>::
Interpolate_Face_Velocities_Valid_To_Direct_Children(const OCTREE_CELL<T>* cell,ARRAY<bool>& face_values)
{
    assert(cell->Has_Children());
    const OCTREE_CHILDREN<T>* children=cell->children;
    face_values(children->Face(0,0))=face_values(children->Face(2,0))=face_values(children->Face(4,0))=face_values(children->Face(6,0))=face_values(cell->Face(0)); // LX external
    face_values(children->Face(1,1))=face_values(children->Face(3,1))=face_values(children->Face(5,1))=face_values(children->Face(7,1))=face_values(cell->Face(1)); // HX external
    face_values(children->Face(0,2))=face_values(children->Face(1,2))=face_values(children->Face(4,2))=face_values(children->Face(5,2))=face_values(cell->Face(2)); // LY external
    face_values(children->Face(2,3))=face_values(children->Face(3,3))=face_values(children->Face(6,3))=face_values(children->Face(7,3))=face_values(cell->Face(3)); // HY external
    face_values(children->Face(0,4))=face_values(children->Face(1,4))=face_values(children->Face(2,4))=face_values(children->Face(3,4))=face_values(cell->Face(4)); // LZ external
    face_values(children->Face(4,5))=face_values(children->Face(5,5))=face_values(children->Face(6,5))=face_values(children->Face(7,5))=face_values(cell->Face(5)); // HZ external
    face_values(children->Face(0,1))=face_values(children->Face(2,1))=face_values(children->Face(4,1))=face_values(children->Face(6,1))=face_values(cell->Face(0)) && face_values(cell->Face(1)); // x internal
    face_values(children->Face(0,3))=face_values(children->Face(1,3))=face_values(children->Face(4,3))=face_values(children->Face(5,3))=face_values(cell->Face(2)) && face_values(cell->Face(3)); // y internal
    face_values(children->Face(0,5))=face_values(children->Face(1,5))=face_values(children->Face(2,5))=face_values(children->Face(3,5))=face_values(cell->Face(4)) && face_values(cell->Face(5)); // z internal
}
//#####################################################################
template class INCOMPRESSIBLE_OCTREE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INCOMPRESSIBLE_OCTREE<double>;
#endif
#endif
