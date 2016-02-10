#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic_Computations/VORTICITY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_OCTREE_HELPER.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_QUADTREE_HELPER.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/EXTRAPOLATION_DYADIC.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_QUADTREE.h>
using namespace PhysBAM;
//#####################################################################
// Function Compute_Vorticity_Confinement_Force
//#####################################################################
template<class T> void INCOMPRESSIBLE_QUADTREE<T>::
Compute_Vorticity_Confinement_Force(QUADTREE_GRID<T>& grid,const ARRAY<VECTOR<T,2> >& V_ghost,ARRAY<VECTOR<T,2> >& F)
{
    int i;ARRAY<T> vorticity(grid.number_of_nodes,false);VORTICITY_DYADIC<T>::Vorticity(grid,V_ghost,vorticity);
    ARRAY<T> vorticity_magnitude(grid.number_of_nodes);for(i=1;i<=grid.number_of_nodes;i++) vorticity_magnitude(i)=abs(vorticity(i));
    T minimum_cell_size=grid.Minimum_Edge_Length();T one_over_two_dx=1/(2*minimum_cell_size);VECTOR<T,2> x_offset(minimum_cell_size,0),y_offset(0,minimum_cell_size);
    for(i=1;i<=grid.number_of_nodes;i++){
        VECTOR<T,2> vortex_normal_vector((LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Nodes(grid,vorticity_magnitude,grid.Node_Location(i)+x_offset)
            -LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Nodes(grid,vorticity_magnitude,grid.Node_Location(i)-x_offset))*one_over_two_dx,
            (LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Nodes(grid,vorticity_magnitude,grid.Node_Location(i)+y_offset)
            -LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Nodes(grid,vorticity_magnitude,grid.Node_Location(i)-y_offset))*one_over_two_dx);
        vortex_normal_vector.Normalize();
        F(i).x=vorticity(i)*vortex_normal_vector.y;F(i).y=-vorticity(i)*vortex_normal_vector.x;}
}
//#####################################################################
// Function Compute_Aggregate_Face_Velocities
//#####################################################################
template<class T> T INCOMPRESSIBLE_QUADTREE<T>::
Compute_Aggregate_Face_Velocities(void* nothing,const T v1,const T v2)
{
    return (T).5*(v1+v2);
}
//#####################################################################
// Function Interpolate_Face_Velocities_Valid_To_Direct_Children
//#####################################################################
template<class T> void INCOMPRESSIBLE_QUADTREE<T>::
Interpolate_Face_Velocities_Valid_To_Direct_Children(const QUADTREE_CELL<T>* cell,ARRAY<bool>& face_values)
{
    assert(cell->Has_Children());
    const QUADTREE_CHILDREN<T>* children=cell->children;
    face_values(children->Face(0,0))=face_values(children->Face(2,0))=face_values(cell->Face(0));
    face_values(children->Face(1,1))=face_values(children->Face(3,1))=face_values(cell->Face(1));
    face_values(children->Face(0,2))=face_values(children->Face(1,2))=face_values(cell->Face(2));
    face_values(children->Face(2,3))=face_values(children->Face(3,3))=face_values(cell->Face(3));
    face_values(children->Face(0,3))=face_values(children->Face(1,3))=face_values(cell->Face(2)) && face_values(cell->Face(3));
    face_values(children->Face(0,1))=face_values(children->Face(2,1))=face_values(cell->Face(0)) && face_values(cell->Face(1));
}
//#####################################################################
template class INCOMPRESSIBLE_QUADTREE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INCOMPRESSIBLE_QUADTREE<double>;
#endif
#endif
