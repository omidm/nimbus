//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_POINT_SIMPLEX_1D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_SEGMENT_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
using namespace PhysBAM;
template<class TV> UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>::
UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO(const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<TV> >& collection)
    :collision_bodies_affecting_fluid(collection),grid(collision_bodies_affecting_fluid.grid),outside_fluid(0),
    use_collision_face_neighbors(collision_bodies_affecting_fluid.use_collision_face_neighbors)
{}
//#####################################################################
// Function Register_Neighbors_As_Collision_Faces
//#####################################################################
template<class TV> void UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>::
Register_Neighbors_As_Collision_Faces()
{
    COLLISION_FACE_INFO<TV> cfi;
    HASHTABLE<FACE_INDEX<TV::m>,int> old_faces;
    HASHTABLE<FACE_INDEX<TV::m>,int> new_faces;
    VECTOR<FACE_INDEX<TV::m>,TV::m*2> faces;
    VECTOR<ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> >,TV::m> simplices,merged;
    for(int i=1;i<=collision_face_info.m;i++) old_faces.Set(FACE_INDEX<TV::m>(collision_face_info(i).axis,collision_face_info(i).index),i);
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next()){
        if((*outside_fluid)(it.index)) continue;
        grid.Neighboring_Faces(faces,it.index);
        for(int i=1;i<=TV::m;i++){simplices(i).Remove_All();merged(i).Remove_All();}
        for(int i=1;i<=faces.m;i++) if(int* a=old_faces.Get_Pointer(faces(i))) simplices((i+1)/2).Append_Unique_Elements(collision_face_info(*a).simplices);
        for(int i=1;i<=TV::m;i++) if(simplices(i).m) for(int j=1;j<=TV::m;j++) if(j!=i) merged(j).Append_Unique_Elements(simplices(i));
        for(int i=1;i<=faces.m;i++) if(merged((i+1)/2).m && !old_faces.Contains(faces(i))){
            if(int* a=new_faces.Get_Pointer(faces(i))) collision_face_info(*a).simplices.Append_Unique_Elements(merged((i+1)/2));
            else{
                cfi.axis=faces(i).axis;
                cfi.index=faces(i).index;
                cfi.side=1;
                collision_face_info.Append(cfi);
                collision_face_info.Last().simplices=merged((i+1)/2);
                new_faces.Set(faces(i),collision_face_info.m);}}}

    Sort(collision_face_info);
}
//#####################################################################
// Function Initialize_Collision_Aware_Face_Iterator
//#####################################################################
template<class TV> void UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>::
Initialize_Collision_Aware_Face_Iterator(const ARRAY<bool,TV_INT>& outside_fluid_input,const ARRAY<bool,FACE_INDEX<d> >& kinematic_faces,int ghost_cells,const bool disable_thinshell)
{
    outside_fluid=&outside_fluid_input;
    ARRAY<ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> >,TV_INT> simplices_in_cell;
    iterator_rasterization_thickness=collision_bodies_affecting_fluid.collision_thickness;
    collision_bodies_affecting_fluid.Compute_Simplices_In_Cell(simplices_in_cell,coupling_bodies,ghost_cells,iterator_rasterization_thickness,true);
    const T collision_body_thickness_over_two=(T).5*iterator_rasterization_thickness;

    COLLISION_FACE_INFO<TV> cfi;
    collision_face_info.Remove_All();
    for(UNIFORM_GRID_ITERATOR_FACE<TV> iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
        // don't couple ghost cells (for mpi assigns face to the correct proc)
        if((!grid.Inside_Domain(first_cell_index) && (*outside_fluid)(second_cell_index))
           || (!grid.Inside_Domain(second_cell_index) && (*outside_fluid)(first_cell_index))) continue;
        // skip Neumann faces or faces with both cells outside fluid
        if(!use_collision_face_neighbors && (kinematic_faces(iterator.Full_Index()) ||
            ((*outside_fluid)(first_cell_index) && (*outside_fluid)(second_cell_index)))) continue;
        // skip if both cells inside fluid for non thinshell cases
        if(!use_collision_face_neighbors && disable_thinshell && !(*outside_fluid)(first_cell_index) && !(*outside_fluid)(second_cell_index)) continue;

        ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> > &list1=simplices_in_cell(first_cell_index);
        ARRAY<PAIR<COLLISION_GEOMETRY_ID,int> > &list2=simplices_in_cell(second_cell_index);
        if(!list1.m && !list2.m) continue;
        cfi.simplices=list1;
        cfi.simplices.Append_Unique_Elements(list2);
        cfi.axis=iterator.Axis();
        cfi.index=iterator.Face_Index();

        bool ray_intersection=false;
        RAY<TV> ray(grid.X(second_cell_index),-TV::Axis_Vector(cfi.axis),true);ray.t_max=grid.dX(cfi.axis);ray.semi_infinite=false;
        for(int i=cfi.simplices.m;i>=1;i--){
            COLLISION_GEOMETRY<TV>* body=coupling_bodies(cfi.simplices(i).x);
            typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX simplex=body->World_Space_Simplex(cfi.simplices(i).y);
            if(!ray_intersection &&
                (simplex.Inside(grid.X(first_cell_index),collision_body_thickness_over_two) || simplex.Inside(grid.X(second_cell_index),collision_body_thickness_over_two) || INTERSECTION::Intersects(ray,simplex,collision_body_thickness_over_two))) ray_intersection=true;}
        if(ray_intersection){
            if(!(*outside_fluid)(first_cell_index) && grid.Inside_Domain(first_cell_index)){
                cfi.side=1;collision_face_info.Append(cfi);}
            if(!(*outside_fluid)(second_cell_index) && grid.Inside_Domain(second_cell_index)){
                cfi.side=2;collision_face_info.Append(cfi);}}}
    if(use_collision_face_neighbors) Register_Neighbors_As_Collision_Faces();
}
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<VECTOR<float,1> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<VECTOR<float,2> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<VECTOR<double,1> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<VECTOR<double,2> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<VECTOR<double,3> >;
#endif
