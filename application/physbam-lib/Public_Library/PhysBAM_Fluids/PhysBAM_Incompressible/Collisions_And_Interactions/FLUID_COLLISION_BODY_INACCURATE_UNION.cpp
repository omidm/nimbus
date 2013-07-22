//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_COLLISION_BODY_INACCURATE_UNION
//##################################################################### 
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/DYADIC_IMPLICIT_OBJECT_ON_A_RAY.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_POLICY_DYADIC.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_POLICY_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/RLE_LEVELSET_ON_A_RAY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/IMPLICIT_OBJECT_ON_A_RAY.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/FLUID_COLLISION_BODY_INACCURATE_UNION.h>
using namespace PhysBAM;
//#####################################################################
// Function Constructor
//#####################################################################
template<class T_GRID> FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
FLUID_COLLISION_BODY_INACCURATE_UNION(T_GRID& grid_input)
    :collision_bodies(grid_input),contour_value(0),grid(grid_input),levelset(grid_input,phi)
{
    collision_geometries_for_rasterization=&collision_bodies.collision_geometry_collection.bodies;
}
template<class T_GRID> FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
FLUID_COLLISION_BODY_INACCURATE_UNION(T_GRID& grid_input,T contour_value_input)
    :contour_value(contour_value_input),grid(grid_input),levelset(grid_input,phi)
{
}    
//#####################################################################
// Function Destructor
//#####################################################################
template<class T_GRID> FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
~FLUID_COLLISION_BODY_INACCURATE_UNION()
{
}
//#####################################################################
// Function Implicit_Geometry_Extended_Value
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
Implicit_Geometry_Extended_Value_Helper(const TV& location,UNIFORM_TAG<TV>) const
{return interpolation.Clamped_To_Array(grid,phi,location);}
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template<class T_GRID> typename T_GRID::SCALAR FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
Implicit_Geometry_Extended_Value_Helper(const TV& location,DYADIC_TAG<TV>) const {PHYSBAM_NOT_IMPLEMENTED();}
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template<class T_GRID> typename T_GRID::SCALAR FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
Implicit_Geometry_Extended_Value_Helper(const TV& location,RLE_TAG<TV>) const {PHYSBAM_NOT_IMPLEMENTED();}
#endif
template<class T_GRID> typename T_GRID::SCALAR FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
Implicit_Geometry_Extended_Value(const TV& location) const
{
    return Implicit_Geometry_Extended_Value_Helper(location,typename T_GRID::GRID_TAG());
}
//#####################################################################
// Function Update_Intersection_Acceleration_Structures
//#####################################################################
template<class T_GRID> void FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
Update_Intersection_Acceleration_Structures(const bool use_swept_simplex_hierarchy,const int state1,const int state2)
{
    for(COLLISION_GEOMETRY_ID i(1);i<=collision_bodies.collision_geometry_collection.bodies.m;i++) if(collision_bodies.Is_Active(i)) collision_bodies.collision_geometry_collection.bodies(i)->Update_Intersection_Acceleration_Structures(use_swept_simplex_hierarchy,state1,state2);
}
//#####################################################################
// Function Restore_State
//#####################################################################
template<class T_GRID> void FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
Restore_State(const int state_index)
{
    for(COLLISION_GEOMETRY_ID i(1);i<=collision_bodies.collision_geometry_collection.bodies.m;i++) if(collision_bodies.Is_Active(i)) collision_bodies.collision_geometry_collection.bodies(i)->Restore_State(state_index);
}
//#####################################################################
// Function Save_State
//#####################################################################
template<class T_GRID> void FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
Save_State(const int state_index,const T time)
{
    for(COLLISION_GEOMETRY_ID i(1);i<=collision_bodies.collision_geometry_collection.bodies.m;i++) if(collision_bodies.Is_Active(i)) collision_bodies.collision_geometry_collection.bodies(i)->Save_State(state_index,time);
}
//#####################################################################
// Function Read_State
//#####################################################################
template<class T_GRID> void FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
Read_State(TYPED_ISTREAM& input,const int state_index)
{
    for(COLLISION_GEOMETRY_ID i(1);i<=collision_bodies.collision_geometry_collection.bodies.m;i++) if(collision_bodies.Is_Active(i)) collision_bodies.collision_geometry_collection.bodies(i)->Read_State(input,state_index);
}
//#####################################################################
// Function Write_State
//#####################################################################
template<class T_GRID> void FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
Write_State(TYPED_OSTREAM& output,const int state_index) const
{
    for(COLLISION_GEOMETRY_ID i(1);i<=collision_bodies.collision_geometry_collection.bodies.m;i++) if(collision_bodies.Is_Active(i)) collision_bodies.collision_geometry_collection.bodies(i)->Write_State(output,state_index);
}
//#####################################################################
// Function Initialize_Grid_Structures
//#####################################################################
template<class T_GRID> void FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
Initialize_Grid_Structures(const T_GRID& grid_input,OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID>& objects_in_cell,const COLLISION_GEOMETRY_ID id) const
{
    if(&grid!=&grid_input) PHYSBAM_FATAL_ERROR();
    const_cast<FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>&>(*this).Initialize_Grid_Structures_Helper(objects_in_cell,id,typename T_GRID::GRID_TAG());
}
//#####################################################################
// Function Initialize_Grid_Structures_Helper
//#####################################################################
template<class T_GRID> void FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
Initialize_Grid_Structures_Helper(OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID>& objects_in_cell,const COLLISION_GEOMETRY_ID id,UNIFORM_TAG<TV>)
{
    collision_bodies.collision_geometry_collection.Update_Bounding_Boxes();
    // phi and velocity
    phi.Resize(grid.Cell_Indices(3),false,false);phi.Fill(10*grid.Minimum_Edge_Length());
    face_velocities.Resize(grid,3,false,false);face_velocities.Fill((T)0);
    face_velocities_set.Resize(grid,3,false,false);face_velocities_set.Fill(false);
    T_FACE_ARRAYS_INT face_velocities_count(grid,3);
    T_FACE_ARRAYS_COLLISION_GEOMETRY_ID face_operations(grid,3);
    for(COLLISION_GEOMETRY_ID i(1);i<=collision_bodies.collision_geometry_collection.bodies.m;i++)if(collision_bodies.Is_Active(i) && collision_bodies.collision_geometry_collection.bodies(i)->active)
        Initialize_Grid_Structures_Subobject(face_velocities_count,face_operations,i,typename T_GRID::GRID_TAG());
    for(typename T_GRID::FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) if(face_velocities_count.Component(iterator.Axis())(iterator.Face_Index())){
        face_velocities.Component(iterator.Axis())(iterator.Face_Index())/=face_velocities_count.Component(iterator.Axis())(iterator.Face_Index());
        face_velocities_set.Component(iterator.Axis())(iterator.Face_Index())=true;}
}
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template<class T_GRID> void FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
Initialize_Grid_Structures_Helper(OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID>& objects_in_cell,const COLLISION_GEOMETRY_ID id,DYADIC_TAG<TV>)
{
    collision_bodies.collision_geometry_collection.Update_Bounding_Boxes();
    // phi and velocity
    phi.Resize(grid.number_of_cells,false,false);phi.Fill(10*grid.Minimum_Edge_Length());
    face_velocities.Resize(grid.number_of_faces,false,false);face_velocities.Fill((T)0);
    face_velocities_set.Resize(grid.number_of_faces,false,false);face_velocities_set.Fill(false);
    ARRAY<int> face_velocities_count(grid.number_of_faces);
    OPERATION_HASH<> face_operations(grid.number_of_faces);
    for(COLLISION_GEOMETRY_ID i(1);i<=collision_bodies.collision_geometry_collection.bodies.m;i++) if(collision_bodies.Is_Active(i) && collision_bodies.collision_geometry_collection.bodies(i)->active)
        Initialize_Grid_Structures_Subobject(face_velocities_count,face_operations,i,typename T_GRID::GRID_TAG());
    for(typename T_GRID::FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()) if(face_velocities_count(iterator.Face_Index())){
        face_velocities(iterator.Face_Index())/=face_velocities_count(iterator.Face_Index());
        face_velocities_set(iterator.Face_Index())=true;}
}
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template<class T_GRID> void FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
Initialize_Grid_Structures_Helper(OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID>& objects_in_cell,const COLLISION_GEOMETRY_ID id,RLE_TAG<TV>)
{
    collision_bodies.collision_geometry_collection.Update_Bounding_Boxes();
    // phi and velocity
    phi.Resize(grid.number_of_cells,false,false);phi.Fill(10*grid.Minimum_Edge_Length());
    face_velocities.Resize(grid.number_of_faces,false,false);face_velocities.Fill((T)0);
    face_velocities_set.Resize(grid.number_of_faces,false,false);face_velocities_set.Fill(false);
    ARRAY<int> face_velocities_count(grid.number_of_faces);
    OPERATION_HASH<> face_operations(grid.number_of_faces);
    for(COLLISION_GEOMETRY_ID i(1);i<=collision_bodies.collision_geometry_collection.bodies.m;i++) if(collision_bodies.Is_Active(i) && collision_bodies.collision_geometry_collection.bodies(i)->active)
        Initialize_Grid_Structures_Subobject(face_velocities_count,face_operations,i,typename T_GRID::GRID_TAG());
    for(int f=1;f<=grid.number_of_faces;f++) if(face_velocities_count(f)){
        face_velocities(f)/=face_velocities_count(f);
        face_velocities_set(f)=true;}
}
#endif
//#####################################################################
// Function Initialize_Grid_Structures_Subobject
//#####################################################################
template<class T_GRID> void FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
Initialize_Grid_Structures_Subobject(T_FACE_ARRAYS_INT& face_velocities_count,T_FACE_ARRAYS_COLLISION_GEOMETRY_ID& face_operations,const COLLISION_GEOMETRY_ID subobject,UNIFORM_TAG<TV>)
{
    COLLISION_GEOMETRY<TV>& collision_body=*collision_bodies.collision_geometry_collection.bodies(subobject);
    RANGE<TV> bounding_box=collision_body.Axis_Aligned_Bounding_Box();
    RANGE<TV_INT> box=grid.Clamp_To_Cell(bounding_box,2).Thickened(1);
    for(CELL_ITERATOR iterator(grid,box);iterator.Valid();iterator.Next()){
        T phi_value=collision_body.Implicit_Geometry_Extended_Value(iterator.Location());
        phi(iterator.Cell_Index())=min(phi_value,phi(iterator.Cell_Index()));
        if(phi_value<0) for(int axis=1;axis<=T_GRID::dimension;axis++){
            TV_INT face1=iterator.First_Face_Index(axis),face2=iterator.Second_Face_Index(axis);
            if(face_operations.Component(axis)(face1)!=subobject){face_operations.Component(axis)(face1)=subobject;
                face_velocities.Component(axis)(face1)=collision_body.Pointwise_Object_Velocity(grid.Face(axis,face1))[axis];face_velocities_count.Component(axis)(face1)++;}
            if(face_operations.Component(axis)(face2)!=subobject){face_operations.Component(axis)(face2)=subobject;
                face_velocities.Component(axis)(face2)=collision_body.Pointwise_Object_Velocity(grid.Face(axis,face2))[axis];face_velocities_count.Component(axis)(face2)++;}}}
}
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template<class T_GRID> void FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
Initialize_Grid_Structures_Subobject(T_FACE_ARRAYS_INT& face_velocities_count,OPERATION_HASH<>& face_operations,const COLLISION_GEOMETRY_ID subobject,DYADIC_TAG<TV>)
{
    COLLISION_GEOMETRY<TV>& collision_body=*collision_bodies.collision_geometry_collection.bodies(subobject);
    RANGE<TV> bounding_box=collision_body.Axis_Aligned_Bounding_Box();
    RANGE<TV_INT> coarse_box=grid.uniform_grid.Clamp_To_Cell(bounding_box,grid.number_of_ghost_cells-1).Thickened(1);
    face_operations.Next_Operation();
    for(CELL_ITERATOR iterator(grid,coarse_box);iterator.Valid();iterator.Next()){
        T phi_value=collision_body.Implicit_Geometry_Extended_Value(iterator.Location());
        phi(iterator.Cell_Index())=min(phi_value,phi(iterator.Cell_Index()));
        if(phi_value<0){
            const typename T_GRID::CELL* cell=iterator.Cell_Pointer();
            for(int f=0;f<T_GRID::number_of_faces_per_cell;f++){int face=cell->Face(f);
                if(!face_operations.Is_Marked_Current(face)){face_operations.Mark(face);
                    face_velocities(face)=collision_body.Pointwise_Object_Velocity(grid.Face_Location(f,cell))[f/2+1];face_velocities_count(face)++;}}}}
}
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template<class T_GRID> void FLUID_COLLISION_BODY_INACCURATE_UNION<T_GRID>::
Initialize_Grid_Structures_Subobject(T_FACE_ARRAYS_INT& face_velocities_count,OPERATION_HASH<>& face_operations,const COLLISION_GEOMETRY_ID subobject,RLE_TAG<TV>)
{
    COLLISION_GEOMETRY<TV>& collision_body=*collision_bodies.collision_geometry_collection.bodies(subobject);
    RANGE<TV> bounding_box=collision_body.Axis_Aligned_Bounding_Box();
    typename T_GRID::BOX_HORIZONTAL_INT box=grid.horizontal_grid.Clamp_To_Cell(bounding_box.Get_Horizontal_Box(),grid.number_of_ghost_cells-1).Thickened(1); // TODO: localize vertically as well
    face_operations.Next_Operation();
    for(CELL_ITERATOR cell(grid,box);cell;cell++)if(cell.Short()){int c=cell.Cell();
        T phi_value=collision_body.Implicit_Geometry_Extended_Value(cell.X());
        phi(c)=min(phi_value,phi(c));
        if(phi_value<0){
            for(int f=1;f<=T_GRID::number_of_faces_per_cell;f++){int face=cell.Face(f);
                if(!face_operations.Is_Marked_Current(face)){face_operations.Mark(face);
                    face_velocities(face)=collision_body.Pointwise_Object_Velocity(cell.Face_Location(f))[(f-1)/2+1];face_velocities_count(face)++;}}}}
}
#endif
//##################################################################### 
#define P(...) __VA_ARGS__
#define INSTANTIATION_HELPER_GRID(T_GRID) \
    template FLUID_COLLISION_BODY_INACCURATE_UNION<P(T_GRID) >::FLUID_COLLISION_BODY_INACCURATE_UNION(P(T_GRID)&); \
    template void FLUID_COLLISION_BODY_INACCURATE_UNION<P(T_GRID) >::Initialize_Grid_Structures(P(T_GRID) const&,OBJECTS_IN_CELL<P(T_GRID),COLLISION_GEOMETRY_ID>&,COLLISION_GEOMETRY_ID) const;
#define INSTANTIATION_HELPER_T(T) \
    INSTANTIATION_HELPER_GRID(P(GRID<VECTOR<T,1> >)) INSTANTIATION_HELPER_GRID(P(GRID<VECTOR<T,2> >)) INSTANTIATION_HELPER_GRID(P(GRID<VECTOR<T,3> >))

INSTANTIATION_HELPER_T(float);
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
INSTANTIATION_HELPER_GRID(RLE_GRID_2D<float>);
INSTANTIATION_HELPER_GRID(RLE_GRID_3D<float>);
#endif
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
INSTANTIATION_HELPER_GRID(OCTREE_GRID<float>);
INSTANTIATION_HELPER_GRID(QUADTREE_GRID<float>);
#endif
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER_T(double);
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
INSTANTIATION_HELPER_GRID(RLE_GRID_2D<double>);
INSTANTIATION_HELPER_GRID(RLE_GRID_3D<double>);
#endif
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
INSTANTIATION_HELPER_GRID(OCTREE_GRID<double>);
INSTANTIATION_HELPER_GRID(QUADTREE_GRID<double>);
#endif
#endif
