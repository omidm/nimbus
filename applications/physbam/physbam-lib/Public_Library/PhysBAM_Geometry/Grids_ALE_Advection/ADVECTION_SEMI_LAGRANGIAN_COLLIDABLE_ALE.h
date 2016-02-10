//#####################################################################
// Copyright 2011, Linhai Qiu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_ALE
//#####################################################################
#ifndef __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_ALE__
#define __ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_ALE__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GRID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/AVERAGING_COLLIDABLE_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_NESTED_LOOKUP,class T_FACE_LOOKUP_COLLIDABLE>
class ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_ALE:public ADVECTION<T_GRID,T2>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;
    typedef typename T_ARRAYS_BOOL::template REBIND<VECTOR<bool,T_GRID::dimension> >::TYPE T_ARRAYS_BOOL_DIMENSION;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE::template REBIND<bool>::TYPE T_ARRAYS_BASE_BOOL;
    typedef typename T_ARRAYS_BASE_BOOL::template REBIND<VECTOR<bool,T_GRID::dimension> >::TYPE T_ARRAYS_BASE_BOOL_DIMENSION;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_BOOL::template REBIND<VECTOR<bool,T_GRID::dimension> >::TYPE T_FACE_ARRAYS_BOOL_DIMENSION;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;typedef typename REBIND<T_BOUNDARY,T2>::TYPE T_BOUNDARY_T2;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename T_FACE_LOOKUP_COLLIDABLE::template REBIND_NESTED_LOOKUP<T_NESTED_LOOKUP>::TYPE T_FACE_LOOKUP_COLLIDABLE_NESTED_LOOKUP;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
public:
    RIGID_GRID<T_GRID>& rigid_grid;
    const T_GRID_BASED_COLLISION_GEOMETRY& body_list;
    T_FACE_ARRAYS_BOOL& face_velocities_valid_mask;
    T_ARRAYS_BOOL &cell_valid_points_current,&cell_valid_points_next;
    T2 cell_crossover_replacement_value;
    bool extrapolate_to_revalidate_interpolation;
private:
    LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM<T_GRID,T2> linear_interpolation_collidable_cell;
    LINEAR_INTERPOLATION_COLLIDABLE_FACE_UNIFORM<T_GRID,T,T_FACE_LOOKUP_COLLIDABLE> linear_interpolation_collidable_face;
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2,T_NESTED_LOOKUP> linear_interpolation;
    AVERAGING_UNIFORM<T_GRID,T_NESTED_LOOKUP> averaging;
    AVERAGING_COLLIDABLE_UNIFORM<T_GRID,T_FACE_LOOKUP_COLLIDABLE> averaging_collidable;
public:
    ADVECTION_SEMI_LAGRANGIAN_COLLIDABLE_ALE(RIGID_GRID<T_GRID>& rigid_grid_input,const T_GRID_BASED_COLLISION_GEOMETRY& body_list_input,T_FACE_ARRAYS_BOOL& face_velocities_valid_mask_input,T_ARRAYS_BOOL& cell_valid_points_current_input,T_ARRAYS_BOOL& cell_valid_points_next_input,const T2& default_cell_replacement_value_input,const bool extrapolate_to_revalidate_interpolation_input)
        :rigid_grid(rigid_grid_input),body_list(body_list_input),face_velocities_valid_mask(face_velocities_valid_mask_input),cell_valid_points_current(cell_valid_points_current_input),cell_valid_points_next(cell_valid_points_next_input),cell_crossover_replacement_value(default_cell_replacement_value_input),extrapolate_to_revalidate_interpolation(extrapolate_to_revalidate_interpolation_input),linear_interpolation_collidable_cell(body_list,&cell_valid_points_current,cell_crossover_replacement_value,extrapolate_to_revalidate_interpolation),averaging_collidable(body_list,0)
    {}

    void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
        const T_NESTED_LOOKUP& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max)
    {
        T_FACE_LOOKUP_COLLIDABLE_NESTED_LOOKUP V_lookup(face_velocities,body_list,&face_velocities_valid_mask);
        T_ARRAYS_VECTOR V(grid.Domain_Indices());
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT cell=iterator.Cell_Index();TV current_X=iterator.Location();
            TV old_X=rigid_grid.Previous_Object_Space_Point(current_X);TV_INT old_cell=rigid_grid.grid.Clamp_To_Cell(old_X);
            TV length_and_direction=-dt*linear_interpolation_collidable_face.Clamped_To_Array_Face(grid,V_lookup.Starting_Point_Cell(cell),old_X),interpolation_point=old_X+length_and_direction;
            if(!body_list.Swept_Occupied_Cell_Center(old_cell)){
                cell_valid_points_next(cell)=true;Z(cell)=linear_interpolation.Clamped_To_Array(grid,Z_ghost,interpolation_point);}
            else{
                COLLISION_GEOMETRY_ID body_id;int aggregate_id;TV initial_hit_point;
                if(body_list.Latest_Crossover(old_X,old_X,dt,body_id,aggregate_id,initial_hit_point)){
                    cell_valid_points_next(cell)=false;
                    Z(cell)=cell_crossover_replacement_value;}
                else{
                    RAY<TV> backtrace_ray;COLLISION_GEOMETRY_ID body_id;
                    if(RAY<TV>::Create_Non_Degenerate_Ray(old_X,length_and_direction,backtrace_ray) && body_list.Closest_Non_Intersecting_Point_Of_Any_Body(backtrace_ray,body_id))
                        interpolation_point=backtrace_ray.Point(backtrace_ray.t_max);
                    Z(cell)=linear_interpolation_collidable_cell.Clamped_To_Array(grid,Z_ghost,interpolation_point,cell_valid_points_next(cell));}}}
        T_ARRAYS_BOOL::Exchange_Arrays(cell_valid_points_current,cell_valid_points_next);
    }

    void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_NESTED_LOOKUP& Z_ghost,
        const T_NESTED_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
        const T_NESTED_LOOKUP* Z_min_ghost,const T_NESTED_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max)
    {
        T_FACE_LOOKUP_COLLIDABLE_NESTED_LOOKUP V_lookup(face_velocities,body_list,&face_velocities_valid_mask);
        T_FACE_LOOKUP_COLLIDABLE_NESTED_LOOKUP Z_ghost_lookup(Z_ghost,body_list,&face_velocities_valid_mask);
        T_FACE_ARRAYS_BOOL face_velocities_valid_mask_next(grid,3,false);
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
            TV_INT face=iterator.Face_Index();int axis=iterator.Axis();
            TV current_X=iterator.Location();TV old_X=rigid_grid.Previous_Object_Space_Point(current_X);TV_INT old_cell=rigid_grid.grid.Clamp_To_Cell(old_X);
            TV length_and_direction=-dt*linear_interpolation_collidable_face.Clamped_To_Array_Face(grid,V_lookup.Starting_Point_Face(axis,face),old_X),interpolation_point=old_X+length_and_direction;
            TV interpolated_velocity;
            if(!body_list.Swept_Occupied_Cell_Center(old_cell)){
                face_velocities_valid_mask_next.Component(axis)(face)=true;
                interpolated_velocity=linear_interpolation.Clamped_To_Array_Face(grid,Z_ghost,interpolation_point);}
            else{
                COLLISION_GEOMETRY_ID body_id;int aggregate_id;TV initial_hit_point;
                if(body_list.Latest_Crossover(old_X,old_X,dt,body_id,aggregate_id,initial_hit_point)){
                    face_velocities_valid_mask_next.Component(axis)(face)=false;
                    interpolated_velocity=body_list.collision_geometry_collection(body_id).Pointwise_Object_Velocity(aggregate_id,initial_hit_point);}
                else{
                    RAY<TV> backtrace_ray;COLLISION_GEOMETRY_ID body_id;
                    if(RAY<TV>::Create_Non_Degenerate_Ray(old_X,length_and_direction,backtrace_ray) && body_list.Closest_Non_Intersecting_Point_Of_Any_Body(backtrace_ray,body_id))
                        interpolation_point=backtrace_ray.Point(backtrace_ray.t_max);
                    const typename T_FACE_LOOKUP_COLLIDABLE::LOOKUP& lookup=Z_ghost_lookup.Starting_Point_Face(axis,face);
                    interpolated_velocity=linear_interpolation_collidable_face.Clamped_To_Array_Face(grid,lookup,interpolation_point);
                    face_velocities_valid_mask_next.Component(axis)(face)=lookup.found_valid_point;}}
            Z.Component(axis)(face)=rigid_grid.Current_Object_Space_Vector(interpolated_velocity)(axis);}
        T_FACE_ARRAYS_BOOL::Exchange_Arrays(face_velocities_valid_mask,face_velocities_valid_mask_next);
        // ghost values should always be valid
        for(int axis=1;axis<=T_GRID::dimension;axis++) grid.Get_Face_Grid(axis).Put_Ghost(true,face_velocities_valid_mask.Component(axis),3);
    }
    
    void Average_To_Invalidated_Cells(const T_GRID& grid,const T2 default_value,T_ARRAYS_T2& values)
    {// average values collision aware in Gauss-Jacobi fashion
        const T_ARRAYS_BOOL_DIMENSION& cell_neighbors_visible=body_list.cell_neighbors_visible;
        bool done=false;ARRAY<PAIR<TV_INT,bool> > invalid_indices; // index and bool true if entry has been validated on iteration
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next())if(!cell_valid_points_current(iterator.Cell_Index()))
            invalid_indices.Append(PAIR<TV_INT,bool>(iterator.Cell_Index(),false));
        grid.Put_Ghost(false,cell_valid_points_current,3); // don't average from boundaries
        
        while(!done){done=true;
            for(int k=1;k<=invalid_indices.m;k++){
                T2 sum=T2();int count=0;
                for(int axis=1;axis<=T_GRID::dimension;axis++){
                    TV_INT min_cell=invalid_indices(k).x-TV_INT::Axis_Vector(axis),max_cell=invalid_indices(k).x+TV_INT::Axis_Vector(axis);
                    if(cell_neighbors_visible(min_cell)(axis) && cell_valid_points_current(min_cell)){sum+=values(min_cell);count++;}
                    if(cell_neighbors_visible(invalid_indices(k).x)(axis) && cell_valid_points_current(max_cell)){sum+=values(max_cell);count++;}}
                if(count){values(invalid_indices(k).x)=sum/(T)count;invalid_indices(k).y=true;done=false;}}
            if(!done) for(int k=invalid_indices.m;k>=1;k--) if(invalid_indices(k).y){cell_valid_points_current(invalid_indices(k).x)=true;invalid_indices.Remove_Index_Lazy(k);}}
        
        // keep a copy of currently valid cells (used for phi so we can revalidate the remaining cells again after collision aware fast marching)
        // but important to initialize ghost cells to true since currently cell_valid_points_current has them set to false
        cell_valid_points_next=cell_valid_points_current;
        grid.Put_Ghost(true,cell_valid_points_next,3);
        
        // average values collision aware in Gauss-Jacobi fashion (here we replace non-visible values with special values defined by Compute_Revalidation_Value())
        done=false;
        while(!done){
            done=true;
            for(int k=1;k<=invalid_indices.m;k++){
                T2 sum=T2();int count=0;
                for(int axis=1;axis<=T_GRID::dimension;axis++){
                    TV_INT min_cell=invalid_indices(k).x-TV_INT::Axis_Vector(axis),max_cell=invalid_indices(k).x+TV_INT::Axis_Vector(axis);
                    if(cell_neighbors_visible(min_cell)(axis)){if(cell_valid_points_current(min_cell)){sum+=values(min_cell);count++;}}
                    else{sum+=Compute_Revalidation_Value(grid.X(invalid_indices(k).x),grid.X(min_cell),values(invalid_indices(k).x),default_value);count++;}
                    if(cell_neighbors_visible(invalid_indices(k).x)(axis)){if(cell_valid_points_current(max_cell)){sum+=values(max_cell);count++;}}
                    else{sum+=Compute_Revalidation_Value(grid.X(invalid_indices(k).x),grid.X(max_cell),values(invalid_indices(k).x),default_value);count++;}}
                if(count){values(invalid_indices(k).x)=sum/(T)count;invalid_indices(k).y=true;done=false;}
                else values(invalid_indices(k).x)=default_value;}
            if(!done) for(int k=invalid_indices.m;k>=1;k--) if(invalid_indices(k).y){cell_valid_points_current(invalid_indices(k).x)=true;invalid_indices.Remove_Index_Lazy(k);}}
        grid.Put_Ghost(true,cell_valid_points_current,3); // set valid for future advection
    }
    
    void Average_To_Invalidated_Face(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& face_values)
    {// average values collision aware in Gauss-Jacobi fashion
        typename TV::template REBIND<ARRAY<PAIR<TV_INT,bool> > >::TYPE face_invalid_indices; // index and bool true if entry has been validated on iteration
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) if(!face_velocities_valid_mask.Component(iterator.Axis())(iterator.Face_Index())) 
            face_invalid_indices[iterator.Axis()].Append(PAIR<TV_INT,bool>(iterator.Face_Index(),false));
        
        for(int arrays_axis=1;arrays_axis<=T_GRID::dimension;arrays_axis++){
            ARRAY<PAIR<TV_INT,bool> >& invalid_indices=face_invalid_indices[arrays_axis];
            T_ARRAYS_BASE_BOOL& valid_points=face_velocities_valid_mask.Component(arrays_axis);T_ARRAYS_BASE& values=face_values.Component(arrays_axis);
            
            bool done=false;
            grid.Put_Ghost(false,valid_points,3); // don't average from boundaries
            
            while(!done){
                done=true;
                for(int k=1;k<=invalid_indices.m;k++){ 
                    T sum=0;int count=0;
                    for(int axis=1;axis<=T_GRID::dimension;axis++){
                        TV_INT min_face=invalid_indices(k).x-TV_INT::Axis_Vector(axis),max_face=invalid_indices(k).x+TV_INT::Axis_Vector(axis);
                        if(body_list.face_neighbors_visible.Component(arrays_axis)(min_face)(axis) && valid_points(min_face)){sum+=values(min_face);count++;}
                        if(body_list.face_neighbors_visible.Component(arrays_axis)(invalid_indices(k).x)(axis) && valid_points(max_face)){sum+=values(max_face);count++;}}
                    if(count){values(invalid_indices(k).x)=sum/(T)count;invalid_indices(k).y=true;done=false;}}
                if(!done) for(int k=invalid_indices.m;k>=1;k--) if(invalid_indices(k).y){valid_points(invalid_indices(k).x)=true;invalid_indices.Remove_Index_Lazy(k);}}
            
            // average values collision aware in Gauss-Jacobi fashion (here we replace non-visible values with special values defined by Compute_Revalidation_Value())
            done=false;
            while(!done){
                done=true;
                for(int k=1;k<=invalid_indices.m;k++){ 
                    T sum=0;int count=0;
                    for(int axis=1;axis<=T_GRID::dimension;axis++){
                        TV_INT min_face=invalid_indices(k).x-TV_INT::Axis_Vector(axis),max_face=invalid_indices(k).x+TV_INT::Axis_Vector(axis);
                        if(body_list.face_neighbors_visible.Component(arrays_axis)(min_face)(axis)){if(valid_points(min_face)){sum+=values(min_face);count++;}}
                        else{sum+=Compute_Revalidation_Value(arrays_axis,grid.X(invalid_indices(k).x),grid.X(min_face),values(invalid_indices(k).x),T());count++;}
                        if(body_list.face_neighbors_visible.Component(arrays_axis)(invalid_indices(k).x)(axis)){if(valid_points(max_face)){sum+=values(max_face);count++;}}
                        else{sum+=Compute_Revalidation_Value(arrays_axis,grid.X(invalid_indices(k).x),grid.X(max_face),values(invalid_indices(k).x),T());count++;}}
                    if(count){values(invalid_indices(k).x)=sum/(T)count;invalid_indices(k).y=true;done=false;}
                    else values(invalid_indices(k).x)=T();}
                if(!done) for(int k=invalid_indices.m;k>=1;k--) if(invalid_indices(k).y){valid_points(invalid_indices(k).x)=true;invalid_indices.Remove_Index_Lazy(k);}}
            grid.Put_Ghost(true,valid_points,3);} // set valid for future advection
    }

    T Compute_Revalidation_Value(const int axis,const TV& from,const TV& to,const T& current_invalid_value,const T& default_value)
    {
        TV point;COLLISION_GEOMETRY_ID body_id;int aggregate_id;
        if(body_list.collision_geometry_collection.Intersection_Between_Points(from,to,body_id,aggregate_id,point)) return body_list.Object_Velocity(body_id,aggregate_id,point)[axis];
        else return default_value;
    }

    T2 Compute_Revalidation_Value(const TV& from,const TV& to,const T2& current_invalid_value,const T2& default_value)
    {return default_value;}
//#####################################################################
};
}
#endif
