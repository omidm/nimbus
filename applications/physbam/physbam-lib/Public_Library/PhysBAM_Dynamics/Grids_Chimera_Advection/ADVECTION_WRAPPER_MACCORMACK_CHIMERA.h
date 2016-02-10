//#####################################################################
// Copyright 2012, Linhai Qiu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_WRAPPER_MACCORMACK_CHIMERA
//#####################################################################
#ifndef __ADVECTION_WRAPPER_MACCORMACK_CHIMERA__
#define __ADVECTION_WRAPPER_MACCORMACK_CHIMERA__

#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_CHIMERA.h>
namespace PhysBAM{

template<class T_GRID> class SOLIDS_FLUIDS_EXAMPLE_CHIMERA;

template<class T_GRID,class T2,class T_NESTED_ADVECTION,class T_NESTED_LOOKUP>
class ADVECTION_WRAPPER_MACCORMACK_CHIMERA:public ADVECTION<T_GRID,T2>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;typedef typename REBIND<T_BOUNDARY,T2>::TYPE T_BOUNDARY_T2;
public:
    T_NESTED_ADVECTION& nested_advection;
    SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>& example;
    const RIGID_GRID_COLLECTION<T_GRID> &rigid_grid_collection;
    const int& n_local_grids;
    const int& number_of_ghost_cells;
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T2,T_NESTED_LOOKUP> interpolation;
    bool clamp_extrema;
    bool enforce_second_order;
    
    ADVECTION_WRAPPER_MACCORMACK_CHIMERA(T_NESTED_ADVECTION& nested_advection_input,SOLIDS_FLUIDS_EXAMPLE_CHIMERA<T_GRID>& example_input,bool use_clamp_extrema,bool enforce_second_order_input)
        :nested_advection(nested_advection_input),example(example_input),rigid_grid_collection(example_input.rigid_grid_collection),n_local_grids(rigid_grid_collection.particles.array_collection->number),number_of_ghost_cells(example.fluids_parameters.number_of_ghost_cells),clamp_extrema(use_clamp_extrema),enforce_second_order(enforce_second_order_input)
    {}
    
    void Update_Advection_Equation_Cell_Chimera(ARRAY<T_ARRAYS_T2*>& Z_array,const ARRAY<T_ARRAYS_T2*>& Z_ghost_array,
        const ARRAY<T_FACE_ARRAYS_SCALAR*>& face_velocities_array,T_BOUNDARY_T2& boundary,const T dt,const T time)
    {
        ARRAY<T_ARRAYS_T2> Z_max_array(n_local_grids),Z_min_array(n_local_grids),Z_forward_array(n_local_grids),Z_backward_array(n_local_grids);
        //forward advection
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
            const RIGID_GRID<T_GRID>& rigid_grid=rigid_grid_collection.Rigid_Grid(grid_index);const T_GRID& grid=rigid_grid.grid;
            const T_ARRAYS_T2& Z_ghost=*Z_ghost_array(grid_index);
            T_NESTED_LOOKUP face_velocities=T_NESTED_LOOKUP(*face_velocities_array(grid_index));
            T_ARRAYS_T2 Z_max_ghost=Z_ghost,Z_min_ghost=Z_ghost;
            T_ARRAYS_T2& Z_max=Z_max_array(grid_index);T_ARRAYS_T2& Z_min=Z_min_array(grid_index);
            Z_max.Resize(grid.Domain_Indices());Z_min.Resize(grid.Domain_Indices());
            T_ARRAYS_T2& Z_forward=Z_forward_array(grid_index);Z_forward.Resize(grid.Domain_Indices(number_of_ghost_cells));
            T_GRID cell_center_grid=T_GRID(grid.numbers_of_cells,RANGE<TV>(grid.domain.min_corner+(T).5*grid.dX,grid.domain.max_corner-(T).5*grid.dX));
            T_ARRAYS_VECTOR V_forward(grid.Domain_Indices());
            for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                TV_INT cell=iterator.Cell_Index();TV current_X=iterator.Location();TV X_new_in_old=rigid_grid.Previous_Object_Space_Point(current_X);
                TV grid_velocity=(X_new_in_old-current_X)/dt;TV velocity_in_old_object_space;
                velocity_in_old_object_space=interpolation.Clamped_To_Array_Face(grid,face_velocities,X_new_in_old);
                V_forward(cell)=velocity_in_old_object_space-grid_velocity;}
            nested_advection.Update_Advection_Equation_Node(cell_center_grid,Z_forward,Z_ghost,V_forward,boundary,dt,time,&Z_min_ghost,&Z_max_ghost,&Z_min,&Z_max);
        }
        example.chimera_grid->Set_To_Use_Current_Grid_Frame_Cell();
        ARRAY<T_ARRAYS_T2*> Z_forward_pointer_array(n_local_grids);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++) Z_forward_pointer_array(grid_index)=&Z_forward_array(grid_index);
        example.Fill_Ghost_Cells_Chimera(number_of_ghost_cells,Z_forward_pointer_array,Z_forward_pointer_array);
        example.chimera_grid->Set_To_Use_Previous_Grid_Frame_Cell();
        //backward advection
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
            const RIGID_GRID<T_GRID>& rigid_grid=rigid_grid_collection.Rigid_Grid(grid_index);const T_GRID& grid=rigid_grid.grid;
            const T_ARRAYS_T2& Z_forward_ghost=Z_forward_array(grid_index);
            T_NESTED_LOOKUP face_velocities=T_NESTED_LOOKUP(*face_velocities_array(grid_index));
            T_ARRAYS_T2& Z_backward=Z_backward_array(grid_index);Z_backward.Resize(grid.Domain_Indices());
            T_GRID cell_center_grid=T_GRID(grid.numbers_of_cells,RANGE<TV>(grid.domain.min_corner+(T).5*grid.dX,grid.domain.max_corner-(T).5*grid.dX));
            T_ARRAYS_VECTOR V_backward(grid.Domain_Indices());
            for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                TV_INT cell=iterator.Cell_Index();TV current_X=iterator.Location();
                TV X_old_in_new=rigid_grid.Object_Space_Point(rigid_grid.previous_state.World_Space_Point(current_X));
                TV grid_velocity_backward=(current_X-X_old_in_new)/dt;
                TV negative_velocity_in_old_object_space;
                negative_velocity_in_old_object_space=-(interpolation.Clamped_To_Array_Face(grid,face_velocities,current_X));
                V_backward(cell)=rigid_grid.Current_Object_Space_Vector(negative_velocity_in_old_object_space)+grid_velocity_backward;}
            nested_advection.Update_Advection_Equation_Node(cell_center_grid,Z_backward,Z_forward_ghost,V_backward,boundary,dt,time+dt,0,0,0,0);
        }
        //error correction
        if(clamp_extrema || enforce_second_order) Apply_Clamped_Extrema_Limiter_Cell_Chimera(Z_array,Z_forward_array,Z_backward_array,Z_min_array,Z_max_array);
        else Apply_Reversion_Limiter_Cell_Chimera(Z_array,Z_forward_array,Z_backward_array,Z_min_array,Z_max_array);
    }

    void Update_Advection_Equation_Face_Chimera(ARRAY<T_FACE_ARRAYS_SCALAR*>& Z_array,const ARRAY<T_FACE_ARRAYS_SCALAR*>& Z_ghost_array,
        const ARRAY<T_FACE_ARRAYS_SCALAR*>& face_velocities_array,T_BOUNDARY& boundary,const T dt,const T time)
    {
        ARRAY<T_FACE_ARRAYS_SCALAR> Z_max_array(n_local_grids),Z_min_array(n_local_grids),Z_forward_array(n_local_grids),Z_backward_array(n_local_grids);
        ARRAY<T_ARRAYS_SCALAR> V_advected(TV::dimension);ARRAY<T_ARRAYS_SCALAR> Z_ghost_components(TV::dimension);
        //forward advection
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
            const RIGID_GRID<T_GRID>& rigid_grid=rigid_grid_collection.Rigid_Grid(grid_index);const T_GRID& grid=rigid_grid.grid;
            const T_FACE_ARRAYS_SCALAR& Z_ghost=*Z_ghost_array(grid_index);
            T_NESTED_LOOKUP face_velocities=T_NESTED_LOOKUP(*face_velocities_array(grid_index));
            T_FACE_ARRAYS_SCALAR Z_max_ghost=Z_ghost,Z_min_ghost=Z_ghost;
            T_FACE_ARRAYS_SCALAR& Z_max=Z_max_array(grid_index);T_FACE_ARRAYS_SCALAR& Z_min=Z_min_array(grid_index);
            Z_max.Resize(grid.Domain_Indices());Z_min.Resize(grid.Domain_Indices());
            T_FACE_ARRAYS_SCALAR& Z_forward=Z_forward_array(grid_index);Z_forward.Resize(grid.Domain_Indices(number_of_ghost_cells));
            ARRAY<T_ARRAYS_VECTOR> V_forward(TV::dimension);
            ARRAY<T_ARRAYS_SCALAR> Z_min_components(TV::dimension);ARRAY<T_ARRAYS_SCALAR> Z_max_components(TV::dimension);
            ARRAY<TV> min_corner_face_grids(TV::dimension);
            for(int axis=1;axis<=TV::dimension;axis++){
                min_corner_face_grids(axis)=grid.Get_Face_Grid(axis).Domain().Minimum_Corner();}
            for(int i=1;i<=TV::dimension;i++){//i is time n+1 grid dimension
                T_GRID face_grid=grid.Get_Face_Grid(i);
                for(int axis=1;axis<=TV::dimension;axis++){
                    V_advected(axis).Resize(face_grid.Domain_Indices());
                    V_forward(axis).Resize(face_grid.Domain_Indices());}
                for(NODE_ITERATOR iterator(face_grid);iterator.Valid();iterator.Next()){
                    TV_INT node=iterator.Node_Index();TV current_X=iterator.Location();
                    TV X_new_in_old=rigid_grid.Previous_Object_Space_Point(current_X);
                    TV grid_velocity=(X_new_in_old-current_X)/dt;
                    TV V_advecting_in_old_object_space=interpolation.Clamped_To_Array_Face(grid,face_velocities,X_new_in_old);
                    for(int axis=1;axis<=TV::dimension;axis++){
                        V_forward(axis)(node)=V_advecting_in_old_object_space-grid_velocity+(min_corner_face_grids(axis)-min_corner_face_grids(i))/dt;}}
                for(int axis=1;axis<=TV::dimension;axis++){ 
                    Z_ghost_components(axis)=Z_ghost.Component(axis);
                    Z_min_components(axis).Resize(face_grid.Domain_Indices());
                    Z_max_components(axis).Resize(face_grid.Domain_Indices());
                    nested_advection.Update_Advection_Equation_Node(face_grid,V_advected(axis),Z_ghost_components(axis),V_forward(axis),boundary,dt,time,&Z_ghost_components(axis),&Z_ghost_components(axis),&Z_min_components(axis),&Z_max_components(axis));}
                for(NODE_ITERATOR iterator(face_grid);iterator.Valid();iterator.Next()){
                    TV V_advected_forward,Z_min_vector,Z_max_vector;
                    for(int axis=1;axis<=TV::dimension;axis++){
                        V_advected_forward(axis)=V_advected(axis)(iterator.Node_Index());
                        Z_min_vector(axis)=Z_min_components(axis)(iterator.Node_Index());
                        Z_max_vector(axis)=Z_max_components(axis)(iterator.Node_Index());}
                    Z_forward.Component(i)(iterator.Node_Index())=rigid_grid.Current_Object_Space_Vector(V_advected_forward)(i);
                    Z_min.Component(i)(iterator.Node_Index())=rigid_grid.Current_Object_Space_Vector(Z_min_vector)(i);
                    Z_max.Component(i)(iterator.Node_Index())=rigid_grid.Current_Object_Space_Vector(Z_max_vector)(i);}}
        }
        example.chimera_grid->Set_To_Use_Current_Grid_Frame_Face();
        ARRAY<T_FACE_ARRAYS_SCALAR*> Z_forward_pointer_array(n_local_grids);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++) Z_forward_pointer_array(grid_index)=&Z_forward_array(grid_index);
        example.Fill_Ghost_Cells_Face_Chimera(number_of_ghost_cells,Z_forward_pointer_array,Z_forward_pointer_array);
        example.chimera_grid->Set_To_Use_Previous_Grid_Frame_Face();
        //backward advection
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
            const RIGID_GRID<T_GRID>& rigid_grid=rigid_grid_collection.Rigid_Grid(grid_index);const T_GRID& grid=rigid_grid.grid;
            const T_FACE_ARRAYS_SCALAR& Z_forward_ghost=Z_forward_array(grid_index);
            T_NESTED_LOOKUP face_velocities=T_NESTED_LOOKUP(*face_velocities_array(grid_index));
            T_FACE_ARRAYS_SCALAR& Z_backward=Z_backward_array(grid_index);Z_backward.Resize(grid.Domain_Indices());
            ARRAY<T_ARRAYS_VECTOR> V_backward(TV::dimension);
            ARRAY<TV> min_corner_face_grids(TV::dimension);
            for(int axis=1;axis<=TV::dimension;axis++){
                min_corner_face_grids(axis)=grid.Get_Face_Grid(axis).Domain().Minimum_Corner();}
            for(int i=1;i<=TV::dimension;i++){//i is time n grid dimension
                T_GRID face_grid=grid.Get_Face_Grid(i);
                for(int axis=1;axis<=TV::dimension;axis++){
                    V_advected(axis).Resize(face_grid.Domain_Indices());
                    V_backward(axis).Resize(face_grid.Domain_Indices());}
                for(NODE_ITERATOR iterator(face_grid);iterator.Valid();iterator.Next()){
                    TV_INT node=iterator.Node_Index();TV current_X=iterator.Location();
                    TV X_old_in_new=rigid_grid.Object_Space_Point(rigid_grid.previous_state.World_Space_Point(current_X));
                    TV grid_velocity_backward=(current_X-X_old_in_new)/dt;
                    TV negative_V_advecting_in_old_object_space=-(interpolation.Clamped_To_Array_Face(grid,face_velocities,current_X));
                    for(int axis=1;axis<=TV::dimension;axis++){
                        V_backward(axis)(node)=rigid_grid.Current_Object_Space_Vector(negative_V_advecting_in_old_object_space)+grid_velocity_backward+(min_corner_face_grids(axis)-min_corner_face_grids(i))/dt;}}
                for(int axis=1;axis<=TV::dimension;axis++){ 
                    Z_ghost_components(axis)=Z_forward_ghost.Component(axis);
                    nested_advection.Update_Advection_Equation_Node(face_grid,V_advected(axis),Z_ghost_components(axis),V_backward(axis),boundary,dt,time+dt,0,0,0,0);}
                for(NODE_ITERATOR iterator(face_grid);iterator.Valid();iterator.Next()){
                    TV V_advected_backward;
                    for(int axis=1;axis<=TV::dimension;axis++){
                        V_advected_backward(axis)=V_advected(axis)(iterator.Node_Index());}
                    Z_backward.Component(i)(iterator.Node_Index())=rigid_grid.Previous_Object_Space_Vector(V_advected_backward)(i);}}
        }
        //error correction
        if(clamp_extrema || enforce_second_order) Apply_Clamped_Extrema_Limiter_Face_Chimera(Z_array,Z_forward_array,Z_backward_array,Z_min_array,Z_max_array);
        else Apply_Reversion_Limiter_Face_Chimera(Z_array,Z_forward_array,Z_backward_array,Z_min_array,Z_max_array);
    }

    void Apply_Clamped_Extrema_Limiter_Cell_Chimera(ARRAY<T_ARRAYS_T2*>& Z_array,const ARRAY<T_ARRAYS_T2>& Z_forward_array,const ARRAY<T_ARRAYS_T2>& Z_backward_array,const ARRAY<T_ARRAYS_T2>& Z_min_array,const ARRAY<T_ARRAYS_T2>& Z_max_array)
    {
        ARRAY<T_ARRAYS_T2> errors(n_local_grids);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
            const T_GRID& grid=rigid_grid_collection.Rigid_Grid(grid_index).grid;
            errors(grid_index).Resize(grid.Domain_Indices(number_of_ghost_cells));
            for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                TV_INT index=iterator.Cell_Index();
                errors(grid_index)(index)=(T).5*((*Z_array(grid_index))(index)-Z_backward_array(grid_index)(index));
            }
        }
        ARRAY<T_ARRAYS_T2*> errors_pointer_array(n_local_grids);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++) errors_pointer_array(grid_index)=&errors(grid_index);
        example.Fill_Ghost_Cells_Chimera(number_of_ghost_cells,errors_pointer_array,errors_pointer_array);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
            const RIGID_GRID<T_GRID>& rigid_grid=rigid_grid_collection.Rigid_Grid(grid_index);const T_GRID& grid=rigid_grid.grid;
            for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                TV_INT index=iterator.Cell_Index();TV location=iterator.Location();TV location_new_in_old=rigid_grid.Previous_Object_Space_Point(location);
                T2 error=interpolation.Clamped_To_Array(grid,errors(grid_index),location_new_in_old);
                if(!enforce_second_order) 
                    (*Z_array(grid_index))(index)=clamp(Z_forward_array(grid_index)(index)+error,Z_min_array(grid_index)(index),Z_max_array(grid_index)(index));
                else 
                    (*Z_array(grid_index))(index)=Z_forward_array(grid_index)(index)+error;
            }
        }
    }

    void Apply_Reversion_Limiter_Cell_Chimera(ARRAY<T_ARRAYS_T2*>& Z_array,const ARRAY<T_ARRAYS_T2>& Z_forward_array,const ARRAY<T_ARRAYS_T2>& Z_backward_array,const ARRAY<T_ARRAYS_T2>& Z_min_array,const ARRAY<T_ARRAYS_T2>& Z_max_array)
    {
        ARRAY<T_ARRAYS_T2> errors(n_local_grids);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
            const T_GRID& grid=rigid_grid_collection.Rigid_Grid(grid_index).grid;
            errors(grid_index).Resize(grid.Domain_Indices(number_of_ghost_cells));
            for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                TV_INT index=iterator.Cell_Index();
                errors(grid_index)(index)=(T).5*((*Z_array(grid_index))(index)-Z_backward_array(grid_index)(index));
            }
        }
        ARRAY<T_ARRAYS_T2*> errors_pointer_array(n_local_grids);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++) errors_pointer_array(grid_index)=&errors(grid_index);
        example.Fill_Ghost_Cells_Chimera(number_of_ghost_cells,errors_pointer_array,errors_pointer_array);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
            const RIGID_GRID<T_GRID>& rigid_grid=rigid_grid_collection.Rigid_Grid(grid_index);const T_GRID& grid=rigid_grid.grid;
            for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                TV_INT index=iterator.Cell_Index();TV location=iterator.Location();TV location_new_in_old=rigid_grid.Previous_Object_Space_Point(location);
                T2 error=interpolation.Clamped_To_Array(grid,errors(grid_index),location_new_in_old);
                (*Z_array(grid_index))(index)=Z_forward_array(grid_index)(index)+error;
                if(!in_bounds((*Z_array(grid_index))(index),Z_min_array(grid_index)(index),Z_max_array(grid_index)(index)))
                    (*Z_array(grid_index))(index)=Z_forward_array(grid_index)(index);
            }
        }
    }

    void Apply_Clamped_Extrema_Limiter_Face_Chimera(ARRAY<T_FACE_ARRAYS_SCALAR*>& Z_array,const ARRAY<T_FACE_ARRAYS_SCALAR>& Z_forward_array,const ARRAY<T_FACE_ARRAYS_SCALAR>& Z_backward_array,const ARRAY<T_FACE_ARRAYS_SCALAR>& Z_min_array,const ARRAY<T_FACE_ARRAYS_SCALAR>& Z_max_array)
    {
        ARRAY<T_FACE_ARRAYS_SCALAR> errors(n_local_grids);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
            const T_GRID& grid=rigid_grid_collection.Rigid_Grid(grid_index).grid;
            errors(grid_index).Resize(grid.Domain_Indices(number_of_ghost_cells));
            for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                errors(grid_index)(iterator.Full_Index())=(T).5*((*Z_array(grid_index))(iterator.Full_Index())-Z_backward_array(grid_index)(iterator.Full_Index()));}
        }
        ARRAY<T_FACE_ARRAYS_SCALAR*> errors_pointer_array(n_local_grids);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++) errors_pointer_array(grid_index)=&errors(grid_index);
        example.Fill_Ghost_Cells_Face_Chimera(number_of_ghost_cells,errors_pointer_array,errors_pointer_array);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
            const RIGID_GRID<T_GRID>& rigid_grid=rigid_grid_collection.Rigid_Grid(grid_index);const T_GRID& grid=rigid_grid.grid;
            for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                TV_INT index=iterator.Face_Index();int axis=iterator.Axis();TV location=iterator.Location();
                TV X_new_in_old=rigid_grid.Previous_Object_Space_Point(location);
                TV error_vector=interpolation.Clamped_To_Array_Face(grid,errors(grid_index),X_new_in_old);
                T error_face=rigid_grid.Current_Object_Space_Vector(error_vector)(axis);
                if(!enforce_second_order) 
                    (*Z_array(grid_index))(axis,index)=clamp(Z_forward_array(grid_index)(axis,index)+error_face,Z_min_array(grid_index)(axis,index),Z_max_array(grid_index)(axis,index));
                else
                    (*Z_array(grid_index))(axis,index)=Z_forward_array(grid_index)(axis,index)+error_face;
            }
        }
    }

    void Apply_Reversion_Limiter_Face_Chimera(ARRAY<T_FACE_ARRAYS_SCALAR*>& Z_array,const ARRAY<T_FACE_ARRAYS_SCALAR>& Z_forward_array,const ARRAY<T_FACE_ARRAYS_SCALAR>& Z_backward_array,const ARRAY<T_FACE_ARRAYS_SCALAR>& Z_min_array,const ARRAY<T_FACE_ARRAYS_SCALAR>& Z_max_array)
    {
        ARRAY<T_FACE_ARRAYS_SCALAR> errors(n_local_grids);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
            const T_GRID& grid=rigid_grid_collection.Rigid_Grid(grid_index).grid;
            errors(grid_index).Resize(grid.Domain_Indices(number_of_ghost_cells));
            for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                errors(grid_index)(iterator.Full_Index())=(T).5*((*Z_array(grid_index))(iterator.Full_Index())-Z_backward_array(grid_index)(iterator.Full_Index()));}
        }
        ARRAY<T_FACE_ARRAYS_SCALAR*> errors_pointer_array(n_local_grids);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++) errors_pointer_array(grid_index)=&errors(grid_index);
        example.Fill_Ghost_Cells_Face_Chimera(number_of_ghost_cells,errors_pointer_array,errors_pointer_array);
        for(int grid_index=1;grid_index<=n_local_grids;grid_index++){
            const RIGID_GRID<T_GRID>& rigid_grid=rigid_grid_collection.Rigid_Grid(grid_index);const T_GRID& grid=rigid_grid.grid;
            for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
                TV_INT index=iterator.Face_Index();int axis=iterator.Axis();TV location=iterator.Location();
                TV X_new_in_old=rigid_grid.Previous_Object_Space_Point(location);
                TV error_vector=interpolation.Clamped_To_Array_Face(grid,errors(grid_index),X_new_in_old);
                T error_face=rigid_grid.Current_Object_Space_Vector(error_vector)(axis);
                (*Z_array(grid_index))(axis,index)=Z_forward_array(grid_index)(axis,index)+error_face;
                if(!in_bounds((*Z_array(grid_index))(axis,index),Z_min_array(grid_index)(axis,index),Z_max_array(grid_index)(axis,index)))
                    (*Z_array(grid_index))(axis,index)=Z_forward_array(grid_index)(axis,index);
            }
        }
    }
};
}
#endif
