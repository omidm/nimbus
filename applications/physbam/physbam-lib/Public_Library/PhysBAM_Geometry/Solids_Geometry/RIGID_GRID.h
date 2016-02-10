//#####################################################################
// Copyright 2011, Mridul Aanjaneya, Linhai Qiu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_GRID
//#####################################################################
#ifndef __RIGID_GRID__
#define __RIGID_GRID__

#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Utilities/Find_Type.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GRID_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_STATE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>

namespace PhysBAM{

template<class T_GRID>
class RIGID_GRID
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SPIN T_SPIN;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::ORIENTED_BOX T_ORIENTED_BOX;

public:
    typedef TV VECTOR_T;

    T_GRID grid;
    RIGID_GRID_COLLECTION<T_GRID>& rigid_grid_collection;
    int index;
    bool is_static; // not saved to file - indicates whether this object is static in the scene

    T_ORIENTED_BOX oriented_box;
    RANGE<TV> axis_aligned_bounding_box;
    bool bounding_box_up_to_date;
    RIGID_GEOMETRY_STATE<TV> previous_state;

    RIGID_GRID(RIGID_GRID_COLLECTION<T_GRID>& rigid_grid_collection_input=*new RIGID_GRID_COLLECTION<T_GRID>(),int input_index=0, RIGID_GEOMETRY_STATE<TV> previous_state_input=RIGID_GEOMETRY_STATE<TV>());
    virtual ~RIGID_GRID();

    RAY<TV> Object_Space_Ray(const RAY<TV>& world_space_ray) const
    {RAY<TV> transformed_ray(Object_Space_Point(world_space_ray.endpoint),Object_Space_Vector(world_space_ray.direction));
    transformed_ray.semi_infinite=world_space_ray.semi_infinite;transformed_ray.t_max=world_space_ray.t_max;
    transformed_ray.aggregate_id=world_space_ray.aggregate_id;
    return transformed_ray;}

    TV Object_Space_Point(const TV& world_space_point) const
    {return Frame().Inverse_Times(world_space_point);}

    TV Object_Space_Vector(const TV& world_space_vector) const
    {return Rotation().Inverse_Rotate(world_space_vector);}

    TV World_Space_Point(const TV& object_space_point) const
    {return Frame()*object_space_point;}

    TV World_Space_Vector(const TV& object_space_vector) const
    {return Rotation().Rotate(object_space_vector);}

    TV& X() PHYSBAM_ALWAYS_INLINE
    {return rigid_grid_collection.particles.X(index);}
    
    const TV& X() const PHYSBAM_ALWAYS_INLINE
    {return rigid_grid_collection.particles.X(index);}

    ROTATION<TV>& Rotation() PHYSBAM_ALWAYS_INLINE
    {return rigid_grid_collection.particles.rotation(index);}
    
    const ROTATION<TV>& Rotation() const PHYSBAM_ALWAYS_INLINE
    {return rigid_grid_collection.particles.rotation(index);}

    const FRAME<TV> Frame() const PHYSBAM_ALWAYS_INLINE
    {return FRAME<TV>(X(),Rotation());}

    void Set_Frame(const FRAME<TV>& frame) PHYSBAM_ALWAYS_INLINE
    {X()=frame.t;Rotation()=frame.r;}

    const TWIST<TV> Twist() const PHYSBAM_ALWAYS_INLINE
    {return TWIST<TV>(V(),Angular());}

    void Set_Twist(const TWIST<TV>& twist) PHYSBAM_ALWAYS_INLINE
    {V()=twist.linear;Angular()=twist.angular;}

    TV& V() PHYSBAM_ALWAYS_INLINE
    {return rigid_grid_collection.particles.V(index);}

    const TV& V() const PHYSBAM_ALWAYS_INLINE
    {return rigid_grid_collection.particles.V(index);}

    T_SPIN& Angular() PHYSBAM_ALWAYS_INLINE
    {return rigid_grid_collection.particles.angular(index);}

    const T_SPIN& Angular() const PHYSBAM_ALWAYS_INLINE
    {return rigid_grid_collection.particles.angular(index);}

    static TV Pointwise_Object_Velocity(const TWIST<TV>& twist,const TV& t,const TV& X)
    {return twist.linear+TV::Cross_Product(twist.angular,X-t);}

    static TV Relative_Velocity(const RIGID_GRID<T_GRID>& grid1,const RIGID_GRID<T_GRID>& grid2,const TV& world_point)
    {return grid1.Pointwise_Object_Velocity(world_point)-grid2.Pointwise_Object_Velocity(world_point);}

    static TV Relative_Velocity_At_grid1_Particle(const RIGID_GRID<T_GRID>& grid1,const RIGID_GRID<T_GRID>& grid2,const TV& world_point,const int index)
    {return grid1.Pointwise_Object_Velocity_At_Particle(world_point,index)-grid2.Pointwise_Object_Velocity(world_point);}

    static T_SPIN Relative_Angular_Velocity(const RIGID_GRID<T_GRID>& grid1,const RIGID_GRID<T_GRID>& grid2) // make sure the angular velocities are updated before calling this!
    {return grid1.Twist().angular-grid2.Twist().angular;}

    static TWIST<TV> Relative_Twist(const RIGID_GRID<T_GRID>& grid1,const RIGID_GRID<T_GRID>& grid2,const TV& world_point)
    {return TWIST<TV>(grid1.Pointwise_Object_Velocity(world_point)-grid2.Pointwise_Object_Velocity(world_point),grid1.Twist().angular-grid2.Twist().angular);}

    const T_ORIENTED_BOX& Oriented_Bounding_Box() const
    {assert(bounding_box_up_to_date);return oriented_box;}

    const RANGE<TV>& Axis_Aligned_Bounding_Box() const
    {assert(bounding_box_up_to_date);return axis_aligned_bounding_box;}

    bool Bounding_Boxes_Intersect(const RIGID_GRID<T_GRID>& rigid_grid,const T thickness=0) const
    {if(!thickness) return Axis_Aligned_Bounding_Box().Lazy_Intersection(rigid_grid.Axis_Aligned_Bounding_Box()) // check axis aligned first for speed
        && Oriented_Bounding_Box().Intersection(rigid_grid.Oriented_Bounding_Box());
    return Axis_Aligned_Bounding_Box().Intersection(rigid_grid.Axis_Aligned_Bounding_Box(),thickness)
        && Oriented_Bounding_Box().Thickened(thickness).Intersection(rigid_grid.Oriented_Bounding_Box());} // Thickened is rather expensive; avoid it where possible

    TV Get_Interpolation_Position_In_Object_Space(TV& world_space_point) const
    {return previous_state.Object_Space_Point(world_space_point);}

    TV Previous_Object_Space_Point(TV& current_object_space_point) const
    {return previous_state.Object_Space_Point(World_Space_Point(current_object_space_point));}

    TV Current_Object_Space_Point(TV& previous_object_space_point) const
    {return Object_Space_Point(previous_state.World_Space_Point(previous_object_space_point));}

    TV Previous_Object_Space_Vector(TV& current_object_space_vector) const
    {return previous_state.Object_Space_Vector(World_Space_Vector(current_object_space_vector));}

    TV Current_Object_Space_Vector(TV& previous_object_space_vector) const
    {return Object_Space_Vector(previous_state.World_Space_Vector(previous_object_space_vector));}

    void Move_Object_Space_Grid_Center_To_Origin()
    {X()=(T)0.5*(grid.domain.min_corner+grid.domain.max_corner);grid.domain.max_corner=(T)0.5*(grid.domain.max_corner-grid.domain.min_corner);grid.domain.min_corner=-grid.domain.max_corner;}

    void Update_Frame_From_Twist(const T dt, const T time);
    void Update_Bounding_Box();
    virtual TV Pointwise_Object_Velocity(const TV& X) const;
    virtual TV Pointwise_Object_Velocity_At_Particle(const TV& X,const int index) const;
    void Interpolate_Between_States(const RIGID_GEOMETRY_STATE<TV>& state1,const RIGID_GEOMETRY_STATE<TV>& state2,const T time,RIGID_GEOMETRY_STATE<TV>& interpolated_state);
    void Compute_Velocity_Between_States(const RIGID_GEOMETRY_STATE<TV>& state1,const RIGID_GEOMETRY_STATE<TV>& state2,RIGID_GEOMETRY_STATE<TV>& result_state); 
};
template<class TV,class SCALAR> struct REPLACE_FLOATING_POINT<RIGID_GRID<GRID<TV> >,SCALAR>{typedef RIGID_GRID<GRID<typename REPLACE_FLOATING_POINT<TV,SCALAR>::TYPE> > TYPE;};
//#####################################################################
}
#endif
