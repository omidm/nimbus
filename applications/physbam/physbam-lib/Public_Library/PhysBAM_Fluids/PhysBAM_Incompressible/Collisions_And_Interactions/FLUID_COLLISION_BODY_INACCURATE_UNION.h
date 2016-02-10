//#####################################################################
// Copyright 2006, Geoffrey Irving, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_COLLISION_BODY_INACCURATE_UNION
//#####################################################################
#ifndef __FLUID_COLLISION_BODY_INACCURATE_UNION__
#define __FLUID_COLLISION_BODY_INACCURATE_UNION__

#include <PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/GRID_BASED_COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Implicit_Objects_Dyadic/DYADIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects_RLE/RLE_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Level_Sets/IMPLICIT_OBJECT_ON_A_RAY.h>
#include <PhysBAM_Geometry/Level_Sets/IMPLICIT_OBJECT_ON_A_RAY_POLICY.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
namespace PhysBAM{

template<class T_GRID>
class FLUID_COLLISION_BODY_INACCURATE_UNION:public COLLISION_GEOMETRY<typename T_GRID::VECTOR_T>
{
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename TV::SCALAR T;
    typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET_IMPLICIT_OBJECT T_LEVELSET_IMPLICIT_OBJECT;
    typedef typename IF<TV::dimension==2,T,typename IF<TV::dimension==1,ONE,TV>::TYPE>::TYPE T_WEIGHTS;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_T;
    typedef typename T_ARRAYS_T::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_T;
    typedef typename T_FACE_ARRAYS_T::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_T::template REBIND<int>::TYPE T_FACE_ARRAYS_INT;
    typedef typename T_FACE_ARRAYS_T::template REBIND<COLLISION_GEOMETRY_ID>::TYPE T_FACE_ARRAYS_COLLISION_GEOMETRY_ID;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX T_SIMPLEX;
    typedef typename IMPLICIT_OBJECT_ON_A_RAY_POLICY<T_GRID>::TYPE T_IMPLICIT_OBJECT_ON_A_RAY;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_MAC_HELPER T_LINEAR_INTERPOLATION_MAC_HELPER;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::BLOCK T_BLOCK;
    typedef COLLISION_GEOMETRY<TV> BASE;

public:
    GRID_BASED_COLLISION_GEOMETRY<T_GRID> collision_bodies;
    using BASE::collision_geometries_for_rasterization;
    T contour_value;
    T_GRID& grid;
    T_ARRAYS_T phi;
    T_LEVELSET_IMPLICIT_OBJECT levelset;
private:
    T_FACE_ARRAYS_T face_velocities;
    T_FACE_ARRAYS_BOOL face_velocities_set;
    T_LINEAR_INTERPOLATION_SCALAR interpolation;
public:

    FLUID_COLLISION_BODY_INACCURATE_UNION(T_GRID& grid_input);
    FLUID_COLLISION_BODY_INACCURATE_UNION(T_GRID& grid_input,T contour_value_input);
    
    ~FLUID_COLLISION_BODY_INACCURATE_UNION();

    TV Pointwise_Object_Pseudo_Velocity(const int aggregate_id,const TV& X,const int state1,const int state2) const PHYSBAM_OVERRIDE
    {return Pointwise_Object_Velocity(0,X);}

    // TODO: check whether this is still valid
    //      in the low accuracy case, the points inside the object are set to invalid in the driver, so this does not need to invalidate anything
    //      (i.e., it ignored crossovers)
    // If it isn't the case, we need to generalize COLLISION_GEOMETRY_COLLECTION::Latest_Crossover to work for this.
    bool Latest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,T_WEIGHTS& weights,int& simplex_id,POINT_SIMPLEX_COLLISION_TYPE& returned_collision_type) const PHYSBAM_OVERRIDE
    {return false;}

    bool Any_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt) const PHYSBAM_OVERRIDE
    {return levelset(end_X)<=0;}

    bool Get_Body_Penetration(const TV& start_X,const TV& end_X,const T contour_value,const T dt,T& hit_time,int& simplex_id,T& start_phi,T& end_phi,TV& end_body_normal,
        TV& body_velocity) const PHYSBAM_OVERRIDE
    {end_phi=levelset(end_X);if(end_phi>contour_value) return false;
    start_phi=levelset(start_X);end_body_normal=levelset.Normal(end_X);
    body_velocity=Pointwise_Object_Velocity(0,end_X);
    hit_time=dt;return true;}

    int Number_Of_Simplices() const PHYSBAM_OVERRIDE
    {return 0;}

    bool Push_Out_Point(TV& X,const T collision_distance,T& distance) const PHYSBAM_OVERRIDE
    {T current_distance=levelset(X);
    if(current_distance<collision_distance){
        X+=(collision_distance-current_distance)*levelset.Normal(X);
        distance=current_distance;return true;}
    return false;}

    bool Inside_Any_Simplex(const TV& location,int& simplex_id) const PHYSBAM_OVERRIDE
    {return levelset(location)<=contour_value;}

    bool Simplex_Closest_Non_Intersecting_Point(RAY<TV>& ray) const PHYSBAM_OVERRIDE
    {T_IMPLICIT_OBJECT_ON_A_RAY implicit_object_on_a_ray(levelset,ray);
    if(implicit_object_on_a_ray(0)<=0){ray.t_max=0;return true;}
    if(implicit_object_on_a_ray(ray.t_max)>0) return false;
    ITERATIVE_SOLVER<T> solver;solver.tolerance=(T).001*grid.Minimum_Edge_Length();
    ray.t_max=solver.Bisection_Secant_Root(implicit_object_on_a_ray,0,ray.t_max);
    ray.t_max=max((T)0,ray.t_max-(T).01*grid.Minimum_Edge_Length()); // TODO: probably make this shift a parameter, or find an entirely different cleaner way
    return true;}

    bool Simplex_Intersection(RAY<TV>& ray) const PHYSBAM_OVERRIDE
    {T_IMPLICIT_OBJECT_ON_A_RAY implicit_object_on_a_ray(levelset,ray);
    if(implicit_object_on_a_ray(0)<=0){ray.t_max=0;return true;}
    if(implicit_object_on_a_ray(ray.t_max)>0) return false;
    ITERATIVE_SOLVER<T> solver;solver.tolerance=(T).001*grid.Minimum_Edge_Length();
    ray.t_max=solver.Bisection_Secant_Root(implicit_object_on_a_ray,0,ray.t_max);
    return true;}

    TV Simplex_Closest_Point_On_Boundary(const TV& location,const T max_distance,const T thickness_over_2=0,int* simplex_id=0,T* returned_distance=0) const PHYSBAM_OVERRIDE
    {T distance=levelset(location);if(returned_distance) *returned_distance=distance;return location-distance*levelset.Normal(location);}

private:
    template<class T_TAG> static bool Block_Valid(const T_BLOCK& block,T_TAG) // TODO: can't work for dyadic, since BLOCK_DYADIC asserts validity on construction
    {return true;}

#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    static bool Block_Valid(const T_BLOCK& block,RLE_TAG<TV>)
    {return block;}
#endif
public:

    TV Pointwise_Object_Velocity(const int aggregate_id,const TV& X) const PHYSBAM_OVERRIDE
    {T_BLOCK block(grid,X);if(!Block_Valid(block,typename T_GRID::GRID_TAG())) return TV();
    return T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Face_Normalized(block,face_velocities,face_velocities_set,X);}

    TV Pointwise_Object_Velocity(const TV& X) const PHYSBAM_OVERRIDE
    {return Pointwise_Object_Velocity(0,X);}

//#####################################################################
    void Update_Intersection_Acceleration_Structures(const bool use_swept_simplex_hierarchy,const int state1=0,const int state2=0) PHYSBAM_OVERRIDE;
    void Save_State(const int state_index,const T time=0) PHYSBAM_OVERRIDE;
    void Restore_State(const int state_index) PHYSBAM_OVERRIDE;
    void Read_State(TYPED_ISTREAM& input,const int state_index) PHYSBAM_OVERRIDE;
    void Write_State(TYPED_OSTREAM& output,const int state_index) const PHYSBAM_OVERRIDE;
    typename T_GRID::SCALAR Implicit_Geometry_Extended_Value(const TV& location) const PHYSBAM_OVERRIDE;
    void Initialize_Grid_Structures(const T_GRID& grid,OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID>& objects_in_cell,const COLLISION_GEOMETRY_ID id) const;
private:
    typename T_GRID::SCALAR Implicit_Geometry_Extended_Value_Helper(const TV& location,UNIFORM_TAG<TV>) const;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    typename T_GRID::SCALAR Implicit_Geometry_Extended_Value_Helper(const TV& location,DYADIC_TAG<TV>) const;
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    typename T_GRID::SCALAR Implicit_Geometry_Extended_Value_Helper(const TV& location,RLE_TAG<TV>) const;
#endif
    void Initialize_Grid_Structures_Helper(OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID>& objects_in_cell,const COLLISION_GEOMETRY_ID id,UNIFORM_TAG<TV>);
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    void Initialize_Grid_Structures_Helper(OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID>& objects_in_cell,const COLLISION_GEOMETRY_ID id,DYADIC_TAG<TV>);
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    void Initialize_Grid_Structures_Helper(OBJECTS_IN_CELL<T_GRID,COLLISION_GEOMETRY_ID>& objects_in_cell,const COLLISION_GEOMETRY_ID id,RLE_TAG<TV>);
#endif
    void Initialize_Grid_Structures_Subobject(T_FACE_ARRAYS_INT& face_velocities_count,T_FACE_ARRAYS_COLLISION_GEOMETRY_ID& face_operations,const COLLISION_GEOMETRY_ID subobject,UNIFORM_TAG<TV>);
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    void Initialize_Grid_Structures_Subobject(T_FACE_ARRAYS_INT& face_velocities_count,OPERATION_HASH<>& face_operations,const COLLISION_GEOMETRY_ID subobject,DYADIC_TAG<TV>);
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    void Initialize_Grid_Structures_Subobject(T_FACE_ARRAYS_INT& face_velocities_count,OPERATION_HASH<>& face_operations,const COLLISION_GEOMETRY_ID subobject,RLE_TAG<TV>);
#endif
//#####################################################################
};
}
#endif
