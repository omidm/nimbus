//#####################################################################
// Copyright 2006-2007, Ron Fedkiw, Geoffrey Irving, Avi Robinson-Mosher, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_OBJECT_FLUID_COLLISIONS
//#####################################################################
#ifndef __DEFORMABLE_OBJECT_FLUID_COLLISIONS__
#define __DEFORMABLE_OBJECT_FLUID_COLLISIONS__

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{

template<class TV>
class DEFORMABLE_OBJECT_FLUID_COLLISIONS:public COLLISION_GEOMETRY<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_BOUNDARY_OBJECT;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT T_VOLUME_OBJECT;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m-1>::SIMPLEX T_SIMPLEX;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX T_VOLUME_SIMPLEX;
    typedef typename IF<TV::m==2,T,typename IF<TV::m==1,ONE,VECTOR<T,3> >::TYPE>::TYPE T_WEIGHTS;
    typedef typename MESH_POLICY<TV::m-1>::MESH T_MESH;

public:
    typedef COLLISION_GEOMETRY<TV> BASE;
    using BASE::collision_thickness;

    T_BOUNDARY_OBJECT& object;
    T_VOLUME_OBJECT* volume_object;
    GEOMETRY_PARTICLES<TV>& particles;

    ARRAY<PAIR<GEOMETRY_PARTICLES<TV>*,T> > saved_states; // TODO: move this into DEFORMABLE_OBJECT and hold a reference here instead

    DEFORMABLE_OBJECT_FLUID_COLLISIONS(T_BOUNDARY_OBJECT& object_input)
        :object(object_input),volume_object(0),particles(dynamic_cast<GEOMETRY_PARTICLES<TV>&>(object.particles))
    {}

    DEFORMABLE_OBJECT_FLUID_COLLISIONS(T_VOLUME_OBJECT& volume_object_input)
        :object(volume_object_input.Get_Boundary_Object()),volume_object(&volume_object_input),particles(dynamic_cast<GEOMETRY_PARTICLES<TV>&>(object.particles))
    {}

    ~DEFORMABLE_OBJECT_FLUID_COLLISIONS()
    {for(int i=1;i<=saved_states.m;i++) delete saved_states(i).x;}

    TV Pointwise_Node_Pseudo_Velocity(const int node_id,const int state1,const int state2)
    {TV dX=saved_states(state2).x->X(node_id)-saved_states(state1).x->X(node_id);
    return dX/(saved_states(state2).y-saved_states(state1).y);}

    bool Has_Volumetric_Geometry() const
    {return volume_object!=0;}

    virtual void Save_State(PAIR<GEOMETRY_PARTICLES<TV>*,T>& state,const T time=0) const
    {state.x->X=object.particles.X;state.x->V=object.particles.V;state.y=time;}

    virtual void Restore_State(const PAIR<GEOMETRY_PARTICLES<TV>*,T>& state)
    {object.particles.X=state.x->X;
    object.particles.V=state.x->V;}

    void Save_State(const int state_index,const T time=0) PHYSBAM_OVERRIDE
    {if(saved_states.m<state_index) saved_states.Resize(state_index);delete saved_states(state_index).x;
    saved_states(state_index).x=new GEOMETRY_PARTICLES<TV>;saved_states(state_index).x->Store_Velocity();
    saved_states(state_index).x->array_collection->Add_Elements(object.particles.array_collection->Size());Save_State(saved_states(state_index),time);}

    void Restore_State(const int state_index) PHYSBAM_OVERRIDE
    {assert(saved_states(state_index).x);Restore_State(saved_states(state_index));}

    void Average_States(const int state1,const int state2,const int result_state,const T interpolation_distance) PHYSBAM_OVERRIDE
    {if(saved_states.m<result_state) saved_states.Resize(result_state);delete saved_states(result_state).x;
    saved_states(result_state).x=new GEOMETRY_PARTICLES<TV>;saved_states(result_state).x->Store_Velocity();
    saved_states(result_state).x->array_collection->Add_Elements(object.particles.array_collection->Size());
    saved_states(result_state).x->X=((T)1-interpolation_distance)*saved_states(state1).x->X+interpolation_distance*saved_states(state2).x->X;
    saved_states(result_state).x->V=((T)1-interpolation_distance)*saved_states(state1).x->V+interpolation_distance*saved_states(state2).x->V;}

    void Delete_State(const int state_index) PHYSBAM_OVERRIDE
    {delete saved_states(state_index).x;saved_states(state_index).x=0;}

    void Read_State(TYPED_ISTREAM& input,const int state_index) PHYSBAM_OVERRIDE
    {if(saved_states.m<state_index) saved_states.Resize(state_index);
    Read_Binary(input,saved_states(state_index));}

    void Write_State(TYPED_OSTREAM& output,const int state_index) const PHYSBAM_OVERRIDE
    {Write_Binary(output,saved_states(state_index));}

private:
    TV Pointwise_Object_Velocity(const TV& location) const PHYSBAM_OVERRIDE
    {PHYSBAM_FATAL_ERROR("need to call version taking simplex_index");}
public:

//#####################################################################
    void Update_Bounding_Box() PHYSBAM_OVERRIDE;
    const RANGE<TV>& Axis_Aligned_Bounding_Box() const PHYSBAM_OVERRIDE;
    TV Pointwise_Object_Velocity(const int simplex_index,const TV& location) const PHYSBAM_OVERRIDE;
    bool Simplex_Intersection(RAY<TV>& ray) const PHYSBAM_OVERRIDE;
    bool Simplex_Closest_Non_Intersecting_Point(RAY<TV>& ray) const PHYSBAM_OVERRIDE;
    bool Inside_Any_Simplex(const TV& location,int& simplex_id) const PHYSBAM_OVERRIDE;
    bool Inside(const TV& location,const T thickness_over_two) const;
    bool Implicit_Geometry_Lazy_Inside(const TV& location,T contour_value=0) const PHYSBAM_OVERRIDE;
    T Implicit_Geometry_Extended_Value(const TV& location) const;
    TV Simplex_Closest_Point_On_Boundary(const TV& location,const T max_distance,const T thickness_over_2=(T)-1,int* simplex_id=0,T* distance=0) const PHYSBAM_OVERRIDE;
    TV Simplex_World_Space_Point_From_Barycentric_Coordinates(const int simplex_id,const T_WEIGHTS& weights) const PHYSBAM_OVERRIDE;
    int Number_Of_Simplices() const PHYSBAM_OVERRIDE;
    void Initialize_For_Thin_Shells_Fluid_Coupling();
    TV Pointwise_Object_Pseudo_Velocity(const int simplex_index,const TV& location,const int state1,const int state2) const PHYSBAM_OVERRIDE;
    POINT_SIMPLEX_COLLISION_TYPE Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,T_WEIGHTS& weights,const int simplex) const;
    T_SIMPLEX World_Space_Simplex(const int simplex_id,const bool use_saved_state=false) const PHYSBAM_OVERRIDE;
    T_SIMPLEX World_Space_Simplex(const int simplex_id,const GEOMETRY_PARTICLES<TV>& state) const;
    bool Earliest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,T_WEIGHTS& weights,int& simplex_id) const PHYSBAM_OVERRIDE;
    bool Latest_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt,T& hit_time,T_WEIGHTS& weights,int& simplex_id,POINT_SIMPLEX_COLLISION_TYPE& returned_collision_type) const PHYSBAM_OVERRIDE;
    bool Any_Simplex_Crossover(const TV& start_X,const TV& end_X,const T dt) const PHYSBAM_OVERRIDE;
    void Get_Simplex_Bounding_Boxes(ARRAY<RANGE<TV> >& bounding_boxes,const bool with_body_motion,const T extra_thickness,const T body_thickness_factor) const PHYSBAM_OVERRIDE;
    void Update_Intersection_Acceleration_Structures(const bool use_swept_simplex_hierarchy,const int state1=0,const int state2=0) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
