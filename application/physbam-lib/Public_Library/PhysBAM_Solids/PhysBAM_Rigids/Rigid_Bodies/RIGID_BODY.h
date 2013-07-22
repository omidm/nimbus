//#####################################################################
// Copyright 2002-2008, Zhaosheng Bao, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Frank Losasso, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY
//#####################################################################
#ifndef __RIGID_BODY__
#define __RIGID_BODY__

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_MASS.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY_PARTICLES;
template<class TV>
class RIGID_BODY:public RIGID_GEOMETRY<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;
    typedef typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR T_INERTIA_TENSOR;
    typedef typename RIGID_BODY_POLICY<TV>::WORLD_SPACE_INERTIA_TENSOR T_WORLD_SPACE_INERTIA_TENSOR;
public:
    typedef RIGID_GEOMETRY<TV> BASE;
    typedef int HAS_TYPED_READ_WRITE;
    using BASE::Rotation;using BASE::X;using BASE::V;using BASE::Angular_Velocity;using BASE::Twist;using BASE::particle_index;using BASE::bounding_box_up_to_date;using BASE::is_static;
    using BASE::simplicial_object;using BASE::implicit_object;using BASE::World_Space_Simplex_Bounding_Box;using BASE::coefficient_of_friction;
    enum WORKAROUND {dimension=TV::m+T_SPIN::m};

    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
protected:
    friend class RIGID_BODY_PARTICLES<TV>; // for structures
    friend class RIGID_BODY_COLLECTION<TV>;
public:
    T coefficient_of_restitution; // not saved to file
    T coefficient_of_rolling_friction; // not saved to file
    bool is_temporarily_static; // not saved to file
    T fracture_threshold;
    bool thin_shell;
    bool CFL_initialized; // true if precalculation is done, false if out of date
    T bounding_box_radius; // needed for CFL calculation

    RIGID_BODY(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,bool create_collision_geometry=false);
    RIGID_BODY(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,bool create_collision_geometry,int index);
    virtual ~RIGID_BODY();

    void Set_Mass(const T mass_input) // rescales moments of inertia
    {T& mass=Mass();Inertia_Tensor()*=mass_input/mass;mass=mass_input;}

    void Set_Rigid_Mass(const RIGID_BODY_MASS<TV>& rigid_mass_input)
    {Inertia_Tensor()=rigid_mass_input.inertia_tensor;Mass()=rigid_mass_input.mass;}

    void Set_Inertia_Tensor(const T_INERTIA_TENSOR& inertia_tensor_input) // already scaled by the mass of the object
    {Inertia_Tensor()=inertia_tensor_input;}

    void Rescale(const T scaling_factor,const bool rescale_mass=true)
    {Inertia_Tensor()*=sqr(scaling_factor);if(rescale_mass) Set_Mass(Mass()*pow(scaling_factor,TV::dimension-thin_shell));}

    T Length_Scale_Squared() const
    {return Inertia_Tensor().Max()/Mass();}

    void Set_Coefficient_Of_Restitution(const T coefficient_input=.5)
    {coefficient_of_restitution=coefficient_input;}

    static T Coefficient_Of_Restitution(const RIGID_BODY<TV>& body1,const RIGID_BODY<TV>& body2)
    {return min(body1.coefficient_of_restitution,body2.coefficient_of_restitution);}

    void Set_Coefficient_Of_Rolling_Friction(const T coefficient_input=.5)
    {coefficient_of_rolling_friction=coefficient_input;}

    static T Coefficient_Of_Rolling_Friction(const RIGID_BODY<TV>& body1,const RIGID_BODY<TV>& body2)
    {return min(body1.coefficient_of_rolling_friction,body2.coefficient_of_rolling_friction);}

    bool& Is_Kinematic() PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particle.kinematic(particle_index);}
    
    const bool& Is_Kinematic() const PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particle.kinematic(particle_index);}

    bool Has_Infinite_Inertia() const
    {return is_static || is_temporarily_static || rigid_body_collection.rigid_body_particle.kinematic(particle_index);}

    bool Is_Simulated() const
    {return !is_static && !rigid_body_collection.rigid_body_particle.kinematic(particle_index);}

    T_WORLD_SPACE_INERTIA_TENSOR World_Space_Inertia_Tensor() const // relative to the center of mass
    {return Rigid_Mass().World_Space_Inertia_Tensor(Rotation());}

    T_WORLD_SPACE_INERTIA_TENSOR World_Space_Inertia_Tensor_Inverse() const // relative to the center of mass
    {return Rigid_Mass().World_Space_Inertia_Tensor_Inverse(Rotation());};

    T_WORLD_SPACE_INERTIA_TENSOR World_Space_Inertia_Tensor(const TV& reference_point) const // relative to a reference point
    {return Rigid_Mass().World_Space_Inertia_Tensor(FRAME<TV>(X(),Rotation()),reference_point);}

    T_SPIN World_Space_Inertia_Tensor_Times(const T_SPIN& angular_velocity) const
    {return Rigid_Mass().World_Space_Inertia_Tensor_Times(Rotation(),angular_velocity);}

    T_SPIN World_Space_Inertia_Tensor_Inverse_Times(const T_SPIN& angular_momentum) const
    {return Rigid_Mass().World_Space_Inertia_Tensor_Inverse_Times(Rotation(),angular_momentum);}

    void Update_Angular_Velocity() // needs to be called to keep the angular velocity valid
    {Angular_Velocity()=World_Space_Inertia_Tensor_Inverse_Times(Angular_Momentum());}

    void Update_Angular_Velocity(RIGID_BODY_STATE<TV>& state) // needs to be called to keep the angular velocity valid
    {state.Update_Angular_Velocity(Inertia_Tensor());}

    void Update_Angular_Momentum() // assumes a valid angular_velocity
    {Angular_Momentum()=World_Space_Inertia_Tensor_Times(Twist().angular);}

    void Update_Angular_Momentum(RIGID_BODY_STATE<TV>& state) // assumes a valid angular_velocity
    {state.Update_Angular_Momentum(Inertia_Tensor());}

    T_SPIN& Angular_Momentum() PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particle.angular_momentum(particle_index);}

    const T_SPIN& Angular_Momentum() const PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particle.angular_momentum(particle_index);}

    T& Mass() PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particle.mass(particle_index);}

    const T& Mass() const PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particle.mass(particle_index);}

    T_INERTIA_TENSOR& Inertia_Tensor() PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particle.inertia_tensor(particle_index);}

    const T_INERTIA_TENSOR& Inertia_Tensor() const PHYSBAM_ALWAYS_INLINE
    {return rigid_body_collection.rigid_body_particle.inertia_tensor(particle_index);}

    const RIGID_BODY_MASS<TV> Rigid_Mass() const PHYSBAM_ALWAYS_INLINE
    {return RIGID_BODY_MASS<TV>(Mass(),Inertia_Tensor());}

    T Kinetic_Energy() const
    {return Translational_Kinetic_Energy()+Rotational_Kinetic_Energy();}

    T Translational_Kinetic_Energy() const
    {return (T).5*Mass()*Twist().linear.Magnitude_Squared();}

    T Rotational_Kinetic_Energy() const
    {return Rigid_Mass().Rotational_Kinetic_Energy(Rotation(),Angular_Momentum());}

    T Kinetic_Energy(const TWIST<TV>& twist) const
    {if(Has_Infinite_Inertia()) return 0;return (T).5*Mass()*twist.linear.Magnitude_Squared()+(T).5*TV::SPIN::Dot_Product(World_Space_Inertia_Tensor_Times(twist.angular),twist.angular);}

    void Update_Bounding_Box()
    {if(bounding_box_up_to_date && is_static) return;
    BASE::Update_Bounding_Box();}

private:
    static MATRIX<T,1> Conjugate_With_Cross_Product_Matrix(const VECTOR<T,1>& vector,const MATRIX<T,0>& tensor)
    {return MATRIX<T,1>();}

    template<class T_TENSOR,int d>
    static SYMMETRIC_MATRIX<T,d> Conjugate_With_Cross_Product_Matrix(const VECTOR<T,d>& vector,const T_TENSOR& tensor)
    {return tensor.Conjugate_With_Cross_Product_Matrix(vector);}
public:

    T_SYMMETRIC_MATRIX Impulse_Factor(const TV& location) const
    {if(Has_Infinite_Inertia()) return T_SYMMETRIC_MATRIX(); // return zero matrix
    return Conjugate_With_Cross_Product_Matrix(location-X(),World_Space_Inertia_Tensor_Inverse())+1/Mass();}

    T_SYMMETRIC_MATRIX Object_Space_Impulse_Factor(const TV& object_space_location) const
    {if(Has_Infinite_Inertia()) return T_SYMMETRIC_MATRIX(); // return zero matrix
    return Conjugate_With_Cross_Product_Matrix(object_space_location,Inertia_Tensor().Inverse())+1/Mass();}

    static T_SYMMETRIC_MATRIX Impulse_Factor(const RIGID_BODY<TV>& body1,const RIGID_BODY<TV>& body2,const TV& location)
    {assert(!body1.Has_Infinite_Inertia() || !body2.Has_Infinite_Inertia());
    return body1.Impulse_Factor(location)+body2.Impulse_Factor(location);}

    TV Simplex_Normal(const int id,const RIGID_BODY_STATE<TV>& state) const
    {return state.World_Space_Vector(simplicial_object->Normal(id));}

    void Save_State(RIGID_BODY_STATE<TV>& state,const T time=0) const
    {state.time=time;state.frame.t=X();state.frame.r=Rotation();state.twist=Twist();state.angular_momentum=Angular_Momentum();}

    void Restore_State(const RIGID_BODY_STATE<TV>& state)
    {X()=state.frame.t;Rotation()=state.frame.r;V()=state.twist.linear;Angular_Velocity()=state.twist.angular;Angular_Momentum()=state.angular_momentum;}

    TWIST<TV> Gather(const TWIST<TV>& wrench,const TV& location) const
    {return TWIST<TV>(wrench.linear,wrench.angular+TV::Cross_Product(location-X(),wrench.linear));}

    TWIST<TV> Scatter(const TWIST<TV>& twist,const TV& location) const
    {return TWIST<TV>(twist.linear+TV::Cross_Product(twist.angular,location-X()),twist.angular);}

    TWIST<TV> Inertia_Inverse_Times(const TWIST<TV>& wrench) const
    {if(Has_Infinite_Inertia()) return TWIST<TV>();
    return TWIST<TV>(wrench.linear/Mass(),World_Space_Inertia_Tensor_Inverse_Times(wrench.angular));}

    TWIST<TV> Inertia_Times(const TWIST<TV>& twist) const
    {PHYSBAM_ASSERT(!Has_Infinite_Inertia());
    return TWIST<TV>(twist.linear*Mass(),World_Space_Inertia_Tensor_Times(twist.angular));}

//#####################################################################
    void Effective_Inertia_Inverse(MATRIX<T,dimension>& extended_mass_inverse,const TV& location) const;
    void Effective_Inertia(MATRIX<T,dimension>& extended_mass_inverse,const TV& location) const;
    void Gather_Matrix(MATRIX<T,dimension>& gather,const TV& location) const;
    void Scatter_Matrix(MATRIX<T,dimension>& scatter,const TV& location) const;
    TWIST<TV> Effective_Inertia_Inverse_Times(const TWIST<TV>& wrench,const TV& location) const;
    TWIST<TV> Effective_Inertia_Times(const TWIST<TV>& twist,const TV& location) const;
    T Volume() const;
    static void Print_Pairs(const RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const ARRAY<VECTOR<int,2> >& pairs);
    void Initialize_CFL();
    virtual T CFL(const T max_distance_per_time_step,const T max_rotation_per_time_step,const bool verbose=false);
    void Interpolate_Between_States(const RIGID_BODY_STATE<TV>& state1,const RIGID_BODY_STATE<TV>& state2,const T time,RIGID_BODY_STATE<TV>& interpolated_state);
    void Compute_Velocity_Between_States(const RIGID_BODY_STATE<TV>& state1,const RIGID_BODY_STATE<TV>& state2,RIGID_BODY_STATE<TV>& result_state);
    void Apply_Impulse_To_Body(const TV& location,const TV& impulse,const T_SPIN& angular_impulse=T_SPIN(),const bool half_impulse_for_accumulator=false);
    static void Apply_Impulse(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const TV& location,const TV& impulse,const T_SPIN& angular_impulse=T_SPIN(),
        const bool half_impulse_for_accumulator=false);
    static void Compute_Clamped_Impulse(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const TV& location,TWIST<TV>& impulse,const ROTATION<TV>& saved_rotation_1,
        const ROTATION<TV>& saved_rotation_2);
    static TWIST<TV> Compute_Collision_Impulse(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const ROTATION<TV>& saved_rotation_1,const ROTATION<TV>& saved_rotation_2,const TV& location,
        const TV& normal,const TV& relative_velocity,const T coefficient_of_restitution,const T coefficient_of_friction=0,const bool clamp_friction_magnitude=true,
        const bool rolling_friction=false,const bool clamp_energy=false);
    static void Apply_Collision_Impulse(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const ROTATION<TV>& saved_rotation_1,const ROTATION<TV>& saved_rotation_2,const TV& location,
        const TV& normal,const TV& relative_velocity,const T coefficient_of_restitution,const T coefficient_of_friction=0,const bool clamp_friction_magnitude=true,
        const bool rolling_friction=false,const bool clamp_energy=false,const bool half_impulse_for_accumulator=false);
    static void Apply_Sticking_And_Angular_Sticking_Impulse(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const TV& location,const TWIST<TV>& delta_relative_twist,
        const MATRIX_MXN<T>& angular_constraint_matrix,const MATRIX_MXN<T>& prismatic_constraint_matrix);
    static TWIST<TV> Apply_Rolling_Friction(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const TV& location,const TV& normal,const T normal_impulse);
    static TWIST<TV> Find_Impulse_And_Angular_Impulse(const RIGID_BODY<TV>& body1,const RIGID_BODY<TV>& body2,const TV& location,const TWIST<TV>& delta_rel_twist_at_location);
    static TWIST<TV> Find_Impulse_And_Angular_Impulse(const RIGID_BODY<TV>& body1,const RIGID_BODY<TV>& body2,const TV& location,const TWIST<TV>& delta_rel_twist_at_location,
        const MATRIX_MXN<T>& angular_constraint_matrix,const MATRIX_MXN<T>& prismatic_constraint_matrix);
    static void Apply_Push(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const TV& location,const TV& normal,const T distance,const bool half_impulse_for_accumulator);
    T Volumetric_Density() const;
    void Diagonalize_Inertia_Tensor(const T_WORLD_SPACE_INERTIA_TENSOR& inertia_tensor_at_center_of_mass);
    template<class T2> void Initialize_From_Tetrahedralized_Volume_And_Triangulated_Surface(TETRAHEDRALIZED_VOLUME<T2>& tetrahedralized_volume,
        TRIANGULATED_SURFACE<T>& triangulated_surface,const T cell_size,const int subdivision_loops=0,const bool (*create_levelset_test)(TETRAHEDRALIZED_VOLUME<T>&)=0,
        const bool use_implicit_surface_maker=true,const int levels_of_octree=0,const T shrink_levelset_amount=0);
//#####################################################################
};
}
#endif
