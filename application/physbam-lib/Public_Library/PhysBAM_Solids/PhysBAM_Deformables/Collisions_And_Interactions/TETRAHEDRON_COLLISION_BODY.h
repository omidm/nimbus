//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Andrew Selle, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TETRAHEDRON_COLLISION_BODY
//#####################################################################
#ifndef __TETRAHEDRON_COLLISION_BODY__
#define __TETRAHEDRON_COLLISION_BODY__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
namespace PhysBAM{

template<class T> class TETRAHEDRALIZED_VOLUME;
template<class T> class TRIANGULATED_SURFACE;
template<class TV> class IMPLICIT_OBJECT;

template<class T_input>
class TETRAHEDRON_COLLISION_BODY:public COLLISION_GEOMETRY<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    using COLLISION_GEOMETRY<TV>::Set_Collision_Thickness;using COLLISION_GEOMETRY<TV>::collision_thickness;

    PARTICLES<TV> &particles,&undeformed_particles;
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume;
    TRIANGULATED_SURFACE<T>& undeformed_triangulated_surface;
    TRIANGULATED_SURFACE<T>& triangulated_surface;
    IMPLICIT_OBJECT<TV>& implicit_surface;
    T max_min_barycentric_weight_tolerance;
    T friction_coefficient;
    T relaxation_factor;
    T normal_interpolation_scale_factor;
    T self_collision_normal_angle_tolerance;
    T min_tet_volume_tolerance;

    TETRAHEDRON_COLLISION_BODY(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume_input,TRIANGULATED_SURFACE<T>& undeformed_triangulated_surface_input,
        IMPLICIT_OBJECT<TV>& implicit_object_input,TRIANGULATED_SURFACE<T>* triangulated_surface_input=0);

    TV Depth_Interpolated_Normal(const T unsigned_depth,const TV& outward_direction,const TV& surface_normal) const
    {T lambda=min(unsigned_depth/normal_interpolation_scale_factor,(T)1);return (lambda*outward_direction+(1-lambda)*surface_normal).Normalized();}

    void Set_Min_Tet_Volume_Tolerance(const T min_tet_volume_tolerance_input=-(T)1e-6)
    {min_tet_volume_tolerance=min_tet_volume_tolerance_input;}

    void Set_Normal_Interpolation_Scale_Factor(const T collision_thickness_scaling=(T)2)
    {normal_interpolation_scale_factor=collision_thickness_scaling*collision_thickness;}

    void Set_Self_Collision_Normal_Angle_Tolerance(const T self_collision_normal_angle_tolerance_input=-(T)1e-5)
    {self_collision_normal_angle_tolerance=self_collision_normal_angle_tolerance_input;}

    void Set_Friction_Coefficient(const T friction_coefficient_input=(T).1)
    {friction_coefficient=friction_coefficient_input;}

    void Set_Max_Min_Barycentric_Weight_Tolerance(const T max_min_barycentric_weight_tolerance_input=(T)-1e-6)
    {max_min_barycentric_weight_tolerance=max_min_barycentric_weight_tolerance_input;}

    void Set_Relaxation_Factor(const T relaxation_factor_input=(T).25)
    {relaxation_factor=relaxation_factor_input;}

//#####################################################################
    TV Surface_Normal(const int triangle,const TV& weights) const;
    bool Implicit_Geometry_Lazy_Inside(const TV& location,T contour_value=0) const PHYSBAM_OVERRIDE;
    bool Implicit_Geometry_Lazy_Inside_And_Value(const TV& location,T& phi,T contour_value=0) const PHYSBAM_OVERRIDE;
    bool Implicit_Geometry_Lazy_Outside_Extended_Levelset_And_Value(const TV& location,T& phi_value,T contour_value=0) const;
    TV Implicit_Geometry_Normal(const TV& location,const int aggregate=-1) const PHYSBAM_OVERRIDE;
    TV Implicit_Geometry_Normal(const TV& location,T& phi_value,const int aggregate=-1,const int location_particle_index=0) const PHYSBAM_OVERRIDE;
    TV Implicit_Geometry_Extended_Normal(const TV& location,T& phi_value,const int aggregate=-1,const int location_particle_index=0) const PHYSBAM_OVERRIDE;
    int Get_Tetrahedron_Near_Point(const TV& point,TV& weights,const ARRAY<int>& particles_to_ignore=ARRAY<int>()) const;
    int Get_Surface_Triangle(const int tetrahedron_index,const TV& tetrahedron_weights,TV& surface_weights,const bool omit_outside_points=false,
        const bool omit_inside_points=false,bool* inside=0) const;
private:
    void Adjust_Point_Face_Collision_Position_And_Velocity(const int triangle_index,TV& X,TV& V,SOFT_BINDINGS<TV>& soft_bindings,const T one_over_mass,const T dt,const TV& weights,
        TV& position_change);
public:
    int Adjust_Nodes_For_Collisions(ARRAY_VIEW<const TV> X_old,PARTICLES<TV>& collision_particles,SOFT_BINDINGS<TV>& soft_bindings,const ARRAY<int>& nodes_to_check,
        const ARRAY<bool>& particle_on_surface,const T collision_tolerance,ARRAY<COLLISION_PARTICLE_STATE<TV> >& collision_particle_state,
        ARRAY<COLLISION_GEOMETRY_ID>& particle_to_collision_body_id,const T max_relative_velocity,const T dt,const HASHTABLE<int,T> *friction_table,const HASHTABLE<int,T> *thickness_table);
    const RANGE<TV>& Axis_Aligned_Bounding_Box() const PHYSBAM_OVERRIDE;
    void Update_Bounding_Box() PHYSBAM_OVERRIDE;
    void Read_State(TYPED_ISTREAM& input,const int state_index) PHYSBAM_OVERRIDE;
    void Write_State(TYPED_OSTREAM& output,const int state_index) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
