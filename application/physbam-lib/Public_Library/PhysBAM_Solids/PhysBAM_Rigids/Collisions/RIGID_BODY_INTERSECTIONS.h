//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Craig Schroeder, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_INTERSECTIONS
//##################################################################### 
#ifndef __RIGID_BODY_INTERSECTIONS__
#define __RIGID_BODY_INTERSECTIONS__

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY;
template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> struct RIGID_BODY_PARTICLE_INTERSECTION;

template<class TV>
class RIGID_BODY_INTERSECTIONS
{
    typedef typename TV::SCALAR T;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::ORIENTED_BOX T_ORIENTED_BOX;
    enum WORKAROUND {d=TV::m};
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d-1>::HIERARCHY T_SIMPLEX_HIERARCHY;
public:
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    bool use_particle_partition,use_particle_partition_center_phi_test;
    VECTOR<int,d> particle_partition_size;
    bool use_particle_hierarchy;
    bool use_triangle_hierarchy,use_triangle_hierarchy_center_phi_test;
    bool use_edge_intersection;
private:
    mutable MATRIX<T,d> rotation,rotation_reverse; // this is not optimal for 1D and 2D, but simpler to be consistent with 3D
    mutable VECTOR<T,d> translation,translation_reverse;
    mutable HASHTABLE<const GEOMETRY_PARTICLES<TV>*,const PARTICLE_HIERARCHY<TV>*> particle_hierarchies;
public:

    RIGID_BODY_INTERSECTIONS(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
        :rigid_body_collection(rigid_body_collection_input),use_particle_partition(false),use_particle_partition_center_phi_test(false),
        use_particle_hierarchy(false),use_triangle_hierarchy(false),use_triangle_hierarchy_center_phi_test(false),
        use_edge_intersection(false)
    {}

    ~RIGID_BODY_INTERSECTIONS();

    void Use_Particle_Partition(const bool use_particle_partition_input=true,const VECTOR<int,d>& particle_partition_size_input=(VECTOR<int,d>()))
    {use_particle_partition=use_particle_partition_input;particle_partition_size=particle_partition_size_input;}

    void Use_Particle_Partition_Center_Phi_Test(const bool use_particle_partition_center_phi_test_input=true)
    {use_particle_partition_center_phi_test=use_particle_partition_center_phi_test_input;}

    void Use_Particle_Hierarchy(const bool use_particle_hierarchy_input=true)
    {use_particle_hierarchy=use_particle_hierarchy_input;}

    void Use_Triangle_Hierarchy(const bool use_triangle_hierarchy_input=true)
    {use_triangle_hierarchy=use_triangle_hierarchy_input;}
    
    void Use_Triangle_Hierarchy_Center_Phi_Test(const bool use_triangle_hierarchy_center_phi_test_input=true)
    {use_triangle_hierarchy_center_phi_test=use_triangle_hierarchy_center_phi_test_input;}
    
    void Use_Edge_Intersection(const bool use_edge_intersection_input=true)
    {use_edge_intersection=use_edge_intersection_input;}

private:
    void Flip_Transformation() const // then from body2 to body1
    {exchange(rotation,rotation_reverse);exchange(translation,translation_reverse);}

public:
    TV Transform_From_Body1_To_Body2_Coordinates(const TV& body1_location) const
    {return rotation*body1_location+translation;}

//#####################################################################
    bool Intersection_Check(const int id_1,const int id_2,int& particle_body,int& levelset_body,const T thickness=0) const;
    bool Bounding_Boxes_Intersect(const int id_1,const int id_2,const T thickness=0) const;
    bool Find_Any_Intersection(const int id_1,const int id_2,int& particle_body,int& levelset_body) const;
    void Append_All_Intersections(const int id_1,const int id_2,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const T contour_value=0) const;
private:
    const PARTICLE_PARTITION<TV>& Particle_Partition(const RIGID_BODY<TV>& rigid_body) const;
    const PARTICLE_HIERARCHY<TV>& Particle_Hierarchy(const RIGID_BODY<TV>& rigid_body) const;
    const T_SIMPLEX_HIERARCHY& Simplex_Hierarchy(const RIGID_BODY<TV>& rigid_body) const;
    void Particles_In_Levelset(const int particle_body_id,const int levelset_body_id,
        ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const T contour_value=0,const bool exit_early=false) const;
    T_ORIENTED_BOX Oriented_Box2_In_Body1_Coordinates(const RIGID_BODY<TV>& body1,const RIGID_BODY<TV>& body2) const;
    void Initialize_Transformation_From_Body1_To_Body2_Coordinates(const RIGID_BODY<TV>& body1,const RIGID_BODY<TV>& body2) const;
//#####################################################################
};
}
#endif
