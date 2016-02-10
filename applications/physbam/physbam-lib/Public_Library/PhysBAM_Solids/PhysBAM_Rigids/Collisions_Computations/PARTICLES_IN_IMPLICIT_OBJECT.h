//#####################################################################
// Copyright 2010, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace PARTICLES_IN_IMPLICIT_OBJECT
//##################################################################### 
#ifndef __PARTICLES_IN_IMPLICIT_OBJECT__
#define __PARTICLES_IN_IMPLICIT_OBJECT__

#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_POLICY.h>

namespace PhysBAM{

template<class TV> class RIGID_BODY;
template<class TV> struct RIGID_BODY_PARTICLE_INTERSECTION;
template<class T,class ID> class ARRAY_VIEW;
template<class T,class ID> class ARRAY;
template<class TV> class IMPLICIT_OBJECT;
template<class T,int d> class VECTOR;
template<class T,int m_input,int n_input> class MATRIX;
template<class K,class T> class HASHTABLE;
template<class TV> class GEOMETRY_PARTICLES;
template<class TV,class T_ARRAY> class PARTICLE_HIERARCHY;

namespace PARTICLES_IN_IMPLICIT_OBJECT
{
template<class TV>
void Append_All_Intersections(RIGID_BODY<TV>& body1, RIGID_BODY<TV>& body2,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool use_triangle_hierarchy,const bool use_edge_intersection,const bool use_triangle_hierarchy_center_phi_test);
template<class TV>
void Append_All_Intersections_Points(RIGID_BODY<TV>& body1, RIGID_BODY<TV>& body2,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value);
template<class TV>
void Append_All_Intersections_Triangles(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool use_triangle_hierarchy_center_phi_test);
template<class TV>
void Append_All_Intersections_Edges(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool use_triangle_hierarchy_center_phi_test);
template<class TV>
void Get_Interfering_Simplices(const RIGID_BODY<TV>& body1,const RIGID_BODY<TV>& body2,ARRAY<int>& simplex_list,MATRIX<typename TV::SCALAR,TV::dimension>& rotation,TV& translation,const bool use_triangle_hierarchy_center_phi_test);
template<class TV>
const typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,VECTOR_POLICY<TV>::DIMMINUSONE>::HIERARCHY& Simplex_Hierarchy(const RIGID_BODY<TV>& rigid_body);
template<class TV>
void Intersections_Using_Hierarchy(RIGID_BODY<TV>& particle_body,RIGID_BODY<TV>& levelset_body,ARRAY<int>& simplex_list,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool exit_early,MATRIX<typename TV::SCALAR,TV::dimension>& rotation,TV& translation);
template<class TV>
void Intersections_Using_Hierarchy_And_Edges(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,ARRAY<int>& simplex_list1,ARRAY<int>& simplex_list2,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool exit_early,MATRIX<typename TV::SCALAR,TV::dimension>& rotation,TV& translation);
template<class TV>
void Particles_In_Implicit_Object(RIGID_BODY<TV>& particle_body,RIGID_BODY<TV>& levelset_body,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool exit_early);
template<class TV>
void Particles_In_Implicit_Object_Hierarchy(RIGID_BODY<TV>& particle_body,RIGID_BODY<TV>& levelset_body,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,HASHTABLE<const GEOMETRY_PARTICLES<TV>*,const PARTICLE_HIERARCHY<TV>*>& particle_hierarchies);
template<class TV>
void Particles_In_Implicit_Object_Partition(RIGID_BODY<TV>& particle_body,RIGID_BODY<TV>& levelset_body,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const typename TV::SCALAR contour_value,const bool use_particle_partition_center_phi_test,const VECTOR<int,TV::dimension>& particle_partition_size,const bool exit_early);
}
}
#endif
