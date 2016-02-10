//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR
//##################################################################### 
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_STRUCTURE_INTERACTION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_TRIANGLE_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Parallel_Computation/MPI_RIGIDS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<TV>::
RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR(ARRAY<VECTOR<int,d+1> >& pairs_internal,ARRAY<VECTOR<int,d+1> >& pairs_external,
    const RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& particle_structure,const RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& face_structure,
    const RIGID_TRIANGLE_COLLISIONS_GEOMETRY<TV>& geometry,const T collision_thickness,MPI_RIGIDS<TV>* mpi_solids)
    :pairs_internal(pairs_internal),pairs_external(pairs_external),particle_active_indices(particle_structure.collision_particles.active_indices),
    faces(face_structure.Face_Mesh_Object()->mesh.elements),particle_structure(particle_structure),face_structure(face_structure),
    collision_thickness(collision_thickness),mpi_solids(mpi_solids)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<TV>::
~RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR()
{}
//#####################################################################
// Function Store
//#####################################################################
template<class TV> void  RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<TV>::
Store(const int point_index,const int face_index)
{
    const VECTOR<int,d>& face_nodes=faces(face_index);int p=particle_active_indices(point_index);
    TV dX_average=(T).5*(particle_structure.X(p)-particle_structure.X_self_collision_free(p)+ARRAYS_COMPUTATIONS::Average(face_structure.X.Subset(face_nodes))-ARRAYS_COMPUTATIONS::Average(face_structure.X_self_collision_free.Subset(face_nodes)));
    RANGE<TV> box1=RANGE<TV>::Bounding_Box(particle_structure.X_self_collision_free(p)+dX_average,particle_structure.X(p));
    RANGE<TV> box2=RANGE<TV>::Bounding_Box(face_structure.X_self_collision_free.Subset(face_nodes))+dX_average;box2.Enlarge_Nonempty_Box_To_Include_Points(face_structure.X.Subset(face_nodes));
    if(!box1.Intersection(box2,collision_thickness)) return;
    VECTOR<int,d+1> nodes=face_nodes.Insert(p,1);
    if (mpi_solids){
        VECTOR<PARTITION_ID,d+1> processors(mpi_solids->partition_id_from_particle_index.Subset(nodes));
        int i; for(i=1;i<=d;i++) if (processors(i)!=processors(d+1)) {pairs_external.Append(nodes);return;}
        pairs_internal.Append(nodes);}
    else pairs_internal.Append(nodes);
}
//####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template struct RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<T,d> >; \
    template void BOX_HIERARCHY<VECTOR<T,d> >::Intersection_List<RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<T,d> >,ZERO>(BOX_HIERARCHY<VECTOR<T,d> > const&, \
        RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<T,d> >&,ZERO) const; \
    template void BOX_HIERARCHY<VECTOR<T,d> >::Intersection_List<BOX_VISITOR_MPI<RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<T,d> > >,ZERO>( \
        BOX_HIERARCHY<VECTOR<T,d> > const&,BOX_VISITOR_MPI<RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<T,d> > >&,ZERO) const; \
    template void BOX_HIERARCHY<VECTOR<T,d> >::Intersection_List<BOX_VISITOR_MPI<RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<T,d> > >,T>(BOX_HIERARCHY<VECTOR<T,d> > const&, \
        BOX_VISITOR_MPI<RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<T,d> > >&,T) const; \
    template void BOX_HIERARCHY<VECTOR<T,d> >::Intersection_List<RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<T,d> >,T>(BOX_HIERARCHY<VECTOR<T,d> > const&, \
        RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<T,d> >&,T) const; \
    template void BOX_HIERARCHY<VECTOR<T,d> >::Swept_Intersection_List<RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<T,d> >,T>(VECTOR<FRAME<VECTOR<T,d> >,2> const&,VECTOR<FRAME<VECTOR<T,d> >,2> const&, \
        BOX_HIERARCHY<VECTOR<T,d> > const&, RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<T,d> >&,T) const; \
    template void BOX_HIERARCHY<VECTOR<T,d> >::Intersection_List<RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<T,d> >&,T>(const FRAME<VECTOR<T,d> >&,const FRAME<VECTOR<T,d> >&, \
        const BOX_HIERARCHY<VECTOR<T,d> >&, RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<T,d> >&,const T) const; \
    template void BOX_HIERARCHY<VECTOR<T,d> >::Intersection_List<RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<T,d> >,T>(const FRAME<VECTOR<T,d> >&,const FRAME<VECTOR<T,d> >&, \
    const BOX_HIERARCHY<VECTOR<T,d> >& other_hierarchy,RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<T,d> >&,const T) const;

template struct RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<float,1> >;
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template struct RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<VECTOR<double,1> >;
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif
