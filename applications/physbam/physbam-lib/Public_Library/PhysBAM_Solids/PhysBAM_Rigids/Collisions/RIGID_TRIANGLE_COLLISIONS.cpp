//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_TRIANGLE_COLLISIONS  
//##################################################################### 
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Matrices/MATRIX_0X0.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/ROTATION.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_STRUCTURE_INTERACTION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_TRIANGLE_COLLISIONS<TV>::
RIGID_TRIANGLE_COLLISIONS(RIGID_TRIANGLE_COLLISIONS_GEOMETRY<TV>& geometry)
    :geometry(geometry),swap_index(0),mpi_solids(0)
{
    // set parameters 
    Set_Collision_Thickness();
    // set checking
    Compute_Point_Face_Collisions();Compute_Edge_Edge_Collisions();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_TRIANGLE_COLLISIONS<TV>::
~RIGID_TRIANGLE_COLLISIONS()
{}
//#####################################################################
// Function Update_Swept_Hierachies
//#####################################################################
template<class TV> void RIGID_TRIANGLE_COLLISIONS<TV>::
Update_Swept_Hierachies()
{
    // Update swept hierarchies
    LOG::SCOPE scope("updating swept hierarchies","updating swept hierarchies");
    for(int structure_index=1;structure_index<=geometry.structure_geometries.m;structure_index++){
        RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure=*geometry.structure_geometries(structure_index);
        if(d==3 && structure.triangulated_surface){
            BOX_HIERARCHY<TV>& hierarchy=structure.Face_Hierarchy();
            for(int kk=1;kk<=structure.triangulated_surface->mesh.elements.m;kk++){
                const VECTOR<int,3>& nodes=structure.triangulated_surface->mesh.elements(kk);
                hierarchy.box_hierarchy(kk)=RANGE<TV>::Combine(RANGE<TV>::Bounding_Box(structure.X.Subset(nodes)),RANGE<TV>::Bounding_Box(structure.X_self_collision_free.Subset(nodes)));}
            hierarchy.Update_Nonleaf_Boxes();}

        if(structure.segmented_curve){
            SEGMENT_HIERARCHY<TV>& hierarchy=*structure.segmented_curve->hierarchy;
            for(int kk=1;kk<=structure.segmented_curve->mesh.elements.m;kk++){
                const VECTOR<int,2>& nodes=structure.segmented_curve->mesh.elements(kk);
                hierarchy.box_hierarchy(kk)=RANGE<TV>::Combine(RANGE<TV>::Bounding_Box(structure.X.Subset(nodes)),RANGE<TV>::Bounding_Box(structure.X_self_collision_free.Subset(nodes)));}
            hierarchy.Update_Nonleaf_Boxes();}

        {PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > >& hierarchy=structure.particle_hierarchy;
            for(int kk=1;kk<=hierarchy.leaves;kk++){
                int p=structure.collision_particles.active_indices(kk);
                hierarchy.box_hierarchy(kk)=RANGE<TV>::Bounding_Box(structure.X(p),structure.X_self_collision_free(p));}
            hierarchy.Update_Nonleaf_Boxes();}}
}
//#####################################################################
// Function Update_Local_Hierachies
//#####################################################################
template<class TV> void RIGID_TRIANGLE_COLLISIONS<TV>::
Update_Local_Hierachies()
{
    LOG::SCOPE scope("updating local hierarchies","updating local hierarchies");
    for(int structure_index=1;structure_index<=geometry.structure_geometries.m;structure_index++){
        RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure=*geometry.structure_geometries(structure_index);
        if(d==3 && structure.triangulated_surface){
            BOX_HIERARCHY<TV>& hierarchy=structure.Face_Hierarchy();
            for(int kk=1;kk<=structure.triangulated_surface->mesh.elements.m;kk++){
                const VECTOR<int,3>& nodes=structure.triangulated_surface->mesh.elements(kk);
                hierarchy.box_hierarchy(kk)=RANGE<TV>::Bounding_Box(structure.full_particles.X.Subset(nodes));}
            hierarchy.Update_Nonleaf_Boxes();}

        if(structure.segmented_curve){
            SEGMENT_HIERARCHY<TV>& hierarchy=*structure.segmented_curve->hierarchy;
            for(int kk=1;kk<=structure.segmented_curve->mesh.elements.m;kk++){
                const VECTOR<int,2>& nodes=structure.segmented_curve->mesh.elements(kk);
                hierarchy.box_hierarchy(kk)=RANGE<TV>::Bounding_Box(structure.full_particles.X.Subset(nodes));}
            hierarchy.Update_Nonleaf_Boxes();}

        {PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > >& hierarchy=structure.particle_hierarchy;
            for(int kk=1;kk<=hierarchy.leaves;kk++){
                int p=structure.collision_particles.active_indices(kk);
                hierarchy.box_hierarchy(kk)=RANGE<TV>(structure.full_particles.X(p));}
            hierarchy.Update_Nonleaf_Boxes();}}
}
//#####################################################################
// Function Compute_Pairs
//#####################################################################
template<class TV> void RIGID_TRIANGLE_COLLISIONS<TV>::
Compute_Pairs(const T detection_thickness,VECTOR<FRAME<TV>,2>* frame1,VECTOR<FRAME<TV>,2>* frame2)
{
    Compute_Pairs(detection_thickness,FRAME<TV>(),FRAME<TV>(),frame1,frame2);
}
//#####################################################################
// Function Compute_Pairs
//#####################################################################
template<class TV> void RIGID_TRIANGLE_COLLISIONS<TV>::
Compute_Pairs(const T detection_thickness,const FRAME<TV>& gframe1,const FRAME<TV>& gframe2,VECTOR<FRAME<TV>,2>* frame1,VECTOR<FRAME<TV>,2>* frame2)
{
    // Compute pairs
    assert((frame1 && frame2) || (!frame1 && !frame2));
    point_face_pairs_internal.Remove_All();edge_edge_pairs_internal.Remove_All();
    point_face_pairs_external.Remove_All();edge_edge_pairs_external.Remove_All();
    for(int pair_i=1;pair_i<=geometry.interacting_structure_pairs.m;pair_i++){VECTOR<int,2>& pair=geometry.interacting_structure_pairs(pair_i);
        if(compute_point_face_collisions){
            for(int i=1;i<=2;i++){if(i==2 && pair[1]==pair[2]) break;
                RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1=*geometry.structure_geometries(pair[i]);
                RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2=*geometry.structure_geometries(pair[3-i]);
                VECTOR<FRAME<TV>,2>* frame_1=frame1,*frame_2=frame2;if(i>1){frame_1=frame2;frame_2=frame1;}
                Get_Moving_Faces_Near_Moving_Points(structure_1,structure_2,point_face_pairs_internal,point_face_pairs_external,detection_thickness,gframe1,gframe2,frame_1,frame_2);
                if(i==1) swap_index=point_face_pairs_internal.m;}}
        if(compute_edge_edge_collisions){
            RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1=*geometry.structure_geometries(pair[1]);
            RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2=*geometry.structure_geometries(pair[2]);
            Get_Moving_Edges_Near_Moving_Edges(structure_1,structure_2,edge_edge_pairs_internal,edge_edge_pairs_external,detection_thickness,gframe1,gframe2,frame1,frame2);}}
}
//#####################################################################
// Function Get_Moving_Faces_Near_Moving_Points
//#####################################################################
template<class TV> void RIGID_TRIANGLE_COLLISIONS<TV>::
Get_Moving_Faces_Near_Moving_Points(RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1,RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2,ARRAY<VECTOR<int,d+1> >& pairs_internal,ARRAY<VECTOR<int,d+1> >& pairs_external,const T detection_thickness,const FRAME<TV>& gframe1,const FRAME<TV>& gframe2,VECTOR<FRAME<TV>,2>* frame1,VECTOR<FRAME<TV>,2>* frame2)
{
    if(!structure_2.Face_Mesh_Object()) return;

    int old_total=pairs_internal.m;
    RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<TV> visitor(pairs_internal,pairs_external,structure_1,structure_2,geometry,detection_thickness,mpi_solids);
    if(mpi_solids){
        BOX_VISITOR_MPI<RIGID_TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<TV> > mpi_visitor(visitor,structure_1.point_processor_masks,structure_2.Face_Processor_Masks());
        structure_1.particle_hierarchy.Intersection_List(structure_2.Face_Hierarchy(),mpi_visitor,detection_thickness);}
    else if(frame1 && frame2) structure_1.particle_hierarchy.Swept_Intersection_List(*frame1,*frame2,structure_2.Face_Hierarchy(),visitor,detection_thickness);
    else structure_1.particle_hierarchy.Intersection_List(gframe1,gframe2,structure_2.Face_Hierarchy(),visitor,detection_thickness);

    if(geometry.output_number_checked && pairs_internal.m-old_total>0) LOG::Stat("checked point face collisions",pairs_internal.m-old_total);
}
template<> void RIGID_TRIANGLE_COLLISIONS<VECTOR<float,1> >::Get_Moving_Faces_Near_Moving_Points(RIGID_STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,RIGID_STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY<VECTOR<int,d+1> >&,ARRAY<VECTOR<int,d+1> >&,const T,const FRAME<VECTOR<float,1> >&,const FRAME<VECTOR<float,1> >&,VECTOR<FRAME<VECTOR<float,1> >,2>*,VECTOR<FRAME<VECTOR<float,1> >,2>*){PHYSBAM_NOT_IMPLEMENTED();}
template<> void RIGID_TRIANGLE_COLLISIONS<VECTOR<double,1> >::Get_Moving_Faces_Near_Moving_Points(RIGID_STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,RIGID_STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY<VECTOR<int,d+1> >&,ARRAY<VECTOR<int,d+1> >&,const T,const FRAME<VECTOR<double,1> >&,const FRAME<VECTOR<double,1> >&,VECTOR<FRAME<VECTOR<double,1> >,2>*,VECTOR<FRAME<VECTOR<double,1> >,2>*){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Get_Moving_Edges_Near_Moving_Edges
//#####################################################################
template<class TV> void RIGID_TRIANGLE_COLLISIONS<TV>::
Get_Moving_Edges_Near_Moving_Edges(RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1,RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2,ARRAY<VECTOR<int,2*d-2> >& pairs_internal,ARRAY<VECTOR<int,2*d-2> >& pairs_external,const T detection_thickness,const FRAME<TV>& gframe1,const FRAME<TV>& gframe2,VECTOR<FRAME<TV>,2>* frame1,VECTOR<FRAME<TV>,2>* frame2)
{
//    LOG::SCOPE scope("computing edge edge collision pairs", "computing edge edge collision pairs");
    if(!structure_1.Has_Edges() || !structure_2.Has_Edges()) return;

    int old_total=pairs_internal.m;
    RIGID_TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<TV> visitor(pairs_internal,pairs_external,structure_1,structure_2,geometry,detection_thickness,mpi_solids);
    if(mpi_solids){
        BOX_VISITOR_MPI<RIGID_TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<TV> > mpi_visitor(visitor,structure_1.Edge_Processor_Masks(),structure_2.Edge_Processor_Masks());
        structure_1.Edge_Hierarchy().Intersection_List(structure_2.Edge_Hierarchy(),mpi_visitor,detection_thickness);}
    else if(frame1 && frame2) structure_1.Edge_Hierarchy().Swept_Intersection_List(*frame1,*frame2,structure_2.Edge_Hierarchy(),visitor,detection_thickness);
    else structure_1.Edge_Hierarchy().Intersection_List(gframe1,gframe2,structure_2.Edge_Hierarchy(),visitor,detection_thickness);

    if(geometry.output_number_checked && pairs_internal.m-old_total>0) LOG::Stat("checked edge edge collisions",pairs_internal.m-old_total);
}
template<> void RIGID_TRIANGLE_COLLISIONS<VECTOR<float,1> >::Get_Moving_Edges_Near_Moving_Edges(RIGID_STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,RIGID_STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY<VECTOR<int,2*d-2> >&,ARRAY<VECTOR<int,2*d-2> >&,const T,const FRAME<VECTOR<float,1> >&,const FRAME<VECTOR<float,1> >&,VECTOR<FRAME<VECTOR<float,1> >,2>*,VECTOR<FRAME<VECTOR<float,1> >,2>*){PHYSBAM_NOT_IMPLEMENTED();}
template<> void RIGID_TRIANGLE_COLLISIONS<VECTOR<double,1> >::Get_Moving_Edges_Near_Moving_Edges(RIGID_STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,RIGID_STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY<VECTOR<int,2*d-2> >&,ARRAY<VECTOR<int,2*d-2> >&,const T,const FRAME<VECTOR<double,1> >&,const FRAME<VECTOR<double,1> >&,VECTOR<FRAME<VECTOR<double,1> >,2>*,VECTOR<FRAME<VECTOR<double,1> >,2>*){PHYSBAM_NOT_IMPLEMENTED();}
//####################################################################
template class RIGID_TRIANGLE_COLLISIONS<VECTOR<float,1> >;
template class RIGID_TRIANGLE_COLLISIONS<VECTOR<float,2> >;
template class RIGID_TRIANGLE_COLLISIONS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_TRIANGLE_COLLISIONS<VECTOR<double,1> >;
template class RIGID_TRIANGLE_COLLISIONS<VECTOR<double,2> >;
template class RIGID_TRIANGLE_COLLISIONS<VECTOR<double,3> >;
#endif
}
