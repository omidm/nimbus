//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Nipun Kwatra, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STRUCTURE_INTERACTION_GEOMETRY
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/STRUCTURE_INTERACTION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDING.h>
using namespace PhysBAM;
//#####################################################################
// Function Segmented_Curve
//#####################################################################
template<class TV> typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,1>::OBJECT* STRUCTURE_INTERACTION_GEOMETRY<TV>::
Segmented_Curve(STRUCTURE<TV>* structure)
{
    if(T_SEGMENTED_CURVE* segmented_curve=dynamic_cast<T_SEGMENTED_CURVE*>(structure)) return segmented_curve;
    else return 0;
}
//#####################################################################
// Function Triangulated_Surface
//#####################################################################
template<class TV> TRIANGULATED_SURFACE<typename TV::SCALAR>* STRUCTURE_INTERACTION_GEOMETRY<TV>::
Triangulated_Surface(STRUCTURE<TV>* structure)
{
    if(TV::m==2) return 0;
    if(TRIANGULATED_SURFACE<T>* triangulated_surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(structure))
        return triangulated_surface;
    else if(TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(structure)){
        if(!tetrahedralized_volume->triangulated_surface) tetrahedralized_volume->Initialize_Triangulated_Surface();
        return tetrahedralized_volume->triangulated_surface;} // TODO: This is broken; long term shallow copy of a temporary auxiliary structure
    else if(HEXAHEDRALIZED_VOLUME<T>* hexahedralized_volume=dynamic_cast<HEXAHEDRALIZED_VOLUME<T>*>(structure)){
        if(!hexahedralized_volume->triangulated_surface) hexahedralized_volume->Initialize_Triangulated_Surface();
        return hexahedralized_volume->triangulated_surface;}
    else if(EMBEDDING<VECTOR<T,3> >* embedding=dynamic_cast<EMBEDDING<VECTOR<T,3> >*>(structure))
        return &embedding->material_surface;
    else return 0;
}
//#####################################################################
// Function Build_Collision_Geometry
//#####################################################################
template<class TV> void STRUCTURE_INTERACTION_GEOMETRY<TV>::
Build_Collision_Geometry(STRUCTURE<TV>& structure)
{
    Clean_Memory();
    if((segmented_curve=Segmented_Curve(&structure)))
        segmented_curve->mesh.elements.Flattened().Get_Unique(collision_particles.active_indices);
    else if((triangulated_surface=Triangulated_Surface(&structure))){
        triangulated_surface->mesh.elements.Flattened().Get_Unique(collision_particles.active_indices);
        triangulated_surface->Update_Number_Nodes();
        if(!triangulated_surface->mesh.segment_mesh) triangulated_surface->mesh.Initialize_Segment_Mesh();
        segmented_curve=new T_SEGMENTED_CURVE(*triangulated_surface->mesh.segment_mesh,full_particles); // TODO: This is broken; long term shallow copy of a temporary auxiliary structure
        need_destroy_segmented_curve=true;}
    else if(FREE_PARTICLES<TV>* free_particles=dynamic_cast<FREE_PARTICLES<TV>*>(&structure))
        collision_particles.active_indices=free_particles->nodes;
    else PHYSBAM_FATAL_ERROR();

    if(segmented_curve && !segmented_curve->mesh.incident_elements) segmented_curve->mesh.Initialize_Incident_Elements();
    collision_particles.Update_Subset_Index_From_Element_Index();
    Build_Topological_Structure_Of_Hierarchies();
}
//#####################################################################
// Function Update_Faces_And_Hierarchies_With_Collision_Free_Positions
//#####################################################################
template<class T,class TV> void Update_Faces_And_Hierarchies_With_Collision_Free_Positions_Helper(TRIANGULATED_SURFACE<T>& triangulated_surface,ARRAY_VIEW<const T> node_thickness,
    const T node_thickness_multiplier,ARRAY_VIEW<const TV> X_old_full)
{}
template<class T> void Update_Faces_And_Hierarchies_With_Collision_Free_Positions_Helper(TRIANGULATED_SURFACE<T>& triangulated_surface,ARRAY_VIEW<const T> node_thickness,
    const T node_thickness_multiplier,ARRAY_VIEW<const VECTOR<T,3> > X_old_full)
{
    triangulated_surface.Update_Triangle_List(X_old_full);triangulated_surface.hierarchy->Update_Leaf_Boxes(X_old_full);
    for(int t=1;t<=triangulated_surface.mesh.elements.m;t++) // increase box size for node thickness
        triangulated_surface.hierarchy->box_hierarchy(t).Change_Size(node_thickness_multiplier*ARRAYS_COMPUTATIONS::Min(node_thickness.Subset(triangulated_surface.mesh.elements(t))));
    triangulated_surface.hierarchy->Update_Nonleaf_Boxes();
}
template<class TV> void STRUCTURE_INTERACTION_GEOMETRY<TV>::
Update_Faces_And_Hierarchies_With_Collision_Free_Positions(ARRAY_VIEW<const T> node_thickness,const T node_thickness_multiplier,ARRAY_VIEW<const TV> X_old_full)
{
    PHYSBAM_ASSERT(node_thickness.Size()==full_particles.array_collection->Size());
    PHYSBAM_ASSERT(node_thickness_multiplier>=1);
    if(segmented_curve){
        segmented_curve->hierarchy->Update_Leaf_Boxes(X_old_full);
        for(int s=1;s<=segmented_curve->mesh.elements.m;s++) // increase box size for node thickness
            segmented_curve->hierarchy->box_hierarchy(s).Change_Size(node_thickness_multiplier*ARRAYS_COMPUTATIONS::Min(node_thickness.Subset(segmented_curve->mesh.elements(s))));
        segmented_curve->hierarchy->Update_Nonleaf_Boxes();}
    if(d==3 && triangulated_surface)
        Update_Faces_And_Hierarchies_With_Collision_Free_Positions_Helper(*triangulated_surface,node_thickness,node_thickness_multiplier,X_old_full);
    particle_hierarchy.Update_Leaf_Boxes(X_old_full.Subset(collision_particles.active_indices));
    for(int k=1;k<=particle_hierarchy.leaves;k++) particle_hierarchy.box_hierarchy(k).Change_Size(node_thickness_multiplier*node_thickness(collision_particles.active_indices(k)));
    particle_hierarchy.Update_Nonleaf_Boxes();
}
//#####################################################################
// Function Update_Processor_Masks_Helper
//#####################################################################
template<class TV> template<class T_OBJECT,class T_HIERARCHY> void STRUCTURE_INTERACTION_GEOMETRY<TV>::
Update_Processor_Masks_Helper(T_OBJECT& object,T_HIERARCHY& hierarchy,const PARTITION_ID processor,const ARRAY<PARTITION_ID>& partition_id_from_particle_index,ARRAY<char>& masks)
{
    typedef typename T_OBJECT::MESH T_MESH;
    // build the leaf masks
    masks.Resize(hierarchy.box_hierarchy.m);
    for(int e=1;e<=object.mesh.elements.m;e++){
        VECTOR<PARTITION_ID,T_MESH::dimension+1> processors(partition_id_from_particle_index.Subset(object.mesh.elements(e)));
        PARTITION_ID simplex_owner=processors.Min();
        masks(e)=0;if(simplex_owner==processor) masks(e)|=2;else if(simplex_owner>processor) masks(e)|=1;}
    // update up the tree
    for(int k=hierarchy.leaves+1;k<=masks.m;k++) masks(k)=masks(hierarchy.children(k-hierarchy.leaves)(1)) | masks(hierarchy.children(k-hierarchy.leaves)(2));
}
//#####################################################################
// Function Update_Processor_Masks
//#####################################################################
template<class TV> void STRUCTURE_INTERACTION_GEOMETRY<TV>::
Update_Processor_Masks(const PARTITION_ID processor,const ARRAY<PARTITION_ID>& partition_id_from_particle_index)
{
    if(segmented_curve)
        Update_Processor_Masks_Helper(*segmented_curve,*segmented_curve->hierarchy,processor,partition_id_from_particle_index,segmented_curve_processor_masks);
    if(triangulated_surface)
        Update_Processor_Masks_Helper(*triangulated_surface,*triangulated_surface->hierarchy,processor,partition_id_from_particle_index,triangulated_surface_processor_masks);
    point_processor_masks.Resize(particle_hierarchy.box_hierarchy.m);
    for(int e=1;e<=collision_particles.active_indices.m;e++){
        int p=collision_particles.active_indices(e);PARTITION_ID particle_processor=partition_id_from_particle_index(p);
        point_processor_masks(e)=0;if(particle_processor==processor) point_processor_masks(e)|=2;else if(particle_processor>processor) point_processor_masks(e)|=1;}
    for(int k=particle_hierarchy.leaves+1;k<=point_processor_masks.m;k++) point_processor_masks(k)=point_processor_masks(particle_hierarchy.children(k-particle_hierarchy.leaves)(1))
        | point_processor_masks(particle_hierarchy.children(k-particle_hierarchy.leaves)(2));
}
//####################################################################
template class STRUCTURE_INTERACTION_GEOMETRY<VECTOR<float,1> >;
template class STRUCTURE_INTERACTION_GEOMETRY<VECTOR<float,2> >;
template class STRUCTURE_INTERACTION_GEOMETRY<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class STRUCTURE_INTERACTION_GEOMETRY<VECTOR<double,1> >;
template class STRUCTURE_INTERACTION_GEOMETRY<VECTOR<double,2> >;
template class STRUCTURE_INTERACTION_GEOMETRY<VECTOR<double,3> >;
#endif
