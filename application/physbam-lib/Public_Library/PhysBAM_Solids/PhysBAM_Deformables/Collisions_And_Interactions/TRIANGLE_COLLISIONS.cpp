//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Neil Molino, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_COLLISIONS  
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
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/STRUCTURE_INTERACTION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS_POINT_FACE_VISITOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> TRIANGLE_COLLISIONS<TV>::
TRIANGLE_COLLISIONS(TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry,const ARRAY<T>& repulsion_thickness)
    :geometry(geometry),repulsion_thickness(repulsion_thickness),final_repulsion_youngs_modulus((T)30),final_repulsion_limiter_fraction((T).1),mpi_solids(0),use_gauss_jacobi(false)
{
    // set parameters 
    Set_Collision_Thickness();Set_Restitution_Coefficient();Set_Gauss_Jacobi();
    // set checking
    Compute_Point_Face_Collisions();Compute_Edge_Edge_Collisions();
    Set_Attempts_For_Rigid_Collisions(false);
    // output
    Output_Collision_Results(false);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> TRIANGLE_COLLISIONS<TV>::
~TRIANGLE_COLLISIONS()
{}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void TRIANGLE_COLLISIONS<TV>::
Initialize(TRIANGLE_COLLISION_PARAMETERS<TV>& triangle_collision_parameters)
{
    Set_Collision_Thickness(triangle_collision_parameters.collisions_collision_thickness);
    if(triangle_collision_parameters.collisions_output_collision_results) Output_Collision_Results();
    Set_Attempts_For_Nonrigid_Collisions(triangle_collision_parameters.collisions_nonrigid_collision_attempts);
    final_repulsion_youngs_modulus=triangle_collision_parameters.collisions_final_repulsion_youngs_modulus;
    final_repulsion_limiter_fraction=triangle_collision_parameters.collisions_final_repulsion_limiter_fraction;
    use_gauss_jacobi=triangle_collision_parameters.use_gauss_jacobi;
}
//#####################################################################
// Function Update_Swept_Hierachies
//#####################################################################
template<class TV> void TRIANGLE_COLLISIONS<TV>::
Update_Swept_Hierachies_And_Compute_Pairs(ARRAY_VIEW<TV> X,ARRAY_VIEW<TV> X_self_collision_free,ARRAY_VIEW<bool> recently_modified,const T detection_thickness)
{
        // Update swept hierarchies
    LOG::SCOPE scope("updating swept hierarchies","updating swept hierarchies");
    for(int structure_index=1;structure_index<=geometry.structure_geometries.m;structure_index++){
        STRUCTURE_INTERACTION_GEOMETRY<TV>& structure=*geometry.structure_geometries(structure_index);
        if(d==3 && structure.triangulated_surface){
            BOX_HIERARCHY<TV>& hierarchy=structure.Face_Hierarchy();
            structure.triangulated_surface_modified.Resize(hierarchy.box_hierarchy.m,false,false);
            for(int kk=1;kk<=structure.triangulated_surface->mesh.elements.m;kk++){
                const VECTOR<int,3>& nodes=structure.triangulated_surface->mesh.elements(kk);
                structure.triangulated_surface_modified(kk)=VECTOR<bool,3>(recently_modified.Subset(nodes)).Contains(true); // TODO: hacking around compiler bug
                if(structure.triangulated_surface_modified(kk))
                    hierarchy.box_hierarchy(kk)=RANGE<TV>::Combine(RANGE<TV>::Bounding_Box(X.Subset(nodes)),RANGE<TV>::Bounding_Box(X_self_collision_free.Subset(nodes)));}
            hierarchy.Update_Modified_Nonleaf_Boxes(structure.triangulated_surface_modified);}

        if(structure.segmented_curve){
            SEGMENT_HIERARCHY<TV>& hierarchy=*structure.segmented_curve->hierarchy;
            structure.segmented_curve_modified.Resize(hierarchy.box_hierarchy.m,false,false);
            for(int kk=1;kk<=structure.segmented_curve->mesh.elements.m;kk++){
                const VECTOR<int,2>& nodes=structure.segmented_curve->mesh.elements(kk);
                structure.segmented_curve_modified(kk)=VECTOR<bool,2>(recently_modified.Subset(nodes)).Contains(true); // TODO: hacking around compiler bug
                if(structure.segmented_curve_modified(kk))
                    hierarchy.box_hierarchy(kk)=RANGE<TV>::Combine(RANGE<TV>::Bounding_Box(X.Subset(nodes)),RANGE<TV>::Bounding_Box(X_self_collision_free.Subset(nodes)));}
            hierarchy.Update_Modified_Nonleaf_Boxes(structure.segmented_curve_modified);}

        {PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > >& hierarchy=structure.particle_hierarchy;
            structure.point_modified.Resize(hierarchy.box_hierarchy.m,false,false);
            for(int kk=1;kk<=hierarchy.leaves;kk++){
                int p=structure.collision_particles.active_indices(kk);
                structure.point_modified(kk)=recently_modified(p);
                if(structure.point_modified(kk)) hierarchy.box_hierarchy(kk)=RANGE<TV>::Bounding_Box(X(p),X_self_collision_free(p));}
            hierarchy.Update_Modified_Nonleaf_Boxes(structure.point_modified);}}
    ARRAYS_COMPUTATIONS::Fill(recently_modified,false);

    // Compute pairs
    point_face_pairs_internal.Remove_All();edge_edge_pairs_internal.Remove_All();
    point_face_pairs_external.Remove_All();edge_edge_pairs_external.Remove_All();
    for(int pair_i=1;pair_i<=geometry.interacting_structure_pairs.m;pair_i++){VECTOR<int,2>& pair=geometry.interacting_structure_pairs(pair_i);
        if(compute_point_face_collisions){
            for(int i=1;i<=2;i++){if(i==2 && pair[1]==pair[2]) break;
                STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1=*geometry.structure_geometries(pair[i]);
                STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2=*geometry.structure_geometries(pair[3-i]);
                Get_Moving_Faces_Near_Moving_Points(structure_1,structure_2,point_face_pairs_internal,point_face_pairs_external,detection_thickness);}}
        if(compute_edge_edge_collisions){
            STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1=*geometry.structure_geometries(pair[1]);
            STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2=*geometry.structure_geometries(pair[2]);
            Get_Moving_Edges_Near_Moving_Edges(structure_1,structure_2,edge_edge_pairs_internal,edge_edge_pairs_external,detection_thickness);}}
}
//#####################################################################
// Function Adjust_Velocity_For_Self_Collisions
//#####################################################################
template<class TV> int TRIANGLE_COLLISIONS<TV>::
Adjust_Velocity_For_Self_Collisions(const T dt,const T time,const bool exit_early)
{
    LOG::SCOPE scope("collisions","checking collisions");
    PARTICLES<TV>& full_particles=geometry.deformable_body_collection.particles;
    ARRAY_VIEW<TV> X(full_particles.X),X_self_collision_free(geometry.X_self_collision_free);ARRAY<bool>& modified_full=geometry.modified_full;
    int collisions=0,collisions_in_attempt=0,
        point_face_collisions=0,edge_edge_collisions=0;
    ARRAY<ARRAY<int> > rigid_lists;ARRAY<int> list_index(full_particles.array_collection->Size()); // index of the rigid list a local node belongs to

    recently_modified_full.Resize(full_particles.array_collection->Size(),false,false);ARRAYS_COMPUTATIONS::Fill(recently_modified_full,true);
    ARRAY<TV> V_save;
    ARRAY<TV> X_save;
    // input velocities are average V.  Also want original velocities?  Delta may be sufficient.

    int attempts=0;bool rigid=false;
    while(!attempts || (!exit_early && collisions_in_attempt)){
        attempts++;if(attempts > nonrigid_collision_attempts) rigid=true;
        if(limit_rigid_collision_attempts && attempts > nonrigid_collision_attempts+rigid_collision_attempts) break; // quit and allow intersections
        LOG::SCOPE scope("collision attempt",STRING_UTILITIES::string_sprintf("collision attempt %d",attempts));

        Update_Swept_Hierachies_And_Compute_Pairs(X,X_self_collision_free,recently_modified_full,collision_thickness);

        T attempt_ratio=(T)attempts/(T)nonrigid_collision_attempts;
        if(mpi_solids){
            //int culled=Prune_Non_Intersecting_Pairs(dt,point_face_pairs,edge_edge_pairs,attempt_ratio);
            //culled=mpi_solids->Reduce_Add_Global(culled);LOG::Stat("Pre-gather collision pairs culled",culled);
            LOG::Time("gathering interaction pairs");
            mpi_solids->Gather_Interaction_Pairs(point_face_pairs_external,edge_edge_pairs_external);} // move all collision pairs to root

        if(geometry.mass_modifier){
            PHYSBAM_FATAL_ERROR();
            geometry.mass_modifier->Reorder_Pairs(edge_edge_pairs_internal,point_face_pairs_internal);}

        int exited_early=1;

        // Make a copy of the particles
        impulse_velocities.Resize(full_particles.array_collection->Size());for(int i=1;i<=full_particles.array_collection->Size();i++) impulse_velocities(i)=full_particles.V(i);
        pf_target_impulses.Resize(point_face_pairs_internal.Size());ARRAYS_COMPUTATIONS::Fill(pf_target_impulses,TV());
        ee_target_impulses.Resize(edge_edge_pairs_internal.Size());ARRAYS_COMPUTATIONS::Fill(ee_target_impulses,TV());
        pf_target_weights.Resize(point_face_pairs_internal.Size());ARRAYS_COMPUTATIONS::Fill(pf_target_weights,VECTOR<T,d+1>());
        ee_target_weights.Resize(edge_edge_pairs_internal.Size());ARRAYS_COMPUTATIONS::Fill(ee_target_weights,VECTOR<T,2*d-2>());
        pf_normals.Resize(point_face_pairs_internal.Size());ARRAYS_COMPUTATIONS::Fill(pf_normals,TV());
        ee_normals.Resize(edge_edge_pairs_internal.Size());ARRAYS_COMPUTATIONS::Fill(ee_normals,TV());
        pf_old_speeds.Resize(point_face_pairs_internal.Size());ARRAYS_COMPUTATIONS::Fill(pf_old_speeds,T());
        ee_old_speeds.Resize(edge_edge_pairs_internal.Size());ARRAYS_COMPUTATIONS::Fill(ee_old_speeds,T());
        
        // point face first for stability
        point_face_collisions=0;edge_edge_collisions=0;collisions_in_attempt=0;
        if(mpi_solids && mpi_solids->rank==0){
            point_face_collisions+=Adjust_Velocity_For_Point_Face_Collision(dt,rigid,rigid_lists,list_index,point_face_pairs_external,attempt_ratio,false,exit_early);
            PHYSBAM_ASSERT(!exit_early);}
        if(mpi_solids){
            LOG::Time("broadcast");
            mpi_solids->Broadcast_Collision_Modified_Data(modified_full,recently_modified_full,full_particles.X,full_particles.V);}
        point_face_collisions+=Adjust_Velocity_For_Point_Face_Collision(dt,rigid,rigid_lists,list_index,point_face_pairs_internal,attempt_ratio,false,exit_early);
        if(exit_early && point_face_collisions) goto EXIT_EARLY_AND_COMMUNICATE;
        if(mpi_solids){
            LOG::Time("gather");
            mpi_solids->Gather_Collision_Modified_Data(modified_full,recently_modified_full,full_particles.X,full_particles.V);}
        // edge edge pairs
        if(mpi_solids && mpi_solids->rank==0){
            edge_edge_collisions+=Adjust_Velocity_For_Edge_Edge_Collision(dt,rigid,rigid_lists,list_index,edge_edge_pairs_external,attempt_ratio,false,exit_early);
            PHYSBAM_ASSERT(!exit_early);}
        if(mpi_solids){
            LOG::Time("broadcast");
            mpi_solids->Broadcast_Collision_Modified_Data(modified_full,recently_modified_full,full_particles.X,full_particles.V);}
        edge_edge_collisions+=Adjust_Velocity_For_Edge_Edge_Collision(dt,rigid,rigid_lists,list_index,edge_edge_pairs_internal,attempt_ratio,false,exit_early);
        if(exit_early && edge_edge_collisions) goto EXIT_EARLY_AND_COMMUNICATE;
        collisions_in_attempt=edge_edge_collisions+point_face_collisions;
        if(mpi_solids) mpi_solids->Gather_Collision_Modified_Data(modified_full,recently_modified_full,full_particles.X,full_particles.V);
        if(mpi_solids){
            LOG::Time("gather");
            mpi_solids->Gather_Collision_Modified_Data(modified_full,recently_modified_full,full_particles.X,full_particles.V);
            LOG::Time("reduce");
            collisions_in_attempt=mpi_solids->Reduce_Add_Global(collisions_in_attempt);}
        collisions+=collisions_in_attempt;

        if(use_gauss_jacobi){
            Scale_And_Apply_Impulses();

            ARRAYS_COMPUTATIONS::Fill(pf_target_impulses,TV());ARRAYS_COMPUTATIONS::Fill(ee_target_impulses,TV());
            ARRAYS_COMPUTATIONS::Fill(pf_target_weights,VECTOR<T,d+1>());ARRAYS_COMPUTATIONS::Fill(ee_target_weights,VECTOR<T,2*d-2>());
            ARRAYS_COMPUTATIONS::Fill(pf_normals,TV());ARRAYS_COMPUTATIONS::Fill(ee_normals,TV());
            ARRAYS_COMPUTATIONS::Fill(pf_old_speeds,T());ARRAYS_COMPUTATIONS::Fill(ee_old_speeds,T());
            Adjust_Velocity_For_Point_Face_Collision(dt,rigid,rigid_lists,list_index,point_face_pairs_internal,attempt_ratio,true,exit_early);
            Adjust_Velocity_For_Edge_Edge_Collision(dt,rigid,rigid_lists,list_index,edge_edge_pairs_internal,attempt_ratio,true,exit_early);
            
            Scale_And_Apply_Impulses();}

        // Apply rigid motions
        if(rigid && collisions_in_attempt && (!mpi_solids || mpi_solids->rank==0)) Apply_Rigid_Body_Motions(dt,rigid_lists);
        // Update positions
        if(collisions_in_attempt){
            for(int p=1;p<=full_particles.array_collection->Size();p++) if(modified_full(p)) full_particles.X(p)=X_self_collision_free(p)+dt*full_particles.V(p);}

        exited_early=0;
        EXIT_EARLY_AND_COMMUNICATE:;
        if(mpi_solids){
            LOG::Time("broadcasting modified data");
            // communicate all data that has modified set to true
            mpi_solids->Broadcast_Collision_Modified_Data(modified_full,recently_modified_full,full_particles.X,full_particles.V);
            LOG::Stop_Time();}

        if(exited_early) break;

        LOG::Stat("processed collisions",collisions_in_attempt);
    }
    if(exit_early && collisions_in_attempt) collisions*=-1; // flag indicating that the collisions were not resolved

    return collisions;
}
template<> int TRIANGLE_COLLISIONS<VECTOR<float,1> >::Adjust_Velocity_For_Self_Collisions(const T,const T time,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> int TRIANGLE_COLLISIONS<VECTOR<double,1> >::Adjust_Velocity_For_Self_Collisions(const T,const T time,const bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Scale_And_Apply_Impulses
//#####################################################################
template<class TV> void TRIANGLE_COLLISIONS<TV>::
Scale_And_Apply_Impulses()
{
    // Go through the computed impulses, and compare them to the actual change in velocity seen. Scale back impulses accordingly
    for(int i=1;i<=pf_target_impulses.m;i++){
        if(pf_target_impulses(i) == TV()) continue;
        const VECTOR<int,d+1>& nodes=point_face_pairs_internal(i);
        // Compute actual new relative_speed
        T relative_speed=TV::Dot_Product(impulse_velocities(nodes(1))-
            (pf_target_weights(i)(2)*impulse_velocities(nodes(2))+pf_target_weights(i)(3)*impulse_velocities(nodes(3))+pf_target_weights(i)(4)*impulse_velocities(nodes(4))),
            pf_normals(i));
        if(relative_speed*pf_old_speeds(i)<0){
            T new_scale=-(relative_speed-pf_old_speeds(i))/pf_old_speeds(i);
            pf_target_impulses(i)/=new_scale;}}
    for(int i=1;i<=ee_target_impulses.m;i++){
        if(ee_target_impulses(i) == TV()) continue;
        const VECTOR<int,2*d-2>& nodes=edge_edge_pairs_internal(i);
        // Compute actual new relative_speed
        T relative_speed=TV::Dot_Product(-(ee_target_weights(i)(1)*impulse_velocities(nodes(1))+ee_target_weights(i)(2)*impulse_velocities(nodes(2))+
                ee_target_weights(i)(3)*impulse_velocities(nodes(3))+ee_target_weights(i)(4)*impulse_velocities(nodes(4))),ee_normals(i));
        if(relative_speed*ee_old_speeds(i)<0){
            T new_scale=-(relative_speed-ee_old_speeds(i))/ee_old_speeds(i);
            ee_target_impulses(i)/=new_scale;}}
    // Apply the newly scaled impulses
    ARRAY_VIEW<T>& one_over_mass=geometry.deformable_body_collection.particles.one_over_mass;
    ARRAY_VIEW<TV> V(geometry.deformable_body_collection.particles.V);
    for(int i=1;i<=pf_target_impulses.m;i++){
        if(pf_target_impulses(i) == TV()) continue;
        const VECTOR<int,d+1>& nodes=point_face_pairs_internal(i);
        VECTOR<T,d+1> one_over_m(one_over_mass.Subset(nodes));
        for(int j=1;j<=d+1;j++) V(nodes(j))+=pf_target_weights(i)(j)*one_over_m[j]*pf_target_impulses(i);}
    for(int i=1;i<=ee_target_impulses.m;i++){
        if(ee_target_impulses(i) == TV()) continue;
        const VECTOR<int,2*d-2>& nodes=edge_edge_pairs_internal(i);
        VECTOR<T,2*d-2> one_over_m(one_over_mass.Subset(nodes));
        for(int j=1;j<=2*d-2;j++) V(nodes(j))+=ee_target_weights(i)(j)*one_over_m[j]*ee_target_impulses(i);}
}
template<> void TRIANGLE_COLLISIONS<VECTOR<float,1> >::Scale_And_Apply_Impulses(){PHYSBAM_NOT_IMPLEMENTED();}
template<> void TRIANGLE_COLLISIONS<VECTOR<double,1> >::Scale_And_Apply_Impulses(){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Get_Moving_Faces_Near_Moving_Points
//#####################################################################
template<class TV> void TRIANGLE_COLLISIONS<TV>::
Get_Moving_Faces_Near_Moving_Points(STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1,STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2,ARRAY<VECTOR<int,d+1> >& pairs_internal,ARRAY<VECTOR<int,d+1> >& pairs_external,const T detection_thickness)
{
    if(!structure_2.Face_Mesh_Object()) return;

    int old_total=pairs_internal.m;
    TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<TV> visitor(pairs_internal,pairs_external,structure_1,structure_2,geometry,detection_thickness,mpi_solids);
    if(mpi_solids){
        BOX_VISITOR_MPI<TRIANGLE_COLLISIONS_POINT_FACE_VISITOR<TV> > mpi_visitor(visitor,structure_1.point_processor_masks,structure_2.Face_Processor_Masks());
        structure_1.particle_hierarchy.Intersection_List(structure_2.Face_Hierarchy(),mpi_visitor,detection_thickness);}
    else structure_1.particle_hierarchy.Intersection_List(structure_2.Face_Hierarchy(),visitor,detection_thickness);

    if(geometry.output_number_checked && pairs_internal.m-old_total>0) LOG::Stat("checked point face collisions",pairs_internal.m-old_total);
}
template<> void TRIANGLE_COLLISIONS<VECTOR<float,1> >::Get_Moving_Faces_Near_Moving_Points(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY<VECTOR<int,d+1> >&,ARRAY<VECTOR<int,d+1> >&,const T){PHYSBAM_NOT_IMPLEMENTED();}
template<> void TRIANGLE_COLLISIONS<VECTOR<double,1> >::Get_Moving_Faces_Near_Moving_Points(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY<VECTOR<int,d+1> >&,ARRAY<VECTOR<int,d+1> >&,const T){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Get_Moving_Edges_Near_Moving_Edges
//#####################################################################
template<class TV> void TRIANGLE_COLLISIONS<TV>::
Get_Moving_Edges_Near_Moving_Edges(STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1,STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2,ARRAY<VECTOR<int,2*d-2> >& pairs_internal,ARRAY<VECTOR<int,2*d-2> >& pairs_external,const T detection_thickness)
{
//    LOG::SCOPE scope("computing edge edge collision pairs", "computing edge edge collision pairs");
    if(!structure_1.Has_Edges() || !structure_2.Has_Edges()) return;

    int old_total=pairs_internal.m;
    TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<TV> visitor(pairs_internal,pairs_external,structure_1,structure_2,geometry,detection_thickness,mpi_solids);
    if(mpi_solids){
        BOX_VISITOR_MPI<TRIANGLE_COLLISIONS_EDGE_EDGE_VISITOR<TV> > mpi_visitor(visitor,structure_1.Edge_Processor_Masks(),structure_2.Edge_Processor_Masks());
        structure_1.Edge_Hierarchy().Intersection_List(structure_2.Edge_Hierarchy(),mpi_visitor,detection_thickness);}
    else structure_1.Edge_Hierarchy().Intersection_List(structure_2.Edge_Hierarchy(),visitor,detection_thickness);

    if(geometry.output_number_checked && pairs_internal.m-old_total>0) LOG::Stat("checked edge edge collisions",pairs_internal.m-old_total);
}
template<> void TRIANGLE_COLLISIONS<VECTOR<float,1> >::Get_Moving_Edges_Near_Moving_Edges(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY<VECTOR<int,2*d-2> >&,ARRAY<VECTOR<int,2*d-2> >&,const T){PHYSBAM_NOT_IMPLEMENTED();}
template<> void TRIANGLE_COLLISIONS<VECTOR<double,1> >::Get_Moving_Edges_Near_Moving_Edges(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY<VECTOR<int,2*d-2> >&,ARRAY<VECTOR<int,2*d-2> >&,const T){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Adjust_Velocity_For_Point_Face_Collision
//#####################################################################
template<class TV> int TRIANGLE_COLLISIONS<TV>::
Adjust_Velocity_For_Point_Face_Collision(const T dt,const bool rigid,ARRAY<ARRAY<int> >& rigid_lists,ARRAY<int>& list_index,const ARRAY<VECTOR<int,d+1> >& pairs,
    const T attempt_ratio,const bool final_repulsion_only,const bool exit_early)
{
    final_point_face_repulsions=final_point_face_collisions=0;
    ARRAY<bool>& modified_full=geometry.modified_full;
    int collisions=0,skipping_already_rigid=0;T collision_time;
    for(int i=1;i<=pairs.m;i++){const VECTOR<int,d+1>& nodes=pairs(i);
        GAUSS_JACOBI_PF_DATA pf_data(pf_target_impulses(i),pf_target_weights(i),pf_normals(i),pf_old_speeds(i));
        if(rigid){VECTOR<int,d+1> node_rigid_indices(list_index.Subset(nodes));if(node_rigid_indices(1) && node_rigid_indices.Elements_Equal()){skipping_already_rigid++;continue;}}
        bool collided;
        if(final_repulsion_only)
            collided=Point_Face_Final_Repulsion(pf_data,nodes,dt,POINT_FACE_REPULSION_PAIR<TV>::Total_Repulsion_Thickness(repulsion_thickness,nodes),collision_time,attempt_ratio,
                exit_early||rigid);
        else collided=Point_Face_Collision(pf_data,nodes,dt,POINT_FACE_REPULSION_PAIR<TV>::Total_Repulsion_Thickness(repulsion_thickness,nodes),collision_time,attempt_ratio,exit_early||rigid);
        if(collided){
            collisions++;
            INDIRECT_ARRAY<ARRAY<bool>,VECTOR<int,d+1>&> modified_subset=modified_full.Subset(nodes);ARRAYS_COMPUTATIONS::Fill(modified_subset,true);
            INDIRECT_ARRAY<ARRAY<bool>,VECTOR<int,d+1>&> recently_modified_subset=recently_modified_full.Subset(nodes);ARRAYS_COMPUTATIONS::Fill(recently_modified_subset,true);
            if(exit_early){if(output_collision_results) {std::stringstream ss;ss<<"exiting collision checking early - point face collision"<<std::endl;LOG::filecout(ss.str());}return collisions;}
            if(rigid) Add_To_Rigid_Lists(rigid_lists,list_index,nodes);}}
    if(output_collision_results && !final_repulsion_only){
        if(final_point_face_repulsions) LOG::Stat("final point face repulsions",final_point_face_repulsions);
        if(final_point_face_collisions) LOG::Stat("final point face collisions",final_point_face_collisions);
        if(skipping_already_rigid) LOG::Stat("rigid collisions where checking was skipped",skipping_already_rigid);
        if(collisions) {std::stringstream ss;ss<<"repeating position update "<<" - point face collisions = "<<collisions<<std::endl;LOG::filecout(ss.str());}}
    return collisions;
}
//#####################################################################
// Function Prune_Non_Intersecting_Pairs
//#####################################################################
template<class T,class T_ARRAY> T Create_Edges(const T_ARRAY& X,const VECTOR<int,2>& nodes,const ARRAY<T>& repulsion_thickness,POINT_2D<T>& point1,POINT_2D<T>& point2)
{
    point1=POINT_2D<T>(X(nodes[1]));
    point2=POINT_2D<T>(X(nodes[2]));
    return EDGE_EDGE_REPULSION_PAIR<typename T_ARRAY::ELEMENT>::Total_Repulsion_Thickness(repulsion_thickness,nodes);
}
template<class T,class T_ARRAY> T Create_Edges(const T_ARRAY& X,const VECTOR<int,4>& nodes,const ARRAY<T>& repulsion_thickness,SEGMENT_3D<T>& segment1,SEGMENT_3D<T>& segment2)
{
    segment1=SEGMENT_3D<T>(X(nodes[1]),X(nodes[2]));
    segment2=SEGMENT_3D<T>(X(nodes[3]),X(nodes[4]));
    return EDGE_EDGE_REPULSION_PAIR<typename T_ARRAY::ELEMENT>::Total_Repulsion_Thickness(repulsion_thickness,nodes);
}
template<class TV> int TRIANGLE_COLLISIONS<TV>::
Prune_Non_Intersecting_Pairs(const T dt,ARRAY<VECTOR<int,d+1> >& point_face_pairs,ARRAY<VECTOR<int,2*d-2> >& edge_edge_pairs,const T attempt_ratio)
{
    int culled=0;T collision_time;
    for(int i=point_face_pairs.m;i>=1;i--){const VECTOR<int,d+1>& nodes=point_face_pairs(i);
        TV temporary_vector;VECTOR<T,d+1> temporary_weights;T temp_old_speed;
        GAUSS_JACOBI_PF_DATA temp_data(temporary_vector,temporary_weights,temporary_vector,temp_old_speed);
        if(!Point_Face_Collision(temp_data,nodes,dt,POINT_FACE_REPULSION_PAIR<TV>::Total_Repulsion_Thickness(repulsion_thickness,nodes),collision_time,attempt_ratio,true)){
            culled++;point_face_pairs.Remove_Index_Lazy(i);}}
    for(int i=edge_edge_pairs.m;i>=1;i--){const VECTOR<int,2*d-2>& nodes=edge_edge_pairs(i);
        TV temporary_vector;VECTOR<T,2*d-2> temporary_weights;T temp_old_speed;
        GAUSS_JACOBI_EE_DATA temp_data(temporary_vector,temporary_weights,temporary_vector,temp_old_speed);
        if(!Edge_Edge_Collision(temp_data,nodes,dt,collision_time,attempt_ratio,true)){
            culled++;edge_edge_pairs.Remove_Index_Lazy(i);}}
    return culled;
}
template<> int TRIANGLE_COLLISIONS<VECTOR<float,1> >::Prune_Non_Intersecting_Pairs(const T,ARRAY<VECTOR<int,d+1> >&,ARRAY<VECTOR<int,2*d-2> >&,const T attempt_ratio)
{PHYSBAM_NOT_IMPLEMENTED();}
template<> int TRIANGLE_COLLISIONS<VECTOR<double,1> >::Prune_Non_Intersecting_Pairs(const T,ARRAY<VECTOR<int,d+1> >&,ARRAY<VECTOR<int,2*d-2> >&,const T attempt_ratio)
{PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Point_Face_Collision
//#####################################################################
namespace{
template<class T,class TV> inline void Update_Velocity_Helper(const T scalar_impulse,const VECTOR<T,2> weights,const TV& normal,const VECTOR<T,3>& one_over_m,
    INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,3>&> V,typename TRIANGLE_COLLISIONS<VECTOR<T,2> >::GAUSS_JACOBI_PF_DATA* pf_data=0)
{
    TV impulse=Pseudo_Divide(scalar_impulse,one_over_m[1]+sqr(weights.x)*one_over_m[2]+sqr(weights.y)*one_over_m[3])*normal;
    V(1)-=one_over_m[1]*impulse;V(2)+=weights.x*one_over_m[2]*impulse;V(3)+=weights.y*one_over_m[3]*impulse;
    if(pf_data){
        pf_data->target_impulse=impulse;
        pf_data->target_weight(1)=-1;pf_data->target_weight(2)=weights.x;pf_data->target_weight(3)=weights.y;
        pf_data->target_normal=normal;}
}
template<class T,class TV> inline void Update_Velocity_Helper(const T scalar_impulse,const VECTOR<T,3> weights,const TV& normal,const VECTOR<T,4>& one_over_m,
    INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,4>&> V,typename TRIANGLE_COLLISIONS<VECTOR<T,3> >::GAUSS_JACOBI_PF_DATA* pf_data=0)
{
    TV impulse=Pseudo_Divide(scalar_impulse,one_over_m[1]+sqr(weights.x)*one_over_m[2]+sqr(weights.y)*one_over_m[3]+sqr(weights.z)*one_over_m[4])*normal;
    V(1)-=one_over_m[1]*impulse;V(2)+=weights.x*one_over_m[2]*impulse;V(3)+=weights.y*one_over_m[3]*impulse;V(4)+=weights.z*one_over_m[4]*impulse;
    if(pf_data){
        pf_data->target_impulse=impulse;
        pf_data->target_weight(1)=-1;pf_data->target_weight(2)=weights.x;pf_data->target_weight(3)=weights.y;pf_data->target_weight(4)=weights.z;
        pf_data->target_normal=normal;}
}
template<class T,class TV> inline SEGMENT_2D<T> Create_Final_Face(const SEGMENT_2D<T>& face,const INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,2>&> V_face,const T dt)
{
    return SEGMENT_2D<T>(face.x1+dt*V_face(1),face.x2+dt*V_face(2));
}
template<class T,class TV> inline TRIANGLE_3D<T> Create_Final_Face(const TRIANGLE_3D<T>& face,const INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,3>&> V_face,const T dt)
{
    return TRIANGLE_3D<T>(face.x1+dt*V_face(1),face.x2+dt*V_face(2),face.x3+dt*V_face(3));
}
}
template<class TV> bool TRIANGLE_COLLISIONS<TV>::
Point_Face_Collision(GAUSS_JACOBI_PF_DATA& pf_data,const VECTOR<int,d+1>& nodes,const T dt,const T repulsion_thickness,T& collision_time,const T attempt_ratio,const bool exit_early)
{
    bool return_type=false;

    ARRAY_VIEW<TV> X(geometry.X_self_collision_free); // TODO: this should be a parameter as we usually want X_self_collision_free but not in Stop_All_Nodes
    ARRAY_VIEW<T>& one_over_mass=geometry.deformable_body_collection.particles.one_over_mass;
    ARRAY_VIEW<TV> V(geometry.deformable_body_collection.particles.V);
    ARRAY_VIEW<TV> V_save(impulse_velocities);
    VECTOR<int,d> face_nodes=nodes.Remove_Index(1);

    T relative_speed;TV normal,weights;T_FACE face(X.Subset(nodes.Remove_Index(1)));
    if(face.Point_Face_Collision(X(nodes[1]),V(nodes[1]),V.Subset(nodes.Remove_Index(1)),dt,collision_thickness,collision_time,normal,weights,relative_speed,exit_early)){
        if(exit_early) return true;
        VECTOR<T,d+1> one_over_m(one_over_mass.Subset(nodes));
        if(geometry.mass_modifier) geometry.mass_modifier->Point_Face_Mass(attempt_ratio,nodes,weights,one_over_m);
        if(use_gauss_jacobi){
            pf_data.old_speed=relative_speed;
            Update_Velocity_Helper((1+restitution_coefficient)*relative_speed,weights,normal,one_over_m,V_save.Subset(nodes),&pf_data);}
        else Update_Velocity_Helper((1+restitution_coefficient)*relative_speed,weights,normal,one_over_m,V.Subset(nodes));
        return_type=true;}

    if(!use_gauss_jacobi) return_type|=Point_Face_Final_Repulsion(pf_data,nodes,dt,repulsion_thickness,collision_time,attempt_ratio,exit_early);

    return return_type;
}
template<> bool TRIANGLE_COLLISIONS<VECTOR<float,1> >::Point_Face_Collision(GAUSS_JACOBI_PF_DATA&,const VECTOR<int,2>&,const T,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> bool TRIANGLE_COLLISIONS<VECTOR<double,1> >::Point_Face_Collision(GAUSS_JACOBI_PF_DATA&,const VECTOR<int,2>&,const T,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Point_Face_Final_Repulsion
//#####################################################################
template<class TV> bool TRIANGLE_COLLISIONS<TV>::
Point_Face_Final_Repulsion(GAUSS_JACOBI_PF_DATA& pf_data,const VECTOR<int,d+1>& nodes,const T dt,const T repulsion_thickness,T& collision_time,const T attempt_ratio,const bool exit_early)
{
    bool return_type=false;

    ARRAY_VIEW<TV> X(geometry.X_self_collision_free); // TODO: this should be a parameter as we usually want X_self_collision_free but not in Stop_All_Nodes
    ARRAY_VIEW<T>& one_over_mass=geometry.deformable_body_collection.particles.one_over_mass;
    ARRAY_VIEW<TV> V(geometry.deformable_body_collection.particles.V);
    ARRAY_VIEW<TV> V_save(impulse_velocities);
    VECTOR<int,d> face_nodes=nodes.Remove_Index(1);INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,d>&> V_face(V,face_nodes);

    // check to see if the final position is too close
    T relative_speed;TV normal,weights;T_FACE face(X.Subset(nodes.Remove_Index(1)));
    T distance;T_FACE face2=Create_Final_Face(face,V_face,dt);TV point(X(nodes[1])+dt*V(nodes[1]));
    if(face2.Point_Face_Interaction(point,V(nodes[1]),V_face,collision_thickness,distance,normal,weights,relative_speed,true,exit_early)){
        collision_time=dt;if(exit_early) return true;
        VECTOR<T,d+1> one_over_m(one_over_mass.Subset(nodes));
        if(geometry.mass_modifier) geometry.mass_modifier->Point_Face_Mass(attempt_ratio,nodes,weights,one_over_m);
        if(!use_gauss_jacobi && relative_speed<0){
            final_point_face_collisions++;
            Update_Velocity_Helper((1+restitution_coefficient)*relative_speed,weights,normal,one_over_m,V.Subset(nodes));
            face2=Create_Final_Face(face,V_face,dt);point=X(nodes[1])+dt*V(nodes[1]); // update point and face and see if repulsion is still necessary
            if(!face2.Point_Face_Interaction(point,V(nodes[1]),V_face,collision_thickness,distance,normal,weights,relative_speed,false,exit_early)) return true;}
        final_point_face_repulsions++;
        T final_relative_speed=final_repulsion_limiter_fraction*(repulsion_thickness-distance)/dt;
        if(relative_speed >= final_relative_speed) return true;
        T ym_over_mass_times_length=final_repulsion_youngs_modulus*one_over_m.Average()/repulsion_thickness;
        T scalar_impulse=min(final_relative_speed-relative_speed,dt*ym_over_mass_times_length*(repulsion_thickness-distance));
        if(use_gauss_jacobi){
            pf_data.old_speed=relative_speed;
            Update_Velocity_Helper(-scalar_impulse,weights,normal,one_over_m,V_save.Subset(nodes),&pf_data);}
        else Update_Velocity_Helper(-scalar_impulse,weights,normal,one_over_m,V.Subset(nodes));
        return_type=true;}

    return return_type;
}
template<> bool TRIANGLE_COLLISIONS<VECTOR<float,1> >::Point_Face_Final_Repulsion(GAUSS_JACOBI_PF_DATA&,const VECTOR<int,2>&,const T,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> bool TRIANGLE_COLLISIONS<VECTOR<double,1> >::Point_Face_Final_Repulsion(GAUSS_JACOBI_PF_DATA&,const VECTOR<int,2>&,const T,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Adjust_Velocity_For_Edge_Edge_Collision
//#####################################################################
template<class TV> int TRIANGLE_COLLISIONS<TV>::
Adjust_Velocity_For_Edge_Edge_Collision(const T dt,const bool rigid,ARRAY<ARRAY<int> >& rigid_lists,ARRAY<int>& list_index,const ARRAY<VECTOR<int,2*d-2> >& pairs,
    const T attempt_ratio,const bool final_repulsion_only,const bool exit_early)
{
    final_edge_edge_repulsions=final_edge_edge_collisions=0;
    ARRAY<bool>& modified_full=geometry.modified_full;
    int collisions=0,skipping_already_rigid=0;T collision_time;

    for(int i=1;i<=pairs.m;i++){const VECTOR<int,2*d-2>& nodes=pairs(i);
        GAUSS_JACOBI_EE_DATA ee_data(ee_target_impulses(i),ee_target_weights(i),ee_normals(i),ee_old_speeds(i));
        if(rigid){VECTOR<int,2*d-2> node_rigid_indices(list_index.Subset(nodes));if(node_rigid_indices(1) && node_rigid_indices.Elements_Equal()){skipping_already_rigid++;continue;}}
        bool collided;
        if(final_repulsion_only)
            collided=Edge_Edge_Final_Repulsion(ee_data,nodes,dt,collision_time,attempt_ratio,(exit_early||rigid));
        else collided=Edge_Edge_Collision(ee_data,nodes,dt,collision_time,attempt_ratio,(exit_early||rigid));
        if(collided){collisions++;
            INDIRECT_ARRAY<ARRAY<bool>,VECTOR<int,2*d-2>&> modified_subset=modified_full.Subset(nodes);ARRAYS_COMPUTATIONS::Fill(modified_subset,true);
            INDIRECT_ARRAY<ARRAY<bool>,VECTOR<int,2*d-2>&> recently_modified_subset=recently_modified_full.Subset(nodes);ARRAYS_COMPUTATIONS::Fill(recently_modified_subset,true);
            if(exit_early){if(output_collision_results) {std::stringstream ss;ss<<"exiting collision checking early - edge collision"<<std::endl;LOG::filecout(ss.str());}return collisions;}
            if(rigid) Add_To_Rigid_Lists(rigid_lists,list_index,nodes);}}
    if(output_collision_results && !final_repulsion_only){
        if(final_edge_edge_repulsions) LOG::Stat("final edge edge repulsions",final_edge_edge_repulsions);
        if(final_edge_edge_collisions) LOG::Stat("final edge edge collisions",final_edge_edge_collisions);
        if(skipping_already_rigid) LOG::Stat("rigid collisions where checking was skipped",skipping_already_rigid);
        if(collisions) {std::stringstream ss;ss<<"repeating position update "<<" - edge edge collisions = "<<collisions<<std::endl;LOG::filecout(ss.str());}}
    return collisions;
}
template<> int TRIANGLE_COLLISIONS<VECTOR<float,1> >::Adjust_Velocity_For_Edge_Edge_Collision(const T,const bool,ARRAY<ARRAY<int> >&,ARRAY<int>&,
    const ARRAY<VECTOR<int,2*d-2> >&,const T,const bool,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> int TRIANGLE_COLLISIONS<VECTOR<double,1> >::Adjust_Velocity_For_Edge_Edge_Collision(const T,const bool,ARRAY<ARRAY<int> >&,ARRAY<int>&,
    const ARRAY<VECTOR<int,2*d-2> >&,const T,const bool,const bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Edge_Edge_Collision
//#####################################################################
namespace{
template<class T,class TV> inline void Update_Velocity_Helper(const T scalar_impulse,const VECTOR<T,2> weights,const VECTOR<T,2>& one_over_m_edges,const TV& normal,
    INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,2>&> V_edges,typename TRIANGLE_COLLISIONS<VECTOR<T,2> >::GAUSS_JACOBI_EE_DATA* ee_data=0)
{
    TV impulse=Pseudo_Divide(scalar_impulse,one_over_m_edges(1)+one_over_m_edges(2))*normal;
    V_edges(1)-=one_over_m_edges(1)*impulse;V_edges(2)+=one_over_m_edges(2)*impulse;
    if(ee_data){
        ee_data->target_impulse=impulse;
        ee_data->target_weight(1)=(T)-1;ee_data->target_weight(2)=(T)1;
        ee_data->target_normal=normal;}
}
template<class T,class TV> inline void Update_Velocity_Helper(const T scalar_impulse,const VECTOR<T,2> weights,const VECTOR<T,4>& one_over_m_edges,const TV& normal,
    INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,4>&> V_edges,typename TRIANGLE_COLLISIONS<VECTOR<T,3> >::GAUSS_JACOBI_EE_DATA* ee_data=0)
{
    TV impulse=Pseudo_Divide(scalar_impulse,sqr(1-weights.x)*one_over_m_edges(1)+sqr(weights.x)*one_over_m_edges(2)+sqr(1-weights.y)*one_over_m_edges(3)+sqr(weights.y)*one_over_m_edges(4))*normal;
    V_edges(1)-=(1-weights.x)*one_over_m_edges(1)*impulse;V_edges(2)-=weights.x*one_over_m_edges(2)*impulse;V_edges(3)+=(1-weights.y)*one_over_m_edges(3)*impulse;V_edges(4)+=weights.y*one_over_m_edges(4)*impulse;
    if(ee_data){
        ee_data->target_impulse=impulse;
        ee_data->target_weight(1)=weights.x-1;ee_data->target_weight(2)=-weights.x;ee_data->target_weight(3)=1-weights.y;ee_data->target_weight(4)=weights.y;
        ee_data->target_normal=normal;}
}
template<class T,class TV> inline void Create_Final_Edges(const POINT_2D<T>& edge1_input,const POINT_2D<T>& edge2_input,const INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,2>&> V_edges,
    const T dt,POINT_2D<T>& edge1_final,POINT_2D<T>& edge2_final)
{
    edge1_final=POINT_2D<T>(edge1_input+dt*V_edges(1));
    edge2_final=POINT_2D<T>(edge2_input+dt*V_edges(2));
}
template<class T,class TV> inline void Create_Final_Edges(const SEGMENT_3D<T>& edge1_input,const SEGMENT_3D<T>& edge2_input,const INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,4>&> V_edges,
    const T dt,SEGMENT_3D<T>& edge1_final,SEGMENT_3D<T>& edge2_final)
{
    edge1_final=SEGMENT_3D<T>(edge1_input.x1+dt*V_edges(1),edge1_input.x2+dt*V_edges(2));
    edge2_final=SEGMENT_3D<T>(edge2_input.x1+dt*V_edges(3),edge2_input.x2+dt*V_edges(4));
}
}
template<class TV> bool TRIANGLE_COLLISIONS<TV>::
Edge_Edge_Collision(GAUSS_JACOBI_EE_DATA& ee_data,const VECTOR<int,2*d-2>& nodes,const T dt,T& collision_time,const T attempt_ratio,const bool exit_early)
{
    ARRAY_VIEW<TV> X(geometry.X_self_collision_free); // TODO: this should be a parameter as we usually want X_self_collision_free but not in Stop_All_Nodes
    ARRAY_VIEW<T>& one_over_mass=geometry.deformable_body_collection.particles.one_over_mass;
    ARRAY_VIEW<TV> V(geometry.deformable_body_collection.particles.V);
    ARRAY_VIEW<TV> V_save(impulse_velocities);
    INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,2*d-2>&> V_edges(V,nodes);
    INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,2*d-2>&> V_save_edges(V_save,nodes);

    bool return_type=false;
    T_EDGE edge1,edge2;Create_Edges(X,nodes,repulsion_thickness,edge1,edge2);

    T relative_speed=0;TV normal;VECTOR<T,2> weights;
    if(edge1.Edge_Edge_Collision(edge2,V_edges,dt,collision_thickness,collision_time,normal,weights,relative_speed,false,geometry.small_number,exit_early)){if(exit_early) return true;
        VECTOR<T,2*d-2> one_over_m_edges(one_over_mass.Subset(nodes));
        if(geometry.mass_modifier) geometry.mass_modifier->Edge_Edge_Mass(attempt_ratio,nodes,weights,one_over_m_edges);
        if(use_gauss_jacobi){
            ee_data.old_speed=relative_speed;
            Update_Velocity_Helper((1+restitution_coefficient)*relative_speed,weights,one_over_m_edges,normal,V_save_edges,&ee_data);}
        else Update_Velocity_Helper((1+restitution_coefficient)*relative_speed,weights,one_over_m_edges,normal,V_edges);
        return_type=true;}

    if(!use_gauss_jacobi) return_type|=Edge_Edge_Final_Repulsion(ee_data,nodes,dt,collision_time,attempt_ratio,exit_early);
    return return_type;
}
template<> bool TRIANGLE_COLLISIONS<VECTOR<float,1> >::Edge_Edge_Collision(GAUSS_JACOBI_EE_DATA&,const VECTOR<int,2*d-2>&,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> bool TRIANGLE_COLLISIONS<VECTOR<double,1> >::Edge_Edge_Collision(GAUSS_JACOBI_EE_DATA&,const VECTOR<int,2*d-2>&,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Edge_Edge_Final_Repulsion
//#####################################################################
template<class TV> bool TRIANGLE_COLLISIONS<TV>::
Edge_Edge_Final_Repulsion(GAUSS_JACOBI_EE_DATA& ee_data,const VECTOR<int,2*d-2>& nodes,const T dt,T& collision_time,const T attempt_ratio,const bool exit_early)
{
    ARRAY_VIEW<TV> X(geometry.X_self_collision_free); // TODO: this should be a parameter as we usually want X_self_collision_free but not in Stop_All_Nodes
    ARRAY_VIEW<T>& one_over_mass=geometry.deformable_body_collection.particles.one_over_mass;
    ARRAY_VIEW<TV> V(geometry.deformable_body_collection.particles.V);
    ARRAY_VIEW<TV> V_save(impulse_velocities);
    INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,2*d-2>&> V_edges(V,nodes);
    INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,2*d-2>&> V_save_edges(V_save,nodes);

    bool return_type=false;
    T_EDGE edge1,edge2;T total_repulsion_thickness=Create_Edges(X,nodes,repulsion_thickness,edge1,edge2);
    T relative_speed=0;TV normal;VECTOR<T,2> weights;

    // check to see if the final position is too close - see if the edge x3-x4 intersects the cylinder around x1-x2
    T distance;T_EDGE edge1_final,edge2_final;Create_Final_Edges(edge1,edge2,V_edges,dt,edge1_final,edge2_final);
    if(edge1_final.Edge_Edge_Interaction(edge2_final,V_edges,collision_thickness,distance,normal,weights,relative_speed,false,geometry.small_number,exit_early)){
        collision_time=dt;if(exit_early) return true;
        VECTOR<T,2*d-2> one_over_m_edges(one_over_mass.Subset(nodes));
        if(geometry.mass_modifier) geometry.mass_modifier->Edge_Edge_Mass(attempt_ratio,nodes,weights,one_over_m_edges);
        if(!use_gauss_jacobi && relative_speed<0){
            final_edge_edge_collisions++;
            Update_Velocity_Helper((1+restitution_coefficient)*relative_speed,weights,one_over_m_edges,normal,V_edges);
            Create_Final_Edges(edge1,edge2,V_edges,dt,edge1_final,edge2_final); // update segments and see if repulsion is still necessary
            if(!edge1_final.Edge_Edge_Interaction(edge2_final,V_edges,collision_thickness,distance,normal,weights,relative_speed,false,geometry.small_number)) return true;}
        final_edge_edge_repulsions++;
        T final_relative_speed=final_repulsion_limiter_fraction*(total_repulsion_thickness-distance)/dt;
        if(relative_speed >= final_relative_speed) return true;
        T ym_over_mass_times_length=final_repulsion_youngs_modulus*one_over_m_edges.Average()/total_repulsion_thickness;
        T scalar_impulse=min(final_relative_speed-relative_speed,dt*ym_over_mass_times_length*(total_repulsion_thickness-distance));
        if(use_gauss_jacobi){
            ee_data.old_speed=relative_speed;
            Update_Velocity_Helper(-scalar_impulse,weights,one_over_m_edges,normal,V_save_edges,&ee_data);}
        else Update_Velocity_Helper(-scalar_impulse,weights,one_over_m_edges,normal,V_edges);
        return_type=true;}
    return return_type;
}
template<> bool TRIANGLE_COLLISIONS<VECTOR<float,1> >::Edge_Edge_Final_Repulsion(GAUSS_JACOBI_EE_DATA&,const VECTOR<int,2*d-2>&,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> bool TRIANGLE_COLLISIONS<VECTOR<double,1> >::Edge_Edge_Final_Repulsion(GAUSS_JACOBI_EE_DATA&,const VECTOR<int,2*d-2>&,const T,T&,const T,const bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Add_To_Rigid_Lists
//#####################################################################
template<class TV> template<int d2> void TRIANGLE_COLLISIONS<TV>::
Add_To_Rigid_Lists(ARRAY<ARRAY<int> >& rigid_lists,ARRAY<int>& list_index,const VECTOR<int,d2>& nodes)
{
    // TODO: make this into union find...

    // make a new list and add the new nodes
    rigid_lists.Resize(rigid_lists.m+1);rigid_lists(rigid_lists.m)=nodes;

    // figure out which list it should be combined with
    int add_list=rigid_lists.m;for(int i=1;i<=nodes.m;i++){int j=list_index(rigid_lists(rigid_lists.m)(i));if(j) add_list=min(add_list,j);}

    // set up a new list or combine with another
    if(add_list == rigid_lists.m) for(int i=1;i<=nodes.m;i++) list_index(rigid_lists(rigid_lists.m)(i))=rigid_lists.m; // label as a new list
    else{ // combine with a pre-existing list
        for(int i=1;i<=nodes.m;i++){
            int node=rigid_lists(rigid_lists.m)(i),current_list=list_index(node);
            if(!current_list){rigid_lists(add_list).Append(node);list_index(node)=add_list;} // add to the add_list
            else if(current_list != add_list){ // not already in the add_list, but in another list - combine current_list with the add_list
                int new_nodes=rigid_lists(current_list).m;rigid_lists(add_list).Resize(rigid_lists(add_list).m+new_nodes);
                for(int j=1;j<=new_nodes;j++){rigid_lists(add_list)(rigid_lists(add_list).m-new_nodes+j)=rigid_lists(current_list)(j);list_index(rigid_lists(current_list)(j))=add_list;}
                rigid_lists(current_list).Resize(0);}} // kill off the current list
        rigid_lists.Resize(rigid_lists.m-1);} // remove the new list since we combined it with another
}
//#####################################################################
// Function Apply_Rigid_Body_Motions
//#####################################################################
namespace{
template<class T> MATRIX<T,0,0> Inertia_Tensor(const ARRAY<int>& rigid_list,const ARRAY_VIEW<T>& mass,const ARRAY<VECTOR<T,1> >& X_self_collision_free,const VECTOR<T,1>& center_of_mass)
{
    MATRIX<T,0,0> I=MATRIX<T,0,0>();
    for(int kk=1;kk<=rigid_list.m;kk++){int p=rigid_list(kk);
        VECTOR<T,1> radial_vector=X_self_collision_free(p)-center_of_mass;
        I+=mass(p)*MATRIX<T,0>();}
    return I;
}
template<class T> MATRIX<T,1,1> Inertia_Tensor(const ARRAY<int>& rigid_list,const ARRAY_VIEW<T>& mass,const ARRAY<VECTOR<T,2> >& X_self_collision_free,const VECTOR<T,2>& center_of_mass)
{
    MATRIX<T,1,1> I=MATRIX<T,1,1>();
    for(int kk=1;kk<=rigid_list.m;kk++){int p=rigid_list(kk);
        VECTOR<T,2> radial_vector=X_self_collision_free(p)-center_of_mass;
        I+=mass(p)*MATRIX<T,1>(radial_vector.Magnitude_Squared());}
    return I;
}
template<class T> SYMMETRIC_MATRIX<T,3> Inertia_Tensor(const ARRAY<int>& rigid_list,const ARRAY_VIEW<T>& mass,const ARRAY<VECTOR<T,3> >& X_self_collision_free,const VECTOR<T,3>& center_of_mass)
{
    SYMMETRIC_MATRIX<T,3> I=SYMMETRIC_MATRIX<T,3>();
    for(int kk=1;kk<=rigid_list.m;kk++){int p=rigid_list(kk);
        VECTOR<T,3> radial_vector=X_self_collision_free(p)-center_of_mass;
        I+=mass(p)*(radial_vector.Magnitude_Squared()-SYMMETRIC_MATRIX<T,3>::Outer_Product(radial_vector));}
    return I;
}
}
template<class TV> void TRIANGLE_COLLISIONS<TV>::
Apply_Rigid_Body_Motions(const T dt,ARRAY<ARRAY<int> >& rigid_lists)
{
    typedef typename TV::SPIN T_SPIN;

    PARTICLES<TV>& full_particles=geometry.deformable_body_collection.particles;ARRAY<TV>& X_self_collision_free=geometry.X_self_collision_free;
    if(output_collision_results){
        {std::stringstream ss;ss<<"TOTAL RIGID GROUPS = "<<rigid_lists.m<<std::endl;LOG::filecout(ss.str());}
        for(int list=1;list<=rigid_lists.m;list++) {std::stringstream ss;ss<<"LIST "<<list<<" = "<<rigid_lists(list).m<<" POINTS"<<std::endl<<rigid_lists(list);LOG::filecout(ss.str());}}
    
    T one_over_dt=1/dt;
    for(int list=1;list<=rigid_lists.m;list++) if(rigid_lists(list).m){
        TV average_velocity,center_of_mass;T total_mass=0;
        for(int kk=1;kk<=rigid_lists(list).m;kk++){int p=rigid_lists(list)(kk);
            total_mass+=full_particles.mass(p);center_of_mass+=full_particles.mass(p)*X_self_collision_free(p);
            average_velocity+=full_particles.mass(p)*(full_particles.X(p)-X_self_collision_free(p));}
        T one_over_total_mass=1/total_mass;center_of_mass*=one_over_total_mass;average_velocity*=one_over_dt*one_over_total_mass;  
        T_SPIN L=T_SPIN(); // moment of inertia & angular momentum
        for(int kk=1;kk<=rigid_lists(list).m;kk++){int p=rigid_lists(list)(kk);
            TV radial_vector=X_self_collision_free(p)-center_of_mass;
            L+=TV::Cross_Product(radial_vector,full_particles.mass(p)*(full_particles.X(p)-X_self_collision_free(p)));} // TODO: figure out why it could be bad to hoist mass out of Cross_Product
        L*=one_over_dt;T_SPIN omega=Inverse(Inertia_Tensor(rigid_lists(list),full_particles.mass,X_self_collision_free,center_of_mass))*L; // TODO: use pseudoinverse
        ROTATION<TV> R=ROTATION<TV>::From_Rotation_Vector(dt*omega);
        for(int kk=1;kk<=rigid_lists(list).m;kk++){int p=rigid_lists(list)(kk);
            TV new_position=center_of_mass+dt*average_velocity+R.Rotate(X_self_collision_free(p)-center_of_mass);
            full_particles.V(p)=one_over_dt*(new_position-X_self_collision_free(p));
            assert(geometry.modified_full(p));recently_modified_full(p)=true;}}
}
//#####################################################################
// Function Stop_Nodes_Before_Self_Collision
//#####################################################################
// stops each node dead in its tracks when it hits something else
template<class TV> void TRIANGLE_COLLISIONS<TV>::
Stop_Nodes_Before_Self_Collision(const T dt)
{
    PHYSBAM_FATAL_ERROR(); // This should be fixed to use X instead of X_self_collision_free by adding the argument to the Point_collision and Edge_edge_collision
    int attempts=0,point_face_collisions=0,edge_edge_collisions=0;T collision_time;bool full_stop=false;
    ARRAY<bool> already_stopped_full(geometry.deformable_body_collection.particles.array_collection->Size());   
    ARRAY_VIEW<TV> X(geometry.deformable_body_collection.particles.X),V(geometry.deformable_body_collection.particles.V);
    ARRAY<bool>& modified_full=geometry.modified_full;
    
    modified_full.Resize(geometry.deformable_body_collection.particles.array_collection->Size(),false,false);ARRAYS_COMPUTATIONS::Fill(modified_full,false);
    while(!attempts || point_face_collisions || edge_edge_collisions){
        LOG::SCOPE scope("Stop Attempt","Stop Attempt");
        if(++attempts > 3) full_stop=true;point_face_collisions=0;edge_edge_collisions=0;
        T attempt_ratio=0; // we don't want to alter mass on subdivision case;

        // update swept hierarchies
        for(int k=1;k<=geometry.structure_geometries.m;k++){ // set up swept hierarchies
            STRUCTURE_INTERACTION_GEOMETRY<TV>& structure=*geometry.structure_geometries(k);
            if(d==3 && structure.triangulated_surface){
                BOX_HIERARCHY<TV>& hierarchy=structure.Face_Hierarchy();
                structure.triangulated_surface_modified.Resize(hierarchy.box_hierarchy.m,false,false);
                for(int kk=1;kk<=structure.triangulated_surface->mesh.elements.m;kk++){
                    const VECTOR<int,3>& nodes=structure.triangulated_surface->mesh.elements(kk);
                    structure.triangulated_surface_modified(kk)=attempts==1 || modified_full.Subset(nodes).Contains(true);
                    if(structure.triangulated_surface_modified(kk))
                        hierarchy.box_hierarchy(kk)=RANGE<TV>::Combine(RANGE<TV>::Bounding_Box(X.Subset(nodes)),
                            RANGE<TV>::Bounding_Box(X(nodes[1])+dt*V(nodes[1]),X(nodes[2])+dt*V(nodes[2]),X(nodes[3])+dt*V(nodes[3])));}
                hierarchy.Update_Modified_Nonleaf_Boxes(structure.triangulated_surface_modified);}
            if(structure.segmented_curve){
                SEGMENT_HIERARCHY<TV>& hierarchy=*structure.segmented_curve->hierarchy;
                structure.segmented_curve_modified.Resize(hierarchy.box_hierarchy.m,false,false);
                for(int kk=1;kk<=structure.segmented_curve->mesh.elements.m;kk++){
                    const VECTOR<int,2>& nodes=structure.segmented_curve->mesh.elements(kk);
                    structure.segmented_curve_modified(kk)=attempts==1 || modified_full.Subset(nodes).Contains(true);
                    if(structure.segmented_curve_modified(kk))
                        hierarchy.box_hierarchy(kk)=RANGE<TV>::Combine(RANGE<TV>::Bounding_Box(X.Subset(nodes)),RANGE<TV>::Bounding_Box(X(nodes[1])+dt*V(nodes[1]),X(nodes[2])+dt*V(nodes[2])));}
                hierarchy.Update_Modified_Nonleaf_Boxes(structure.segmented_curve_modified);}
            
            {PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > >& hierarchy=structure.particle_hierarchy;
            structure.point_modified.Resize(hierarchy.box_hierarchy.m,false,false);
            for(int kk=1;kk<=hierarchy.leaves;kk++){
                int p=structure.collision_particles.active_indices(kk);
                structure.point_modified(kk)=attempts==1 || modified_full(p);
                if(structure.point_modified(kk)) hierarchy.box_hierarchy(kk)=RANGE<TV>::Bounding_Box(X(p),X(p)+dt*V(p));}
            hierarchy.Update_Modified_Nonleaf_Boxes(structure.point_modified);}}

        // compute pairs
        point_face_pairs_internal.Remove_All();edge_edge_pairs_internal.Remove_All();
        for(int pair_i=1;pair_i<=geometry.interacting_structure_pairs.m;pair_i++){VECTOR<int,2>& pair=geometry.interacting_structure_pairs(pair_i);
            if(compute_point_face_collisions){
                for(int i=1;i<=2;i++){if(i==2 && pair[1]==pair[2]) break;
                    STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1=*geometry.structure_geometries(pair[i]);
                    STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2=*geometry.structure_geometries(pair[3-i]);
                    Get_Moving_Faces_Near_Moving_Points(structure_1,structure_2,point_face_pairs_internal,point_face_pairs_external,collision_thickness);}}
            if(compute_edge_edge_collisions){
                STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1=*geometry.structure_geometries(pair[1]);
                STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2=*geometry.structure_geometries(pair[2]);
                Get_Moving_Edges_Near_Moving_Edges(structure_1,structure_2,edge_edge_pairs_internal,edge_edge_pairs_external,collision_thickness);}}
        LOG::Stat("point triangle collision pairs",point_face_pairs_internal.m);
        LOG::Stat("edge edge collision pairs",edge_edge_pairs_internal.m);

        for(int pair_index=1;pair_index<=point_face_pairs_internal.m;pair_index++){
            const VECTOR<int,d+1>& nodes=point_face_pairs_internal(pair_index);
            if(already_stopped_full.Subset(nodes).Number_True()==d+1) continue;
            // TODO: this should not be X self collision free
            TV temporary_vector;VECTOR<T,d+1> temporary_weights;T temp_old_speed;
            GAUSS_JACOBI_PF_DATA temp_data(temporary_vector,temporary_weights,temporary_vector,temp_old_speed);
            if(Point_Face_Collision(temp_data,nodes,dt,0,collision_time,attempt_ratio,true)){
                point_face_collisions++;T scale=collision_time/dt;
                if(full_stop || scale < (T).01){scale=0;
                    INDIRECT_ARRAY<ARRAY<bool>,VECTOR<int,d+1>&> already_stopped_subset=already_stopped_full.Subset(nodes);
                    ARRAYS_COMPUTATIONS::Fill(already_stopped_subset,true);}
                else scale*=(T).95;
                V.Subset(nodes)*=scale;
                INDIRECT_ARRAY<ARRAY<bool>,VECTOR<int,d+1>&> modified_subset=modified_full.Subset(nodes);
                ARRAYS_COMPUTATIONS::Fill(modified_subset,true);}}

        for(int pair_index=1;pair_index<=edge_edge_pairs_internal.m;pair_index++){
            const VECTOR<int,2*d-2>& nodes=edge_edge_pairs_internal(pair_index);
            if(already_stopped_full.Subset(nodes).Number_True()==2*d-2) continue;
            TV temporary_vector;VECTOR<T,2*d-2> temporary_weights;T temp_old_speed;
            GAUSS_JACOBI_EE_DATA temp_data(temporary_vector,temporary_weights,temporary_vector,temp_old_speed);
            if(Edge_Edge_Collision(temp_data,nodes,dt,collision_time,attempt_ratio,true)){
                edge_edge_collisions++;T scale=collision_time/dt;
                if(full_stop || scale < (T).1){scale=0;
                    INDIRECT_ARRAY<ARRAY<bool>,VECTOR<int,2*d-2>&> already_stopped_subset=already_stopped_full.Subset(nodes);
                    ARRAYS_COMPUTATIONS::Fill(already_stopped_subset,true);}
                else scale*=(T).9;
                V.Subset(nodes)*=scale;
                INDIRECT_ARRAY<ARRAY<bool>,VECTOR<int,2*d-2>&> modified_subset=modified_full.Subset(nodes);
                ARRAYS_COMPUTATIONS::Fill(modified_subset,true);}}
        LOG::Stat("point triangle collisions",point_face_collisions);
        LOG::Stat("edge edge collisions",edge_edge_collisions);}
    if(point_face_collisions || edge_edge_collisions) LOG::Stat("attempts at stopping nodes",attempts);

    geometry.deformable_body_collection.particles.Euler_Step_Position(dt);
}
template<> void TRIANGLE_COLLISIONS<VECTOR<float,1> >::Stop_Nodes_Before_Self_Collision(const T){PHYSBAM_NOT_IMPLEMENTED();}
template<> void TRIANGLE_COLLISIONS<VECTOR<double,1> >::Stop_Nodes_Before_Self_Collision(const T){PHYSBAM_NOT_IMPLEMENTED();}
//####################################################################
template class TRIANGLE_COLLISIONS<VECTOR<float,1> >;
template class TRIANGLE_COLLISIONS<VECTOR<float,2> >;
template class TRIANGLE_COLLISIONS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TRIANGLE_COLLISIONS<VECTOR<double,1> >;
template class TRIANGLE_COLLISIONS<VECTOR<double,2> >;
template class TRIANGLE_COLLISIONS<VECTOR<double,3> >;
#endif
}
