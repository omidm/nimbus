//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Neil Molino, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_REPULSIONS
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/STRUCTURE_INTERACTION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_POINT_FACE_VISITOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> TRIANGLE_REPULSIONS<TV>::
TRIANGLE_REPULSIONS(TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry)
    :geometry(geometry),youngs_modulus(30),spring_limiter_fraction((T).1),perform_attractions(true),attractions_threshold(-2),
    hierarchy_repulsion_thickness_multiplier(1),repulsion_thickness_detection_multiplier((T)1.1),mpi_solids(0),use_gauss_jacobi(false)
{
    // set parameters
    Set_Repulsion_Thickness();Set_Friction_Coefficient();
    // set checking
    Compute_Point_Face_Friction();Compute_Edge_Edge_Friction();
    Compute_Point_Face_Inelastic_Collision_Repulsion();Compute_Edge_Edge_Inelastic_Collision_Repulsion();
    Compute_Point_Face_Repulsion();Compute_Edge_Edge_Repulsion();
    // output
    Output_Repulsion_Results(false);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> TRIANGLE_REPULSIONS<TV>::
~TRIANGLE_REPULSIONS()
{}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS<TV>::
Clean_Memory()
{
    repulsion_thickness.Clean_Memory();
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS<TV>::
Initialize(TRIANGLE_COLLISION_PARAMETERS<TV>& triangle_collision_parameters)
{
    Clean_Memory();
    Set_Repulsion_Thickness(triangle_collision_parameters.collisions_repulsion_thickness);
    perform_attractions=triangle_collision_parameters.perform_repulsion_pair_attractions;
    {std::stringstream ss;ss<<"Solids Parameters set perform_attractions "<<perform_attractions<<std::endl;LOG::filecout(ss.str());}
    attractions_threshold=triangle_collision_parameters.repulsion_pair_attractions_threshold;
    {std::stringstream ss;ss<<"Solids Parameters set attractions_threshold "<<attractions_threshold<<std::endl;LOG::filecout(ss.str());}
    if(triangle_collision_parameters.clamp_repulsion_thickness) Clamp_Repulsion_Thickness_With_Meshes(triangle_collision_parameters.collisions_repulsion_clamp_fraction);
    Output_Repulsion_Results(triangle_collision_parameters.collisions_output_repulsion_results);
    Set_Friction_Coefficient(triangle_collision_parameters.self_collision_friction_coefficient);
    if(triangle_collision_parameters.collisions_disable_repulsions_based_on_proximity_factor){
        Turn_Off_Repulsions_Based_On_Current_Proximity(triangle_collision_parameters.collisions_disable_repulsions_based_on_proximity_factor);}
    youngs_modulus=triangle_collision_parameters.repulsions_youngs_modulus;
    spring_limiter_fraction=triangle_collision_parameters.repulsions_limiter_fraction;
    Set_Gauss_Jacobi(triangle_collision_parameters.use_gauss_jacobi);
}
//#####################################################################
// Function Clamp_Repulsion_Thickness_With_Meshes
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS<TV>::
Clamp_Repulsion_Thickness_With_Meshes(ARRAY_VIEW<const TV> X,const T scale)
{
    for(int k=1;k<=geometry.structure_geometries.m;k++){
        STRUCTURE_INTERACTION_GEOMETRY<TV>& structure=*geometry.structure_geometries(k);
        if(structure.segmented_curve) for(int s=1;s<=structure.segmented_curve->mesh.elements.m;s++){
            int i,j;structure.segmented_curve->mesh.elements(s).Get(i,j);T d=scale*(X(i)-X(j)).Magnitude();
            repulsion_thickness(i)=min(repulsion_thickness(i),d);repulsion_thickness(j)=min(repulsion_thickness(j),d);}}
}
template<> void TRIANGLE_REPULSIONS<VECTOR<float,1> >::Clamp_Repulsion_Thickness_With_Meshes(ARRAY_VIEW<const VECTOR<T,1> >,const T){PHYSBAM_NOT_IMPLEMENTED();}
template<> void TRIANGLE_REPULSIONS<VECTOR<double,1> >::Clamp_Repulsion_Thickness_With_Meshes(ARRAY_VIEW<const VECTOR<T,1> >,const T){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Turn_Off_Repulsions_Based_On_Current_Proximity
//#####################################################################
static inline VECTOR<int,2> Flip_Edges(const VECTOR<int,2>& nodes)
{
    return VECTOR<int,2>(nodes[2],nodes[1]);
}
static inline VECTOR<int,4> Flip_Edges(const VECTOR<int,4>& nodes)
{
    return VECTOR<int,4>(nodes[3],nodes[4],nodes[1],nodes[2]);
}
template<class TV> void TRIANGLE_REPULSIONS<TV>::
Turn_Off_Repulsions_Based_On_Current_Proximity(const T extra_factor_on_distance)
{
    if(geometry.intersecting_point_face_pairs.Size() || geometry.intersecting_edge_edge_pairs.Size()) PHYSBAM_FATAL_ERROR();
    const ARRAY<TV>& X_self_collision_free=geometry.X_self_collision_free;
    bool output_number_checked_save=geometry.output_number_checked;geometry.output_number_checked=false;
    repulsion_thickness_detection_multiplier*=extra_factor_on_distance; // enlarge to catch repulsions at a further distance away
    Update_Faces_And_Hierarchies_With_Collision_Free_Positions(0);
    // find interacting pairs
    point_face_interaction_pairs.Remove_All();edge_edge_interaction_pairs.Remove_All();
    // TODO: consider checking cross-structure interactions as well
    for(int a=1;a<=geometry.interacting_structure_pairs.m;a++){VECTOR<int,2>& pair=geometry.interacting_structure_pairs(a);
        if(pair[1]!=pair[2]) continue;
        STRUCTURE_INTERACTION_GEOMETRY<TV>& structure=*geometry.structure_geometries(pair[1]);
        Get_Faces_Near_Points(structure,structure,X_self_collision_free,false);
        Get_Edges_Near_Edges(structure,structure,X_self_collision_free,false);}
    // construct omit hashtables
    omit_point_face_repulsion_pairs.Remove_All();
    for(int k=1;k<=point_face_interaction_pairs.m;k++) omit_point_face_repulsion_pairs.Insert(point_face_interaction_pairs(k).nodes);
    omit_edge_edge_repulsion_pairs.Remove_All();
    for(int k=1;k<=edge_edge_interaction_pairs.m;k++){const VECTOR<int,2*d-2>& nodes=edge_edge_interaction_pairs(k).nodes;
        omit_edge_edge_repulsion_pairs.Insert(nodes); // add edge edge pairs in both orders since the hierarchy can change after this is called
        omit_edge_edge_repulsion_pairs.Insert(Flip_Edges(nodes));}
    // cleanup
    repulsion_thickness_detection_multiplier/=extra_factor_on_distance;
    geometry.output_number_checked=output_number_checked_save;
    LOG::Stat("omitted point face repulsion pairs",point_face_interaction_pairs.m);
    LOG::Stat("omitted edge edge repulsion pairs",edge_edge_interaction_pairs.m);
    point_face_interaction_pairs.Remove_All();edge_edge_interaction_pairs.Remove_All();
}
template<> void TRIANGLE_REPULSIONS<VECTOR<float,1> >::Turn_Off_Repulsions_Based_On_Current_Proximity(const T){PHYSBAM_NOT_IMPLEMENTED();}
template<> void TRIANGLE_REPULSIONS<VECTOR<double,1> >::Turn_Off_Repulsions_Based_On_Current_Proximity(const T){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Update_Faces_And_Hierarchies_With_Collision_Free_Positions
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS<TV>::
Update_Faces_And_Hierarchies_With_Collision_Free_Positions(const ARRAY_VIEW<TV>* X_other)
{
    LOG::SCOPE scope("Update Triangles and Hierarchies","Update Triangles and Hierarchies");
    const ARRAY_VIEW<TV>& X=X_other?*X_other:geometry.X_self_collision_free;
    T multiplier=repulsion_thickness_detection_multiplier*hierarchy_repulsion_thickness_multiplier;
    for(int k=1;k<=geometry.structure_geometries.m;k++)
        geometry.structure_geometries(k)->Update_Faces_And_Hierarchies_With_Collision_Free_Positions(repulsion_thickness,multiplier,X);
}
//#####################################################################
// Function Compute_Interaction_Pairs
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS<TV>::
Compute_Interaction_Pairs(ARRAY_VIEW<const TV> X_other)
{
    LOG::SCOPE scope("computing repulsion pairs", "computing repulsion pairs");
    point_face_interaction_pairs.Remove_All();edge_edge_interaction_pairs.Remove_All();
    for(int pair_i=1;pair_i<=geometry.interacting_structure_pairs.m;pair_i++){VECTOR<int,2>& pair=geometry.interacting_structure_pairs(pair_i);
        for(int i=1;i<=2;i++){if(i==2 && pair[1]==pair[2]) break;
            if(compute_point_face_friction || compute_point_face_inelastic_collision_repulsion || compute_point_face_repulsion){
                Get_Faces_Near_Points(*geometry.structure_geometries(pair[i]),*geometry.structure_geometries(pair[3-i]),X_other,true);}}
        if(compute_edge_edge_friction || compute_edge_edge_inelastic_collision_repulsion || compute_edge_edge_repulsion){
            Get_Edges_Near_Edges(*geometry.structure_geometries(pair[1]),*geometry.structure_geometries(pair[2]),X_other,true);}}

    if(mpi_solids){
        mpi_solids->Gather_Interaction_Pairs(point_face_interaction_pairs,edge_edge_interaction_pairs);

        mpi_solids->Distribute_Repulsion_Pairs(point_face_interaction_pairs,point_face_send_particles,point_face_receive_particles,point_face_boundary_pairs,point_face_internal_pairs);
        
        mpi_solids->Distribute_Repulsion_Pairs(edge_edge_interaction_pairs,edge_edge_send_particles,edge_edge_receive_particles,edge_edge_boundary_pairs,edge_edge_internal_pairs);}
}
//#####################################################################
// Function Adjust_Velocity_For_Self_Repulsion
//#####################################################################
template<class T,class TV,class T_ARRAY> void
Edge_Edge_Interaction_Data_Helper(ARRAY_VIEW<const VECTOR<T,2> > X,EDGE_EDGE_REPULSION_PAIR<TV>& pair,const INDIRECT_ARRAY<T_ARRAY,VECTOR<int,2>&> V_edges,const T& small_number)
{
    pair.normal=(X(pair.nodes[1])-X(pair.nodes[2])).Normalized();
}
template<class T,class TV,class T_ARRAY> void
Edge_Edge_Interaction_Data_Helper(ARRAY_VIEW<const VECTOR<T,3> > X,EDGE_EDGE_REPULSION_PAIR<TV>& pair,const INDIRECT_ARRAY<T_ARRAY,VECTOR<int,4>&> V_edges,const T& small_number)
{
    SEGMENT_3D<T> segment1(X.Subset(VECTOR<int,2>(pair.nodes[1],pair.nodes[2]))),segment2(X.Subset(VECTOR<int,2>(pair.nodes[3],pair.nodes[4])));
    segment1.Edge_Edge_Interaction_Data(segment2,V_edges,pair.distance,pair.normal,pair.weights,small_number);
}
template<class TV> int TRIANGLE_REPULSIONS<TV>::
Adjust_Velocity_For_Self_Repulsion(const T dt,bool use_saved_pairs)
{
    PHYSBAM_ASSERT(!mpi_solids); // Only per time step repulsions are supported in MPI

    LOG::SCOPE scope("repulsions","checking repulsions");
    ARRAY_VIEW<const TV> X_self_collision_free(geometry.X_self_collision_free);ARRAY<bool>& modified_full=geometry.modified_full;ARRAY_VIEW<const TV> V(geometry.deformable_body_collection.particles.V);
    ARRAY<POINT_FACE_REPULSION_PAIR<TV> > point_face_pairs(point_face_interaction_pairs);
    ARRAY<EDGE_EDGE_REPULSION_PAIR<TV> > edge_edge_pairs(edge_edge_interaction_pairs);

    for(int pair_index=1;pair_index<=point_face_pairs.m;pair_index++){
        POINT_FACE_REPULSION_PAIR<TV>& pair=point_face_pairs(pair_index);
        T_FACE face(X_self_collision_free.Subset(pair.nodes.Remove_Index(1)));
        face.Point_Face_Interaction_Data(X_self_collision_free(pair.nodes[1]),pair.distance,pair.normal,pair.weights,perform_attractions);
        INDIRECT_ARRAY<ARRAY<bool>,VECTOR<int,TV::m+1>&> modified_subset=modified_full.Subset(pair.nodes);ARRAYS_COMPUTATIONS::Fill(modified_subset,true);}

    // TODO: do we need update binding here?
    // TODO: MPI check fragments in super fragment for whether this fragment is in the processor?
    geometry.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities();

    for(int pair_index=1;pair_index<=edge_edge_pairs.m;pair_index++){
        EDGE_EDGE_REPULSION_PAIR<TV>& pair=edge_edge_pairs(pair_index);
        // Note: don't call this, all it does is mess up the normal and has already been called when the pairs are created
        //Edge_Edge_Interaction_Data_Helper(X_self_collision_free,pair,V.Subset(pair.nodes),geometry.small_number);
        INDIRECT_ARRAY<ARRAY<bool>,VECTOR<int,2*TV::m-2>&> modified_subset=modified_full.Subset(pair.nodes);
        ARRAYS_COMPUTATIONS::Fill(modified_subset,true);}

    int repulsions=Apply_Repulsions_To_Velocities(dt,point_face_pairs,edge_edge_pairs,true,use_saved_pairs);
    LOG::Stat("adjusting velocity for repulsions",repulsions);
    if(repulsions) for(int p=1;p<=geometry.deformable_body_collection.particles.array_collection->Size();p++) if(modified_full(p)) geometry.deformable_body_collection.particles.X(p)=X_self_collision_free(p)+dt*V(p);
// TODO: this is broken
//    for(int i=1;i<=geometry.deformable_body_collection.rigid_body_particles.array_collection->Size();i++){
//        geometry.deformable_body_collection.rigid_body_particles.Frame(i)=geometry.rigid_body_particle_state_collision_free(i).frame;
//        geometry.deformable_body_collection.rigid_body_particles.Euler_Step_Position(VECTOR<int,1>(i),dt);}
    return repulsions;
}
template<> int TRIANGLE_REPULSIONS<VECTOR<float,1> >::Adjust_Velocity_For_Self_Repulsion(const T,bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> int TRIANGLE_REPULSIONS<VECTOR<double,1> >::Adjust_Velocity_For_Self_Repulsion(const T,bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Pair_Is_Separating
//#####################################################################
template<class TV> bool Pair_Is_Separating(POINT_FACE_REPULSION_PAIR<TV>& pair,ARRAY_VIEW<const TV> V)
{
    TV relative_velocity=V(pair.nodes(1));for(int i=1;i<=TV::dimension;i++) relative_velocity-=V(pair.nodes(i+1))*pair.weights(i);
    return TV::Dot_Product(relative_velocity,pair.normal)>=0;
}
template<class T> bool Pair_Is_Separating(EDGE_EDGE_REPULSION_PAIR<VECTOR<T,1> >& pair,ARRAY_VIEW<const VECTOR<T,1> > V)
{
    return true;
}
template<class T> bool Pair_Is_Separating(EDGE_EDGE_REPULSION_PAIR<VECTOR<T,2> >& pair,ARRAY_VIEW<const VECTOR<T,2> > V)
{
    VECTOR<T,2> relative_velocity=V(pair.nodes(1))-V(pair.nodes(2));
    return VECTOR<T,2>::Dot_Product(relative_velocity,pair.normal)>=0;
}
template<class T> bool Pair_Is_Separating(EDGE_EDGE_REPULSION_PAIR<VECTOR<T,3> >& pair,ARRAY_VIEW<const VECTOR<T,3> > V)
{
    VECTOR<T,3> relative_velocity=V(pair.nodes(1))*(1-pair.weights(1))+V(pair.nodes(2))*pair.weights(1)-V(pair.nodes(3))*(1-pair.weights(2))-V(pair.nodes(4))*pair.weights(2);
    return VECTOR<T,3>::Dot_Product(relative_velocity,pair.normal)>=0;
}
//#####################################################################
// Function Update_Repulsion_Pairs_Using_History
//#####################################################################
template<class T,class TV> bool Edge_Edge_Interaction_Helper(ARRAY_VIEW<const VECTOR<T,2> > X,EDGE_EDGE_REPULSION_PAIR<TV>& pair,const ARRAY<T>& repulsion_thickness,
    const T repulsion_thickness_detection_multiplier)
{
    pair.normal=X(pair.nodes[1])-X(pair.nodes[2]);
    pair.distance=pair.normal.Magnitude();
    T total_repulsion_thickness=repulsion_thickness_detection_multiplier*pair.Total_Repulsion_Thickness(repulsion_thickness);
    return pair.distance<=total_repulsion_thickness;
}
template<class T,class TV> bool Edge_Edge_Interaction_Helper(ARRAY_VIEW<const VECTOR<T,3> > X,EDGE_EDGE_REPULSION_PAIR<TV>& pair,const ARRAY<T>& repulsion_thickness,
    const T repulsion_thickness_detection_multiplier)
{
    SEGMENT_3D<T> segment(X.Subset(VECTOR<int,2>(pair.nodes[1],pair.nodes[2])));
    T total_repulsion_thickness=repulsion_thickness_detection_multiplier*pair.Total_Repulsion_Thickness(repulsion_thickness);
    return segment.Edge_Edge_Interaction(SEGMENT_3D<T>(X(pair.nodes[3]),X(pair.nodes[4])),total_repulsion_thickness,pair.distance,pair.normal,pair.weights,false);
}
template<class TV> template<class T_ARRAY1,class T_ARRAY2> void TRIANGLE_REPULSIONS<TV>::
Update_Repulsion_Pairs_Using_History(T_ARRAY1& point_face_pairs,T_ARRAY2& edge_edge_pairs,bool prune_separating)
{
    ARRAY_VIEW<const TV> X(geometry.deformable_body_collection.particles.X),V(geometry.deformable_body_collection.particles.V);
    for(int pair_index=point_face_pairs.Size();pair_index>=1;pair_index--){
        POINT_FACE_REPULSION_PAIR<TV>& pair=point_face_pairs(pair_index);VECTOR<int,d> face_nodes=pair.nodes.Remove_Index(1);
        T_FACE face(X.Subset(face_nodes));
        if(!face.Point_Face_Interaction(X(pair.nodes(1)),repulsion_thickness_detection_multiplier*pair.Total_Repulsion_Thickness(repulsion_thickness),false,pair.distance) ||
            (prune_separating && Pair_Is_Separating(pair,V))){
            point_face_pairs.Remove_Index_Lazy(pair_index);}
        else face.Point_Face_Interaction_Data(X(pair.nodes[1]),pair.distance,pair.normal,pair.weights,perform_attractions);}

    // TODO: do we need update binding here? (what about the other fragments?)
    geometry.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities();

    for(int pair_index=edge_edge_pairs.Size();pair_index>=1;pair_index--){
        EDGE_EDGE_REPULSION_PAIR<TV>& pair=edge_edge_pairs(pair_index);
        if(!Edge_Edge_Interaction_Helper(X,pair,repulsion_thickness,repulsion_thickness_detection_multiplier) ||
            (prune_separating && Pair_Is_Separating(pair,V))){
            edge_edge_pairs.Remove_Index_Lazy(pair_index);}
        else{
            Edge_Edge_Interaction_Data_Helper(X,pair,V.Subset(pair.nodes),geometry.small_number);
            if(perform_attractions && TV::Dot_Product(pair.normal,pair.collision_free_normal)<attractions_threshold){pair.distance*=-1;pair.normal*=-1;}}}
}
//#####################################################################
// Function Adjust_Velocity_For_Self_Repulsion_Using_History
//#####################################################################
template<class TV> int TRIANGLE_REPULSIONS<TV>::
Adjust_Velocity_For_Self_Repulsion_Using_History(const T dt,const bool use_repulsions,bool use_saved_pairs)
{
   // TODO: MPI
    LOG::SCOPE scope("history repulsions","checking history repulsions");
    if(use_repulsions) geometry.deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Positions(true); // otherwise assume we're in the velocity update, and positions are correct
    geometry.deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Velocities(true);
    int repulsions=0;
    if(mpi_solids){
        ARRAY_VIEW<TV> X(geometry.deformable_body_collection.particles.X),V(geometry.deformable_body_collection.particles.V),X_self_collision_free(geometry.X_self_collision_free);
        mpi_solids->Gather_Repulsion_Inputs(X_self_collision_free,X,V,point_face_send_particles,point_face_receive_particles);
        mpi_solids->Gather_Repulsion_Inputs(X_self_collision_free,X,V,edge_edge_send_particles,edge_edge_receive_particles);
        Update_Repulsion_Pairs_Using_History(point_face_boundary_pairs,edge_edge_boundary_pairs,false);
        Update_Repulsion_Pairs_Using_History(point_face_internal_pairs,edge_edge_internal_pairs,false);
        repulsions+=Apply_Repulsions_To_Velocities(dt,point_face_boundary_pairs,edge_edge_boundary_pairs,point_face_internal_pairs,edge_edge_internal_pairs,use_repulsions);
    }
    else{
        Update_Repulsion_Pairs_Using_History(point_face_interaction_pairs,edge_edge_interaction_pairs,false);
        repulsions+=Apply_Repulsions_To_Velocities(dt,point_face_interaction_pairs,edge_edge_interaction_pairs,use_repulsions,use_saved_pairs);}

    LOG::Stat("adjusting velocity for repulsions using history",repulsions);
    if(repulsions) geometry.deformable_body_collection.Adjust_Mesh_For_Self_Repulsion();
    return repulsions;
}
template<> int TRIANGLE_REPULSIONS<VECTOR<float,1> >::Adjust_Velocity_For_Self_Repulsion_Using_History(const T,const bool,bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> int TRIANGLE_REPULSIONS<VECTOR<double,1> >::Adjust_Velocity_For_Self_Repulsion_Using_History(const T,const bool,bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Adjust_Velocity_For_Self_Repulsion
//#####################################################################
template<class TV> template<class T_ARRAY1,class T_ARRAY2> int TRIANGLE_REPULSIONS<TV>::
Apply_Repulsions_To_Velocities(const T dt,T_ARRAY1& point_face_interaction_pairs,T_ARRAY2& edge_edge_interaction_pairs,const bool use_repulsions,bool use_saved_pairs)
{
    int repulsions=0;

    if(compute_point_face_friction){
        Adjust_Velocity_For_Point_Face_Repulsion(dt,point_face_interaction_pairs,use_repulsions,true,use_repulsions); // if not using repulsions, only do friction for inelastic component
        repulsions+=point_face_interaction_pairs.Size();
        if(output_repulsion_results) LOG::Stat("total point face friction",point_face_interaction_pairs.Size());}
    if(compute_edge_edge_friction){
        Adjust_Velocity_For_Edge_Edge_Repulsion(dt,edge_edge_interaction_pairs,use_repulsions,true,use_repulsions); // if not using repulsions, only do friction for inelastic component
        repulsions+=edge_edge_interaction_pairs.Size();
        if(output_repulsion_results) LOG::Stat("total edge edge friction",edge_edge_interaction_pairs.Size());}
    if(compute_point_face_inelastic_collision_repulsion){
        Adjust_Velocity_For_Point_Face_Repulsion(dt,point_face_interaction_pairs,false,false,use_repulsions);
        repulsions+=point_face_inelastic_collision_repulsion_attempts*point_face_interaction_pairs.Size();
        if(output_repulsion_results) LOG::Stat("total point face inelastic collision repulsion",point_face_inelastic_collision_repulsion_attempts*point_face_interaction_pairs.Size());}
    if(compute_edge_edge_inelastic_collision_repulsion){
        Adjust_Velocity_For_Edge_Edge_Repulsion(dt,edge_edge_interaction_pairs,false,false,use_repulsions);
        repulsions+=edge_edge_inelastic_collision_repulsion_attempts*edge_edge_interaction_pairs.Size();
        if(output_repulsion_results) LOG::Stat("total edge edge inelastic collision repulsion",edge_edge_inelastic_collision_repulsion_attempts*edge_edge_interaction_pairs.Size());}
    if(use_repulsions){
        if(compute_point_face_repulsion){
            Adjust_Velocity_For_Point_Face_Repulsion(dt,point_face_interaction_pairs,true,false,use_repulsions);
            repulsions+=point_face_interaction_pairs.Size();
            if(output_repulsion_results) LOG::Stat("total point face repulsion",point_face_interaction_pairs.Size());}
        if(compute_edge_edge_repulsion){
            Adjust_Velocity_For_Edge_Edge_Repulsion(dt,edge_edge_interaction_pairs,true,false,use_repulsions);
            repulsions+=edge_edge_interaction_pairs.Size();
            if(output_repulsion_results) LOG::Stat("total edge edge repulsion",edge_edge_interaction_pairs.Size());}}
    return repulsions;
}
template<class TV> template<class T_ARRAY1,class T_ARRAY2> int TRIANGLE_REPULSIONS<TV>::
Apply_Repulsions_To_Velocities(const T dt,T_ARRAY1& point_face_boundary_pairs,T_ARRAY2& edge_edge_boundary_pairs,T_ARRAY1& point_face_internal_pairs,T_ARRAY2& edge_edge_internal_pairs,
    const bool use_repulsions)
{
    ARRAY_VIEW<TV> X(geometry.deformable_body_collection.particles.X),V(geometry.deformable_body_collection.particles.V),X_self_collision_free(geometry.X_self_collision_free);
    int repulsions=0;bool used=false;
    if(compute_point_face_friction){
        used=true;
        Adjust_Velocity_For_Point_Face_Repulsion(dt,point_face_boundary_pairs,use_repulsions,true,use_repulsions); // if not using repulsions, only do friction for inelastic component
        mpi_solids->Scatter_Repulsion_Outputs(X_self_collision_free,X,V,point_face_send_particles,point_face_receive_particles);
        Adjust_Velocity_For_Point_Face_Repulsion(dt,point_face_internal_pairs,use_repulsions,true,use_repulsions); // if not using repulsions, only do friction for inelastic component
        int applied_repulsions=point_face_boundary_pairs.Size()+point_face_internal_pairs.Size();
        repulsions+=applied_repulsions;
        if(output_repulsion_results) LOG::Stat("total point face friction",applied_repulsions);}
    if(compute_edge_edge_friction){
        if(used) mpi_solids->Gather_Repulsion_Inputs(X_self_collision_free,X,V,edge_edge_send_particles,edge_edge_receive_particles);used=true;
        Adjust_Velocity_For_Edge_Edge_Repulsion(dt,edge_edge_boundary_pairs,use_repulsions,true,use_repulsions); // if not using repulsions, only do friction for inelastic component
        mpi_solids->Scatter_Repulsion_Outputs(X_self_collision_free,X,V,edge_edge_send_particles,edge_edge_receive_particles);
        Adjust_Velocity_For_Edge_Edge_Repulsion(dt,edge_edge_internal_pairs,use_repulsions,true,use_repulsions); // if not using repulsions, only do friction for inelastic component
        int applied_repulsions=edge_edge_boundary_pairs.Size()+edge_edge_internal_pairs.Size();
        repulsions+=applied_repulsions;
        if(output_repulsion_results) LOG::Stat("total edge edge friction",applied_repulsions);}
    if(compute_point_face_inelastic_collision_repulsion){
        if(used) mpi_solids->Gather_Repulsion_Inputs(X_self_collision_free,X,V,point_face_send_particles,point_face_receive_particles);used=true;
        Adjust_Velocity_For_Point_Face_Repulsion(dt,point_face_boundary_pairs,false,false,use_repulsions);
        mpi_solids->Scatter_Repulsion_Outputs(X_self_collision_free,X,V,point_face_send_particles,point_face_receive_particles);
        Adjust_Velocity_For_Point_Face_Repulsion(dt,point_face_internal_pairs,false,false,use_repulsions);
        int applied_repulsions=point_face_inelastic_collision_repulsion_attempts*(point_face_boundary_pairs.Size()+point_face_internal_pairs.Size());
        repulsions+=applied_repulsions;
        if(output_repulsion_results) LOG::Stat("total point face inelastic collision repulsion",applied_repulsions);}
    if(compute_edge_edge_inelastic_collision_repulsion){
        if(used) mpi_solids->Gather_Repulsion_Inputs(X_self_collision_free,X,V,edge_edge_send_particles,edge_edge_receive_particles);used=true;
        Adjust_Velocity_For_Edge_Edge_Repulsion(dt,edge_edge_boundary_pairs,false,false,use_repulsions);
        mpi_solids->Scatter_Repulsion_Outputs(X_self_collision_free,X,V,edge_edge_send_particles,edge_edge_receive_particles);
        Adjust_Velocity_For_Edge_Edge_Repulsion(dt,edge_edge_internal_pairs,false,false,use_repulsions);
        int applied_repulsions=edge_edge_inelastic_collision_repulsion_attempts*(edge_edge_boundary_pairs.Size()+edge_edge_internal_pairs.Size());
        repulsions+=applied_repulsions;
        if(output_repulsion_results) LOG::Stat("total edge edge inelastic collision repulsion",applied_repulsions);}
    if(use_repulsions){
        if(compute_point_face_repulsion){
            if(used) mpi_solids->Gather_Repulsion_Inputs(X_self_collision_free,X,V,point_face_send_particles,point_face_receive_particles);used=true;
            Adjust_Velocity_For_Point_Face_Repulsion(dt,point_face_boundary_pairs,true,false,use_repulsions);
            mpi_solids->Scatter_Repulsion_Outputs(X_self_collision_free,X,V,point_face_send_particles,point_face_receive_particles);
            Adjust_Velocity_For_Point_Face_Repulsion(dt,point_face_internal_pairs,true,false,use_repulsions);
            int applied_repulsions=point_face_boundary_pairs.Size()+point_face_internal_pairs.Size();
            repulsions+=applied_repulsions;
            if(output_repulsion_results) LOG::Stat("total point face repulsion",applied_repulsions);}
        if(compute_edge_edge_repulsion){
            if(used) mpi_solids->Gather_Repulsion_Inputs(X_self_collision_free,X,V,edge_edge_send_particles,edge_edge_receive_particles);used=true;
            Adjust_Velocity_For_Edge_Edge_Repulsion(dt,edge_edge_boundary_pairs,true,false,use_repulsions);
            mpi_solids->Scatter_Repulsion_Outputs(X_self_collision_free,X,V,edge_edge_send_particles,edge_edge_receive_particles);
            Adjust_Velocity_For_Edge_Edge_Repulsion(dt,edge_edge_internal_pairs,true,false,use_repulsions);
            int applied_repulsions=edge_edge_boundary_pairs.Size()+edge_edge_internal_pairs.Size();
            repulsions+=applied_repulsions;
            if(output_repulsion_results) LOG::Stat("total edge edge repulsion",applied_repulsions);}}
    return repulsions;
}
//#####################################################################
// Function Get_Faces_Near_Points
//#####################################################################
template<class TV> int TRIANGLE_REPULSIONS<TV>::
Get_Faces_Near_Points(STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1,STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2,ARRAY_VIEW<const TV> X_other,const bool use_processor_cull)
{
    if(!structure_2.Face_Mesh_Object()) return 0;

    int start_count=point_face_interaction_pairs.m;int pruned=0;
    TRIANGLE_REPULSIONS_POINT_FACE_VISITOR<TV> visitor(point_face_interaction_pairs,structure_1,structure_2,X_other,*this,pruned);
    if(use_processor_cull && mpi_solids){
        BOX_VISITOR_MPI<TRIANGLE_REPULSIONS_POINT_FACE_VISITOR<TV> > mpi_visitor(visitor,structure_1.point_processor_masks,structure_2.Face_Processor_Masks());
        structure_1.particle_hierarchy.Intersection_List(structure_2.Face_Hierarchy(),mpi_visitor,ZERO());}
    else structure_1.particle_hierarchy.Intersection_List(structure_2.Face_Hierarchy(),visitor,ZERO());

    int checked=point_face_interaction_pairs.m-start_count;
    if(checked){
//        LOG::cout<<"pruned "<<pruned<<" out of "<<checked+pruned<<std::endl;
        if(geometry.output_number_checked) LOG::Stat("checked point face repulsions",checked);}
    return checked;
}
template<> int TRIANGLE_REPULSIONS<VECTOR<float,1> >::Get_Faces_Near_Points(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY_VIEW<const VECTOR<T,1> >,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> int TRIANGLE_REPULSIONS<VECTOR<double,1> >::Get_Faces_Near_Points(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY_VIEW<const VECTOR<T,1> >,const bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Output_Interaction_Pairs
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS<TV>::
Output_Interaction_Pairs(const STREAM_TYPE stream_type,const std::string& filename) const
{
    FILE_UTILITIES::Write_To_File(stream_type,filename,point_face_interaction_pairs,edge_edge_interaction_pairs);
}
template<> void TRIANGLE_REPULSIONS<VECTOR<float,1> >::Output_Interaction_Pairs(const STREAM_TYPE,const std::string&) const {PHYSBAM_NOT_IMPLEMENTED();}
template<> void TRIANGLE_REPULSIONS<VECTOR<double,1> >::Output_Interaction_Pairs(const STREAM_TYPE,const std::string&) const {PHYSBAM_NOT_IMPLEMENTED();}

//#####################################################################
// Function Get_Edges_Near_Edges
//#####################################################################
template<class TV> int TRIANGLE_REPULSIONS<TV>::
Get_Edges_Near_Edges(STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1,STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2,ARRAY_VIEW<const TV> X_other,const bool use_processor_cull)
{
    if(!structure_1.Has_Edges() || !structure_2.Has_Edges()) return 0;

    int start_count=edge_edge_interaction_pairs.m;int pruned=0;
    TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<TV> visitor(edge_edge_interaction_pairs,structure_1,structure_2,X_other,*this,pruned);
    if(use_processor_cull && mpi_solids){
        BOX_VISITOR_MPI<TRIANGLE_REPULSIONS_EDGE_EDGE_VISITOR<TV> > mpi_visitor(visitor,structure_1.Edge_Processor_Masks(),structure_2.Edge_Processor_Masks());
        structure_1.Edge_Hierarchy().Intersection_List(structure_2.Edge_Hierarchy(),mpi_visitor,ZERO());}
    else structure_1.Edge_Hierarchy().Intersection_List(structure_2.Edge_Hierarchy(),visitor,ZERO());

    int checked=edge_edge_interaction_pairs.m-start_count;
    if(checked){
        //LOG::cout<<"pruned "<<pruned<<" out of "<<checked+pruned<<std::endl;
        if(geometry.output_number_checked) LOG::Stat("checked edge edge repulsions",checked);}
    return checked;
}
template<> int TRIANGLE_REPULSIONS<VECTOR<float,1> >::Get_Edges_Near_Edges(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY_VIEW<const VECTOR<T,1> >,const bool){PHYSBAM_NOT_IMPLEMENTED();}
template<> int TRIANGLE_REPULSIONS<VECTOR<double,1> >::Get_Edges_Near_Edges(STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,STRUCTURE_INTERACTION_GEOMETRY<VECTOR<T,1> >&,
    ARRAY_VIEW<const VECTOR<T,1> >,const bool){PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Repulsion_Impulse
//#####################################################################
template<class TV> template<class T_PAIR> inline typename TV::SCALAR TRIANGLE_REPULSIONS<TV>::
Repulsion_Impulse(TV& direction,const T dt,const T_PAIR& pair,const TV& relative_velocity,const bool elastic_repulsion,const bool friction)
{
    const T spring_limiter_fraction_over_dt=spring_limiter_fraction/dt;
    const T relative_speed=TV::Dot_Product(relative_velocity,pair.normal);

    // compute scalar impulse
    T scalar_impulse;
    if(elastic_repulsion){
        T total_thickness=pair.Total_Repulsion_Thickness(repulsion_thickness);
        T repulsion_thickness_minus_distance=total_thickness-pair.distance;
        T final_relative_speed=spring_limiter_fraction_over_dt*repulsion_thickness_minus_distance;
        if(relative_speed>=final_relative_speed) return 0;
        T ym_over_mass_times_length=youngs_modulus*VECTOR<T,T_PAIR::count>(
            geometry.deformable_body_collection.particles.one_over_effective_mass.Subset(pair.nodes)).Average()/total_thickness;
        T spring_impulse=dt*ym_over_mass_times_length*max((T)0,repulsion_thickness_minus_distance);
        if(friction) scalar_impulse=min(final_relative_speed-relative_speed,spring_impulse-min((T)0,relative_speed));
        else scalar_impulse=min(final_relative_speed-relative_speed,spring_impulse);}
    else if(relative_speed>=0) return 0;
    else{
        T total_thickness=pair.Total_Repulsion_Thickness(repulsion_thickness);
        T final_relative_speed=total_thickness>pair.distance?0:spring_limiter_fraction_over_dt*(total_thickness-pair.distance);
        if(relative_speed>=final_relative_speed) return 0;
        scalar_impulse=final_relative_speed-relative_speed;}

    // compute friction if necessary and set direction
    if(friction){
        TV tangent=relative_velocity.Projected_Orthogonal_To_Unit_Direction(pair.normal);
        T relative_tangent_velocity_magnitude=tangent.Magnitude(),friction_based_velocity_change=friction_coefficient*scalar_impulse;
        scalar_impulse=-1;if(friction_based_velocity_change<relative_tangent_velocity_magnitude) scalar_impulse=-friction_based_velocity_change/relative_tangent_velocity_magnitude;
        direction=tangent;}
    else direction=pair.normal;

    return scalar_impulse;
}
//#####################################################################
// Function Adjust_Velocity_For_Point_Face_Repulsion
//#####################################################################
template<class TV,class T_ARRAY,class T_MASS_ARRAY> inline void Update_Velocity_Helper(const TV& impulse,const TV& weights,INDIRECT_ARRAY<T_MASS_ARRAY,VECTOR<int,TV::m+1>&> one_over_m,INDIRECT_ARRAY<T_ARRAY,VECTOR<int,TV::m+1>&> V)
{
    V(1)+=one_over_m(1)*impulse;
    for(int i=1;i<=TV::m;i++) V(i+1)-=weights(i)*one_over_m(i+1)*impulse;
}
template<class TV> template<class T_ARRAY> void TRIANGLE_REPULSIONS<TV>::
Adjust_Velocity_For_Point_Face_Repulsion(const T dt,const T_ARRAY& pairs,const bool elastic_repulsion,const bool friction,const bool use_repulsions)
{
    PARTICLES<TV>& particles=geometry.deformable_body_collection.particles;
    int attempts=0,total_attempts=1;if(!elastic_repulsion && !friction) total_attempts=point_face_inelastic_collision_repulsion_attempts;
    int inverted_pairs=0,applied_impulses=0;
    ARRAY_VIEW<TV> V(particles.V);
    ARRAY_VIEW<const T> one_over_effective_mass(particles.one_over_effective_mass);
    while(++attempts<=total_attempts){
        impulse_velocities.Resize(particles.array_collection->Size());impulse_velocities=V;
        pf_target_impulses.Resize(pairs.Size());ARRAYS_COMPUTATIONS::Fill(pf_target_impulses,TV());
        pf_normals.Resize(pairs.Size());ARRAYS_COMPUTATIONS::Fill(pf_normals,TV());
        pf_old_speeds.Resize(pairs.Size());ARRAYS_COMPUTATIONS::Fill(pf_old_speeds,T());
        {std::stringstream ss;ss<<"Repulsion application step "<<attempts<<std::endl;LOG::filecout(ss.str());}

        for(int pair_index=1;pair_index<=pairs.Size();pair_index++){
            const POINT_FACE_REPULSION_PAIR<TV>& pair=pairs(pair_index);
            int p=pair.nodes[1];VECTOR<int,d> face_nodes=pair.nodes.Remove_Index(1);
            if(pair.distance<0) inverted_pairs++;

            INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,d>&> V_subset=V.Subset(face_nodes);
            TV relative_velocity=V(p)-ARRAYS_COMPUTATIONS::Weighted_Sum(pair.weights,V_subset);
            TV direction;T scalar_impulse=Repulsion_Impulse(direction,dt,pair,relative_velocity,elastic_repulsion,friction);
            
            if(scalar_impulse){
                applied_impulses++;
                if(0 && geometry.mass_modifier&&use_repulsions) geometry.mass_modifier->Point_Face_Mass((T)1.,pair.nodes,pair.weights,particles.one_over_mass);
                T one_over_mass=one_over_effective_mass(p);
                for(int i=1;i<=face_nodes.m;i++) one_over_mass+=sqr(pair.weights[i])*one_over_effective_mass(face_nodes[i]);
                TV impulse=Pseudo_Divide(scalar_impulse*direction,one_over_mass);
                if(use_gauss_jacobi && !friction && !elastic_repulsion){
                    Update_Velocity_Helper(impulse,pair.weights,one_over_effective_mass.Subset(pair.nodes),impulse_velocities.Subset(pair.nodes));
                    pf_old_speeds(pair_index)=TV::Dot_Product(relative_velocity,direction);
                    pf_target_impulses(pair_index)=impulse;
                    pf_normals(pair_index)=direction; // tangential for friction.
                }
                else
                    Update_Velocity_Helper(impulse,pair.weights,one_over_effective_mass.Subset(pair.nodes),V.Subset(pair.nodes));
                if(0 && geometry.mass_modifier&&use_repulsions) geometry.mass_modifier->Point_Face_Mass_Revert(pair.nodes,particles.one_over_mass);}}

        if(use_gauss_jacobi) Scale_And_Apply_Point_Face_Impulses(pairs);}
    if(inverted_pairs) LOG::Stat("inverted point face repulsion pairs",inverted_pairs);
    if(applied_impulses) LOG::Stat("applied point face repulsion impulses",applied_impulses);
}
//#####################################################################
// Function Adjust_Velocity_For_Edge_Edge_Repulsion
//#####################################################################
template<class T,class T_ARRAY,class T_MASS_ARRAY> inline void Edge_Edge_Update_Velocity_Helper(const VECTOR<T,3>& impulse,const VECTOR<T,2>& weights,INDIRECT_ARRAY<T_MASS_ARRAY,VECTOR<int,4>&> one_over_effective_mass,
    INDIRECT_ARRAY<T_ARRAY,VECTOR<int,4>&> V)
{
    V(1)+=(1-weights[1])*one_over_effective_mass(1)*impulse;
    V(2)+=weights[1]*one_over_effective_mass(2)*impulse;
    V(3)-=(1-weights[2])*one_over_effective_mass(3)*impulse;
    V(4)-=weights[2]*one_over_effective_mass(4)*impulse;
}
template<class T,class T_ARRAY,class T_MASS_ARRAY> inline void Edge_Edge_Update_Velocity_Helper(const VECTOR<T,2>& impulse,const VECTOR<T,2>& weights,INDIRECT_ARRAY<T_MASS_ARRAY,VECTOR<int,2>&> one_over_effective_mass,
    INDIRECT_ARRAY<T_ARRAY,VECTOR<int,2>&> V)
{
    V(1)+=one_over_effective_mass(1)*impulse;
    V(2)-=one_over_effective_mass(2)*impulse;
}
template<class T,class T_ARRAY> VECTOR<T,3> Edge_Edge_Relative_Velocity_Helper(const VECTOR<T,2>& weights,INDIRECT_ARRAY<T_ARRAY,VECTOR<int,4>&> V)
{
    return (1-weights[1])*V(1)+weights[1]*V(2)-(1-weights[2])*V(3)-weights[2]*V(4);
}
template<class T,class T_ARRAY> VECTOR<T,2> Edge_Edge_Relative_Velocity_Helper(const VECTOR<T,2>& weights,INDIRECT_ARRAY<T_ARRAY,VECTOR<int,2>&> V)
{
    return V(1)-V(2);
}
template<class T,class T_MASS_ARRAY> T Edge_Edge_One_Over_Mass_Helper(const VECTOR<T,2>& weights,INDIRECT_ARRAY<T_MASS_ARRAY,VECTOR<int,4>&> one_over_effective_mass)
{
    return sqr(1-weights[1])*one_over_effective_mass(1)+sqr(weights[1])*one_over_effective_mass(2)
        +sqr(1-weights[2])*one_over_effective_mass(3)+sqr(weights[2])*one_over_effective_mass(4);
}
template<class T,class T_MASS_ARRAY> T Edge_Edge_One_Over_Mass_Helper(const VECTOR<T,2>& weights,INDIRECT_ARRAY<T_MASS_ARRAY,VECTOR<int,2>&> one_over_effective_mass)
{
    return one_over_effective_mass(1)+one_over_effective_mass(2);
}
template<class TV> template<class T_ARRAY,class S> void TRIANGLE_REPULSIONS<TV>::
Adjust_Velocity_For_Edge_Edge_Repulsion_Helper(const T dt,const T_ARRAY& pairs,const bool elastic_repulsion,const bool friction,const VECTOR<S,2>&,const bool use_repulsions)
{
    Adjust_Velocity_For_Edge_Edge_Repulsion_Helper(dt,pairs,elastic_repulsion,friction,VECTOR<T,3>(),use_repulsions);
}
template<class TV> template<class T_ARRAY,class S> void TRIANGLE_REPULSIONS<TV>::
Adjust_Velocity_For_Edge_Edge_Repulsion_Helper(const T dt,const T_ARRAY& pairs,const bool elastic_repulsion,const bool friction,const VECTOR<S,3>&,const bool use_repulsions)
{
    PARTICLES<TV>& particles=geometry.deformable_body_collection.particles;
    int attempts=0,total_attempts=1;if(!elastic_repulsion && !friction) total_attempts=edge_edge_inelastic_collision_repulsion_attempts;
    int inverted_pairs=0,applied_impulses=0;
    ARRAY_VIEW<TV> V(particles.V);
    ARRAY_VIEW<const T> one_over_effective_mass(particles.one_over_effective_mass);
    while(++attempts<=total_attempts){
        impulse_velocities.Resize(particles.array_collection->Size());impulse_velocities=V;
        ee_target_impulses.Resize(pairs.Size());ARRAYS_COMPUTATIONS::Fill(ee_target_impulses,TV());
        ee_normals.Resize(pairs.Size());ARRAYS_COMPUTATIONS::Fill(ee_normals,TV());
        ee_old_speeds.Resize(pairs.Size());ARRAYS_COMPUTATIONS::Fill(ee_old_speeds,T());

        for(int pair_index=1;pair_index<=pairs.Size();pair_index++){
            const EDGE_EDGE_REPULSION_PAIR<TV>& pair=pairs(pair_index);
            const VECTOR<int,2*d-2>& nodes=pair.nodes;const VECTOR<T,2>& w=pair.weights;
            if(pair.distance<0) inverted_pairs++;

            TV relative_velocity=Edge_Edge_Relative_Velocity_Helper(w,V.Subset(nodes));
            TV direction;T scalar_impulse=Repulsion_Impulse(direction,dt,pair,relative_velocity,elastic_repulsion,friction);
            if(scalar_impulse){
                applied_impulses++;
                if(0 && geometry.mass_modifier && use_repulsions) geometry.mass_modifier->Edge_Edge_Mass((T)1.,pair.nodes,pair.weights,geometry.deformable_body_collection.particles.one_over_mass);
                T one_over_mass=Edge_Edge_One_Over_Mass_Helper(w,one_over_effective_mass.Subset(nodes));
                TV impulse=Pseudo_Divide(scalar_impulse*direction,one_over_mass);
                if(use_gauss_jacobi && !friction && !elastic_repulsion){
                    Edge_Edge_Update_Velocity_Helper(impulse,w,one_over_effective_mass.Subset(nodes),impulse_velocities.Subset(nodes));
                    ee_old_speeds(pair_index)=TV::Dot_Product(relative_velocity,direction);
                    ee_target_impulses(pair_index)=impulse;
                    ee_normals(pair_index)=direction.Normalized();} // something weird going on here (unnormalized directions)
                else
                    Edge_Edge_Update_Velocity_Helper(impulse,w,one_over_effective_mass.Subset(nodes),V.Subset(nodes));
                if(0 && geometry.mass_modifier && use_repulsions) geometry.mass_modifier->Edge_Edge_Mass_Revert(pair.nodes,geometry.deformable_body_collection.particles.one_over_mass);}}

        if(use_gauss_jacobi) Scale_And_Apply_Edge_Edge_Impulses(pairs);
    }

    if(inverted_pairs) LOG::Stat("inverted edge edge repulsion pairs",inverted_pairs);
    if(applied_impulses) LOG::Stat("applied edge edge repulsion impulses",applied_impulses);
}
template<class TV> template<class T_ARRAY> void TRIANGLE_REPULSIONS<TV>::
Adjust_Velocity_For_Edge_Edge_Repulsion(const T dt,const T_ARRAY& pairs,const bool elastic_repulsion,const bool friction,const bool use_repulsions)
{
    Adjust_Velocity_For_Edge_Edge_Repulsion_Helper(dt,pairs,elastic_repulsion,friction,TV(),use_repulsions);
}
//#####################################################################
// Function Project_All_Moving_Constraints
//#####################################################################
template<class TV> void TRIANGLE_REPULSIONS<TV>::
Project_All_Moving_Constraints(const ARRAY<PRECOMPUTE_PROJECT_POINT_FACE<VECTOR<T,3> > >& point_face_precomputed,
    const ARRAY<PRECOMPUTE_PROJECT_EDGE_EDGE<VECTOR<T,3> > >& edge_edge_precomputed,ARRAY_VIEW<VECTOR<T,3> >& field)
{
    for(int i=1;i<=point_face_precomputed.m;i++) point_face_precomputed(i).Project(field.Subset(point_face_precomputed(i).nodes));
    for(int i=1;i<=edge_edge_precomputed.m;i++) edge_edge_precomputed(i).Project(field.Subset(edge_edge_precomputed(i).nodes));
    for(int i=edge_edge_precomputed.m-1;i>=1;i--) edge_edge_precomputed(i).Project(field.Subset(edge_edge_precomputed(i).nodes));
    for(int i=point_face_precomputed.m;i>=1;i--) point_face_precomputed(i).Project(field.Subset(point_face_precomputed(i).nodes));
}
//#####################################################################
// Function Set_Collision_Pairs
//#####################################################################
template<class TV> template<class TV2> void TRIANGLE_REPULSIONS<TV>::
Set_Collision_Pairs(ARRAY<PRECOMPUTE_PROJECT_POINT_FACE<VECTOR<T,3> > >& point_face_precomputed,
    ARRAY<PRECOMPUTE_PROJECT_EDGE_EDGE<VECTOR<T,3> > >& edge_edge_precomputed,ARRAY<POINT_FACE_REPULSION_PAIR<VECTOR<T,3> > >& point_face_pairs,
    ARRAY<EDGE_EDGE_REPULSION_PAIR<TV2> >& edge_edge_pairs,const T repulsion_thickness_multiplier)
{
    STATIC_ASSERT((IS_SAME<VECTOR<T,3>,TV2>::value));
    point_face_pairs.Remove_All();
    edge_edge_pairs.Remove_All();
    point_face_pairs.Append_Elements(point_face_interaction_pairs);
    edge_edge_pairs.Append_Elements(edge_edge_interaction_pairs);
    repulsion_thickness*=repulsion_thickness_multiplier;
    Update_Repulsion_Pairs_Using_History(point_face_pairs,edge_edge_pairs,true); // TODO: Something is wrong here...
    repulsion_thickness/=repulsion_thickness_multiplier;
    point_face_precomputed.Resize(point_face_pairs.m);
    edge_edge_precomputed.Resize(edge_edge_pairs.m);
    for(int i=1;i<=point_face_pairs.m;i++){const POINT_FACE_REPULSION_PAIR<VECTOR<T,3> >& pr=point_face_pairs(i);
        point_face_precomputed(i).Precompute(geometry.deformable_body_collection.particles.mass.Subset(pr.nodes),pr.weights,pr.normal);}
    for(int i=1;i<=edge_edge_pairs.m;i++){const EDGE_EDGE_REPULSION_PAIR<VECTOR<T,3> >& pr=edge_edge_pairs(i);
        edge_edge_precomputed(i).Precompute(geometry.deformable_body_collection.particles.mass.Subset(pr.nodes),pr.weights,pr.normal);}
    internal_point_face_precomputed=point_face_precomputed;
    internal_edge_edge_precomputed=edge_edge_precomputed;
}
//#####################################################################
// Function Scale_And_Apply_Point_Face_Impulses
//#####################################################################
template<class TV> template<class T_ARRAY> void TRIANGLE_REPULSIONS<TV>::
Scale_And_Apply_Point_Face_Impulses(const T_ARRAY& pairs)
{
    // Go through the computed impulses, and compare them to the actual change in velocity seen.  Scale back impulses accordingly
    for(int i=1;i<=pf_target_impulses.m;i++){
        if(pf_target_impulses(i)==TV()) continue;
        const POINT_FACE_REPULSION_PAIR<TV>& pair=pairs(i);
        const VECTOR<int,d+1>& nodes=pair.nodes;
        // Compute actual new relative_speed
        int p=nodes[1];VECTOR<int,d> face_nodes=nodes.Remove_Index(1);
        T relative_speed=TV::Dot_Product(impulse_velocities(p)-ARRAYS_COMPUTATIONS::Weighted_Sum(pair.weights,impulse_velocities.Subset(face_nodes)),pf_normals(i));
        if(relative_speed*pf_old_speeds(i)<0){
            T new_scale=-(relative_speed-pf_old_speeds(i))/pf_old_speeds(i);
            pf_target_impulses(i)/=new_scale;}} // TODO: should we be doing this, or storing a maximum scale factor at each node?
    // Apply the newly scaled impulses
    ARRAY_VIEW<T>& one_over_effective_mass=geometry.deformable_body_collection.particles.one_over_effective_mass;
    ARRAY_VIEW<TV> V(geometry.deformable_body_collection.particles.V);
    for(int i=1;i<=pf_target_impulses.m;i++){
        if(pf_target_impulses(i)==TV()) continue;
        const POINT_FACE_REPULSION_PAIR<TV>& pair=pairs(i);const VECTOR<int,d+1>& nodes=pair.nodes;
        Update_Velocity_Helper(pf_target_impulses(i),pair.weights,one_over_effective_mass.Subset(nodes),V.Subset(nodes));}
}
//#####################################################################
// Function Scale_And_Apply_Edge_Edge_Impulses
//#####################################################################
template<class TV> template<class T_ARRAY> void TRIANGLE_REPULSIONS<TV>::
Scale_And_Apply_Edge_Edge_Impulses(const T_ARRAY& pairs)
{
    // Go through the computed impulses, and compare them to the actual change in velocity seen.  Scale back impulses accordingly
    for(int i=1;i<=ee_target_impulses.m;i++){
        if(ee_target_impulses(i)==TV()) continue;
        const VECTOR<int,2*d-2>& nodes=pairs(i).nodes;
        // Compute actual new relative_speed
        TV relative_velocity;
        T relative_speed=TV::Dot_Product(Edge_Edge_Relative_Velocity_Helper(pairs(i).weights,impulse_velocities.Subset(nodes)),ee_normals(i));
        if(relative_speed*ee_old_speeds(i)<0){
            T new_scale=-(relative_speed-ee_old_speeds(i))/ee_old_speeds(i);
            ee_target_impulses(i)/=new_scale;}}
    // Apply the newly scaled impulses
    ARRAY_VIEW<T>& one_over_effective_mass=geometry.deformable_body_collection.particles.one_over_effective_mass;
    ARRAY_VIEW<TV> V(geometry.deformable_body_collection.particles.V);
    for(int i=1;i<=ee_target_impulses.m;i++){
        if(ee_target_impulses(i)==TV()) continue;
        const VECTOR<int,2*d-2>& nodes=pairs(i).nodes;
        Edge_Edge_Update_Velocity_Helper(ee_target_impulses(i),pairs(i).weights,one_over_effective_mass.Subset(nodes),V.Subset(nodes));}
}
//#####################################################################
// Precompute
//#####################################################################
template<class TV> void PRECOMPUTE_PROJECT_POINT_FACE<TV>::
Precompute(const INDIRECT_ARRAY<ARRAY_VIEW<T>,VECTOR<int,4>&> mass,const TV& weights_input,const TV& normal_input)
{
    weights=weights_input;normal=normal_input;
    T tau=1/(1/mass(1)+weights(1)*weights(1)/mass(2)+weights(2)*weights(2)/mass(3)+weights(3)*weights(3)/mass(4));
    v_scaled_normals(1)=-tau/mass(1)*normal;
    v_scaled_normals(2)=tau*weights(1)/mass(2)*normal;
    v_scaled_normals(3)=tau*weights(2)/mass(3)*normal;
    v_scaled_normals(4)=tau*weights(3)/mass(4)*normal;
    nodes=mass.indices;
}
//#####################################################################
// Precompute
//#####################################################################
template<class TV> void PRECOMPUTE_PROJECT_EDGE_EDGE<TV>::
Precompute(const INDIRECT_ARRAY<ARRAY_VIEW<T>,VECTOR<int,4>&> mass,const VECTOR<T,2>& weights_input,const TV& normal_input)
{
    weights=weights_input;normal=normal_input;
    T tau=1/(sqr(1-weights(1))/mass(1)+sqr(weights(1))/mass(2)+sqr(1-weights(2))/mass(3)+sqr(weights(2))/mass(4));
    v_scaled_normals(1)=-tau/mass(1)*(1-weights(1))*normal;
    v_scaled_normals(2)=-tau/mass(2)*weights(1)*normal;
    v_scaled_normals(3)=tau/mass(3)*(1-weights(2))*normal;
    v_scaled_normals(4)=tau/mass(4)*weights(2)*normal;
    normal_13=normal;normal_12=weights(1)*normal;normal_34=-weights(2)*normal;
    nodes=mass.indices;
}
//####################################################################
template<> VECTOR<float,2> EDGE_EDGE_REPULSION_PAIR<VECTOR<float,2> >::weights=VECTOR<float,2>();
template class TRIANGLE_REPULSIONS<VECTOR<float,1> >;
template class TRIANGLE_REPULSIONS<VECTOR<float,2> >;
template class TRIANGLE_REPULSIONS<VECTOR<float,3> >;
template void TRIANGLE_REPULSIONS<VECTOR<float,3> >::Set_Collision_Pairs<VECTOR<float,3> >(ARRAY<PRECOMPUTE_PROJECT_POINT_FACE<VECTOR<float,3> >,int>&,
    ARRAY<PRECOMPUTE_PROJECT_EDGE_EDGE<VECTOR<float,3> >,int>&,ARRAY<POINT_FACE_REPULSION_PAIR<VECTOR<float,3> >,int>&,
    ARRAY<EDGE_EDGE_REPULSION_PAIR<VECTOR<float,3> >,int>&,float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template<> VECTOR<double,2> EDGE_EDGE_REPULSION_PAIR<VECTOR<double,2> >::weights=VECTOR<double,2>();
template class TRIANGLE_REPULSIONS<VECTOR<double,1> >;
template class TRIANGLE_REPULSIONS<VECTOR<double,2> >;
template class TRIANGLE_REPULSIONS<VECTOR<double,3> >;
template void TRIANGLE_REPULSIONS<VECTOR<double,3> >::Set_Collision_Pairs<VECTOR<double,3> >(ARRAY<PRECOMPUTE_PROJECT_POINT_FACE<VECTOR<double,3> >,int>&,
    ARRAY<PRECOMPUTE_PROJECT_EDGE_EDGE<VECTOR<double,3> >,int>&,ARRAY<POINT_FACE_REPULSION_PAIR<VECTOR<double,3> >,int>&,
    ARRAY<EDGE_EDGE_REPULSION_PAIR<VECTOR<double,3> >,int>&,double);
#endif
}
