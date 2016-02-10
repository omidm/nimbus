//#####################################################################
// Copyright 2005-2006, Zhaosheng Bao, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRACTURE_TRIANGULATED_OBJECT
//##################################################################### 
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_TRIANGULATED_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/HYPOTHETICAL_CUT_TRIANGLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Add_First_Embedded_Segment
//#####################################################################
template<class TV> void FRACTURE_TRIANGULATED_OBJECT<TV>::
Add_First_Embedded_Segment(const int triangle,const TV& fracture_normal)
{
    if(Fracture_Phi_Index(triangle)) Add_First_Cut_Based_On_Phi(triangle);
    else{
        HYPOTHETICAL_CUT_TRIANGLES<TV> hypothetical_cut(embedded_object);
        if(Add_Intersected_Points_To_Embedded_Triangulated_Object(triangle,fracture_normal,hypothetical_cut))
            embedded_object.Add_Embedded_Segment(hypothetical_cut.hypothetical_nodes(1).index_in_embedded_particles,
                hypothetical_cut.hypothetical_nodes(2).index_in_embedded_particles);}
}
//#####################################################################
// Function Add_First_Cut_Based_On_Phi
//#####################################################################
template<class TV> void FRACTURE_TRIANGULATED_OBJECT<TV>::
Add_First_Cut_Based_On_Phi(const int triangle)
{
    assert(embedded_object.Number_Of_Embedded_Cuts(triangle)==0);
    int i,j,k;embedded_object.simplicial_object.mesh.elements(triangle).Get(i,j,k);
    T phi1,phi2,phi3;Get_Phi(triangle).Get(phi1,phi2,phi3);
    int positive_count=(phi1>0)+(phi2>0)+(phi3>0);
    int endpoint1,endpoint2;
    if(positive_count==1){
        if(phi1>0){
            T lambda=LEVELSET_UTILITIES<T>::Theta((T)phi2,phi1);endpoint1=embedded_object.Add_Embedded_Particle_If_Not_Already_There(j,i,lambda);
            lambda=LEVELSET_UTILITIES<T>::Theta((T)phi3,phi1);endpoint2=embedded_object.Add_Embedded_Particle_If_Not_Already_There(k,i,lambda);}
        else if(phi2>0){
            T lambda=LEVELSET_UTILITIES<T>::Theta((T)phi1,phi2);endpoint1=embedded_object.Add_Embedded_Particle_If_Not_Already_There(i,j,lambda);
            lambda=LEVELSET_UTILITIES<T>::Theta((T)phi3,phi2);endpoint2=embedded_object.Add_Embedded_Particle_If_Not_Already_There(k,j,lambda);}
        else{assert(phi3>0);
            T lambda=LEVELSET_UTILITIES<T>::Theta((T)phi1,phi3);endpoint1=embedded_object.Add_Embedded_Particle_If_Not_Already_There(i,k,lambda);
            lambda=LEVELSET_UTILITIES<T>::Theta((T)phi2,phi3);endpoint2=embedded_object.Add_Embedded_Particle_If_Not_Already_There(j,k,lambda);}}
    else{assert(positive_count==2);
        if(phi1<=0){
            T lambda=LEVELSET_UTILITIES<T>::Theta((T)phi1,phi2);endpoint1=embedded_object.Add_Embedded_Particle_If_Not_Already_There(i,j,lambda);
            lambda=LEVELSET_UTILITIES<T>::Theta((T)phi1,phi3);endpoint2=embedded_object.Add_Embedded_Particle_If_Not_Already_There(i,k,lambda);}
        else if(phi2<=0){
            T lambda=LEVELSET_UTILITIES<T>::Theta((T)phi2,phi1);endpoint1=embedded_object.Add_Embedded_Particle_If_Not_Already_There(j,i,lambda);
            lambda=LEVELSET_UTILITIES<T>::Theta((T)phi2,phi3);endpoint2=embedded_object.Add_Embedded_Particle_If_Not_Already_There(j,k,lambda);}
        else{assert(phi3<=0);
            T lambda=LEVELSET_UTILITIES<T>::Theta((T)phi3,phi1);endpoint1=embedded_object.Add_Embedded_Particle_If_Not_Already_There(k,i,lambda);
            lambda=LEVELSET_UTILITIES<T>::Theta((T)phi3,phi2);endpoint2=embedded_object.Add_Embedded_Particle_If_Not_Already_There(k,j,lambda);}}
    embedded_object.Add_Embedded_Segment(endpoint1,endpoint2);
}
//#####################################################################
// Function Add_Second_Embedded_Segment
//#####################################################################
template<class TV> bool
Cut_31_Better_Than_Cut_32(const TV& x1,const TV& x2,const TV& x3,const TV& fracture_normal)
{
    return TV::Cross_Product((x3-x1).Normalized(),fracture_normal).Magnitude_Squared()>TV::Cross_Product((x3-x2).Normalized(),fracture_normal).Magnitude_Squared();
}
template<class TV> void FRACTURE_TRIANGULATED_OBJECT<TV>::
Add_Second_Embedded_Segment(const int triangle,const TV& fracture_normal)
{
    if(Fracture_Phi_Index(triangle)){
        int isolated_node=embedded_object.Isolated_Node(triangle);
        if(Phi_In_Simplex(isolated_node,triangle) > 0) return;} // we can't fracture (the quad is inside the levelset)
   
    ARRAY<int>& active_indices=embedded_object.embedded_particles.active_indices;
    HYPOTHETICAL_CUT_TRIANGLES<TV> hypothetical_cut(embedded_object);
    if(Add_Intersected_Points_To_Embedded_Triangulated_Object(triangle,fracture_normal,hypothetical_cut)){
        int embedded_particle_1=embedded_object.Embedded_Particle_On_Segment(hypothetical_cut.hypothetical_nodes(1).parents),
            embedded_particle_2=embedded_object.Embedded_Particle_On_Segment(hypothetical_cut.hypothetical_nodes(2).parents);
        if(embedded_particle_1 && embedded_particle_2 && embedded_object.embedded_mesh.Segment(active_indices(embedded_particle_1),active_indices(embedded_particle_2))){
            int embedded_particle_3=0;
            if(embedded_object.Is_Parent(hypothetical_cut.hypothetical_nodes(1).parents[1],embedded_particle_2)){
                int node_3=embedded_object.Other_Parent(hypothetical_cut.hypothetical_nodes(1).parents[1],embedded_particle_2);
                embedded_particle_3=embedded_object.Embedded_Particle_On_Segment(hypothetical_cut.hypothetical_nodes(1).parents[2],node_3);}
            else{
                int node_3=embedded_object.Other_Parent(hypothetical_cut.hypothetical_nodes(1).parents[2],embedded_particle_2);
                embedded_particle_3=embedded_object.Embedded_Particle_On_Segment(hypothetical_cut.hypothetical_nodes(1).parents[1],node_3);}
            if(embedded_particle_3){
                assert(!embedded_object.embedded_mesh.Segment(active_indices(embedded_particle_3),active_indices(embedded_particle_1))
                    && !embedded_object.embedded_mesh.Segment(active_indices(embedded_particle_3),active_indices(embedded_particle_2)));
                if(Cut_31_Better_Than_Cut_32(embedded_object.Position_Of_Embedded_Particle(embedded_particle_1),embedded_object.Position_Of_Embedded_Particle(embedded_particle_2),
                    embedded_object.Position_Of_Embedded_Particle(embedded_particle_3),fracture_normal))
                    embedded_object.Add_Embedded_Segment(embedded_particle_3,hypothetical_cut.hypothetical_nodes(1).index_in_embedded_particles);
                else embedded_object.Add_Embedded_Segment(embedded_particle_3,hypothetical_cut.hypothetical_nodes(2).index_in_embedded_particles);}}
        else embedded_object.Add_Embedded_Segment(hypothetical_cut.hypothetical_nodes(1).index_in_embedded_particles,
            hypothetical_cut.hypothetical_nodes(2).index_in_embedded_particles);}
}
//#####################################################################
// Function Add_Intersected_Points_To_Embedded_Triangulated_Object
//#####################################################################
template<class TV> int FRACTURE_TRIANGULATED_OBJECT<TV>::
Add_Intersected_Points_To_Embedded_Triangulated_Object(const int triangle,const TV& fracture_normal,HYPOTHETICAL_CUT_TRIANGLES<TV>& hypothetical_cut)
{
    HYPOTHETICAL_CUT_TRIANGLES<TV> hc_ij(embedded_object),hc_ik(embedded_object),hc_jk(embedded_object),*best_hc=0;
    int i,j,k;embedded_object.simplicial_object.mesh.elements(triangle).Get(i,j,k);
    T best_cut_quality=(T)-FLT_MAX;
    if(embedded_object.Embedded_Particle_On_Segment(i,j) &&
        hc_ij.Initialize_Hypothetical_Cut(T_HYPERPLANE(fracture_normal,embedded_object.Position_Of_Embedded_Particle(i,j)),triangle)){
        T cut_quality=hc_ij.Quality_Of_Cut();if(cut_quality > best_cut_quality){best_hc=&hc_ij;best_cut_quality=cut_quality;}}
    if(embedded_object.Embedded_Particle_On_Segment(i,k) &&
        hc_ik.Initialize_Hypothetical_Cut(T_HYPERPLANE(fracture_normal,embedded_object.Position_Of_Embedded_Particle(i,k)),triangle)){
        T cut_quality=hc_ik.Quality_Of_Cut();if(cut_quality > best_cut_quality){best_hc=&hc_ik;best_cut_quality=cut_quality;}}
    if(embedded_object.Embedded_Particle_On_Segment(j,k) &&
        hc_jk.Initialize_Hypothetical_Cut(T_HYPERPLANE(fracture_normal,embedded_object.Position_Of_Embedded_Particle(j,k)),triangle)){
        T cut_quality=hc_jk.Quality_Of_Cut();if(cut_quality > best_cut_quality){best_hc=&hc_jk;best_cut_quality=cut_quality;}}

    // control -- only add if continues an existing fracture (edge-connected)
    if(force_edge_connected_fracture && best_hc && !best_hc->Number_Of_Nodes_Shared_With_Existing_Embedded_Curve()) best_hc=0;

    // add the centroid one as last resort 
    HYPOTHETICAL_CUT_TRIANGLES<TV> hc_centroid(embedded_object);
    if(!best_hc && Initiation_Point(triangle)){
        best_hc=&hc_centroid;best_hc->Initialize_Hypothetical_Cut(T_HYPERPLANE(fracture_normal,embedded_object.simplicial_object.Centroid(triangle)),triangle);}
    if(!best_hc||!best_hc->Valid_Cut()) return 0;
    best_hc->Add_Hypothetical_Nodes_To_Embedded_Object(embedded_object);
    hypothetical_cut=(*best_hc);
    return hypothetical_cut.hypothetical_nodes.m;
}
//#####################################################################
// Function Add_Cut
//#####################################################################
template<class TV> void FRACTURE_TRIANGULATED_OBJECT<TV>::
Add_Cut(const int triangle,const TV& fracture_normal)
{
    int cut=embedded_object.Number_Of_Embedded_Cuts(triangle)+1;
    if(cut==1) Add_First_Embedded_Segment(triangle,fracture_normal);
    else if(cut==2) Add_Second_Embedded_Segment(triangle,fracture_normal);
    else PHYSBAM_FATAL_ERROR();
}
//#####################################################################
template class FRACTURE_TRIANGULATED_OBJECT<VECTOR<float,2> >;
template class FRACTURE_TRIANGULATED_OBJECT<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FRACTURE_TRIANGULATED_OBJECT<VECTOR<double,2> >;
template class FRACTURE_TRIANGULATED_OBJECT<VECTOR<double,3> >;
#endif
