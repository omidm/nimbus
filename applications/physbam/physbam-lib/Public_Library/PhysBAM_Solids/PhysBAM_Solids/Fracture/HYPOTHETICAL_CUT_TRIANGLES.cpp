//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Neil Molino, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HYPOTHETICAL_CUT_TRIANGLES
//##################################################################### 
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TRIANGULATED_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/HYPOTHETICAL_CUT_TRIANGLES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/HYPOTHETICAL_NODE.h>
using namespace PhysBAM;
//#####################################################################
// Function operator=
//##################################################################### 
template<class TV> HYPOTHETICAL_CUT_TRIANGLES<TV>& HYPOTHETICAL_CUT_TRIANGLES<TV>::
operator=(const HYPOTHETICAL_CUT_TRIANGLES<TV>& old_cut)
{
    hypothetical_nodes=old_cut.hypothetical_nodes;
    cut_quality_metric=old_cut.cut_quality_metric;
    cuts_ij=old_cut.cuts_ij;cuts_ik=old_cut.cuts_ik;cuts_jk=old_cut.cuts_jk;
    interpolation_fraction_ij=old_cut.interpolation_fraction_ij;interpolation_fraction_ik=old_cut.interpolation_fraction_ik;interpolation_fraction_jk=old_cut.interpolation_fraction_jk;
    triangle=old_cut.triangle;hyperplane=old_cut.hyperplane;
    return *this;
}
//#####################################################################
// Function Initialize_Hypothetical_Cut
//##################################################################### 
template<class TV> void HYPOTHETICAL_CUT_TRIANGLES<TV>::
Initialize_Quality_Metric()
{
    if(!Valid_Cut()){cut_quality_metric=0;return;}
    TV hypothetical_cut_direction=(Position(2)-Position(1)).Normalized();
    cut_quality_metric=TV::Cross_Product(hyperplane.normal,hypothetical_cut_direction).Magnitude();
}
template<class TV> bool HYPOTHETICAL_CUT_TRIANGLES<TV>::
Initialize_Hypothetical_Cut(const T_HYPERPLANE& hyperplane_input,const int triangle_input)
{
    hyperplane=hyperplane_input;triangle=triangle_input;
    int i,j,k;embedded_object.simplicial_object.mesh.elements(triangle).Get(i,j,k);
    TV xi=embedded_object.particles.X(i),xj=embedded_object.particles.X(j),xk=embedded_object.particles.X(k);
    Compute_Cuts(xi,xj,xk);
    if(cuts_ij) Add_Hypothetical_Node(i,j,interpolation_fraction_ij);
    if(cuts_ik) Add_Hypothetical_Node(i,k,interpolation_fraction_ik);
    if(cuts_jk) Add_Hypothetical_Node(j,k,interpolation_fraction_jk);
    Initialize_Quality_Metric();
    return Valid_Cut();
}
//#####################################################################
// Function Segment_Cut_Already_In_Embedded_Triangulated_Object
//##################################################################### 
template<class TV> bool HYPOTHETICAL_CUT_TRIANGLES<TV>::
Segment_Cut_Already_In_Embedded_Triangulated_Object()
{
    int emb_node1=embedded_object.Embedded_Particle_On_Segment(hypothetical_nodes(1).parents);
    int emb_node2=embedded_object.Embedded_Particle_On_Segment(hypothetical_nodes(2).parents);  
    const ARRAY<int>& active_indices=embedded_object.embedded_particles.active_indices;
    return embedded_object.embedded_mesh.Segment(active_indices(emb_node1),active_indices(emb_node2))!=0;
}
//#####################################################################
// Function Number_Of_Nodes_Shared_With_Existing_Embedded_Curve
//##################################################################### 
template<class TV> int HYPOTHETICAL_CUT_TRIANGLES<TV>::
Number_Of_Nodes_Shared_With_Existing_Embedded_Curve()
{
    int number_of_nodes_shared=0;
    int i=hypothetical_nodes(1).index_in_embedded_particles,j=hypothetical_nodes(2).index_in_embedded_particles;
    const ARRAY<int>& active_indices=embedded_object.embedded_particles.active_indices;
    if(i && (*embedded_object.embedded_mesh.incident_elements)(active_indices(i)).m) number_of_nodes_shared++;
    if(j && (*embedded_object.embedded_mesh.incident_elements)(active_indices(j)).m) number_of_nodes_shared++;
    return number_of_nodes_shared;
}
//##################################################################### 
template class HYPOTHETICAL_CUT_TRIANGLES<VECTOR<float,2> >;
template class HYPOTHETICAL_CUT_TRIANGLES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class HYPOTHETICAL_CUT_TRIANGLES<VECTOR<double,2> >;
template class HYPOTHETICAL_CUT_TRIANGLES<VECTOR<double,3> >;
#endif
