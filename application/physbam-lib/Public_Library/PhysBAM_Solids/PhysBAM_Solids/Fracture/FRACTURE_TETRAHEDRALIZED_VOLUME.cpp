//#####################################################################
// Copyright 2004-2007, Zhaosheng Bao, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRACTURE_TETRAHEDRALIZED_VOLUME
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Math_Tools/cyclic_shift.h>
#include <PhysBAM_Tools/Math_Tools/permutation.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/HYPOTHETICAL_CUT_TETRAHEDRONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/HYPOTHETICAL_NODE.h>
using namespace PhysBAM;
//#####################################################################
// Function Add_First_Cut
//#####################################################################
template<class T> void FRACTURE_TETRAHEDRALIZED_VOLUME<T>::
Add_First_Cut(const int tetrahedron,const TV& fracture_normal)
{
    assert((embedded_object.Embedded_Subelements_In_Element(tetrahedron)==VECTOR<int,4>()));
    VECTOR<int,4> nodes=embedded_object.simplicial_object.mesh.elements(tetrahedron);
    if(Fracture_Phi_Index(tetrahedron)){Add_First_Cut_Based_On_Phi(tetrahedron,Get_Phi(tetrahedron));return;}
    HYPOTHETICAL_CUT_TETRAHEDRONS<T> hypothetical_cut(embedded_object);
    int segments_intersected=Add_Intersected_Points_To_Embedded_Tetrahedralized_Volume(tetrahedron,fracture_normal,hypothetical_cut);
    if(segments_intersected==3){
        int count=1,i,j,k,l;nodes.Get(i,j,k,l);
        while(!Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(i,j))
            || !Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(i,k))
            || !Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(i,l))){
                cyclic_shift(i,j,k,l);next_cyclic_permutation_of_four_index(count);}
        embedded_object.orientation_index(tetrahedron)=count;
        embedded_object.Add_Embedded_Triangle(i,j,i,k,i,l);}
    else if(segments_intersected==4){
        int count=1,i,j,k,l;nodes.Get(i,j,k,l);
        while(!Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(i,k))
            || !Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(k,j))
            || !Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(j,l))
            || !Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(l,i))
            || !Quad_Diagonal_Is_Correct(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(i,k),embedded_object.Embedded_Particle_On_Segment(j,l)))
            permute_four(nodes,++count).Get(i,j,k,l);
        embedded_object.orientation_index(tetrahedron)=count;
        embedded_object.Add_Embedded_Triangle(i,k,k,j,j,l);
        embedded_object.Add_Embedded_Triangle(j,l,l,i,i,k);}
}
//#####################################################################
// Function Add_First_Cut_Based_On_Phi
//#####################################################################
template<class T> void FRACTURE_TETRAHEDRALIZED_VOLUME<T>::
Add_First_Cut_Based_On_Phi(const int tetrahedron,const VECTOR<T,4>& tetrahedron_phi)
{
    VECTOR<int,4> nodes=embedded_object.simplicial_object.mesh.elements(tetrahedron);
    int positive_count=0;for(int i=1;i<=4;i++)if(tetrahedron_phi(i)>0)positive_count++;
    assert(positive_count==1 || positive_count==2 || positive_count==3);
    int i,j,k,l,pi=1,pj=2,pk=3,pl=4;int count=1;
    if(positive_count==2){
        while(!(tetrahedron_phi(pi)>0) || !(tetrahedron_phi(pj)>0)) permute_four(VECTOR<int,4>(1,2,3,4),++count).Get(pi,pj,pk,pl);
        permute_four(nodes,count).Get(i,j,k,l);embedded_object.orientation_index(tetrahedron)=count;
        T phi1=tetrahedron_phi(pi),phi2=tetrahedron_phi(pj),phi3=tetrahedron_phi(pk),phi4=tetrahedron_phi(pl);
        int i_k=embedded_object.Add_Embedded_Particle_If_Not_Already_There(i,k,LEVELSET_UTILITIES<T>::Theta(phi1,phi3));
        int i_l=embedded_object.Add_Embedded_Particle_If_Not_Already_There(i,l,LEVELSET_UTILITIES<T>::Theta(phi1,phi4));
        int j_k=embedded_object.Add_Embedded_Particle_If_Not_Already_There(j,k,LEVELSET_UTILITIES<T>::Theta(phi2,phi3));
        int j_l=embedded_object.Add_Embedded_Particle_If_Not_Already_There(j,l,LEVELSET_UTILITIES<T>::Theta(phi2,phi4));
        embedded_object.Add_Embedded_Triangle(i_k,j_k,j_l);
        embedded_object.Add_Embedded_Triangle(j_l,i_l,i_k);}
    else{ // count is 1 or 3
        if(positive_count==1)while(!(tetrahedron_phi(pi)>0)) permute_four(VECTOR<int,4>(1,2,3,4),++count).Get(pi,pj,pk,pl);
        else while(!(tetrahedron_phi(pi)<=0)) permute_four(VECTOR<int,4>(1,2,3,4),++count).Get(pi,pj,pk,pl);
        permute_four(nodes,count).Get(i,j,k,l);embedded_object.orientation_index(tetrahedron)=count;
        T phi1=tetrahedron_phi(pi),phi2=tetrahedron_phi(pj),phi3=tetrahedron_phi(pk),phi4=tetrahedron_phi(pl);
        int i_j=embedded_object.Add_Embedded_Particle_If_Not_Already_There(i,j,LEVELSET_UTILITIES<T>::Theta(phi1,phi2));
        int i_k=embedded_object.Add_Embedded_Particle_If_Not_Already_There(i,k,LEVELSET_UTILITIES<T>::Theta(phi1,phi3));
        int i_l=embedded_object.Add_Embedded_Particle_If_Not_Already_There(i,l,LEVELSET_UTILITIES<T>::Theta(phi1,phi4));
        embedded_object.Add_Embedded_Triangle(i_j,i_k,i_l);}
}
//#####################################################################
// Function Add_Second_Cut_Based_On_Phi
//#####################################################################
template<class T> void FRACTURE_TETRAHEDRALIZED_VOLUME<T>::
Add_Next_Cut_Based_On_Phi(const int tetrahedron,const VECTOR<T,4>& tetrahedron_phi,const bool second_cut)
{
    VECTOR<int,4> nodes=embedded_object.simplicial_object.mesh.elements(tetrahedron);
    int positive_count=0;for(int i=1;i<=4;i++)if(tetrahedron_phi(i)>0)positive_count++;
    assert(positive_count==1 || positive_count==2 || positive_count==3);
    int i,j,k,l,pi=1,pj=2,pk=3,pl=4;int count=1;
    if(positive_count==2){
        while(!(tetrahedron_phi(pi)>0) || !(tetrahedron_phi(pj)>0)) permute_four(VECTOR<int,4>(1,2,3,4),++count).Get(pi,pj,pk,pl);
        permute_four(nodes,count).Get(i,j,k,l);T phi1=tetrahedron_phi(pi),phi2=tetrahedron_phi(pj),phi3=tetrahedron_phi(pk),phi4=tetrahedron_phi(pl);
        embedded_object.Add_Embedded_Particle_If_Not_Already_There(i,k,LEVELSET_UTILITIES<T>::Theta(phi1,phi3));
        embedded_object.Add_Embedded_Particle_If_Not_Already_There(i,l,LEVELSET_UTILITIES<T>::Theta(phi1,phi4));
        embedded_object.Add_Embedded_Particle_If_Not_Already_There(j,k,LEVELSET_UTILITIES<T>::Theta(phi2,phi3));
        embedded_object.Add_Embedded_Particle_If_Not_Already_There(j,l,LEVELSET_UTILITIES<T>::Theta(phi2,phi4));
        VECTOR<T,3> xik=embedded_object.Position_Of_Embedded_Particle(i,k),xil=embedded_object.Position_Of_Embedded_Particle(i,l);
        VECTOR<T,3> xjk=embedded_object.Position_Of_Embedded_Particle(j,k),xjl=embedded_object.Position_Of_Embedded_Particle(j,l);
        if(second_cut) Add_Second_Cut(tetrahedron,(TRIANGLE_3D<T>::Normal(xik,xjk,xjl)+TRIANGLE_3D<T>::Normal(xjl,xil,xik))/(T)2);
        else Add_Third_Cut(tetrahedron,(TRIANGLE_3D<T>::Normal(xik,xjk,xjl)+TRIANGLE_3D<T>::Normal(xjl,xil,xik))/(T)2);}
    else{ // count is 1 or 3
        if(positive_count==1)while(!(tetrahedron_phi(pi)>0)) permute_four(VECTOR<int,4>(1,2,3,4),++count).Get(pi,pj,pk,pl);
        else while(!(tetrahedron_phi(pi)<=0)) permute_four(VECTOR<int,4>(1,2,3,4),++count).Get(pi,pj,pk,pl);
        permute_four(nodes,count).Get(i,j,k,l);T phi1=tetrahedron_phi(pi),phi2=tetrahedron_phi(pj),phi3=tetrahedron_phi(pk),phi4=tetrahedron_phi(pl);
        embedded_object.Add_Embedded_Particle_If_Not_Already_There(i,j,LEVELSET_UTILITIES<T>::Theta(phi1,phi2));
        embedded_object.Add_Embedded_Particle_If_Not_Already_There(i,k,LEVELSET_UTILITIES<T>::Theta(phi1,phi3));
        embedded_object.Add_Embedded_Particle_If_Not_Already_There(i,l,LEVELSET_UTILITIES<T>::Theta(phi1,phi4));
        VECTOR<T,3> xij=embedded_object.Position_Of_Embedded_Particle(i,j),xik=embedded_object.Position_Of_Embedded_Particle(i,k),xil=embedded_object.Position_Of_Embedded_Particle(i,l);
        if(second_cut) Add_Second_Cut(tetrahedron,TRIANGLE_3D<T>::Normal(xij,xik,xil));
        else Add_Third_Cut(tetrahedron,TRIANGLE_3D<T>::Normal(xij,xik,xil));}
}
//#####################################################################
// Function Add_Intersected_Points_To_Embedded_Tetrahedralized_Volume
//#####################################################################
template<class T> int FRACTURE_TETRAHEDRALIZED_VOLUME<T>::
Add_Intersected_Points_To_Embedded_Tetrahedralized_Volume(const int tetrahedron,const TV& fracture_normal,HYPOTHETICAL_CUT_TETRAHEDRONS<T>& hypothetical_cut)
{
    T best_cut_quality=(T)-FLT_MAX;
    for(int a=1;a<=7;a++){
        HYPOTHETICAL_CUT_TETRAHEDRONS<T> hc(embedded_object);hc.Initialize_Hypothetical_Cut(fracture_normal,a,tetrahedron);
        if(!hc.hypothetical_nodes.m) continue;
        if(!hc.Edges_Shared_With_Existing_Embedded_Surface()) continue;
        if(hc.Cut_Already_Exists()) continue;
        T cut_quality=hc.Quality_Of_Cut();
        if(cut_quality>best_cut_quality){best_cut_quality=cut_quality;hypothetical_cut=hc;}}
    if(hypothetical_cut.hypothetical_nodes.m==0 && Initiation_Point(tetrahedron)){
        HYPOTHETICAL_CUT_TETRAHEDRONS<T> hc(embedded_object);
        hc.Initialize_Initiation_Point_Cut(PLANE<T>(fracture_normal,embedded_object.simplicial_object.Centroid(tetrahedron)),tetrahedron);
        hypothetical_cut=hc;best_cut_quality=hc.Quality_Of_Cut();}
    if(best_cut_quality<fracture_quality_threshold) return 0;        
    hypothetical_cut.Add_Hypothetical_Nodes_To_Embedded_Object(embedded_object);
    return hypothetical_cut.hypothetical_nodes.m;
}
//#####################################################################
// Function Emb_Node_Added_By_This_Tet
//#####################################################################
template<class T> bool FRACTURE_TETRAHEDRALIZED_VOLUME<T>::
Emb_Node_Added_By_This_Tet(const HYPOTHETICAL_CUT_TETRAHEDRONS<T>& hypothetical_cut,const int emb_node) const
{
    if(emb_node==0) return false;
    else return hypothetical_cut.Contains_Embedded_Particle(embedded_object,emb_node);
}
//#####################################################################
// Function Quad_Diagonal_Is_Correct
//#####################################################################
template<class T> bool FRACTURE_TETRAHEDRALIZED_VOLUME<T>::
Quad_Diagonal_Is_Correct(const HYPOTHETICAL_CUT_TETRAHEDRONS<T>& hypothetical_cut,const int diagonal_embedded_particle1,const int diagonal_embedded_particle2) const
{
    assert(hypothetical_cut.hypothetical_nodes.m==4);
    int node1=hypothetical_cut.hypothetical_nodes(hypothetical_cut.quad_diagonal_indices.x).index_in_embedded_particles;
    int node2=hypothetical_cut.hypothetical_nodes(hypothetical_cut.quad_diagonal_indices.y).index_in_embedded_particles;
    return (diagonal_embedded_particle1==node1 && diagonal_embedded_particle2==node2) || (diagonal_embedded_particle1==node2 && diagonal_embedded_particle2==node1);
}
//#####################################################################
// Function Add_Second_Cut
//#####################################################################
template<class T> void FRACTURE_TETRAHEDRALIZED_VOLUME<T>::
Add_Second_Cut(const int tetrahedron,const TV& fracture_normal,const VECTOR<T,4>* tetrahedron_phi)
{
    VECTOR<int,4> nodes=embedded_object.simplicial_object.mesh.elements(tetrahedron);
    VECTOR<int,4> emb_triangles=embedded_object.Embedded_Subelements_In_Element(tetrahedron);
    int emb_triangle1=emb_triangles[1],emb_triangle2=emb_triangles[2];assert(emb_triangles[3]==0 && emb_triangles[4]==0);
    if(emb_triangle1 && !emb_triangle2){
        int isolated_node=embedded_object.Node_Separated_By_Embedded_Subelement(emb_triangle1);
        int other_node1,other_node2,other_node3;embedded_object.simplicial_object.mesh.Other_Three_Nodes(isolated_node,tetrahedron,other_node1,other_node2,other_node3);
        if(!embedded_object.Node_In_Simplex_Is_Material(other_node1,tetrahedron)) return;
        if(!tetrahedron_phi && Fracture_Phi_Index(tetrahedron) && Phi_In_Simplex(other_node1,tetrahedron)<=0) return;
        HYPOTHETICAL_CUT_TETRAHEDRONS<T> hypothetical_cut(embedded_object);
        int segments_intersected=Add_Intersected_Points_To_Embedded_Tetrahedralized_Volume(tetrahedron,fracture_normal,hypothetical_cut);
        if(segments_intersected==3){
            assert(!(hypothetical_cut.Triangle_Cut_Already_In_Embedded_Tetrahedralized_Volume()));
            int count=1,i,j,k,l;nodes.Get(i,j,k,l);
            while(!Emb_Node_In_Triangle(embedded_object,embedded_object.Embedded_Particle_On_Segment(i,k),emb_triangle1)
                || !Emb_Node_In_Triangle(embedded_object,embedded_object.Embedded_Particle_On_Segment(i,l),emb_triangle1)
                || !Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(j,k))
                || !Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(j,l))
                || !Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(j,i)))
                permute_four(nodes,++count).Get(i,j,k,l);
            embedded_object.orientation_index(tetrahedron)=count;
            embedded_object.Add_Embedded_Triangle(j,k,j,l,j,i);}
        else if(segments_intersected==4){
            assert(!(hypothetical_cut.Quad_Cut_Already_In_Embedded_Tetrahedralized_Volume()));
            int count=1,i,j,k,l;nodes.Get(i,j,k,l);
            while(!Emb_Node_In_Triangle(embedded_object,embedded_object.Embedded_Particle_On_Segment(i,j),emb_triangle1)
                || !Emb_Node_In_Triangle(embedded_object,embedded_object.Embedded_Particle_On_Segment(i,k),emb_triangle1)
                || !Emb_Node_In_Triangle(embedded_object,embedded_object.Embedded_Particle_On_Segment(i,l),emb_triangle1)
                || !Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(i,k))
                || !Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(k,j))
                || !Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(j,l))
                || !Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(l,i))
                || !Quad_Diagonal_Is_Correct(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(i,k),embedded_object.Embedded_Particle_On_Segment(j,l)))
                    permute_four(nodes,++count).Get(i,j,k,l);
            embedded_object.orientation_index(tetrahedron)=count;
            embedded_object.Add_Embedded_Triangle(i,k,k,j,j,l);embedded_object.Add_Embedded_Triangle(j,l,l,i,i,k);}}
    else if(Add_Best_Embedded_Triangle_With_Quad(fracture_normal,tetrahedron,tetrahedron_phi)){
        int count=1,i,j,k,l;nodes.Get(i,j,k,l);
        while(!Emb_Node_In_Quad(embedded_object,embedded_object.Embedded_Particle_On_Segment(i,k),emb_triangle1,emb_triangle2)
              || !Emb_Node_In_Quad(embedded_object,embedded_object.Embedded_Particle_On_Segment(k,j),emb_triangle1,emb_triangle2)
              || !Emb_Node_In_Quad(embedded_object,embedded_object.Embedded_Particle_On_Segment(j,l),emb_triangle1,emb_triangle2)
              || !Emb_Node_In_Quad(embedded_object,embedded_object.Embedded_Particle_On_Segment(l,i),emb_triangle1,emb_triangle2)
              || !embedded_object.Embedded_Particle_On_Segment(i,j) || !embedded_object.Embedded_Particle_On_Segment(i,k)
              || !embedded_object.Embedded_Particle_On_Segment(i,l)
              || !embedded_object.embedded_object.mesh.Triangle(embedded_object.Particle_Embedded_On_Segment(i,j),
                      embedded_object.Particle_Embedded_On_Segment(i,k),embedded_object.Particle_Embedded_On_Segment(i,l)))
            permute_four(nodes,++count).Get(i,j,k,l);
        embedded_object.orientation_index(tetrahedron)=count;}
}
//#####################################################################
// Function Emb_Node_In_Triangle
//#####################################################################
template<class T> bool FRACTURE_TETRAHEDRALIZED_VOLUME<T>::
Emb_Node_In_Triangle(const EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_tetrahedralized_volume,const int emb_node,const int emb_tri)
{
    if(emb_node==0) return false;
    else return embedded_tetrahedralized_volume.embedded_object.mesh.Node_In_Triangle(embedded_tetrahedralized_volume.embedded_particles.active_indices(emb_node),emb_tri);
}
//####################################################################
// Function Add_Best_Embedded_Triangle_With_Quad
//###################################################################
template<class T> bool FRACTURE_TETRAHEDRALIZED_VOLUME<T>::
Add_Best_Embedded_Triangle_With_Quad(const TV& fracture_normal,const int tetrahedron,const VECTOR<T,4>* tetrahedron_phi)
{
    int i,j,k,l;permute_four(embedded_object.simplicial_object.mesh.elements(tetrahedron),embedded_object.orientation_index(tetrahedron)).Get(i,j,k,l);
    int pi,pj,pk,pl;permute_four(VECTOR<int,4>(1,2,3,4),embedded_object.orientation_index(tetrahedron)).Get(pi,pj,pk,pl);
    T best_matched=(T)-FLT_MAX,best_lambda=0;int best_matched_index=0;
    TV ik=embedded_object.Position_Of_Embedded_Particle(i,k),kj=embedded_object.Position_Of_Embedded_Particle(k,j),
        jl=embedded_object.Position_Of_Embedded_Particle(j,l),li=embedded_object.Position_Of_Embedded_Particle(l,i),
        &xi=embedded_object.simplicial_object.particles.X(i),&xj=embedded_object.simplicial_object.particles.X(j),
        &xk=embedded_object.simplicial_object.particles.X(k),&xl=embedded_object.simplicial_object.particles.X(l);
    if(embedded_object.Node_In_Simplex_Is_Material(i,tetrahedron)){assert(embedded_object.Node_In_Simplex_Is_Material(j,tetrahedron));
        if(!tetrahedron_phi || (*tetrahedron_phi)(pi)>0){
            T lambda=Interpolation_Fraction_For_Best_Normal(fracture_normal,tetrahedron,ik,li,i,j),
                candidate=abs(TV::Dot_Product(fracture_normal,TRIANGLE_3D<T>::Normal(ik,li,LINEAR_INTERPOLATION<T,TV>::Linear(xi,xj,lambda))));
            if(candidate>best_matched){best_matched=candidate;best_matched_index=i;best_lambda=lambda;}}
        if(!tetrahedron_phi || (*tetrahedron_phi)(pj)>0){
            T lambda=Interpolation_Fraction_For_Best_Normal(fracture_normal,tetrahedron,jl,kj,j,i),
                candidate=abs(TV::Dot_Product(fracture_normal,TRIANGLE_3D<T>::Normal(jl,kj,LINEAR_INTERPOLATION<T,TV>::Linear(xj,xi,lambda))));
            if(candidate>best_matched){best_matched=candidate;best_matched_index=j;best_lambda=lambda;}}}
    else{assert(embedded_object.Node_In_Simplex_Is_Material(k,tetrahedron) && embedded_object.Node_In_Simplex_Is_Material(k,tetrahedron));
        if(!tetrahedron_phi || (*tetrahedron_phi)(pk)>0){
            T lambda=Interpolation_Fraction_For_Best_Normal(fracture_normal,tetrahedron,ik,kj,k,l),
                candidate=abs(TV::Dot_Product(fracture_normal,TRIANGLE_3D<T>::Normal(ik,kj,LINEAR_INTERPOLATION<T,TV>::Linear(xk,xl,lambda))));
            if(candidate>best_matched){best_matched=candidate;best_matched_index=k;best_lambda=lambda;}}
        if(!tetrahedron_phi || (*tetrahedron_phi)(pl)>0){
            T lambda=Interpolation_Fraction_For_Best_Normal(fracture_normal,tetrahedron,jl,li,l,k),
                candidate=abs(TV::Dot_Product(fracture_normal,TRIANGLE_3D<T>::Normal(jl,li,LINEAR_INTERPOLATION<T,TV>::Linear(xl,xk,lambda))));
            if(candidate>best_matched){best_matched=candidate;best_matched_index=l;best_lambda=lambda;}}}
    if(best_matched_index==i){
        embedded_object.Add_Embedded_Particle_If_Not_Already_There(i,j,best_lambda);embedded_object.Add_Embedded_Triangle(i,j,i,k,i,l);return true;}
    else if(best_matched_index==j){
        embedded_object.Add_Embedded_Particle_If_Not_Already_There(j,i,best_lambda);embedded_object.Add_Embedded_Triangle(j,i,j,l,j,k);return true;}
    else if(best_matched_index==k){
        embedded_object.Add_Embedded_Particle_If_Not_Already_There(k,l,best_lambda);embedded_object.Add_Embedded_Triangle(k,l,k,i,k,j);return true;}
    else if(best_matched_index==l){
        embedded_object.Add_Embedded_Particle_If_Not_Already_There(l,k,best_lambda);embedded_object.Add_Embedded_Triangle(l,k,l,j,l,i);return true;}
    else return false;
}
//#####################################################################
// Function Interpolation_Fraction_For_Best_Normal
//#####################################################################
template<class T> T FRACTURE_TETRAHEDRALIZED_VOLUME<T>::
Interpolation_Fraction_For_Best_Normal(const TV& fracture_normal,const int tetrahedron,const TV& ik,const TV& il,const int i,const int j)
{
    if(int surface_particle_ij=embedded_object.Embedded_Particle_On_Segment(i,j)){
        if(embedded_object.parent_particles(surface_particle_ij)(1)==i) return embedded_object.interpolation_fraction(surface_particle_ij);
        else return (T)1-embedded_object.interpolation_fraction(surface_particle_ij);}
    TV &xi=embedded_object.simplicial_object.particles.X(i),&xj=embedded_object.simplicial_object.particles.X(j);
    TV d=(ik-il).Normalized();
    TV projected_normal=fracture_normal-TV::Dot_Product(fracture_normal,d)*d;
    if(projected_normal.Magnitude_Squared()<(T)1e-8) return (T).5; // arbitrary since orthogonal -- and the dot product with fracture normal will be zero
    T denominator=TV::Dot_Product(xj-xi,projected_normal);
    if(abs(denominator)<(T)1e-8){
        TV ni=TRIANGLE_3D<T>::Normal(xi,ik,il),nj=TRIANGLE_3D<T>::Normal(xj,ik,il);
        T ni_dot_projected_normal=abs(TV::Dot_Product(ni,projected_normal));T nj_dot_projected_normal=abs(TV::Dot_Product(nj,projected_normal));
        if(ni_dot_projected_normal>nj_dot_projected_normal) return embedded_object.Clamp_Interpolation_Fraction((T)0);
        else return embedded_object.Clamp_Interpolation_Fraction((T)1);}
    return embedded_object.Clamp_Interpolation_Fraction(TV::Dot_Product(ik-xi,projected_normal)/denominator);
}
//#####################################################################
// Function Emb_Node_In_Quad
//#####################################################################
template<class T> bool FRACTURE_TETRAHEDRALIZED_VOLUME<T>::
Emb_Node_In_Quad(const EMBEDDED_TETRAHEDRALIZED_VOLUME<T>& embedded_tetrahedralized_volume,const int emb_node,const int emb_tri1,const int emb_tri2)
{
    if(emb_node==0) return false;
    int global_emb_node=embedded_tetrahedralized_volume.embedded_particles.active_indices(emb_node);
    return embedded_tetrahedralized_volume.embedded_object.mesh.Node_In_Triangle(global_emb_node,emb_tri1) ||
        embedded_tetrahedralized_volume.embedded_object.mesh.Node_In_Triangle(global_emb_node,emb_tri2);
}
//#####################################################################
// Function Add_Third_Cut
//#####################################################################
template<class T> void FRACTURE_TETRAHEDRALIZED_VOLUME<T>::
Add_Third_Cut(const int tetrahedron,const TV& fracture_normal,const VECTOR<T,4>* tetrahedron_phi)
{
    VECTOR<int,4> nodes=embedded_object.simplicial_object.mesh.elements(tetrahedron);
    VECTOR<int,4> emb_triangles=embedded_object.Embedded_Subelements_In_Element(tetrahedron);
    int emb_triangle1=emb_triangles[1],emb_triangle3=emb_triangles[3];assert(emb_triangles[4]==0);
    HYPOTHETICAL_CUT_TETRAHEDRONS<T> hypothetical_cut(embedded_object);
    if(emb_triangle3){
        assert(embedded_object.Cut_By_Quad(tetrahedron) 
            && (embedded_object.Node_Separated_By_Embedded_Subelement(emb_triangle1) 
                || embedded_object.Node_Separated_By_Embedded_Subelement(emb_triangles[2])
                || embedded_object.Node_Separated_By_Embedded_Subelement(emb_triangle3)));
        Add_Best_Embedded_Triangle_With_Quad_And_Triangle(fracture_normal,tetrahedron,tetrahedron_phi);}
    else if(int segments_intersected=Add_Best_Embedded_Triangle_Or_Quad_With_Two_Triangles(fracture_normal,tetrahedron,hypothetical_cut)){
        if(segments_intersected==3){// we know we have three elements
            assert(!(hypothetical_cut.Triangle_Cut_Already_In_Embedded_Tetrahedralized_Volume()));
            int i,j,k,l;permute_four(nodes,embedded_object.orientation_index(tetrahedron)).Get(i,j,k,l);
            if(Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(i,l))) 
                embedded_object.Add_Embedded_Triangle(l,k,l,j,l,i);
            else{
                assert(Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(i,k)));
                embedded_object.Add_Embedded_Triangle(k,l,k,i,k,j);}}
            // the two elements mandated by the bdry code are already in the right place (either of the two third elements work)
        else if(segments_intersected==4){ // we know we have two triangles and a quad (we just added the quad)
            assert(!(hypothetical_cut.Quad_Cut_Already_In_Embedded_Tetrahedralized_Volume()));
            int count=1,i,j,k,l;nodes.Get(i,j,k,l);
            while(!Emb_Node_In_Triangle(embedded_object,embedded_object.Embedded_Particle_On_Segment(i,j),emb_triangle1)
                || !Emb_Node_In_Triangle(embedded_object,embedded_object.Embedded_Particle_On_Segment(i,k),emb_triangle1)
                || !Emb_Node_In_Triangle(embedded_object,embedded_object.Embedded_Particle_On_Segment(i,l),emb_triangle1)
                || !Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(i,k))
                || !Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(k,j))
                || !Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(j,l))
                || !Emb_Node_Added_By_This_Tet(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(l,i))
                || !Quad_Diagonal_Is_Correct(hypothetical_cut,embedded_object.Embedded_Particle_On_Segment(i,k),embedded_object.Embedded_Particle_On_Segment(j,l)))
                permute_four(nodes,++count).Get(i,j,k,l);
            embedded_object.orientation_index(tetrahedron)=count;
            embedded_object.Add_Embedded_Triangle(i,k,k,j,j,l);embedded_object.Add_Embedded_Triangle(j,l,l,i,i,k);}}
}
//####################################################################
// Function Add_Best_Embedded_Triangle_With_Quad_And_Triangle
//####################################################################
template<class T> void FRACTURE_TETRAHEDRALIZED_VOLUME<T>::
Add_Best_Embedded_Triangle_With_Quad_And_Triangle(const TV& fracture_normal,const int tetrahedron,const VECTOR<T,4>* tetrahedron_phi)
{
    int i,j,k,l;permute_four(embedded_object.simplicial_object.mesh.elements(tetrahedron),embedded_object.orientation_index(tetrahedron)).Get(i,j,k,l);
    if(!embedded_object.Node_In_Simplex_Is_Material(k,tetrahedron)) return; // only break if it's material (note: k material --> l is material)
    assert(embedded_object.Node_In_Simplex_Is_Material(l,tetrahedron));
    if(!tetrahedron_phi && Fracture_Phi_Index(tetrahedron) && Phi_In_Simplex(k,tetrahedron)<=0) return;
    TV ik=embedded_object.Position_Of_Embedded_Particle(i,k),kj=embedded_object.Position_Of_Embedded_Particle(k,j),
        jl=embedded_object.Position_Of_Embedded_Particle(j,l),li=embedded_object.Position_Of_Embedded_Particle(l,i),
        &xk=embedded_object.simplicial_object.particles.X(k),&xl=embedded_object.simplicial_object.particles.X(l);
    T lambda_k=Interpolation_Fraction_For_Best_Normal(fracture_normal,tetrahedron,ik,kj,k,l);
    T lambda_l=Interpolation_Fraction_For_Best_Normal(fracture_normal,tetrahedron,jl,li,l,k);
    T quality_metric_k=abs(TV::Dot_Product(fracture_normal,TRIANGLE_3D<T>::Normal(ik,kj,LINEAR_INTERPOLATION<T,TV>::Linear(xk,xl,lambda_k))));
    T quality_metric_l=abs(TV::Dot_Product(fracture_normal,TRIANGLE_3D<T>::Normal(jl,li,LINEAR_INTERPOLATION<T,TV>::Linear(xl,xk,lambda_l))));
    if(quality_metric_k>quality_metric_l){
        embedded_object.Add_Embedded_Particle_If_Not_Already_There(k,l,lambda_k);embedded_object.Add_Embedded_Triangle(k,l,k,i,k,j);}
    else{embedded_object.Add_Embedded_Particle_If_Not_Already_There(l,k,lambda_l);embedded_object.Add_Embedded_Triangle(l,k,l,j,l,i);}
}
//####################################################################
// Function Add_Best_Embedded_Triangle_Or_Quad_With_Two_Triangles
//###################################################################
template<class T> int FRACTURE_TETRAHEDRALIZED_VOLUME<T>::
Add_Best_Embedded_Triangle_Or_Quad_With_Two_Triangles(const TV& fracture_normal,const int tetrahedron,HYPOTHETICAL_CUT_TETRAHEDRONS<T>& hypothetical_cut)
{
    T best_cut_quality=(T)-FLT_MAX;
    for(int a=1;a<=7;a++){
        HYPOTHETICAL_CUT_TETRAHEDRONS<T> hc(embedded_object);
        hc.Initialize_Hypothetical_Cut(fracture_normal,a,tetrahedron);
        if(!hc.hypothetical_nodes.m)continue;
        if(!hc.Edges_Shared_With_Existing_Embedded_Surface())continue;
        if(hc.Cut_Already_Exists())continue;
        if(hc.Would_Orphan_Half_Oct())continue;
        T cut_quality=hc.Quality_Of_Cut();
        if(cut_quality>best_cut_quality){best_cut_quality=cut_quality;hypothetical_cut=hc;}}
    hypothetical_cut.Add_Hypothetical_Nodes_To_Embedded_Object(embedded_object);
    return hypothetical_cut.hypothetical_nodes.m;
}
//#####################################################################
// Function Add_Cut
//#####################################################################
template<class T> void FRACTURE_TETRAHEDRALIZED_VOLUME<T>::
Add_Cut(const int tetrahedron,const TV& fracture_normal)
{
    embedded_object.Initialize_Orientation_Index_If_Necessary();

    int cut=embedded_object.Number_Of_Embedded_Cuts(tetrahedron)+1;
    if(cut==1) Add_First_Cut(tetrahedron,fracture_normal);
    else if(cut==2) Add_Second_Cut(tetrahedron,fracture_normal);
    else if(cut==3) Add_Third_Cut(tetrahedron,fracture_normal);
    else PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Cut
//#####################################################################
template<class T> void FRACTURE_TETRAHEDRALIZED_VOLUME<T>::
Add_Cut_Based_On_Phi(const int tetrahedron,const VECTOR<T,4>& tetrahedron_phi)
{
    embedded_object.Initialize_Orientation_Index_If_Necessary();

    int cut=embedded_object.Number_Of_Embedded_Cuts(tetrahedron)+1;
    if(cut==1) Add_First_Cut_Based_On_Phi(tetrahedron,tetrahedron_phi);
    else if(cut==2) Add_Next_Cut_Based_On_Phi(tetrahedron,tetrahedron_phi,true);
    else if(cut==3) Add_Next_Cut_Based_On_Phi(tetrahedron,tetrahedron_phi,false);
    //else PHYSBAM_FATAL_ERROR();
}
//#####################################################################
template class FRACTURE_TETRAHEDRALIZED_VOLUME<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FRACTURE_TETRAHEDRALIZED_VOLUME<double>;
#endif
