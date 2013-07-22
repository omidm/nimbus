//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Neil Molino, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HYPOTHETICAL_CUT_TETRAHEDRONS
//##################################################################### 
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/HYPOTHETICAL_CUT_TETRAHEDRONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/HYPOTHETICAL_NODE.h>
using namespace PhysBAM;
//#####################################################################
// Function operator=
//##################################################################### 
template<class T> HYPOTHETICAL_CUT_TETRAHEDRONS<T>& HYPOTHETICAL_CUT_TETRAHEDRONS<T>::
operator=(const HYPOTHETICAL_CUT_TETRAHEDRONS& old_cut)
{
    hypothetical_nodes=old_cut.hypothetical_nodes;
    cut_quality_metric=old_cut.cut_quality_metric;
    fracture_normal=old_cut.fracture_normal;cut_index=old_cut.cut_index;
    tetrahedron=old_cut.tetrahedron;
    quad_diagonal_indices=old_cut.quad_diagonal_indices;
    return *this;
}
//#####################################################################
// Function Initialize_Hypothetical_Cut
//##################################################################### 
template<class T> void HYPOTHETICAL_CUT_TETRAHEDRONS<T>::
Initialize_Hypothetical_Cut(const VECTOR<T,3>& fracture_normal_input,const int cut_index_input,const int tetrahedron_input)
{
    assert(1<=cut_index_input && cut_index_input <= 7);
    cut_index=cut_index_input;
    fracture_normal=fracture_normal_input;
    tetrahedron=tetrahedron_input;
    int i,j,k,l;embedded_object.simplicial_object.mesh.elements(tetrahedron).Get(i,j,k,l);
    if(cut_index==1) Initialize_Triangle_Cut(i);   // i
    else if(cut_index==2) Initialize_Triangle_Cut(j);  // j
    else if(cut_index==3) Initialize_Triangle_Cut(k);  // k
    else if(cut_index==4) Initialize_Triangle_Cut(l);  // l
    else if(cut_index==5) Initialize_Quad_Cut(i,j); // ij-kl
    else if(cut_index==6) Initialize_Quad_Cut(i,k); // ik-jl
    else{assert(cut_index==7);Initialize_Quad_Cut(i,l);}// il-jk        
    Initialize_Quality_Metric_And_Quad_Diagonal_Indices();
}
//#####################################################################
// Function Initialize_Triangle_Cut
//##################################################################### 
template<class T> void HYPOTHETICAL_CUT_TETRAHEDRONS<T>::
Initialize_Triangle_Cut(const int isolated_node)
{
    int on1,on2,on3;embedded_object.simplicial_object.mesh.Other_Three_Nodes(isolated_node,tetrahedron,on1,on2,on3);
    int ep1=embedded_object.Embedded_Particle_On_Segment(isolated_node,on1),
        ep2=embedded_object.Embedded_Particle_On_Segment(isolated_node,on2),
        ep3=embedded_object.Embedded_Particle_On_Segment(isolated_node,on3);
    VECTOR<T,3> x1,x2,x3,x4;
    x1=embedded_object.simplicial_object.particles.X(isolated_node);
    x2=embedded_object.simplicial_object.particles.X(on1);
    x3=embedded_object.simplicial_object.particles.X(on2);
    x4=embedded_object.simplicial_object.particles.X(on3);
    assert(TETRAHEDRON<T>::Signed_Volume(x1,x2,x3,x4)>0);
    if(ep1&&ep2&&ep3){
        Add_Hypothetical_Node(isolated_node,on1,embedded_object.interpolation_fraction(ep1));
        Add_Hypothetical_Node(isolated_node,on2,embedded_object.interpolation_fraction(ep2));
        Add_Hypothetical_Node(isolated_node,on3,embedded_object.interpolation_fraction(ep3));}
    else if(ep1&&ep2&&!ep3){
        Add_Hypothetical_Node(isolated_node,on1,embedded_object.interpolation_fraction(ep1));
        Add_Hypothetical_Node(isolated_node,on2,embedded_object.interpolation_fraction(ep2));
        T interpolation_fraction3=Interpolation_Fraction_For_Best_Normal(Position(1),Position(2),isolated_node,on3);
        Add_Hypothetical_Node(isolated_node,on3,interpolation_fraction3);}
    else if(ep1&&!ep2&&ep3){
        Add_Hypothetical_Node(isolated_node,on1,embedded_object.interpolation_fraction(ep1));
        Add_Hypothetical_Node(isolated_node,on3,embedded_object.interpolation_fraction(ep3));
        T interpolation_fraction2=Interpolation_Fraction_For_Best_Normal(Position(1),Position(2),isolated_node,on2);
        Add_Hypothetical_Node(isolated_node,on2,interpolation_fraction2);}
    else if(!ep1&&ep2&&ep3){
        Add_Hypothetical_Node(isolated_node,on2,embedded_object.interpolation_fraction(ep2));
        Add_Hypothetical_Node(isolated_node,on3,embedded_object.interpolation_fraction(ep3));
        T interpolation_fraction1=Interpolation_Fraction_For_Best_Normal(Position(1),Position(2),isolated_node,on1);
        Add_Hypothetical_Node(isolated_node,on1,interpolation_fraction1);}
}
//#####################################################################
// Function Initialize_Quad_Cut
//##################################################################### 
template<class T> void HYPOTHETICAL_CUT_TETRAHEDRONS<T>::
Initialize_Quad_Cut(const int i,const int j)
{
    int k,l;embedded_object.simplicial_object.mesh.Other_Two_Nodes(i,j,tetrahedron,k,l);
    int ep1=embedded_object.Embedded_Particle_On_Segment(i,l);
    int ep2=embedded_object.Embedded_Particle_On_Segment(l,j);
    int ep3=embedded_object.Embedded_Particle_On_Segment(j,k);
    int ep4=embedded_object.Embedded_Particle_On_Segment(k,i);
    PARTICLES<TV>& particles=dynamic_cast<PARTICLES<TV>&>(embedded_object.simplicial_object.particles);
    assert(TETRAHEDRON<T>::Signed_Volume(particles.X(i),particles.X(j),particles.X(k),particles.X(l))>0);
    if(ep1&&ep2&&!ep3&&!ep4){
        Add_Hypothetical_Node(i,l,embedded_object.interpolation_fraction(ep1));
        Add_Hypothetical_Node(l,j,embedded_object.interpolation_fraction(ep2));
        T interpolation_fraction_kj,interpolation_fraction_ki;
        Interpolation_Fractions_For_Best_Normal(Position(1),Position(2),k,j,i,interpolation_fraction_kj,interpolation_fraction_ki);
        Add_Hypothetical_Node(k,j,interpolation_fraction_kj);
        Add_Hypothetical_Node(k,i,interpolation_fraction_ki);} 
    if(ep1&&!ep2&&ep3&&!ep4){} // do nothing -- no way it could be edge connected  
    if(ep1&&!ep2&&!ep3&&ep4){
        Add_Hypothetical_Node(i,l,embedded_object.interpolation_fraction(ep1));
        Add_Hypothetical_Node(k,i,embedded_object.interpolation_fraction(ep4));
        T interpolation_fraction_jk,interpolation_fraction_jl;
        Interpolation_Fractions_For_Best_Normal(Position(1),Position(2),j,k,l,interpolation_fraction_jk,interpolation_fraction_jl);
        Add_Hypothetical_Node(j,k,interpolation_fraction_jk);
        Add_Hypothetical_Node(j,l,interpolation_fraction_jl);}
    if(!ep1&&!ep2&&ep3&&ep4){
        Add_Hypothetical_Node(j,k,embedded_object.interpolation_fraction(ep3));
        Add_Hypothetical_Node(k,i,embedded_object.interpolation_fraction(ep4));
        T interpolation_fraction_li,interpolation_fraction_lj;
        Interpolation_Fractions_For_Best_Normal(Position(1),Position(2),l,i,j,interpolation_fraction_li,interpolation_fraction_lj);
        Add_Hypothetical_Node(l,i,interpolation_fraction_li);
        Add_Hypothetical_Node(l,j,interpolation_fraction_lj);}
    if(!ep1&&ep2&&!ep3&&ep4){} // do nothing -- no way it could be edge connected    
    if(!ep1&&ep2&&ep3&&!ep4){
        Add_Hypothetical_Node(l,j,embedded_object.interpolation_fraction(ep2));
        Add_Hypothetical_Node(k,j,embedded_object.interpolation_fraction(ep3));
        T interpolation_fraction_il,interpolation_fraction_ik;
        Interpolation_Fractions_For_Best_Normal(Position(1),Position(2),i,l,k,interpolation_fraction_il,interpolation_fraction_ik);
        Add_Hypothetical_Node(i,k,interpolation_fraction_ik);
        Add_Hypothetical_Node(i,l,interpolation_fraction_il);}
    if(!ep1&&ep2&&ep3&&ep4){    
        Add_Hypothetical_Node(l,j,embedded_object.interpolation_fraction(ep2));
        Add_Hypothetical_Node(k,j,embedded_object.interpolation_fraction(ep3));
        Add_Hypothetical_Node(k,i,embedded_object.interpolation_fraction(ep4));
        VECTOR<T,3> normal=TRIANGLE_3D<T>::Normal(Position(1),Position(2),Position(3));
        PLANE<T> plane(normal,Position(1));
        T interpolation_fraction_il;
        plane.Segment_Plane_Intersection(particles.X(i),particles.X(l),interpolation_fraction_il);
        Add_Hypothetical_Node(i,l,interpolation_fraction_il);}
    if(ep1&&!ep2&&ep3&&ep4){
        Add_Hypothetical_Node(k,j,embedded_object.interpolation_fraction(ep3));
        Add_Hypothetical_Node(k,i,embedded_object.interpolation_fraction(ep4));
        Add_Hypothetical_Node(i,l,embedded_object.interpolation_fraction(ep1));
        VECTOR<T,3> normal=TRIANGLE_3D<T>::Normal(Position(1),Position(2),Position(3));
        PLANE<T> plane(normal,Position(1));
        T interpolation_fraction_jl;
        plane.Segment_Plane_Intersection(particles.X(j),particles.X(l),interpolation_fraction_jl);
        Add_Hypothetical_Node(j,l,interpolation_fraction_jl);}
    if(ep1&&ep2&&!ep3&&ep4){
        Add_Hypothetical_Node(k,i,embedded_object.interpolation_fraction(ep4));
        Add_Hypothetical_Node(i,l,embedded_object.interpolation_fraction(ep1));
        Add_Hypothetical_Node(l,j,embedded_object.interpolation_fraction(ep2));
        VECTOR<T,3> normal=TRIANGLE_3D<T>::Normal(Position(1),Position(2),Position(3));
        PLANE<T> plane(normal,Position(1));
        T interpolation_fraction_jk;
        plane.Segment_Plane_Intersection(particles.X(j),particles.X(k),interpolation_fraction_jk);
        Add_Hypothetical_Node(j,k,interpolation_fraction_jk);}
    if(ep1&&ep2&&ep3&&!ep4){
        Add_Hypothetical_Node(i,l,embedded_object.interpolation_fraction(ep1));
        Add_Hypothetical_Node(l,j,embedded_object.interpolation_fraction(ep2));
        Add_Hypothetical_Node(k,j,embedded_object.interpolation_fraction(ep3));
        VECTOR<T,3> normal=TRIANGLE_3D<T>::Normal(Position(1),Position(2),Position(3));
        PLANE<T> plane(normal,Position(1));
        T interpolation_fraction_ki;
        plane.Segment_Plane_Intersection(particles.X(k),particles.X(i),interpolation_fraction_ki);
        Add_Hypothetical_Node(k,i,interpolation_fraction_ki);}
    if(ep1&&ep2&&ep3&&ep4){
        Add_Hypothetical_Node(i,l,embedded_object.interpolation_fraction(ep1));
        Add_Hypothetical_Node(l,j,embedded_object.interpolation_fraction(ep2));
        Add_Hypothetical_Node(k,j,embedded_object.interpolation_fraction(ep3));
        Add_Hypothetical_Node(k,i,embedded_object.interpolation_fraction(ep4));}
}
//#####################################################################
// Function Initialize_Initiation_Point_Cut
//##################################################################### 
template<class T> void HYPOTHETICAL_CUT_TETRAHEDRONS<T>::
Initialize_Initiation_Point_Cut(const PLANE<T>& plane,const int tetrahedron_input)
{
    tetrahedron=tetrahedron_input;
    fracture_normal=plane.normal;cut_index=0;
    int i,j,k,l;embedded_object.simplicial_object.mesh.elements(tetrahedron).Get(i,j,k,l);
    VECTOR<T,3> xi=embedded_object.particles.X(i),xj=embedded_object.particles.X(j),xk=embedded_object.particles.X(k),xl=embedded_object.particles.X(l);
    T interpolation_fraction_ij,interpolation_fraction_ik,interpolation_fraction_il,interpolation_fraction_jk,interpolation_fraction_jl,interpolation_fraction_kl;
    bool cuts_ij=plane.Segment_Plane_Intersection(xi,xj,interpolation_fraction_ij);bool cuts_ik=plane.Segment_Plane_Intersection(xi,xk,interpolation_fraction_ik);
    bool cuts_il=plane.Segment_Plane_Intersection(xi,xl,interpolation_fraction_il);bool cuts_jk=plane.Segment_Plane_Intersection(xj,xk,interpolation_fraction_jk);
    bool cuts_jl=plane.Segment_Plane_Intersection(xj,xl,interpolation_fraction_jl);bool cuts_kl=plane.Segment_Plane_Intersection(xk,xl,interpolation_fraction_kl);
    if(cuts_ij) Add_Hypothetical_Node(i,j,interpolation_fraction_ij);
    if(cuts_ik) Add_Hypothetical_Node(i,k,interpolation_fraction_ik);
    if(cuts_il) Add_Hypothetical_Node(i,l,interpolation_fraction_il);
    if(cuts_jk) Add_Hypothetical_Node(j,k,interpolation_fraction_jk);
    if(cuts_jl) Add_Hypothetical_Node(j,l,interpolation_fraction_jl);
    if(cuts_kl) Add_Hypothetical_Node(k,l,interpolation_fraction_kl);
    if(!Valid_Cut(cuts_ij,cuts_ik,cuts_il,cuts_jk,cuts_jl,cuts_kl)){hypothetical_nodes.Remove_All();return;}
    Initialize_Quality_Metric_And_Quad_Diagonal_Indices_For_Initiation_Point();
}
//#####################################################################
// Function Valid_Cut
//##################################################################### 
template<class T> bool HYPOTHETICAL_CUT_TETRAHEDRONS<T>::
Valid_Cut(const bool cuts_ij,const bool cuts_ik,const bool cuts_il,const bool cuts_jk,const bool cuts_jl,const bool cuts_kl) const
{
    if(hypothetical_nodes.m==3){ 
        VECTOR<int,2> parents=hypothetical_nodes(1).parents;
        if(hypothetical_nodes(2).Is_Parent(parents[1]) && hypothetical_nodes(3).Is_Parent(parents[1])) return true;
        else if(hypothetical_nodes(2).Is_Parent(parents[2]) && hypothetical_nodes(3).Is_Parent(parents[2])) return true;
        else return false;}
    else if(hypothetical_nodes.m==4){// note in all cases quad diagonal is either 14 or 23
        if(!cuts_ij && !cuts_kl)return true; // 1 = ik, 2 = il, 3 = jk, 4 = jl
        else if(!cuts_ik && !cuts_jl)return true; // 1 = ij, 2 = il, 3 = jk, 4 = kl            
        else if(!cuts_il && !cuts_jk)return true; // 1 = ij, 2 = ik, 3 = jl, 4 = kl 
        else return false;}
    else return false;
}
//#####################################################################
// Function Embedded_Edge_Exists
//##################################################################### 
template<class T> bool HYPOTHETICAL_CUT_TETRAHEDRONS<T>::
Embedded_Edge_Exists(const int emb_node1,const int emb_node2)
{
    if(!emb_node1 || !emb_node2) return false;
    const ARRAY<int>& active_indices=embedded_object.embedded_particles.active_indices;
    return embedded_object.embedded_mesh.Triangle_On_Edge(active_indices(emb_node1),active_indices(emb_node2));
}
//#####################################################################
// Function Embedded_Edge_Exists
//##################################################################### 
template<class T> bool HYPOTHETICAL_CUT_TETRAHEDRONS<T>::
Embedded_Edge_Exists(const int ppa1,const int ppa2,const int ppb1,const int ppb2)
{
    return Embedded_Edge_Exists(embedded_object.Embedded_Particle_On_Segment(ppa1,ppa2),
        embedded_object.Embedded_Particle_On_Segment(ppb1,ppb2));
}
//#####################################################################
// Function Interpolation_Fraction_For_Best_Normal
//##################################################################### 
template<class T> T HYPOTHETICAL_CUT_TETRAHEDRONS<T>::
Interpolation_Fraction_For_Best_Normal(const VECTOR<T,3>& ik,const VECTOR<T,3>& il,const int i,const int j)
{
    assert(!embedded_object.Embedded_Particle_On_Segment(i,j));
    VECTOR<T,3> &xi=embedded_object.simplicial_object.particles.X(i),&xj=embedded_object.simplicial_object.particles.X(j);
    VECTOR<T,3> d=(ik-il).Normalized();
    VECTOR<T,3> projected_normal=fracture_normal-VECTOR<T,3>::Dot_Product(fracture_normal,d)*d;
    if(projected_normal.Magnitude_Squared()<(T)1e-8) return (T).5; // arbitrary since orthogonal -- and the dot product with fracture normal will be zero
    T denominator=VECTOR<T,3>::Dot_Product(xj-xi,projected_normal);
    if(abs(denominator)<(T)1e-8){
        VECTOR<T,3> ni=TRIANGLE_3D<T>::Normal(xi,ik,il),nj=TRIANGLE_3D<T>::Normal(xj,ik,il);
        T ni_dot_projected_normal=abs(VECTOR<T,3>::Dot_Product(ni,projected_normal));T nj_dot_projected_normal=abs(VECTOR<T,3>::Dot_Product(nj,projected_normal));
        if(ni_dot_projected_normal > nj_dot_projected_normal) return embedded_object.Clamp_Interpolation_Fraction((T)0);            
        else return embedded_object.Clamp_Interpolation_Fraction((T)1);}
    return embedded_object.Clamp_Interpolation_Fraction(VECTOR<T,3>::Dot_Product(ik-xi,projected_normal)/denominator);
}
//#####################################################################
// Function Interpolation_Fractions_For_Best_Normal
//##################################################################### 
template<class T> void HYPOTHETICAL_CUT_TETRAHEDRONS<T>::
Interpolation_Fractions_For_Best_Normal(const VECTOR<T,3>& ik,const VECTOR<T,3>& il,const int j,const int k,const int l,T& interpolation_fraction_jk,T& interpolation_fraction_jl)
{
    interpolation_fraction_jk=Interpolation_Fraction_For_Best_Normal(ik,il,j,k);
    interpolation_fraction_jl=Interpolation_Fraction_For_Best_Normal(ik,il,j,l);
}
//#####################################################################
// Function Initialize_Quality_Metric_And_Quad_Diagonal_Indices_For_Initiation_Point
//##################################################################### 
template<class T> void HYPOTHETICAL_CUT_TETRAHEDRONS<T>::
Initialize_Quality_Metric_And_Quad_Diagonal_Indices_For_Initiation_Point()
{
    if(hypothetical_nodes.m==3){
        VECTOR<T,3> position1=Position(1);
        VECTOR<T,3> hypothetical_cut_normal=VECTOR<T,3>::Cross_Product(Position(2)-position1,Position(3)-position1).Normalized();
        cut_quality_metric=abs(VECTOR<T,3>::Dot_Product(fracture_normal,hypothetical_cut_normal));}
    else if(hypothetical_nodes.m==4){ 
        VECTOR<T,3> position1=Position(1),position2=Position(2),position3=Position(3),position4=Position(4);
        VECTOR<T,3> n123=VECTOR<T,3>::Cross_Product(position2-position1,position3-position1).Normalized(),
                     n234=VECTOR<T,3>::Cross_Product(position2-position4,position3-position4).Normalized(),
                     n124=VECTOR<T,3>::Cross_Product(position2-position1,position4-position1).Normalized(),
                     n143=VECTOR<T,3>::Cross_Product(position3-position1,position4-position1).Normalized();
        T metric23=abs(VECTOR<T,3>::Dot_Product(n123,fracture_normal))+abs(VECTOR<T,3>::Dot_Product(n234,fracture_normal));
        T metric14=abs(VECTOR<T,3>::Dot_Product(n124,fracture_normal))+abs(VECTOR<T,3>::Dot_Product(n143,fracture_normal));
        if(metric23 > metric14){quad_diagonal_indices.x=2;quad_diagonal_indices.y=3;cut_quality_metric=(T).5*metric23;}
        else{quad_diagonal_indices.x=1;quad_diagonal_indices.y=4;cut_quality_metric=(T).5*metric14;}}
}
//#####################################################################
// Function Initialize_Quality_Metric_And_Quad_Diagonal_Indices
//##################################################################### 
template<class T> void HYPOTHETICAL_CUT_TETRAHEDRONS<T>::
Initialize_Quality_Metric_And_Quad_Diagonal_Indices()
{
    if(hypothetical_nodes.m==3){
        VECTOR<T,3> position1=Position(1);
        VECTOR<T,3> hypothetical_cut_normal=VECTOR<T,3>::Cross_Product(Position(2)-position1,Position(3)-position1).Normalized();
        cut_quality_metric=abs(VECTOR<T,3>::Dot_Product(fracture_normal,hypothetical_cut_normal));}
    else if(hypothetical_nodes.m==4){ 
        VECTOR<T,3> position1=Position(1),position2=Position(2),position3=Position(3),position4=Position(4);
        VECTOR<T,3> n123=VECTOR<T,3>::Cross_Product(position2-position1,position3-position1).Normalized(),
                     n134=VECTOR<T,3>::Cross_Product(position3-position1,position4-position1).Normalized(),
                     n234=VECTOR<T,3>::Cross_Product(position2-position4,position3-position4).Normalized(),
                     n124=VECTOR<T,3>::Cross_Product(position2-position1,position4-position1).Normalized();
        T metric13=abs(VECTOR<T,3>::Dot_Product(n123,fracture_normal))+abs(VECTOR<T,3>::Dot_Product(n134,fracture_normal));
        T metric24=abs(VECTOR<T,3>::Dot_Product(n124,fracture_normal))+abs(VECTOR<T,3>::Dot_Product(n234,fracture_normal));
        if(metric13 > metric24){quad_diagonal_indices.x=1;quad_diagonal_indices.y=3;cut_quality_metric=(T).5*metric13;}
        else{quad_diagonal_indices.x=2;quad_diagonal_indices.y=4;cut_quality_metric=(T).5*metric24;}}
}
//#####################################################################
// Function Cut_Already_Exists
//##################################################################### 
template<class T> bool HYPOTHETICAL_CUT_TETRAHEDRONS<T>::
Cut_Already_Exists()
{
    if(hypothetical_nodes.m==3) return Triangle_Cut_Already_In_Embedded_Tetrahedralized_Volume();
    else if(hypothetical_nodes.m==4) return Quad_Cut_Already_In_Embedded_Tetrahedralized_Volume();
    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Triangle_Cut_Already_In_Embedded_Tetrahedralized_Volume
//##################################################################### 
template<class T> bool HYPOTHETICAL_CUT_TETRAHEDRONS<T>::
Triangle_Cut_Already_In_Embedded_Tetrahedralized_Volume()
{
    assert(hypothetical_nodes.m==3);
    int emb_node1=embedded_object.Embedded_Particle_On_Segment(hypothetical_nodes(1).parents);
    int emb_node2=embedded_object.Embedded_Particle_On_Segment(hypothetical_nodes(2).parents);
    int emb_node3=embedded_object.Embedded_Particle_On_Segment(hypothetical_nodes(3).parents);
    if(!emb_node1 || !emb_node2 || !emb_node3) return false;
    const ARRAY<int>& active_indices=embedded_object.embedded_particles.active_indices;
    return embedded_object.embedded_mesh.Triangle(active_indices(emb_node1),active_indices(emb_node2),active_indices(emb_node3))!=0;
}
//#####################################################################
// Function Quad_Cut_Already_In_Embedded_Tetrahedralized_Volume
//##################################################################### 
template<class T> bool HYPOTHETICAL_CUT_TETRAHEDRONS<T>::
Quad_Cut_Already_In_Embedded_Tetrahedralized_Volume()
{
    assert(hypothetical_nodes.m==4);
    int emb_node1=embedded_object.Embedded_Particle_On_Segment(hypothetical_nodes(1).parents);
    int emb_node2=embedded_object.Embedded_Particle_On_Segment(hypothetical_nodes(2).parents);
    int emb_node3=embedded_object.Embedded_Particle_On_Segment(hypothetical_nodes(3).parents);  
    int emb_node4=embedded_object.Embedded_Particle_On_Segment(hypothetical_nodes(4).parents);
    if(!emb_node1 || !emb_node2 || !emb_node3 || !emb_node4) return false;
    const ARRAY<int>& active_indices=embedded_object.embedded_particles.active_indices;
    int global_node1=active_indices(emb_node1),global_node2=active_indices(emb_node2),global_node3=active_indices(emb_node3),global_node4=active_indices(emb_node4);
    return (embedded_object.embedded_mesh.Triangle(global_node1,global_node2,global_node3) && embedded_object.embedded_mesh.Triangle(global_node1,global_node3,global_node4))
        || (embedded_object.embedded_mesh.Triangle(global_node1,global_node2,global_node4) && embedded_object.embedded_mesh.Triangle(global_node2,global_node3,global_node4));
}
//#####################################################################
// Function Quality_Of_Cut
//##################################################################### 
template<class T> T HYPOTHETICAL_CUT_TETRAHEDRONS<T>::
Quality_Of_Cut(const T extra_quad_penalty_multiplier)const 
{
    if(hypothetical_nodes.m==3) return cut_quality_metric;
    return extra_quad_penalty_multiplier*cut_quality_metric;
}
//#####################################################################
// Function Number_Of_Edges_Shared_With_Existing_Embedded_Surface
//##################################################################### 
template<class T> bool HYPOTHETICAL_CUT_TETRAHEDRONS<T>::
Edges_Shared_With_Existing_Embedded_Surface()
{
    const TRIANGLE_MESH& embedded_mesh=embedded_object.embedded_mesh;
    const ARRAY<int>& active_indices=embedded_object.embedded_particles.active_indices;  
    if(hypothetical_nodes.m==3){
        int ri=hypothetical_nodes(1).index_in_embedded_particles,rj=hypothetical_nodes(2).index_in_embedded_particles,
            rk=hypothetical_nodes(3).index_in_embedded_particles;
        if(ri&&rj&&embedded_mesh.Triangle_On_Edge(active_indices(ri),active_indices(rj)))return true;
        if(ri&&rk&&embedded_mesh.Triangle_On_Edge(active_indices(ri),active_indices(rk)))return true;
        if(rj&&rk&&embedded_mesh.Triangle_On_Edge(active_indices(rj),active_indices(rk)))return true;}
    else{
        assert(hypothetical_nodes.m==4);
        int ri=hypothetical_nodes(1).index_in_embedded_particles,rj=hypothetical_nodes(2).index_in_embedded_particles,
            rk=hypothetical_nodes(3).index_in_embedded_particles,rl=hypothetical_nodes(4).index_in_embedded_particles;
        if(ri&&rj&&embedded_mesh.Triangle_On_Edge(active_indices(ri),active_indices(rj)))return true; 
        assert(!(ri&&rk&&embedded_mesh.Triangle_On_Edge(active_indices(ri),active_indices(rk)))); 
        if(ri&&rl&&embedded_mesh.Triangle_On_Edge(active_indices(ri),active_indices(rl)))return true;
        if(rj&&rk&&embedded_mesh.Triangle_On_Edge(active_indices(rj),active_indices(rk)))return true;
        assert(!(rj&&rl&&embedded_mesh.Triangle_On_Edge(active_indices(rj),active_indices(rl)))); 
        if(rk&&rl&&embedded_mesh.Triangle_On_Edge(active_indices(rk),active_indices(rl)))return true;}
    return false;
}
//#####################################################################
// Function Would_Orphan_Half_Oct
//##################################################################### 
template<class T> bool HYPOTHETICAL_CUT_TETRAHEDRONS<T>::
Would_Orphan_Half_Oct()
{
    if(cut_index < 5) return false;assert(cut_index<=7);
    assert(!embedded_object.Cut_By_Quad(tetrahedron));
    assert(embedded_object.Number_Of_Embedded_Cuts(tetrahedron)==2);
    VECTOR<int,4> emb_tris=embedded_object.Embedded_Subelements_In_Element(tetrahedron);
    assert(emb_tris[1]&&emb_tris[2]&&!emb_tris[3]&&!emb_tris[4]);
    int isolated1=embedded_object.Node_Separated_By_Embedded_Subelement(emb_tris[1]);
    int isolated2=embedded_object.Node_Separated_By_Embedded_Subelement(emb_tris[2]);
    int i,j,k,l;embedded_object.simplicial_object.mesh.elements(tetrahedron).Get(i,j,k,l);
    if(cut_index==5){
        if(i==isolated1 && j==isolated2) return true;
        if(j==isolated1 && i==isolated2) return true;
        if(k==isolated1 && l==isolated2) return true;
        if(l==isolated1 && k==isolated2) return true;}
    else if(cut_index==6){
        if(i==isolated1 && k==isolated2) return true;
        if(k==isolated1 && i==isolated2) return true;
        if(j==isolated1 && l==isolated2) return true;
        if(l==isolated1 && j==isolated2) return true;}
    else if(cut_index==7){
        if(i==isolated1 && l==isolated2) return true;
        if(l==isolated1 && i==isolated2) return true;
        if(j==isolated1 && k==isolated2) return true;
        if(k==isolated1 && j==isolated2) return true;}
    return false;
}
//##################################################################### 
template class HYPOTHETICAL_CUT_TETRAHEDRONS<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class HYPOTHETICAL_CUT_TETRAHEDRONS<double>;
#endif
