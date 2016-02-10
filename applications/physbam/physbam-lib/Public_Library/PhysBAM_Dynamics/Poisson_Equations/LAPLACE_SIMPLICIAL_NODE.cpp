//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_SIMPLICIAL_NODE
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/FACTORIAL.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SOLIDS_FORCES_POLICY.h>
#include <PhysBAM_Dynamics/Poisson_Equations/LAPLACE_SIMPLICIAL_NODE.h>
using namespace PhysBAM;
//#####################################################################
// Function Initialize_Object
//#####################################################################
template<class T> static inline T Determinant(const MATRIX<T,2>& A,const SYMMETRIC_MATRIX<T,2>& AA){return A.Determinant();}
template<class T> static inline T Determinant(const MATRIX<T,3>& A,const SYMMETRIC_MATRIX<T,3>& AA){return A.Determinant();}
template<class T> static inline T Determinant(const MATRIX<T,3,2>& A,const SYMMETRIC_MATRIX<T,2>& AA){return sqrt(AA.Determinant());}

template<class T_DIAGONAL,class T_OFFDIAGONAL,class T> static inline void Distribute_W(T_DIAGONAL diagonal,T_OFFDIAGONAL offdiagonal,const SYMMETRIC_MATRIX<T,2>& W)
{
    diagonal(2)+=W.x11;offdiagonal(2)+=W.x21;diagonal(3)+=W.x22;
    offdiagonal(1)-=W.x11+W.x21;offdiagonal(3)-=W.x21+W.x22;
    diagonal(1)+=W.x11+W.x21+W.x21+W.x22;
}
template<class T_DIAGONAL,class T_OFFDIAGONAL,class T> static inline void Distribute_W(T_DIAGONAL diagonal,T_OFFDIAGONAL offdiagonal,const SYMMETRIC_MATRIX<T,3>& W)
{
    diagonal(2)+=W.x11;offdiagonal(2)+=W.x21;offdiagonal(5)+=W.x31;
    diagonal(3)+=W.x22;offdiagonal(6)+=W.x32;diagonal(4)+=W.x33;
    offdiagonal(1)-=W.x11+W.x21+W.x31;offdiagonal(3)-=W.x21+W.x22+W.x32;offdiagonal(4)-=W.x31+W.x32+W.x33;
    diagonal(1)+=W.x11+W.x21+W.x31+W.x21+W.x22+W.x32+W.x31+W.x32+W.x33;
}
template<class TV,int d> void LAPLACE_SIMPLICIAL_NODE<TV,d>::
Initialize_Object()
{
    typedef typename MESH_POLICY<d>::MESH T_MESH;
    typedef MATRIX<T,TV::m,d> T_MATRIX;
    static const T volume_over_determinant=(T)1/FACTORIAL<d>::value;

    T_MESH& mesh=object.mesh;
    if(!mesh.segment_mesh) mesh.Initialize_Segment_Mesh();
    bool element_edges_defined=mesh.element_edges!=0;if(!element_edges_defined) mesh.Initialize_Element_Edges();
    psi_D.Resize(particles.array_collection->Size());
    diagonal.Resize(particles.array_collection->Size(),false,false);ARRAYS_COMPUTATIONS::Fill(diagonal,(T)0);
    INDIRECT_ARRAY<ARRAY<T> > diagonal_full(diagonal,particles.subset_index_from_point_cloud_index);
    offdiagonal.Resize(mesh.segment_mesh->elements.m,false);ARRAYS_COMPUTATIONS::Fill(offdiagonal,(T)0);
    for(int t=1;t<=mesh.elements.m;t++){
        VECTOR<int,d+1>& nodes=mesh.elements(t);
        VECTOR<int,d==2?3:6>& edges=(*mesh.element_edges)(t);
        T_MATRIX Dm=STRAIN_MEASURE<TV,d>::Ds(object.particles.X,nodes);
        SYMMETRIC_MATRIX<T,d> Dm_normal_equations=Dm.Normal_Equations_Matrix();
        SYMMETRIC_MATRIX<T,d> W=volume_over_determinant*Determinant(Dm,Dm_normal_equations)*Dm_normal_equations.Inverse();
        Distribute_W(diagonal_full.Subset(nodes),offdiagonal.Subset(edges),W);}
    if(!element_edges_defined){delete mesh.element_edges;mesh.element_edges=0;}
}
//#####################################################################
// Function Initialize_Dirichlet_Boundary
//#####################################################################
template<class T_OBJECT,class T_PSI_D_EDGES,class T> static void Initialize_Dirichlet_Boundary_Helper(T_OBJECT& object,T_PSI_D_EDGES& psi_D_edges,TRIANGULATED_SURFACE<T>& triangulated_surface,
    const T thickness_over_2)
{
    if(!triangulated_surface.triangle_list) triangulated_surface.Update_Triangle_List();
    if(!triangulated_surface.hierarchy) triangulated_surface.Initialize_Hierarchy();
    SEGMENT_MESH& segment_mesh=*object.mesh.segment_mesh;
    psi_D_edges.Resize(segment_mesh.elements.m);
    for(int s=1;s<=segment_mesh.elements.m;s++){
        VECTOR<int,2>& nodes=segment_mesh.elements(s);
        for(int i=1;i<=2;i++){
            RAY<VECTOR<T,3> > ray(SEGMENT_3D<T>(object.particles.X(nodes[i]),object.particles.X(nodes[3-i])));T magnitude=ray.t_max;
            if(triangulated_surface.hierarchy->Intersection(ray,psi_D_edges(s)(i).weights,thickness_over_2)){
                psi_D_edges(s)(i).theta=ray.t_max/magnitude;psi_D_edges(s)(i).simplex=ray.aggregate_id;}}
        if(!psi_D_edges(s)(1).simplex != !psi_D_edges(s)(2).simplex){ // make sure cuts are symmetric
            int cut=psi_D_edges(s)(1).simplex?1:2;    
            psi_D_edges(s)(3-cut)=psi_D_edges(s)(cut);
            psi_D_edges(s)(3-cut).theta=1-psi_D_edges(s)(3-cut).theta;}}
}
template<class T_OBJECT,class T_PSI_D_EDGES,class T> static void Initialize_Dirichlet_Boundary_Helper(T_OBJECT&,T_PSI_D_EDGES&,SEGMENTED_CURVE_2D<T>&,const T)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
template<class TV,int d> void LAPLACE_SIMPLICIAL_NODE<TV,d>::
Initialize_Dirichlet_Boundary(T_BOUNDARY_OBJECT& boundary_object,const T thickness_over_2,const bool verbose)
{
    if(verbose) LOG::Time("initializing Dirichlet surface");
    Initialize_Dirichlet_Boundary_Helper(object,psi_D_edges,boundary_object,thickness_over_2);
    int intersections=0;
    for(int s=1;s<=psi_D_edges.m;s++)for(int i=1;i<=2;i++)if(psi_D_edges(s)(i).simplex){
        intersections++;psi_D_edges(s)(i).theta=max(psi_D_edges(s)(i).theta,second_order_cut_cell_threshold);}
    u_boundary.Resize(boundary_object.particles.array_collection->Size());
    psi_D_mesh=&boundary_object.mesh;
    if(verbose) {std::stringstream ss;ss<<intersections<<" intersections found"<<std::endl;LOG::filecout(ss.str());}
    if(verbose) LOG::Stop_Time();
}
//#####################################################################
// Function Solve
//#####################################################################
template<class TV,int d> void LAPLACE_SIMPLICIAL_NODE<TV,d>::
Solve()
{
    ARRAY<int> matrix_indices;int count=Find_Matrix_Indices(matrix_indices);
    ARRAY<int> row_lengths(count);Find_Row_Lengths(row_lengths,matrix_indices);
    SPARSE_MATRIX_FLAT_NXN<T> A;A.Set_Row_Lengths(row_lengths);Find_A(A,matrix_indices);
    VECTOR_ND<T> b(count);Find_b(b,matrix_indices);
    Solve(matrix_indices,A,b);
}
//#####################################################################
// Function Solve
//#####################################################################
template<class TV,int d> void LAPLACE_SIMPLICIAL_NODE<TV,d>::
Solve(const ARRAY<int>& matrix_indices,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& b,const bool recompute_preconditioner)
{
    if(u.m != particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    Find_Tolerance(b);
    VECTOR_ND<T> x(b.n),q,s,r,k,z;for(int p=1;p<=diagonal.m;p++)if(matrix_indices(p)) x(matrix_indices(p))=u(p);
    pcg.Enforce_Compatibility(b.n==diagonal.m);
    if(pcg.show_results) {std::stringstream ss;ss<<"solving "<<b.n<<" nodes to tolerance "<<tolerance<<std::endl;LOG::filecout(ss.str());}
    pcg.Solve(A,x,b,q,s,r,k,z,tolerance,recompute_preconditioner);
    for(int p=1;p<=diagonal.m;p++)if(matrix_indices(p)) u(p)=x(matrix_indices(p));
}
//#####################################################################
// Function Find_Row_Lengths
//#####################################################################
template<class TV,int d> void LAPLACE_SIMPLICIAL_NODE<TV,d>::
Find_Row_Lengths(ARRAY<int>& row_lengths,const ARRAY<int>& matrix_indices)
{
    row_lengths+=1;
    SEGMENT_MESH& segment_mesh=*object.mesh.segment_mesh;
    for(int s=1;s<=segment_mesh.elements.m;s++){
        VECTOR<int,2> nodes(particles.subset_index_from_point_cloud_index.Subset(segment_mesh.elements(s)));
        int mi=matrix_indices(nodes[1]),mj=matrix_indices(nodes[2]);
        if(mi && mj && !(psi_D_edges.m && psi_D_edges(s)(1).simplex)){row_lengths(mi)++;row_lengths(mj)++;}}
}
//#####################################################################
// Function Find_A
//#####################################################################
template<class TV,int d> void LAPLACE_SIMPLICIAL_NODE<TV,d>::
Find_A(SPARSE_MATRIX_FLAT_NXN<T>& A,const ARRAY<int>& matrix_indices)
{
    SEGMENT_MESH& segment_mesh=*object.mesh.segment_mesh;
    for(int p=1;p<=diagonal.m;p++){int m=matrix_indices(p);
        if(m) A.Set_Element(m,m,diagonal(p));}
    for(int s=1;s<=segment_mesh.elements.m;s++){
        VECTOR<int,2> nodes(particles.subset_index_from_point_cloud_index.Subset(segment_mesh.elements(s)));
        int mi=matrix_indices(nodes[1]),mj=matrix_indices(nodes[2]);
        if(psi_D_edges.m && psi_D_edges(s)(1).simplex){
            T theta1=psi_D_edges(s)(1).theta,theta2=psi_D_edges(s)(2).theta;
            if(mi) A.Add_Element(mi,mi,-offdiagonal(s)*(1-theta1)/theta1);
            if(mj) A.Add_Element(mj,mj,-offdiagonal(s)*(1-theta2)/theta2);}
        else if(mi && mj) A.Set_Symmetric_Elements(mi,mj,offdiagonal(s));}
}
//#####################################################################
// Function Find_b
//#####################################################################
template<class TV,int d> void LAPLACE_SIMPLICIAL_NODE<TV,d>::
Find_b(VECTOR_ND<T>& b,const ARRAY<int>& matrix_indices)
{
    b.Set_Zero();
    if(u.m != particles.array_collection->Size()) PHYSBAM_FATAL_ERROR();
    if(psi_D_mesh && u_boundary.Size() != psi_D_mesh->number_nodes) PHYSBAM_FATAL_ERROR();
    SEGMENT_MESH& segment_mesh=*object.mesh.segment_mesh;
    for(int s=1;s<=offdiagonal.m;s++){
        VECTOR<int,2> nodes(particles.subset_index_from_point_cloud_index.Subset(segment_mesh.elements(s)));
        int mi=matrix_indices(nodes[1]),mj=matrix_indices(nodes[2]);
        if(psi_D_edges.m && psi_D_edges(s)(1).simplex){
            PSI_D_EDGE e1,e2;psi_D_edges(s).Get(e1,e2);
            if(mi){
                VECTOR<int,TV::m>& nodes=psi_D_mesh->elements(e1.simplex);
                T u_point=0;for(int i=1;i<=TV::m;i++) u_point+=e1.weights[i]*u_boundary(nodes[i]);
                b(mi)-=offdiagonal(s)/e1.theta*u_point;}
            if(mj){
                VECTOR<int,TV::m>& nodes=psi_D_mesh->elements(e2.simplex);
                T u_point=0;for(int i=1;i<=TV::m;i++) u_point+=e2.weights[i]*u_boundary(nodes[i]);
                b(mj)-=offdiagonal(s)/e2.theta*u_point;}}
        else if(mi && !mj) b(mi)-=offdiagonal(s)*u(nodes[2]);
        else if(!mi && mj) b(mj)-=offdiagonal(s)*u(nodes[1]);}
}
//#####################################################################
// Function Find_Matrix_Indices
//#####################################################################
template<class TV,int d> int LAPLACE_SIMPLICIAL_NODE<TV,d>::
Find_Matrix_Indices(ARRAY<int>& matrix_indices)
{
    int count=0;
    matrix_indices.Resize(diagonal.m);
    for(int p=1;p<=diagonal.m;p++)if(diagonal(p) && !psi_D(p)) matrix_indices(p)=++count;
    return count;
}
//#####################################################################
template class LAPLACE_SIMPLICIAL_NODE<VECTOR<float,2>,2>;
template class LAPLACE_SIMPLICIAL_NODE<VECTOR<float,3>,2>;
template class LAPLACE_SIMPLICIAL_NODE<VECTOR<float,3>,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LAPLACE_SIMPLICIAL_NODE<VECTOR<double,2>,2>;
template class LAPLACE_SIMPLICIAL_NODE<VECTOR<double,3>,2>;
template class LAPLACE_SIMPLICIAL_NODE<VECTOR<double,3>,3>;
#endif
