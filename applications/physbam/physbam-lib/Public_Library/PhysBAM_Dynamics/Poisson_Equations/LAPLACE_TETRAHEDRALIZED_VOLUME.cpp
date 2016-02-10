//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Dynamics/Poisson_Equations/LAPLACE_TETRAHEDRALIZED_VOLUME.h>
using namespace PhysBAM;
//#####################################################################
// Function Solve
//#####################################################################
template<class T> void LAPLACE_TETRAHEDRALIZED_VOLUME<T>::
Initialize_Tetrahedralized_Volume()
{
    TETRAHEDRON_MESH& mesh=tetrahedralized_volume.mesh;
    if(!mesh.triangle_tetrahedrons) mesh.Initialize_Triangle_Tetrahedrons();
    tetrahedralized_volume.Compute_Tetrahedron_Volumes();
    f.Resize(mesh.elements.m);
}
//#####################################################################
// Function Solve
//#####################################################################
template<class T> void LAPLACE_TETRAHEDRALIZED_VOLUME<T>::
Solve(const T time)
{
    SPARSE_MATRIX_FLAT_NXN<T> A;Find_A(A);
    VECTOR_ND<T> b;Find_b(b);
    Solve(A,b);
}
//#####################################################################
// Function Solve
//#####################################################################
template<class T> void LAPLACE_TETRAHEDRALIZED_VOLUME<T>::
Solve(SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& b,const bool enforce_compatability,const bool recompute_preconditioner)
{
    pcg.Enforce_Compatibility(enforce_compatability);
    Find_Tolerance(b); // needs to happen after b is completely set up
    VECTOR_ND<T> x(u.m),q,s,r,k,z;for(int t=1;t<=u.m;t++)x(t)=u(t);
    pcg.Solve(A,x,b,q,s,r,k,z,tolerance,recompute_preconditioner);
    for(int t=1;t<=u.m;t++)u(t)=x(t);
}
//#####################################################################
// Function Find_A
//#####################################################################
template<class T> void LAPLACE_TETRAHEDRALIZED_VOLUME<T>::
Find_A(SPARSE_MATRIX_FLAT_NXN<T>& A)
{
    TETRAHEDRON_MESH& mesh=tetrahedralized_volume.mesh;
    GEOMETRY_PARTICLES<VECTOR<T,3> >& particles=tetrahedralized_volume.particles;
    TRIANGLE_MESH& face_mesh=*mesh.triangle_mesh;ARRAY<VECTOR<int,2> >& face_tetrahedrons=*mesh.triangle_tetrahedrons;
    ARRAY<VECTOR<T,3> > center(mesh.elements.m);
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        center(t)=(T).25*(particles.X(i)+particles.X(j)+particles.X(k)+particles.X(l));}
    ARRAY<int> row_lengths(mesh.elements.m,false);ARRAYS_COMPUTATIONS::Fill(row_lengths,1);
    for(int f=1;f<=face_tetrahedrons.m;f++){
        int t1,t2;face_tetrahedrons(f).Get(t1,t2);if(!t1||!t2) continue;
        row_lengths(t1)++;row_lengths(t2)++;}
    A.Set_Row_Lengths(row_lengths);
    ARRAY<T> row_sum(A.n);
    for(int f=1;f<=face_tetrahedrons.m;f++){
        int t1,t2;face_tetrahedrons(f).Get(t1,t2);if(!t1||!t2) continue;
        int i,j,k;face_mesh.elements(f).Get(i,j,k);
        VECTOR<T,3> normal=TRIANGLE_3D<T>::Normal_Direction(particles.X(i),particles.X(j),particles.X(k));
        T entry=(T)-.5*normal.Magnitude_Squared()/VECTOR<T,3>::Dot_Product(center(t2)-center(t1),normal); // area/distance
        assert(entry<0);row_sum(t1)+=entry;row_sum(t2)+=entry;
        A.Set_Symmetric_Elements(t1,t2,entry);}
    for(int t=1;t<=mesh.elements.m;t++)A.Set_Element(t,t,-row_sum(t));
}
//#####################################################################
// Function Find_A
//#####################################################################
template<class T> void LAPLACE_TETRAHEDRALIZED_VOLUME<T>::
Find_b(VECTOR_ND<T>& b)
{
    ARRAY<T>& volume=*tetrahedralized_volume.tetrahedron_volumes;
    b.Resize(volume.m);
    for(int t=1;t<=volume.m;t++){assert(volume(t)>0);b(t)=volume(t)*f(t);}
}
//#####################################################################
// Function Add_Scalar_Term
//#####################################################################
template<class T> void LAPLACE_TETRAHEDRALIZED_VOLUME<T>::
Add_Scalar_Term(SPARSE_MATRIX_FLAT_NXN<T>& A,const T s)
{
    ARRAY<T>& volume=*tetrahedralized_volume.tetrahedron_volumes;
    for(int t=1;t<=volume.m;t++)A.Set_Element(t,t,A(t,t)+s*volume(t));
}
//#####################################################################
template class LAPLACE_TETRAHEDRALIZED_VOLUME<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LAPLACE_TETRAHEDRALIZED_VOLUME<double>;
#endif
