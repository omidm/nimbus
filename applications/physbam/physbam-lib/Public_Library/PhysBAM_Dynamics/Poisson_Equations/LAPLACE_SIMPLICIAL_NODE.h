//#####################################################################
// Copyright 2005-2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_SIMPLICIAL_NODE
//#####################################################################
//
// Solves -laplace u = 0, with Neumann boundary conditions at the surface.
// Set psi_D=true at Dirichlet boundary condition points.
//
//#####################################################################  
#ifndef __LAPLACE_SIMPLICIAL_NODE__
#define __LAPLACE_SIMPLICIAL_NODE__

#include <PhysBAM_Tools/Grids_PDE_Linear/LAPLACE.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD_SUBSET.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
namespace PhysBAM{

template<class T> class SPARSE_MATRIX_FLAT_NXN;

template<class TV,int d>
class LAPLACE_SIMPLICIAL_NODE:public LAPLACE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d>::OBJECT T_OBJECT;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_BOUNDARY_OBJECT;
    typedef typename MESH_POLICY<TV::m-1>::MESH T_BOUNDARY_MESH;
public:
    using LAPLACE<T>::tolerance;

    POINT_CLOUD_SUBSET<TV,PARTICLES<TV> >& particles;
    T_OBJECT& object;
    ARRAY<T>& u;
    PCG_SPARSE<T> pcg;
    ARRAY<bool> psi_D;
    ARRAY<T> u_boundary;
private:
    ARRAY<T> diagonal;
    ARRAY<T> offdiagonal;
    T_BOUNDARY_MESH* psi_D_mesh;
    struct PSI_D_EDGE
    {
        T theta;int simplex;TV weights;
        PSI_D_EDGE():theta(0),simplex(0){}
    };
    ARRAY<VECTOR<PSI_D_EDGE,2> > psi_D_edges;
public:
    T second_order_cut_cell_threshold;

    LAPLACE_SIMPLICIAL_NODE(POINT_CLOUD_SUBSET<TV,PARTICLES<TV> >& particles_input,T_OBJECT& object_input,ARRAY<T>& u_input)
        :particles(particles_input),object(object_input),u(u_input),psi_D_mesh(0),second_order_cut_cell_threshold((T)1e-3)
    {}

//#####################################################################
    void Initialize_Object();
    void Initialize_Dirichlet_Boundary(T_BOUNDARY_OBJECT& boundary_object,const T thickness_over_2,const bool verbose);
    virtual void Solve();
    void Solve(const ARRAY<int>& matrix_indices,SPARSE_MATRIX_FLAT_NXN<T>& A,VECTOR_ND<T>& b,const bool recompute_preconditioner=true);
    virtual void Find_Row_Lengths(ARRAY<int>& row_lengths,const ARRAY<int>& matrix_indices);
    virtual void Find_A(SPARSE_MATRIX_FLAT_NXN<T>& A,const ARRAY<int>& matrix_indices);
    virtual void Find_b(VECTOR_ND<T>& b,const ARRAY<int>& matrix_indices);
    int Find_Matrix_Indices(ARRAY<int>& matrix_indices);
//#####################################################################
};
}
#endif
