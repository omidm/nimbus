//#####################################################################
// Copyright 2012
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_CHIMERA_MPI
//#####################################################################
#ifndef __LAPLACE_CHIMERA_MPI__
#define __LAPLACE_CHIMERA_MPI__

#include <PhysBAM_Dynamics/Grids_Chimera/Parallel_Computation/CHIMERA_GRID.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COUPLED_SYSTEM_VECTOR.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology/TETRAHEDRON_MESH.h>
#include <PhysBAM_Geometry/Basic_Geometry/POLYGON.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/CHIMERA_SYSTEM.h>
namespace PhysBAM{

template<class T_GRID>
class LAPLACE_CHIMERA_MPI_CALLBACKS:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
public:
    virtual bool Get_Neumann_Boundary_Condition(const TV& location,const TV& normal,T& normal_gradient,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return false;}
    virtual bool Get_Dirichlet_Boundary_Condition(const TV& location,T& value,const T time){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return false;}
    virtual T Get_Density(const TV& location){PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return (T)1;};
};

//NOTES
//flow is from first to second cell on voronoi faces
//normal=normalize(location(voronoi_face(2))-location(voronoi_face(1)))
//volume weighted divergence is net flow out of cell

template<class T_GRID> class LAPLACE_CHIMERA_GRID_MPI;

template<class T_GRID>
class LAPLACE_CHIMERA_MPI:public NONCOPYABLE
{
    typedef char BYTE;
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    typedef FACE_INDEX<TV::dimension> D_FACE_INDEX;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef typename T_ARRAYS_SCALAR::template REBIND<BYTE>::TYPE T_ARRAYS_BYTE;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_BASE::template REBIND<bool>::TYPE T_ARRAYS_BASE_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<ARRAY<int> >::TYPE T_ARRAYS_ARRAYS_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<int>::TYPE T_FACE_ARRAYS_INT;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_BOOL::template REBIND<VECTOR<bool,T_GRID::dimension> >::TYPE T_FACE_ARRAYS_BOOL_DIMENSION;
public:
    typedef PAIR<int,TV_INT> GRID_CELL_INDEX;
    typedef PAIR<int,TV_INT> GRID_FACE_INDEX;//negative grid index indicates voronoi face
    typedef VECTOR<GRID_CELL_INDEX,2> VORONOI_FACE_INDICES;

    LAPLACE_CHIMERA_MPI(CHIMERA_GRID<T_GRID>& chimera_grid_input,LAPLACE_CHIMERA_GRID_MPI<T_GRID>& laplace_grid_input)
        :chimera_grid(chimera_grid_input),laplace_grid(laplace_grid_input)
    {}
    
    CHIMERA_GRID<T_GRID>& chimera_grid;
    LAPLACE_CHIMERA_GRID_MPI<T_GRID>& laplace_grid;
    
    SPARSE_MATRIX_FLAT_MXN<T> divergence_matrix;
    SPARSE_MATRIX_FLAT_NXN<T> system_matrix;//not diffusion:negative laplacian. diffusion:volume matrix minus dt time negative laplacian

    //STATISTICS
    int n_global_grids;
    int n_local_grids;

    int n_matrix_cells;
    int n_local_matrix_cells;
    int n_matrix_faces;

    //LOCAL GRID CHIMERA CELL INDICES
    //ARRAY<ARRAY<HASHTABLE<TV_INT,bool> > > face_invalid_indices;
    
    ARRAY<ARRAY<PAIR<bool,T> > > boundary_cell_psi_D;
    ARRAY<ARRAY<T> > boundary_cell_phi;
    ARRAY<ARRAY<bool> > boundary_cell_active;

    //MPI structures
    SPARSE_MATRIX_PARTITION partition;

    ARRAY<PAIR<bool,T> > voronoi_face_psi_N_velocities;

    //MATRIX INDICES
    ARRAY<T_ARRAYS_INT> cell_indices_to_matrix_cell_indices;
    ARRAY<HASHTABLE<TV_INT,int> > boundary_cell_indices_to_matrix_cell_indices;
    ARRAY<T_FACE_ARRAYS_INT> face_indices_to_matrix_face_indices;
    ARRAY<int> voronoi_face_indices_to_matrix_face_indices;
    
    VECTOR_ND<T> dual_cell_inverse_mass;

    void Clean_Memory()
    {
        divergence_matrix.Reset(0);
        system_matrix.Reset();
        boundary_cell_psi_D.Clean_Memory();
        boundary_cell_phi.Clean_Memory();
        boundary_cell_active.Clean_Memory();
        voronoi_face_psi_N_velocities.Clean_Memory();
        cell_indices_to_matrix_cell_indices.Clean_Memory();
        boundary_cell_indices_to_matrix_cell_indices.Clean_Memory();
        face_indices_to_matrix_face_indices.Clean_Memory();
        voronoi_face_indices_to_matrix_face_indices.Clean_Memory();
    }

    bool Psi_D(const int grid_index,const TV_INT& cell_index,T& pressure) const;
    T Phi(const int grid_index,const TV_INT& cell_index) const;
    T Dual_Cell_Volume_Fraction(const int grid_index_1,const TV_INT& cell_index_1,const int grid_index_2,const TV_INT& cell_index_2) const;
    bool Neumann_Pocket(const int grid_index,const TV_INT& cell_index) const;
    bool Active_Cell_Compute(const int grid_index,const TV_INT& cell_index) const;
    
    int Matrix_Cell_Index(const int grid_index,const TV_INT& cell_index) const;
    int Matrix_Face_Index(const int grid_index,const D_FACE_INDEX& face_index) const;
    int Matrix_Face_Index(const int voronoi_face_index) const;

    void Construct_Matrix_Indices(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time); //COMPUTE MATRIX INDICES FOR CELLS AND FACES ON LOCAL GRID AND NEIGHBORING GRID BOUNDARY CELLS
    void Construct_Matrix_Partitions(); //BUILD MATRIX PARTITION FOR LOCAL AND NONLOCAL BOUNDARY GRIDS
    void Construct_Laplacian(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks);
    void Construct_Divergence(); //BUILD DIVERGENCE MATRIX
    
    T Inverse_Mass(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const int grid_index,const D_FACE_INDEX& face_index);
    T Inverse_Mass(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const int voronoi_face_index);
    void Construct_Inverse_Mass(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks); //BUILD INVERSE MASS MATRIX

    void Compute_Neumann_Divergence(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time,VECTOR_ND<T>& divergence);
    void Compute_Dirichlet_Divergence(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time,VECTOR_ND<T>& divergence);

    void Construct_System(LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time,const T dt,const bool diffusion);
    void Multiply_By_Volume_Weighted_Divergence(const VECTOR_ND<T>& x,VECTOR_ND<T>& result);
    void Multiply_By_Gradient(const VECTOR_ND<T>& x,VECTOR_ND<T>& result);
    void Multiply_By_Inverse_Mass(const VECTOR_ND<T>& x,VECTOR_ND<T>& result);
    void Multiply_By_Cell_Size(const VECTOR_ND<T>& x,VECTOR_ND<T>& result);
    void Multiply_By_Dual_Cell_Size(const VECTOR_ND<T>& x,VECTOR_ND<T>& result);
};

template<class TV> typename BASIC_GEOMETRY_POLICY<TV>::CONVEX_POLYGON Unclipped_Voronoi_Face_Helper(const TV& normal,const TV& center,const typename TV::SCALAR radius);

}
#endif
