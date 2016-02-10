//#####################################################################
// Copyright 2012
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_CHIMERA_GRID_MPI
//#####################################################################
#ifndef __LAPLACE_CHIMERA_GRID_MPI__
#define __LAPLACE_CHIMERA_GRID_MPI__

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
#include <PhysBAM_Dynamics/Grids_Chimera_Interpolation/LINEAR_INTERPOLATION_CHIMERA.h>

namespace PhysBAM{

template<class TV>
class VORONOI_GEOMETRY
{
public:
    typedef char BYTE;
    typedef GRID<TV> T_GRID;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension>::SIMPLEX_FACE T_SIMPLEX_FACE;
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::dimension> TV_INT;
    ARRAY<TV> cell_vertices;
    typename MESH_POLICY<TV::dimension>::MESH delaunay_mesh;
    typedef PAIR<int,TV_INT> GRID_CELL_INDEX;
    typedef VECTOR<GRID_CELL_INDEX,2> VORONOI_FACE_INDICES;
    typedef PAIR<VORONOI_FACE_INDICES,T> VORONOI_FACE;
    //typedef PAIR<VORONOI_FACE_INDICES,ARRAY<T_SIMPLEX_FACE> > VORONOI_FACE;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef typename T_ARRAYS_SCALAR::template REBIND<BYTE>::TYPE T_ARRAYS_BYTE;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<int>::TYPE T_FACE_ARRAYS_INT;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    
    ARRAY<VORONOI_FACE> voronoi_faces;
    ARRAY<typename BASIC_GEOMETRY_POLICY<TV>::CONVEX_POLYGON> voronoi_faces_geometry;
    ARRAY<PAIR<bool,T> > voronoi_psi_N;
    ARRAY<T_ARRAYS_INT> cell_indices_to_chimera_cell_status;
    ARRAY<T_ARRAYS_INT> cell_indices_to_matrix_cell_indices;
    ARRAY<T_FACE_ARRAYS_INT> face_indices_to_matrix_face_indices;
    ARRAY<int> voronoi_face_indices_to_matrix_face_indices;
    ARRAY<VECTOR<TV,TV::dimension+1> > boundary_delaunay_simplex_vertices;
    ARRAY<T> voronoi_face_velocities;
    ARRAY<VECTOR<TV,2> > vectors;

    VORONOI_GEOMETRY& operator=(const VORONOI_GEOMETRY& other)
    {
        voronoi_faces=other.voronoi_faces;
        voronoi_faces_geometry=other.voronoi_faces_geometry;
        voronoi_psi_N=other.voronoi_psi_N;
        cell_indices_to_chimera_cell_status=other.cell_indices_to_chimera_cell_status;
        cell_indices_to_matrix_cell_indices=other.cell_indices_to_matrix_cell_indices;
        face_indices_to_matrix_face_indices=other.face_indices_to_matrix_face_indices;
        voronoi_face_indices_to_matrix_face_indices=other.voronoi_face_indices_to_matrix_face_indices;
        boundary_delaunay_simplex_vertices=other.boundary_delaunay_simplex_vertices;
        voronoi_face_velocities=other.voronoi_face_velocities;
        vectors=other.vectors;
        return *this;
    }
};

template<class RW,class TV>
class Read_Write<VORONOI_GEOMETRY<TV>,RW>
{
public:
    static void Read(std::istream& input,VORONOI_GEOMETRY<TV>& object)
    {
        Read_Binary<RW>(input,
            object.voronoi_faces,
            object.voronoi_faces_geometry,
            object.voronoi_psi_N,
            object.cell_indices_to_chimera_cell_status,
            object.cell_indices_to_matrix_cell_indices,
            object.voronoi_face_indices_to_matrix_face_indices,
                        object.boundary_delaunay_simplex_vertices,
                        object.voronoi_face_velocities,
                        object.vectors);
    }
    static void Write(std::ostream& output,const VORONOI_GEOMETRY<TV>& object)
    {
        Write_Binary<RW>(output,
            object.voronoi_faces,
            object.voronoi_faces_geometry,
            object.voronoi_psi_N,
            object.cell_indices_to_chimera_cell_status,
            object.cell_indices_to_matrix_cell_indices,
            object.voronoi_face_indices_to_matrix_face_indices,
                         object.boundary_delaunay_simplex_vertices,
                         object.voronoi_face_velocities,
                         object.vectors);
    }
};

//NOTES
//flow is from second to first cell on voronoi faces
//normal=normalize(location(voronoi_face(1))-location(voronoi_face(2)))
template<class T_GRID>
class LAPLACE_CHIMERA_GRID_MPI:public NONCOPYABLE
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
    enum CHIMERA_STATUS{CHIMERA_INACTIVE=0,CHIMERA_ACTIVE=1,CHIMERA_BOUNDARY=2};
    
    typedef PAIR<int,TV_INT> GRID_CELL_INDEX;
    typedef PAIR<int,TV_INT> GRID_FACE_INDEX;//negative grid index indicates voronoi face
    
    //VORONOI FACES
    typedef VECTOR<GRID_CELL_INDEX,2> VORONOI_FACE_INDICES;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::CONVEX_POLYGON CONVEX_POLYGON;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::HYPERPLANE HYPERPLANE;
    typedef PAIR<VORONOI_FACE_INDICES,T> VORONOI_FACE;
    
    LAPLACE_CHIMERA_GRID_MPI(CHIMERA_GRID<T_GRID>& chimera_grid_input)
        :chimera_grid(chimera_grid_input),interpolation(*this){}
    
    CHIMERA_GRID<T_GRID>& chimera_grid;
    int n_local_grids;
    int n_global_grids;
    ARRAY<ARRAY<bool> > is_boundary_grid;
    
    int grid_index_coarsest; //EETODO: we shouldn't need this really
    
    //CHIMERA CELL STATUS 0=inactive 1=active 2=active+boundary
    ARRAY<T_ARRAYS_INT> cell_indices_to_chimera_cell_status;
    int n_chimera_cells; //includes all non-inactive cells

    //BOUNDARY GRID CHIMERA CELL INDICES
    
    //LOCAL BOUNDARY INFORMATION
    ARRAY<ARRAY<TV_INT> > boundary_cell_indices;
    ARRAY<ARRAY<bool> > boundary_cell_chimera_cell_status;
    ARRAY<ARRAY<bool> > boundary_cell_split;
    
    //NONLOCAL BOUNDARY INFORMATION
    ARRAY<HASHTABLE<TV_INT,int> > boundary_cell_indices_to_linear_index;
    ARRAY<HASHTABLE<TV_INT,int> > boundary_cell_indices_to_chimera_cell_status;
    
    ARRAY<VORONOI_FACE> voronoi_faces;
    ARRAY<CONVEX_POLYGON> voronoi_faces_geometry;
    HASHTABLE<GRID_CELL_INDEX,ARRAY<int> > cell_indices_to_incident_voronoi_face_indices;
    ARRAY<T> cell_sizes,dual_cell_sizes;
    
    LINEAR_INTERPOLATION_CHIMERA<T_GRID> interpolation;
    
    ///////////////DEBUG
    VORONOI_GEOMETRY<TV> voronoi;

    void Clean_Memory()
    {
        is_boundary_grid.Clean_Memory();
        cell_indices_to_chimera_cell_status.Clean_Memory();
        boundary_cell_indices.Clean_Memory();
        boundary_cell_chimera_cell_status.Clean_Memory();
        boundary_cell_indices_to_linear_index.Clean_Memory();
        boundary_cell_indices_to_chimera_cell_status.Clean_Memory();
        voronoi_faces.Clean_Memory();
        voronoi_faces_geometry.Clean_Memory();
        cell_indices_to_incident_voronoi_face_indices.Clean_Memory();
        cell_sizes.Clean_Memory();
        dual_cell_sizes.Clean_Memory();
        interpolation.Clean_Memory();
    }
    
    bool Chimera_Cell_Compute(const int grid_index,const TV_INT& cell_index) const;
    bool Chimera_Cell(const int grid_index,const TV_INT& cell_index) const;
    int Chimera_Cell_Status(const int grid_index,const TV_INT& cell_index) const;
    bool Chimera_Face(const int grid_index,const D_FACE_INDEX& face_index) const;
    bool Boundary_Cell(const int grid_index,const TV_INT& cell_index) const;
    bool Boundary_Cell_Incident_To_Local_Grid(const int grid_index,const TV_INT& cell_index) const;

    void Construct_Chimera_Cells(); //COMPUTE CELLS NOT COVERED BY FINER GRIDS
    void Construct_Boundary_Cells(); //CONSTRUCT ARRAYS OF CELLS ON THE BOUNDARY OF THE VORONOI REGION
    typename BASIC_GEOMETRY_POLICY<TV>::CONVEX_POLYGON Unclipped_Voronoi_Face(const TV& normal,const TV& center,const T radius) const;
    typename BASIC_GEOMETRY_POLICY<TV>::CONVEX_POLYGON Cartesian_Face(const int grid_index,const D_FACE_INDEX& face_index) const;
    RANGE<TV_INT> Compute_Candidate_Neighbor_Cells(const int grid_index_1,const TV_INT& cell_index_1,const int grid_index_2,const T& search_radius) const;
    int Find_Containing_Grid(const TV& location);
    bool Cells_Occluded(const int grid_index_1,const TV_INT& cell_index_1,const int grid_index_2,const TV_INT& cell_index_2);
    void Construct_Voronoi_Faces(); //CONSTRUCT FACES OF THE PRESSURE SAMPLE VORONOI DIAGRAM
    void Construct_Mesh();

    T Face_Size(const int grid_index,const D_FACE_INDEX& face_index) const;
    T Face_Size(const int voronoi_face_index) const;
    T Cell_Size(const int grid_index,const TV_INT& cell_index) const;

    //SINGLE GRID ACCESS FUNCTIONS, PROBABLY SHOULD JUST HAVE THIS FUNCTIONALITY GO STRAIGHT TO CHIMERA_GRID
    T_GRID& Grid(const int grid_index) const;
    RIGID_GRID<T_GRID>& Rigid_Grid(const int grid_index) const;
    const FRAME<TV> Frame(const int grid_index) const;
    bool Local_Grid(const int grid_index) const;
    int Local_Grid_Index(const int global_grid_index) const;
    int Global_Grid_Index(const int local_grid_index) const;
    bool Boundary_Grid(const int grid_index) const;
    bool Inside_Domain(const int grid_index,const TV_INT& cell_index) const;

    inline void Voronoi_Face(const int voronoi_face_index,TV& location_1,TV& location_2,TV& location,TV& normal,T& distance)
    {
        const VECTOR<GRID_CELL_INDEX,2>& grid_cell_indices=voronoi_faces(voronoi_face_index).x;
        location_1=Frame(grid_cell_indices(1).x)*Grid(grid_cell_indices(1).x).X(grid_cell_indices(1).y);
        location_2=Frame(grid_cell_indices(2).x)*Grid(grid_cell_indices(2).x).X(grid_cell_indices(2).y);
        location=(T).5*(location_1+location_2);
        normal=location_2-location_1;
        distance=normal.Magnitude();
        normal/=distance;
    }
};
}

#endif

