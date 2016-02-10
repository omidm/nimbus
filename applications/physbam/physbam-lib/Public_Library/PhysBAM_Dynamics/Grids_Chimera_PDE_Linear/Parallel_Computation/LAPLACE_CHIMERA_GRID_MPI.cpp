//#####################################################################
// Copyright 2012
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAPLACE_CHIMERA_GRID_MPI
//#####################################################################
#include <PhysBAM_Dynamics/Grids_Chimera_PDE_Linear/Parallel_Computation/LAPLACE_CHIMERA_GRID_MPI.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Utilities/INTERRUPTS.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Dynamics/Meshing/POLYGONAL_TRIANGULATION.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/POLYGON_HYPERPLANE_INTERSECTION.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/POLYGON.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Utilities/TIMER.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
using namespace PhysBAM;
//#####################################################################
// Grid
//#####################################################################
template<class T_GRID> T_GRID& LAPLACE_CHIMERA_GRID_MPI<T_GRID>::
Grid(const int grid_index) const
{
    return chimera_grid.global_grid.Rigid_Grid(grid_index).grid;
}
//#####################################################################
// Rigid_Grid
//#####################################################################
template<class T_GRID> RIGID_GRID<T_GRID>& LAPLACE_CHIMERA_GRID_MPI<T_GRID>::
Rigid_Grid(const int grid_index) const
{
    return chimera_grid.global_grid.Rigid_Grid(grid_index);
}
//#####################################################################
// Frame
//#####################################################################
template<class T_GRID> const FRAME<typename T_GRID::VECTOR_T> LAPLACE_CHIMERA_GRID_MPI<T_GRID>::
Frame(const int grid_index) const
{
    return chimera_grid.global_grid.Rigid_Grid(grid_index).Frame();
}
//#####################################################################
// Local_Grid
//#####################################################################
template<class T_GRID> bool LAPLACE_CHIMERA_GRID_MPI<T_GRID>::
Local_Grid(const int grid_index) const
{
    return (chimera_grid.global_to_local_grid_index_map(grid_index))==0?false:true;
}
//#####################################################################
// Local_Grid_Index
//#####################################################################
template<class T_GRID> int LAPLACE_CHIMERA_GRID_MPI<T_GRID>::Local_Grid_Index(const int global_grid_index) const
{
    return chimera_grid.global_to_local_grid_index_map(global_grid_index);
}
//#####################################################################
// Global_Grid_Index
//#####################################################################
template<class T_GRID> int LAPLACE_CHIMERA_GRID_MPI<T_GRID>::Global_Grid_Index(const int local_grid_index) const
{
    return chimera_grid.local_to_global_grid_index_map(local_grid_index);
}
//#####################################################################
// Boundary_Grid
//#####################################################################
template<class T_GRID> bool LAPLACE_CHIMERA_GRID_MPI<T_GRID>::
Boundary_Grid(const int grid_index) const
{
    return boundary_cell_indices_to_chimera_cell_status(grid_index).Size()==0?false:true;
}
//#####################################################################
// Inside_Domain
//#####################################################################
template<class T_GRID> bool LAPLACE_CHIMERA_GRID_MPI<T_GRID>::
Inside_Domain(const int grid_index,const TV_INT& cell_index) const
{
    TV location=Rigid_Grid(grid_index).Frame()*Grid(grid_index).X(cell_index);
    for(int grid_index_containing=1;grid_index_containing<=n_global_grids;grid_index_containing++)
        if(Grid(grid_index_containing).Domain().Lazy_Inside(Rigid_Grid(grid_index_containing).Frame().Inverse()*location))
            return true;
    return false;
}
//#####################################################################
// Chimera_Cell_Compute
//#####################################################################
template<class T_GRID> bool LAPLACE_CHIMERA_GRID_MPI<T_GRID>::
Chimera_Cell_Compute(const int grid_index,const TV_INT& cell_index) const
{
    T tolerance=(T)1e-6;

    if(Grid(grid_index).Inside_Domain(cell_index) || (Grid(grid_index).Inside_Domain(cell_index,1) && !Inside_Domain(grid_index,cell_index))){
        //first, check if this cell's center is in a finer cell
        TV location=Rigid_Grid(grid_index).Frame()*Grid(grid_index).X(cell_index);
        T cell_size=Grid(grid_index).DX().Magnitude();
        for(int grid_index_fine=1;grid_index_fine<=n_global_grids;grid_index_fine++){
            T cell_size_fine=Grid(grid_index_fine).DX().Magnitude();
            if(cell_size_fine<cell_size || (cell_size_fine==cell_size && grid_index_fine<grid_index)){
                TV location_fine=Rigid_Grid(grid_index_fine).Frame().Inverse_Times(location);
                if(Grid(grid_index_fine).Domain().Inside(location_fine,-tolerance*cell_size) && Chimera_Cell_Compute(grid_index_fine,Grid(grid_index_fine).Clamp_To_Cell(location_fine)))
                    return false;}}
        
        //second, check if this cell contains a finer cell center
        RANGE<TV> cell_domain=Grid(grid_index).Cell_Domain(cell_index);
        cell_domain.Change_Size(cell_domain.Edge_Lengths()*tolerance);
        for(int grid_index_fine=1;grid_index_fine<=n_global_grids;grid_index_fine++){
            T cell_size_fine=Grid(grid_index_fine).DX().Magnitude();
            if(grid_index_fine!=grid_index && (cell_size_fine<cell_size || (cell_size_fine==cell_size && grid_index>grid_index_fine))){
                RANGE<TV> range_fine_location=ORIENTED_BOX<TV>(cell_domain,Rigid_Grid(grid_index_fine).Frame().Inverse()*Rigid_Grid(grid_index).Frame()).Axis_Aligned_Bounding_Box();
                RANGE<TV_INT> range_fine(Grid(grid_index_fine).Cell(range_fine_location.min_corner),Grid(grid_index_fine).Cell(range_fine_location.max_corner));
                for(CELL_ITERATOR iterator_fine(Grid(grid_index_fine),RANGE<TV_INT>::Intersect(range_fine,Grid(grid_index_fine).Domain_Indices()));iterator_fine.Valid();iterator_fine.Next()){
                    TV_INT cell_index_fine=iterator_fine.Cell_Index();
                    TV location=Rigid_Grid(grid_index).Frame().Inverse()*Rigid_Grid(grid_index_fine).Frame()*iterator_fine.Location();
                    if(cell_domain.Lazy_Inside(location) && Chimera_Cell_Compute(grid_index_fine,cell_index_fine))
                        return false;}}}

        return true;}
    return false;
}
//#####################################################################
// Chimera_Cell
//#####################################################################
template<class T_GRID> bool LAPLACE_CHIMERA_GRID_MPI<T_GRID>::Chimera_Cell(const int grid_index,const TV_INT& cell_index) const
{
    if(Grid(grid_index).Inside_Domain(cell_index,1)){
        if(Local_Grid(grid_index)) return cell_indices_to_chimera_cell_status(grid_index)(cell_index)==0?false:true;
        else return boundary_cell_indices_to_chimera_cell_status(grid_index).Get_Default(cell_index,CHIMERA_INACTIVE)==0?false:true;}
    return false;
}
//#####################################################################
// Chimera_Cell_Status
//#####################################################################
template<class T_GRID> int LAPLACE_CHIMERA_GRID_MPI<T_GRID>::Chimera_Cell_Status(const int grid_index,const TV_INT& cell_index) const
{
    if(Local_Grid(grid_index)) return cell_indices_to_chimera_cell_status(grid_index)(cell_index);
    else return boundary_cell_indices_to_chimera_cell_status(grid_index).Get_Default(cell_index,CHIMERA_INACTIVE);
}
//#####################################################################
// Valid_Face
//#####################################################################
template<class T_GRID> bool LAPLACE_CHIMERA_GRID_MPI<T_GRID>::Chimera_Face(const int grid_index,const FACE_INDEX<TV::dimension>& face_index) const
{
    TV_INT cell_index_1=face_index.First_Cell_Index();
    TV_INT cell_index_2=face_index.Second_Cell_Index();
    return Chimera_Cell(grid_index,cell_index_1) && Chimera_Cell(grid_index,cell_index_2) && (!Boundary_Cell(grid_index,cell_index_1) || !Boundary_Cell(grid_index,cell_index_2));
}
//#####################################################################
// Boundary_Cell
//#####################################################################
template<class T_GRID> bool LAPLACE_CHIMERA_GRID_MPI<T_GRID>::Boundary_Cell(const int grid_index,const TV_INT& cell_index) const
{
    if(Local_Grid(grid_index)) return cell_indices_to_chimera_cell_status(grid_index)(cell_index)==CHIMERA_BOUNDARY;
    else return boundary_cell_indices_to_chimera_cell_status(grid_index).Get_Default(cell_index,CHIMERA_INACTIVE)==CHIMERA_BOUNDARY;
}
//#####################################################################
// Function Boundary_Cell_Incident_To_Local_Grid
//#####################################################################
template<class T_GRID> bool LAPLACE_CHIMERA_GRID_MPI<T_GRID>::Boundary_Cell_Incident_To_Local_Grid(const int grid_index,const TV_INT& cell_index) const
{
    GRID_CELL_INDEX grid_cell_index(grid_index,cell_index);
    if(Boundary_Grid(grid_index) && Boundary_Cell(grid_index,cell_index) && cell_indices_to_incident_voronoi_face_indices.Contains(grid_cell_index)){
        const ARRAY<int>& incident_voronoi_face_indices=cell_indices_to_incident_voronoi_face_indices.Get(grid_cell_index);
        for(int voronoi_face_index=1;voronoi_face_index<=incident_voronoi_face_indices.Size();voronoi_face_index++){
            const VORONOI_FACE_INDICES& indices=voronoi_faces(incident_voronoi_face_indices(voronoi_face_index)).x;
            if((Local_Grid(indices(1).x) || Local_Grid(indices(2).x)))
                return true;}}
    return false;
}
//#####################################################################
// Construct_Chimera_Cells
//#####################################################################
template<class T_GRID> void LAPLACE_CHIMERA_GRID_MPI<T_GRID>::
Construct_Chimera_Cells()
{
    LOG::SCOPE scope("Compute_Chimera_Cells");
    
    int n_cells=0;
    n_chimera_cells=0;

    LOG::SCOPE scope_chimera("Chimera_Cells");
    //DETERMINE CELL STATUS BY CHECKING FOR PROXIMITY TO HIGHER RESOLUTION GRIDS
    //EETODO: FOR EACH CELL CONSIDER EACH OVERLAPPING CELL FROM EACH GRID AND ALLOCATE ONLY IF THAT CELL IS LOWER RESOLUTION
    //EETODO: FIND BETTER DISTANCE METRIC
    cell_indices_to_chimera_cell_status.Resize(n_global_grids);
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(Local_Grid(grid_index)){
            cell_indices_to_chimera_cell_status(grid_index).Resize(Grid(grid_index).Domain_Indices(1));
            n_cells+=Grid(grid_index).Counts().Product();
            for(CELL_ITERATOR iterator(Grid(grid_index),1);iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                if(Chimera_Cell_Compute(grid_index,cell_index)){
                    if(Grid(grid_index).Inside_Domain(cell_index)) n_chimera_cells++;
                    cell_indices_to_chimera_cell_status(grid_index)(cell_index)=CHIMERA_ACTIVE;}
                else
                    cell_indices_to_chimera_cell_status(grid_index)(cell_index)=CHIMERA_INACTIVE;}}
    scope_chimera.Pop();

    LOG::SCOPE scope_boundary("Boundary_Cells");
    //CALCULATE BOUNDARY CELLS TO DETERMINE WHICH CELLS TO CONSIDER FOR VORONOI FACES
    boundary_cell_indices.Resize(n_global_grids);
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
        boundary_cell_indices(grid_index).Remove_All();
        if(Local_Grid(grid_index))
            for(CELL_ITERATOR iterator(Grid(grid_index));iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                if(Chimera_Cell(grid_index,cell_index)){
                    bool all_neighbors_chimera=true;
                    for(CELL_ITERATOR iterator_neighbor(Grid(grid_index),RANGE<TV_INT>(cell_index-1,cell_index+1));iterator_neighbor.Valid();iterator_neighbor.Next())
                        if(!cell_indices_to_chimera_cell_status(grid_index)(iterator_neighbor.Cell_Index())){
                            all_neighbors_chimera=false;break;}
                    if(!all_neighbors_chimera){
                        cell_indices_to_chimera_cell_status(grid_index)(cell_index)=CHIMERA_BOUNDARY;
                        boundary_cell_indices(grid_index).Append(cell_index);}}}}
    scope_boundary.Pop();
    
    LOG::cout << "number of cells " << n_cells << " chimera " << n_chimera_cells << std::endl;
    
    //DEBUG OUTPUT
    voronoi.cell_indices_to_chimera_cell_status=cell_indices_to_chimera_cell_status;
    //voronoi.boundary_cells=cell_indices_to_boundary_cells;
}
//#####################################################################
// Construct_Boundary_Cell_Arrays
//#####################################################################
template<class T_GRID> void LAPLACE_CHIMERA_GRID_MPI<T_GRID>::Construct_Boundary_Cells()
{
    LOG::SCOPE scope("Construct_Boundary_Cell_Arrays");

    is_boundary_grid.Resize(n_local_grids);

    //EETODO: make this robust and prove the distance bounds
    T dx=Grid(grid_index_coarsest).DX().Magnitude();//this should be computed per grid as the coarsest containing grid account for the search radius
    for(int local_grid_index=1;local_grid_index<=n_local_grids;local_grid_index++){
        is_boundary_grid(local_grid_index).Resize(n_global_grids);
        int this_grid_index=chimera_grid.local_to_global_grid_index_map(local_grid_index);
        ORIENTED_BOX<TV> this_grid_box(Grid(this_grid_index).Domain().Thickened(dx),Rigid_Grid(this_grid_index).Frame());
        for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
            if(Local_Grid(grid_index)){is_boundary_grid(local_grid_index)(grid_index)=false;continue;}
            ORIENTED_BOX<TV> grid_box(Grid(grid_index).Domain().Thickened(dx),Rigid_Grid(grid_index).Frame());
            is_boundary_grid(local_grid_index)(grid_index)=(this_grid_index<grid_index)?grid_box.Intersection(this_grid_box):this_grid_box.Intersection(grid_box);}}
   
    chimera_grid.Exchange_Boundary_Cell_Arrays(boundary_cell_indices,is_boundary_grid);
    
    boundary_cell_indices_to_chimera_cell_status.Resize(n_global_grids);
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
        boundary_cell_indices_to_chimera_cell_status(grid_index).Remove_All();
        if(!Local_Grid(grid_index)){
            for(int index=1;index<=boundary_cell_indices(grid_index).Size();index++){
                TV_INT cell_index=boundary_cell_indices(grid_index)(index);
                if(Grid(grid_index).Inside_Domain(cell_index)) n_chimera_cells++;
                boundary_cell_indices_to_chimera_cell_status(grid_index).Insert(cell_index,CHIMERA_BOUNDARY);}
            //we compute and add non-local boundary cell stencil-neighbor chimera cell status here, computed locally
            for(int index=1;index<=boundary_cell_indices(grid_index).Size();index++){
                TV_INT cell_index=boundary_cell_indices(grid_index)(index);
                for(CELL_ITERATOR iterator(Grid(grid_index),RANGE<TV_INT>(cell_index-1,cell_index+1));iterator.Valid();iterator.Next())
                    if(!boundary_cell_indices_to_chimera_cell_status(grid_index).Contains(iterator.Cell_Index()))
                        boundary_cell_indices_to_chimera_cell_status(grid_index).Insert(iterator.Cell_Index(),Chimera_Cell_Compute(grid_index,iterator.Cell_Index())?CHIMERA_ACTIVE:CHIMERA_INACTIVE);}}}
    
    boundary_cell_indices_to_linear_index.Resize(n_global_grids);
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(Local_Grid(grid_index) || Boundary_Grid(grid_index)){
            boundary_cell_indices_to_linear_index(grid_index).Clean_Memory();
            for(int boundary_cell_index=1;boundary_cell_index<=boundary_cell_indices(grid_index).Size();boundary_cell_index++)
                boundary_cell_indices_to_linear_index(grid_index).Insert(boundary_cell_indices(grid_index)(boundary_cell_index),boundary_cell_index);}
    
    //handle split grids explicitly
    boundary_cell_split.Resize(n_global_grids);
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
        boundary_cell_split(grid_index).Resize(boundary_cell_indices(grid_index).Size());
        if(Local_Grid(grid_index)){
            for(int boundary_cell_index=1;boundary_cell_index<=boundary_cell_indices(grid_index).Size();boundary_cell_index++){
                const TV_INT& cell_index=boundary_cell_indices(grid_index)(boundary_cell_index);
                bool all_neighbors_chimera=true;
                for(CELL_ITERATOR iterator_neighbor(Grid(grid_index),RANGE<TV_INT>(cell_index-1,cell_index+1));iterator_neighbor.Valid();iterator_neighbor.Next()){
                    GRID_CELL_INDEX real_index=chimera_grid.Find_Real_Grid_Cell_Index(GRID_CELL_INDEX(grid_index,iterator_neighbor.Cell_Index()));
                    if(!Chimera_Cell(real_index.x,real_index.y)){
                        all_neighbors_chimera=false;break;}}
                boundary_cell_split(grid_index)(boundary_cell_index)=all_neighbors_chimera;
                voronoi.cell_indices_to_chimera_cell_status(grid_index)(cell_index)=all_neighbors_chimera?CHIMERA_ACTIVE:CHIMERA_BOUNDARY;}}}
    
    chimera_grid.Exchange_Boundary_Scalar_Values(boundary_cell_split,is_boundary_grid);

    //voronoi.cell_indices_to_chimera_cell_status=cell_indices_to_chimera_cell_status;
}
//#####################################################################
// Function Unclipped_Voronoi_Face
//#####################################################################
template<class T>
typename BASIC_GEOMETRY_POLICY<VECTOR<T,1> >::CONVEX_POLYGON Unclipped_Voronoi_Face_Helper(const VECTOR<T,1>& normal,const VECTOR<T,1>& center,const T radius)
{
    return typename BASIC_GEOMETRY_POLICY<VECTOR<T,1> >::CONVEX_POLYGON(center,normal(1)>=0?true:false);
}
template<class T>
typename BASIC_GEOMETRY_POLICY<VECTOR<T,2> >::CONVEX_POLYGON Unclipped_Voronoi_Face_Helper(const VECTOR<T,2>& normal,const VECTOR<T,2>& center,const T radius)
{
    VECTOR<T,2> perpendicular=normal.Unit_Orthogonal_Vector()*radius;
    return typename BASIC_GEOMETRY_POLICY<VECTOR<T,2> >::CONVEX_POLYGON(center-perpendicular,center+perpendicular);
}
template<class T>
typename BASIC_GEOMETRY_POLICY<VECTOR<T,3> >::CONVEX_POLYGON Unclipped_Voronoi_Face_Helper(const VECTOR<T,3>& normal,const VECTOR<T,3>& center,const T radius)
{
    VECTOR<T,3> u=normal.Unit_Orthogonal_Vector();
    VECTOR<T,3> v=VECTOR<T,3>::Cross_Product(u,normal);
    typename BASIC_GEOMETRY_POLICY<VECTOR<T,3> >::CONVEX_POLYGON polygon;
    polygon.X.Resize(4);
    polygon.X(1)=center+radius*(-u-v);
    polygon.X(2)=center+radius*(u-v);
    polygon.X(3)=center+radius*(u+v);
    polygon.X(4)=center+radius*(-u+v);
    return polygon;
}
template<class T_GRID>
typename BASIC_GEOMETRY_POLICY<typename T_GRID::VECTOR_T>::CONVEX_POLYGON LAPLACE_CHIMERA_GRID_MPI<T_GRID>::Unclipped_Voronoi_Face(const TV& normal,const TV& center,const T radius) const
{
    return Unclipped_Voronoi_Face_Helper(normal,center,radius);
}
//#####################################################################
// Function Cartesian_Face
//#####################################################################
template<class T>
typename BASIC_GEOMETRY_POLICY<VECTOR<T,1> >::CONVEX_POLYGON Cartesian_Face_Helper(const GRID<VECTOR<T,1> >& grid,const FRAME<VECTOR<T,1> >& frame,const FACE_INDEX<1>& face_index)
{
    return typename BASIC_GEOMETRY_POLICY<VECTOR<T,1> >::CONVEX_POLYGON(frame*grid.Node(grid.Face_Node_Index(face_index.axis,face_index.index,1)),true);
}
template<class T>
typename BASIC_GEOMETRY_POLICY<VECTOR<T,2> >::CONVEX_POLYGON Cartesian_Face_Helper(const GRID<VECTOR<T,2> >& grid,const FRAME<VECTOR<T,2> >& frame,const FACE_INDEX<2>& face_index)
{
    return typename BASIC_GEOMETRY_POLICY<VECTOR<T,2> >::CONVEX_POLYGON(frame*grid.Node(grid.Face_Node_Index(face_index.axis,face_index.index,1)),frame*grid.Node(grid.Face_Node_Index(face_index.axis,face_index.index,2)));
}
template<class T>
typename BASIC_GEOMETRY_POLICY<VECTOR<T,3> >::CONVEX_POLYGON Cartesian_Face_Helper(const GRID<VECTOR<T,3> >& grid,const FRAME<VECTOR<T,3> >& frame,const FACE_INDEX<3>& face_index)
{
    typename BASIC_GEOMETRY_POLICY<VECTOR<T,3> >::CONVEX_POLYGON polygon;
    polygon.X.Resize(4);
    polygon.X(1)=frame*grid.Node(grid.Face_Node_Index(face_index.axis,face_index.index,1));
    polygon.X(2)=frame*grid.Node(grid.Face_Node_Index(face_index.axis,face_index.index,2));
    polygon.X(3)=frame*grid.Node(grid.Face_Node_Index(face_index.axis,face_index.index,4));
    polygon.X(4)=frame*grid.Node(grid.Face_Node_Index(face_index.axis,face_index.index,3));
    return polygon;
}
template<class T_GRID>
typename BASIC_GEOMETRY_POLICY<typename T_GRID::VECTOR_T>::CONVEX_POLYGON LAPLACE_CHIMERA_GRID_MPI<T_GRID>::Cartesian_Face(const int grid_index,const D_FACE_INDEX& face_index) const
{
    return Cartesian_Face_Helper(Grid(grid_index),Frame(grid_index),face_index);
}
//#####################################################################
// _Voronoi_Faces
//#####################################################################
template<class T_GRID> RANGE<typename T_GRID::VECTOR_INT> LAPLACE_CHIMERA_GRID_MPI<T_GRID>::Compute_Candidate_Neighbor_Cells(const int grid_index_1,const TV_INT& cell_index_1,const int grid_index_2,const T& search_radius) const
{
    RANGE<TV_INT> range;
    
    TV location=Rigid_Grid(grid_index_1).Frame()*Grid(grid_index_1).X(cell_index_1);
    TV location_local_2=Rigid_Grid(grid_index_2).Frame().Inverse_Times(location);
    T local_search_radius=(T)2*Grid(grid_index_2).DX().Magnitude();
    if(Grid(grid_index_2).Domain().Lazy_Inside(location_local_2)){
        range=Grid(grid_index_2).Clamp_To_Cell(RANGE<TV>(location_local_2-local_search_radius,location_local_2+local_search_radius),0);
    }else{
        range.min_corner=Grid(grid_index_2).Cell(location_local_2-search_radius);
        range.max_corner=Grid(grid_index_2).Cell(location_local_2+search_radius);
        RANGE<TV_INT> domain_indices=Grid(grid_index_2).Domain_Indices(0);
        range=RANGE<TV_INT>::Intersect(range,domain_indices);}

    return range;
}
template<class T_GRID> int LAPLACE_CHIMERA_GRID_MPI<T_GRID>::Find_Containing_Grid(const TV& location)
{
    int grid_index_containing=0;
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
        TV location_local=Rigid_Grid(grid_index).Frame().Inverse_Times(location);
        T dx=Grid(grid_index).DX().Magnitude();
        if(Grid(grid_index).Domain().Inside(location_local,dx) && (!grid_index_containing || dx<Grid(grid_index_containing).DX().Magnitude()))
            grid_index_containing=grid_index;}
    if(!grid_index_containing) return grid_index_coarsest;
    return grid_index_containing;
}
template<class T_GRID> bool LAPLACE_CHIMERA_GRID_MPI<T_GRID>::Cells_Occluded(const int grid_index_1,const TV_INT& cell_index_1,const int grid_index_2,const TV_INT& cell_index_2)
{
    if(grid_index_1==grid_index_2){
        if((cell_index_1-cell_index_2).L1_Norm()<=1)
            return false;
        for(int axis=1;axis<=TV::dimension;axis++){
            if(cell_index_1(axis)<cell_index_2(axis)){
                if(!Chimera_Cell(grid_index_1,cell_index_1+TV_INT::Axis_Vector(axis)))
                    return false;}
            else if(cell_index_1(axis)>cell_index_2(axis)){
                if(!Chimera_Cell(grid_index_1,cell_index_1-TV_INT::Axis_Vector(axis)))
                    return false;}}
        return true;}
    return false;
}
template<class T_GRID>
class INCIDENT_CELL
{
public:
    typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::HYPERPLANE HYPERPLANE;

    T distance;
    int grid_index;
    TV_INT cell_index;
    TV location;
    HYPERPLANE plane;

    INCIDENT_CELL() {}
    INCIDENT_CELL(const T distance_input,const int grid_index_input,const TV_INT& cell_index_input,const TV& location_input,const HYPERPLANE& plane_input):
        distance(distance_input),grid_index(grid_index_input),cell_index(cell_index_input),location(location_input),plane(plane_input) {}
    bool operator<(const INCIDENT_CELL<T_GRID>& other) const
    {return distance<other.distance;}
};
template<class T_GRID> void LAPLACE_CHIMERA_GRID_MPI<T_GRID>::
Construct_Voronoi_Faces()
{
    TIMER timer;
    int timer_global=timer.Register_Timer();
    int timer_intragrid_local=timer.Register_Timer();
    int timer_intragrid_nonlocal=timer.Register_Timer();
    //int timer_intragrid_clipping=timer.Register_Timer();
    int timer_intergrid_local=timer.Register_Timer();
    int timer_intergrid_nonlocal=timer.Register_Timer();
    //int timer_intergrid_clipping=timer.Register_Timer();
    int timer_clipping=timer.Register_Timer();

    int clips_active=0,clips_boundary=0;
    
    LOG::SCOPE scope("Construct_Voronoi_Faces");
    //DEBUG OUTPUT
    timer.Start(timer_global);
    
    ARRAY<INCIDENT_CELL<T_GRID> > incident_cells;
    
    voronoi_faces.Remove_All();
    voronoi_faces_geometry.Remove_All();
    cell_indices_to_incident_voronoi_face_indices.Remove_All();
    for(int grid_index_1=1;grid_index_1<=n_global_grids;grid_index_1++)
        if(Local_Grid(grid_index_1)){
            T dx_1=Grid(grid_index_1).DX().Magnitude();
            for(int boundary_cell_index=1;boundary_cell_index<=boundary_cell_indices(grid_index_1).Size();boundary_cell_index++){
                if(boundary_cell_split(grid_index_1)(boundary_cell_index)) continue;
                
                TV_INT cell_index_1=boundary_cell_indices(grid_index_1)(boundary_cell_index);
                TV location_1=Frame(grid_index_1)*Grid(grid_index_1).X(cell_index_1);
                int grid_index_containing=Find_Containing_Grid(location_1); //THIS SHOULD BE DONE GENERICALLY WITH PARAMETER TO CONTROL RADIUS
                T search_radius=2*Grid(grid_index_containing).DX().Magnitude();
                    
                incident_cells.Remove_All();
                for(int grid_index_2=1;grid_index_2<=n_global_grids;grid_index_2++){
                    RANGE<TV_INT> range_2=Compute_Candidate_Neighbor_Cells(grid_index_1,cell_index_1,grid_index_2,search_radius);
                    for(CELL_ITERATOR iterator_2(Grid(grid_index_2),range_2);iterator_2.Valid();iterator_2.Next()){
                        TV_INT cell_index_2=iterator_2.Cell_Index();
                        if((grid_index_1!=grid_index_2 || cell_index_1!=cell_index_2) && Boundary_Cell(grid_index_2,cell_index_2) && !boundary_cell_split(grid_index_2)(boundary_cell_indices_to_linear_index(grid_index_2).Get(cell_index_2)) && !Cells_Occluded(grid_index_1,cell_index_1,grid_index_2,cell_index_2)){
                            TV location_2=Rigid_Grid(grid_index_2).Frame()*iterator_2.Location();
                            TV normal=location_2-location_1;
                            T distance=normal.Magnitude();
                            if(distance<search_radius){
                                normal/=distance;
                                incident_cells.Append(INCIDENT_CELL<T_GRID>(distance,grid_index_2,cell_index_2,location_2,HYPERPLANE(normal,(T).5*(location_1+location_2))));}}}}
                Sort(incident_cells);
                
                for(int i=incident_cells.Size();i>=1;i--){
                    int grid_index_2=incident_cells(i).grid_index;
                    if(Local_Grid(grid_index_1) || Local_Grid(grid_index_2)){
                        T dx_2=Grid(grid_index_2).DX().Magnitude();
                        if(dx_1<dx_2 || (dx_1==dx_2 && grid_index_1<grid_index_2) || (grid_index_1==grid_index_2 && incident_cells(i).cell_index>cell_index_1)){
                            bool face_exists=true;
                            CONVEX_POLYGON face_simplex=Unclipped_Voronoi_Face(incident_cells(i).plane.normal,incident_cells(i).plane.x1,search_radius);
                            for(int j=1;j<=incident_cells.Size()&&face_exists;j++)
                                if(i!=j && incident_cells(j).grid_index)
                                    face_exists=INTERSECTION::Cut_Convex_Polygon_With_Hyperplane_And_Discard_Outside_Polygon(face_simplex,incident_cells(j).plane);
                            if(face_exists){
                                TV_INT& cell_index_2=incident_cells(i).cell_index;
                                for(int axis=1;axis<=TV::dimension&&face_exists;axis++)
                                    for(int side=1;side<=2&&face_exists;side++){
                                        GRID_CELL_INDEX real_grid_cell_index=chimera_grid.Find_Real_Grid_Cell_Index(GRID_CELL_INDEX(grid_index_1,cell_index_1+(2*side-3)*TV_INT::Axis_Vector(axis)));
                                        CHIMERA_STATUS status=(CHIMERA_STATUS)Chimera_Cell_Status(real_grid_cell_index.x,real_grid_cell_index.y);
                                        if(status==CHIMERA_ACTIVE || (status==CHIMERA_BOUNDARY && boundary_cell_split(real_grid_cell_index.x)(boundary_cell_indices_to_linear_index(real_grid_cell_index.x).Get(real_grid_cell_index.y)))){
                                            TV location_clipping=Rigid_Grid(real_grid_cell_index.x).Frame()*Grid(real_grid_cell_index.x).X(real_grid_cell_index.y);
                                            face_exists=INTERSECTION::Cut_Convex_Polygon_With_Hyperplane_And_Discard_Outside_Polygon(
                                                face_simplex,HYPERPLANE((location_clipping-location_1).Normalized(),(T).5*(location_clipping+location_1)));}}
                                if(face_exists){
                                    //should change this to accept all faces and then alter usage of faces to handle zero size faces
                                    T size=face_simplex.Size();
                                    face_exists=size>(T)1e-8*min(Grid(grid_index_1).Face_Sizes().Min(),Grid(grid_index_2).Face_Sizes().Min());
                                    if(face_exists){
                                        VORONOI_FACE_INDICES voronoi_face_indices(GRID_CELL_INDEX(grid_index_1,cell_index_1),GRID_CELL_INDEX(grid_index_2,cell_index_2));
                                        voronoi_faces.Append(VORONOI_FACE(voronoi_face_indices,size));
                                        voronoi_faces_geometry.Append(face_simplex);
                                        
                                        int voronoi_face_index=voronoi_faces.Size();
                                        cell_indices_to_incident_voronoi_face_indices.Get_Or_Insert(voronoi_face_indices.x).Append(voronoi_face_index);
                                        cell_indices_to_incident_voronoi_face_indices.Get_Or_Insert(voronoi_face_indices.y).Append(voronoi_face_index);}}}
                            if(!face_exists)
                                incident_cells(i).grid_index=0;}}}}} //MARK FACE AS INVALID SO THAT IT IS NOT USED TO CLIP SINCE IT IS NOT ACTUALLY INCIDENT
    
    for(int grid_index_1=1;grid_index_1<=n_global_grids;grid_index_1++)
        if(Local_Grid(grid_index_1)){
            TV face_sizes=Grid(grid_index_1).Face_Sizes();
            for(int boundary_cell_index_1=1;boundary_cell_index_1<=boundary_cell_indices(grid_index_1).Size();boundary_cell_index_1++){
                if(!boundary_cell_split(grid_index_1)(boundary_cell_index_1)) continue;
                TV_INT cell_index_1=boundary_cell_indices(grid_index_1)(boundary_cell_index_1);
                GRID_CELL_INDEX grid_cell_index_1(grid_index_1,cell_index_1);
                //LOG::cout << "adding split grid faces " << cell_index_1 << std::endl;
                    
                for(int axis=1;axis<=TV::dimension;axis++)
                    for(int side=1;side<=2;side++){
                        GRID_CELL_INDEX grid_cell_index_2=chimera_grid.Find_Real_Grid_Cell_Index(GRID_CELL_INDEX(grid_index_1,cell_index_1+(2*side-3)*TV_INT::Axis_Vector(axis)));
                        if(Boundary_Cell(grid_cell_index_2.x,grid_cell_index_2.y)){
                            int boundary_cell_index_2=boundary_cell_indices_to_linear_index(grid_cell_index_2.x).Get(grid_cell_index_2.y);
                            if(!boundary_cell_split(grid_cell_index_2.x)(boundary_cell_index_2) || grid_cell_index_1.x<grid_cell_index_2.x || (grid_cell_index_1.x==grid_cell_index_2.x && grid_cell_index_2.y>grid_cell_index_1.y)){
                            //LOG::cout << "adding split grid face " << grid_index_1 << " "<< cell_index_1 << " " << grid_cell_index_2.x << " " << grid_cell_index_2.y << std::endl;
                            VORONOI_FACE_INDICES voronoi_face_indices(grid_cell_index_1,grid_cell_index_2);
                            voronoi_faces.Append(VORONOI_FACE(voronoi_face_indices,face_sizes(axis)));
                            voronoi_faces_geometry.Append(Cartesian_Face(grid_index_1,D_FACE_INDEX(axis,cell_index_1+(side-1)*TV_INT::Axis_Vector(axis))));
                            int voronoi_face_index=voronoi_faces.Size();
                            cell_indices_to_incident_voronoi_face_indices.Get_Or_Insert(voronoi_face_indices.x).Append(voronoi_face_index);
                            cell_indices_to_incident_voronoi_face_indices.Get_Or_Insert(voronoi_face_indices.y).Append(voronoi_face_index);}}}}}
                
    //DEBUG OUTPUT
    voronoi.voronoi_faces=voronoi_faces;
    voronoi.voronoi_faces_geometry=voronoi_faces_geometry;

    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("after construct voronoi faces"),2,2);
    
    timer.Stop(timer_global);
    
    LOG::cout << "timing " << timer.Get_Accumulated_Time(timer_global) << std::endl;
    LOG::cout << "intragrid local " << timer.Get_Accumulated_Time(timer_intragrid_local) << " nonlocal " << timer.Get_Accumulated_Time(timer_intragrid_nonlocal) << std::endl;
    LOG::cout << "intergrid local " << timer.Get_Accumulated_Time(timer_intergrid_local) << " nonlocal " << timer.Get_Accumulated_Time(timer_intergrid_nonlocal) << std::endl;
    LOG::cout << " clipping " << timer.Get_Accumulated_Time(timer_clipping) << std::endl;
    LOG::cout << "counts clipping cells " << clips_boundary << " " << clips_active << std::endl;
    LOG::cout << "simplices " << voronoi_faces.Size() << std::endl;
    //Write_Substep("added voronoi faces",1,1);
    
    //EXCHANGING FACES - WE WANT TO CHANGE THIS
    
    //EXCHANGE VORONOI FACES BETWEEN BOUNDARY GRIDS IF CELLS ARE INCIDENT TO CURRENT GRID
    //REQUIRED TO BUILD INTERPOLATION SIMPLICES COVERING ALL LOCAL PRESSURE SAMPLES
    //The outer array is the array corresponds to the grids that need to receive
    if(chimera_grid.use_mpi && n_local_grids){
        ARRAY<ARRAY<VORONOI_FACE> > boundary_cell_voronoi_faces_to_send(n_global_grids),boundary_cell_voronoi_faces_to_receive(n_global_grids);
        for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
            if(Local_Grid(grid_index)){
                for(int boundary_cell_index=1;boundary_cell_index<=boundary_cell_indices(grid_index).Size();boundary_cell_index++){
                    TV_INT cell_index=boundary_cell_indices(grid_index)(boundary_cell_index);
                    if(cell_indices_to_incident_voronoi_face_indices.Contains(GRID_CELL_INDEX(grid_index,cell_index))){
                        ARRAY<int> incident_voronoi_faces=cell_indices_to_incident_voronoi_face_indices.Get(GRID_CELL_INDEX(grid_index,cell_index));
                        for(int voronoi_face_index=1;voronoi_face_index<=incident_voronoi_faces.Size();voronoi_face_index++){
                            VORONOI_FACE& voronoi_face=voronoi_faces(incident_voronoi_faces(voronoi_face_index));
                            for(int neighbor_grid_index=1;neighbor_grid_index<=n_global_grids;neighbor_grid_index++){
                                if(is_boundary_grid(1)(neighbor_grid_index)){
                                    boundary_cell_voronoi_faces_to_send(neighbor_grid_index).Append(voronoi_face);}}}}}}}
        
        chimera_grid.Exchange_Boundary_Voronoi_Faces(boundary_cell_voronoi_faces_to_send,boundary_cell_voronoi_faces_to_receive);
        
        HASHTABLE<VORONOI_FACE> added_boundary_voronoi_faces;
        for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
            for(int boundary_voronoi_face_index=1;boundary_voronoi_face_index<=boundary_cell_voronoi_faces_to_receive(grid_index).Size();boundary_voronoi_face_index++){
                VORONOI_FACE boundary_voronoi_face=boundary_cell_voronoi_faces_to_receive(grid_index)(boundary_voronoi_face_index);
                if(!added_boundary_voronoi_faces.Contains(boundary_voronoi_face)){
                    voronoi_faces.Append(boundary_voronoi_face);
                    added_boundary_voronoi_faces.Set(boundary_voronoi_face);
                    int voronoi_face_index=voronoi_faces.Size();
                    cell_indices_to_incident_voronoi_face_indices.Get_Or_Insert(boundary_voronoi_face.x(1)).Append(voronoi_face_index);
                    cell_indices_to_incident_voronoi_face_indices.Get_Or_Insert(boundary_voronoi_face.x(2)).Append(voronoi_face_index);}}}

    //LOG::cout<<"lqiu debug number of voronoi faces after exchanging voronoi faces is "<<voronoi_faces.Size()<<std::endl;
}
//#####################################################################
// Function Construct_Dual_Cell_And_Cell_Sizes
//#####################################################################
/*template<class T_GRID> void LAPLACE_CHIMERA_GRID_MPI<T_GRID>::Construct_Dual_Cell_And_Cell_Sizes()
{
    LOG::SCOPE scope("Compute_Dual_Cell_And_Cell_Sizes");

    dual_cell_sizes.Resize(n_matrix_faces);
    cell_sizes.Resize(n_matrix_cells);
    ARRAYS_COMPUTATIONS::Fill(cell_sizes,(T)0);
    
    for(int voronoi_face_index=1;voronoi_face_index<=voronoi_faces.Size();voronoi_face_index++){
        int matrix_face_index=voronoi_face_indices_to_matrix_face_indices(voronoi_face_index);
        if(matrix_face_index){
            T face_size=laplace_grid.voronoi_faces(voronoi_face_index).y;
            GRID_CELL_INDEX& grid_cell_index_1=laplace_grid.voronoi_faces(voronoi_face_index).x(1);
            GRID_CELL_INDEX& grid_cell_index_2=laplace_grid.voronoi_faces(voronoi_face_index).x(2);
            TV location_1=laplace_grid.Rigid_Grid(grid_cell_index_1.x).Frame()*laplace_grid.Grid(grid_cell_index_1.x).X(grid_cell_index_1.y);
            TV location_2=laplace_grid.Rigid_Grid(grid_cell_index_2.x).Frame()*laplace_grid.Grid(grid_cell_index_2.x).X(grid_cell_index_2.y);
            T dx=(location_2-location_1).Magnitude();
            dual_cell_sizes(matrix_face_index)=dx*face_size;}}
    
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index))
            for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                int matrix_cell_index=Matrix_Cell_Index(grid_index,cell_index);
                if(matrix_cell_index){
                    GRID_CELL_INDEX grid_cell_index(grid_index,cell_index);
                    T cell_size=0;
                    if(cell_indices_to_incident_voronoi_face_indices.Contains(grid_cell_index)){
                        ARRAY<int>& incident_voronoi_face_indices=cell_indices_to_incident_voronoi_face_indices.Get(grid_cell_index);
                        for(int i=1;i<=incident_voronoi_face_indices.Size();++i){
                            int matrix_face_index=Matrix_Face_Index(incident_voronoi_face_indices(i));
                            if(matrix_face_index) cell_size+=(T)0.5*dual_cell_sizes(matrix_face_index)/TV::dimension;}}
                    for(int axis=1;axis<=TV::dimension;axis++)
                        for(int side=1;side<=2;side++){
                            int matrix_face_index=Matrix_Face_Index(grid_index,D_FACE_INDEX(axis,cell_index+(side-1)*TV_INT::Axis_Vector(axis)));
                            if(matrix_face_index)
                                cell_size+=(T)0.5*dual_cell_sizes(matrix_face_index)/TV::dimension;}
                    cell_sizes(matrix_cell_index)=cell_size;}}
}*/
//#####################################################################
// Function Face_Size
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LAPLACE_CHIMERA_GRID_MPI<T_GRID>::Face_Size(const int grid_index,const D_FACE_INDEX& face_index) const
{
    if(Chimera_Face(grid_index,face_index))
        return Grid(grid_index).Cell_Size();
    return 0;
}
template<class T_GRID> typename T_GRID::SCALAR LAPLACE_CHIMERA_GRID_MPI<T_GRID>::Face_Size(const int voronoi_face_index) const
{
    const GRID_CELL_INDEX& grid_cell_index_1=voronoi_faces(voronoi_face_index).x(1);
    const GRID_CELL_INDEX& grid_cell_index_2=voronoi_faces(voronoi_face_index).x(2);
    TV location_1=Rigid_Grid(grid_cell_index_1.x).Frame()*Grid(grid_cell_index_1.x).X(grid_cell_index_1.y);
    TV location_2=Rigid_Grid(grid_cell_index_2.x).Frame()*Grid(grid_cell_index_2.x).X(grid_cell_index_2.y);
    T dx=(location_2-location_1).Magnitude();
    T size=voronoi_faces(voronoi_face_index).y;
    return dx*size;
}
//#####################################################################
// Function Cell_Size
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR LAPLACE_CHIMERA_GRID_MPI<T_GRID>::Cell_Size(const int grid_index,const TV_INT& cell_index) const
{
    if(Chimera_Cell(grid_index,cell_index)){
        if(Boundary_Cell(grid_index,cell_index)){
            T size=0;
            for(int axis=1;axis<=TV::dimension;axis++)
                for(int side=1;side<=2;side++)
                    size+=(T)0.5/TV::dimension*Face_Size(grid_index,D_FACE_INDEX(axis,cell_index+(side-1)*TV_INT::Axis_Vector(axis)));
            GRID_CELL_INDEX grid_cell_index(grid_index,cell_index);
            if(cell_indices_to_incident_voronoi_face_indices.Contains(grid_cell_index)){
                const ARRAY<int>& incident_voronoi_face_indices=cell_indices_to_incident_voronoi_face_indices.Get(grid_cell_index);
                for(int i=1;i<=incident_voronoi_face_indices.Size();++i)
                    size+=(T)0.5/TV::dimension*Face_Size(incident_voronoi_face_indices(i));}
            return size;
        }
        else return Grid(grid_index).Cell_Size();}
    return 0;
}

//#####################################################################
// Function Construct_Mesh
//#####################################################################
template<class T_GRID> void LAPLACE_CHIMERA_GRID_MPI<T_GRID>::Construct_Mesh()
{
    n_global_grids=chimera_grid.number_of_global_grids;
    n_local_grids=chimera_grid.number_of_local_grids;
    
    grid_index_coarsest=1;
    for(int grid_index=2;grid_index<=n_global_grids;grid_index++)
        if(Grid(grid_index).DX().Magnitude()>Grid(grid_index_coarsest).DX().Magnitude())
            grid_index_coarsest=grid_index;
    
    Construct_Chimera_Cells();
    Construct_Boundary_Cells();
    Construct_Voronoi_Faces();
    {LOG::SCOPE scope("Building Delaunay Mesh");
    interpolation.Construct_Delaunay_Simplices();
    }
    //Construct_Dual_Cell_And_Cell_Sizes();
}
//#####################################################################

#define INSTANTIATION_HELPER(T,D)                               \
    template class LAPLACE_CHIMERA_GRID_MPI<GRID<VECTOR<T,D> > >;

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif
