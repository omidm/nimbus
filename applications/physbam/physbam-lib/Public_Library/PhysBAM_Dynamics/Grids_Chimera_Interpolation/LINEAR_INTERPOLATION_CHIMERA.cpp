//#####################################################################
// Copyright 2011
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INTERPOLATION_CHIMERA
//#####################################################################

#include <PhysBAM_Dynamics/Grids_Chimera/Parallel_Computation/CHIMERA_GRID.h>
#include <PhysBAM_Dynamics/Grids_Chimera/Grids_Chimera_Boundaries/BOUNDARY_CHIMERA.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_PACKAGE.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Tools/Parallel_Computation/PCG_SPARSE_UNSTRUCTURE_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Dynamics/Grids_Chimera_Interpolation/LINEAR_INTERPOLATION_CHIMERA.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Dynamics/Grids_Chimera_Interpolation/CELL_LOOKUP_CHIMERA.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology/TETRAHEDRON_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY_2D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Dynamics/Grids_Chimera_PDE_Linear/Parallel_Computation/LAPLACE_CHIMERA_GRID_MPI.h>

using namespace PhysBAM;

template<int d>
class COMBINATION_ITERATOR
{
private:
    int n;
    VECTOR<int,d> index;

public:
    COMBINATION_ITERATOR(int n_input):n(n_input){
        assert(n>=d);
        for(int i=1;i<=d;i++)
            index(i)=i;}

    bool Valid(){return index.Max()<=n;}
    void Next(){
        for(int i=d;i>=1;i--)
            if(index(i)<(n+i-d)){
                index(i)++;
                for(int j=i+1;j<=d;j++)
                    index(j)=index(j-1)+1;
                return;
            }
        index(1)=n+1;
    }
    VECTOR<int,d> Index(){return index;}
};
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> LINEAR_INTERPOLATION_CHIMERA<T_GRID>::
LINEAR_INTERPOLATION_CHIMERA(T_LAPLACE_GRID& laplace_grid_input)
    :laplace_grid(laplace_grid_input)
{}
//#####################################################################
// Construct_Delaunay_Simplices
//#####################################################################
template<class T_GRID> bool LINEAR_INTERPOLATION_CHIMERA<T_GRID>::
Cell_Intersect_Local_Domain(GRID_CELL_INDEX cell)
{
    RANGE<TV> cell_domain=laplace_grid.Grid(cell.x).Cell_Domain(cell.y);
    T dx_coarsest=laplace_grid.Grid(laplace_grid.grid_index_coarsest).DX().Magnitude();

    ORIENTED_BOX<TV> cell_box(cell_domain.Thickened(dx_coarsest),laplace_grid.Rigid_Grid(cell.x).Frame());
    for(int grid_index=1;grid_index<=laplace_grid.n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index)){
            ORIENTED_BOX<TV> grid_box(laplace_grid.Grid(grid_index).Domain().Thickened(dx_coarsest),laplace_grid.Rigid_Grid(grid_index).Frame());
            if(grid_box.Intersection(cell_box)) return true;}
    return false;
}
template<class T_GRID> void LINEAR_INTERPOLATION_CHIMERA<T_GRID>::
Construct_Delaunay_Simplices()
{
    //laplace_grid.voronoi.boundary_delaunay_simplex_vertices.Remove_All();
    HASHTABLE<VECTOR<GRID_CELL_INDEX,TV::dimension+1> > delaunay_simplex_hashtable;
    LOG::cout<<"building delaunay simplices"<<std::endl;

    for(typename HASHTABLE<GRID_CELL_INDEX,ARRAY<int> >::ITERATOR iterator(laplace_grid.cell_indices_to_incident_voronoi_face_indices);iterator.Valid();iterator.Next()){
        const int grid_index=iterator.Key().x;
        const TV_INT& cell_index=iterator.Key().y;
        
        if(!(laplace_grid.Local_Grid(grid_index) || laplace_grid.Boundary_Grid(grid_index)) || laplace_grid.boundary_cell_split(grid_index)(laplace_grid.boundary_cell_indices_to_linear_index(grid_index).Get(cell_index)) || !Cell_Intersect_Local_Domain(GRID_CELL_INDEX(grid_index,cell_index))) continue;

        TV location=laplace_grid.Rigid_Grid(grid_index).Frame()*laplace_grid.Grid(grid_index).X(cell_index);
        GRID_CELL_INDEX grid_cell_index(grid_index,cell_index);
        if(laplace_grid.cell_indices_to_incident_voronoi_face_indices.Contains(grid_cell_index)){
            ARRAY<GRID_CELL_INDEX> neighbor_cells;
            ARRAY<int>& incident_voronoi_face_indices=laplace_grid.cell_indices_to_incident_voronoi_face_indices.Get(grid_cell_index);
            for(int voronoi_face_index=1;voronoi_face_index<=incident_voronoi_face_indices.Size();voronoi_face_index++){
                const VECTOR<GRID_CELL_INDEX,2>& grid_cell_indices=laplace_grid.voronoi_faces(incident_voronoi_face_indices(voronoi_face_index)).x;
                const GRID_CELL_INDEX& neighbor=grid_cell_indices(1)==grid_cell_index?grid_cell_indices(2):grid_cell_indices(1);

                if(!(laplace_grid.Local_Grid(neighbor.x) || laplace_grid.Boundary_Grid(neighbor.x)) || laplace_grid.boundary_cell_split(neighbor.x)(laplace_grid.boundary_cell_indices_to_linear_index(neighbor.x).Get(neighbor.y)))
                    continue;
                neighbor_cells.Append(neighbor);}
            ARRAY<TV> neighbor_points;
            for(int p=1;p<=neighbor_cells.Size();p++)
                neighbor_points.Append(laplace_grid.Rigid_Grid(neighbor_cells(p).x).Frame()*laplace_grid.Grid(neighbor_cells(p).x).X(neighbor_cells(p).y));

            if(neighbor_cells.Size()<TV::dimension) continue;

            for(COMBINATION_ITERATOR<TV::dimension> iterator(neighbor_points.Size());iterator.Valid();iterator.Next()){
                VECTOR<int,TV::dimension> indices=iterator.Index();
                T_SIMPLEX simplex(VECTOR<TV,TV::dimension>(neighbor_points.Subset(indices)).Append(location));
                if(simplex.Minimum_Altitude()>((T)1e-4*simplex.Maximum_Edge_Length())){ //this should be minimum altitude or aspect ratio
                    bool delaunay=true;
                    TV circumcenter=simplex.Circumcenter();
                    T radius=(circumcenter-location).Magnitude();
                    for(int i=1;i<=neighbor_cells.Size();++i)
                        if(!indices.Contains(i) && (neighbor_points(i)-circumcenter).Magnitude()<((T)0.9999*radius)){
                            delaunay=false;break;}
                    if(delaunay){
                        VECTOR<GRID_CELL_INDEX,TV::dimension+1> grid_cell_indices=VECTOR<GRID_CELL_INDEX,TV::dimension>(neighbor_cells.Subset(indices)).Append(grid_cell_index).Sorted();
                        bool valid_stencil=true;
                        RANGE<TV_INT> cell_centers_cube(cell_index);
                        for(int d=1;d<=TV::dimension+1;d++){
                            if(grid_cell_indices(d).x!=grid_index){
                                valid_stencil=false;break;}
                            cell_centers_cube.Enlarge_To_Include_Point(grid_cell_indices(d).y);}
                        if(valid_stencil)
                            for(CELL_ITERATOR stencil_iterator(laplace_grid.Grid(grid_index),cell_centers_cube);stencil_iterator.Valid();stencil_iterator.Next())
                                    if(!laplace_grid.Chimera_Cell(grid_index,stencil_iterator.Cell_Index())){
                                        valid_stencil=false;break;}
                        if(!valid_stencil)
                            delaunay_simplex_hashtable.Set(grid_cell_indices);}}}}}
    
    LOG::cout << "building interpolation acceleration structures" << std::endl;
    simplex_indices.Resize(delaunay_simplex_hashtable.Size());
    int c=0;
    for(typename HASHTABLE<VECTOR<GRID_CELL_INDEX,TV::dimension+1> >::ITERATOR iterator(delaunay_simplex_hashtable);iterator.Valid();iterator.Next())
        simplex_indices(++c)=iterator.Key();
    
    ARRAY<RANGE<TV> > simplex_ranges(simplex_indices.Size());
    simplex_vertices.Resize(simplex_indices.Size());
    for(int i=1;i<=simplex_indices.Size();i++){
        RANGE<TV> range;
        for(int j=1;j<=TV::dimension+1;j++){
            simplex_vertices(i)(j)=laplace_grid.Frame(simplex_indices(i)(j).x)*laplace_grid.Grid(simplex_indices(i)(j).x).X(simplex_indices(i)(j).y);
            simplex_ranges(i).Enlarge_To_Include_Point(simplex_vertices(i)(j));}
        if(T_SIMPLEX::Signed_Size(simplex_vertices(i))<0){
            exchange(simplex_indices(i)(1),simplex_indices(i)(2));
            exchange(simplex_vertices(i)(1),simplex_vertices(i)(2));}}
    
    hierarchy.Clean_Memory();
    if(simplex_ranges.Size()) hierarchy.Set_Leaf_Boxes(simplex_ranges,true);
    LOG::cout<<"total number of boundary_delaunay_simplices: "<<simplex_indices.Size()<<std::endl;
    
    laplace_grid.voronoi.boundary_delaunay_simplex_vertices=simplex_vertices;
    //DEBUG_UTILITIES::Debug_Breakpoint();
}
//#####################################################################
// Function 
//#####################################################################
template<class T_GRID> typename LINEAR_INTERPOLATION_CHIMERA<T_GRID>::T LINEAR_INTERPOLATION_CHIMERA<T_GRID>::
From_Cell_Centers(const CELL_LOOKUP_CHIMERA<T_GRID>& u,const TV& X/*,bool* interpolated*/) const
{
    if(simplex_indices.Size()){
        ARRAY<int> candidate_simplices;
        T thickness=(T)1e-8*laplace_grid.Grid(laplace_grid.grid_index_coarsest).Minimum_Edge_Length();
        //T thickness=0;
        hierarchy.Intersection_List(X,candidate_simplices,thickness);

        T largest_area=0;
        int largest_simplex_index=0;

        for(int i=1;i<=candidate_simplices.Size();i++){
            int simplex_index=candidate_simplices(i);
            T_SIMPLEX simplex(simplex_vertices(simplex_index));
            if(!simplex.Outside(X,(T)1e-5*simplex.Minimum_Edge_Length())){
                T area=simplex.Size();
                if(area>largest_area){
                    largest_area=area;
                    largest_simplex_index=simplex_index;}}}
        
        if(largest_simplex_index){
            //if(!simplex.Outside(X,(T)1e-8*simplex.Minimum_Edge_Length())){//-thickness)){
            T_SIMPLEX simplex(simplex_vertices(largest_simplex_index));
            VECTOR<T,TV::dimension+1> bary_coordinates=simplex.Barycentric_Coordinates(X);
            T interpolated_value(0);
            for(int p=1;p<=TV::dimension+1;++p)
                interpolated_value+=bary_coordinates(p)*u(simplex_indices(largest_simplex_index)(p));
            return interpolated_value;}}
    
    //NEED TO SUPPORT SPLIT GRIDS IN THIS CASE
    //we do this by handling the overlapped case using the grid coupling code
    /*LINEAR_INTERPOLATION_UNIFORM<T_GRID,T> interpolation;
    int grid_index_containing_finest=0;
    for(int grid_index_containing=1;grid_index_containing<=laplace_grid.n_global_grids;grid_index_containing++)
        if(laplace_grid.Local_Grid(grid_index_containing)){
            TV object_location=laplace_grid.Rigid_Grid(grid_index_containing).Frame().Inverse_Times(X);
            const T_GRID& grid_containing=laplace_grid.Grid(grid_index_containing);
            if(grid_containing.domain.Inside(object_location,(T)0.5*grid_containing.dX.Max())) //THIS SHOULD BE ABSTRACTIFIED TO SUPPORT SPLIT GRIDS, SOME CALLBACK TO DETERMINE WHETHER A POINT IS IN THE INTERPOLATION DOMAIN OF A GRID
                if(!grid_index_containing_finest || grid_containing.dX.Product()<laplace_grid.Grid(grid_index_containing_finest).dX.Product())
                    grid_index_containing_finest=grid_index_containing;
        } 
    if(grid_index_containing_finest)
        return interpolation.Clamped_To_Array_Cell(laplace_grid.Grid(grid_index_containing_finest),u(grid_index_containing_finest),laplace_grid.Rigid_Grid(grid_index_containing_finest).Frame().Inverse_Times(X));*/
    
    return std::numeric_limits<T>::quiet_NaN();
}
//#####################################################################
template class LINEAR_INTERPOLATION_CHIMERA<GRID<VECTOR<float,1> > >;
template class LINEAR_INTERPOLATION_CHIMERA<GRID<VECTOR<float,2> > >;
template class LINEAR_INTERPOLATION_CHIMERA<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LINEAR_INTERPOLATION_CHIMERA<GRID<VECTOR<double,1> > >;
template class LINEAR_INTERPOLATION_CHIMERA<GRID<VECTOR<double,2> > >;
template class LINEAR_INTERPOLATION_CHIMERA<GRID<VECTOR<double,3> > >;
#endif
