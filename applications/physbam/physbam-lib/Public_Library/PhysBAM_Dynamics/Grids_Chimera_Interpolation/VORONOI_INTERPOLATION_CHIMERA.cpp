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
#include <PhysBAM_Dynamics/Grids_Chimera_Interpolation/VORONOI_INTERPOLATION_CHIMERA.h>
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
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>

using namespace PhysBAM;

//#####################################################################
// Compute_Least_Squares_Velocity
//#####################################################################
template<class T_GRID> typename T_GRID::VECTOR_T VORONOI_INTERPOLATION_CHIMERA<T_GRID>::Compute_Least_Squares_Vector(const TV& x,const T search_radius,const ARRAY<T_FACE_ARRAYS_SCALAR*>& face_values,const ARRAY<T>& voronoi_face_values) const
{
    int n_global_grids=laplace_grid.n_global_grids;

    HASHTABLE<GRID_CELL_INDEX> incident_cells;
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++){
        TV x_local=laplace_grid.Frame(grid_index).Inverse_Times(x);
        for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index),laplace_grid.Grid(grid_index).Clamp_To_Cell(RANGE<TV>(x_local-search_radius,x_local+search_radius),0));iterator.Valid();iterator.Next())
            if(laplace_grid.Chimera_Cell(grid_index,iterator.Cell_Index()) && (x_local-iterator.Location()).Magnitude()<=search_radius)
                incident_cells.Set(GRID_CELL_INDEX(grid_index,iterator.Cell_Index()));}
    
    D_MATRIX A;
    TV rhs;
    for(typename HASHTABLE<GRID_CELL_INDEX>::ITERATOR iterator(incident_cells);iterator.Valid();iterator.Next()){
        if(laplace_grid.cell_indices_to_incident_voronoi_face_indices.Contains(iterator.Key())){
            ARRAY<int> incident=laplace_grid.cell_indices_to_incident_voronoi_face_indices.Get(iterator.Key());
            for(int k=1;k<=incident.Size();k++){
                VECTOR<GRID_CELL_INDEX,2>& incident_grid_cell_indices=laplace_grid.voronoi_faces(incident(k)).x;
                if(incident_grid_cell_indices(1)==iterator.Key() && incident_cells.Contains(incident_grid_cell_indices(2))){
                    TV location_1,location_2,location,normal;T distance;
                    laplace_grid.Voronoi_Face(incident(k),location_1,location_2,location,normal,distance);
                    T size=1;//laplace_grid.voronoi_faces(incident(k)).y;
                    A+=size*D_MATRIX::Outer_Product(normal,normal);
                    rhs+=size*normal*voronoi_face_values(incident(k));}}}
        for(int axis=1;axis<=TV::dimension;axis++){
            D_FACE_INDEX face_index(axis,iterator.Key().y);
            if(laplace_grid.Chimera_Face(iterator.Key().x,face_index) && incident_cells.Contains(GRID_CELL_INDEX(iterator.Key().x,iterator.Key().y-TV_INT::Axis_Vector(axis)))){
                TV normal=laplace_grid.Frame(iterator.Key().x).r.Rotate(TV::Axis_Vector(axis));
                T size=1;//laplace_grid.Grid(iterator.Key().x).Face_Size(axis);
                A+=size*D_MATRIX::Outer_Product(normal,normal);
                rhs+=size*normal*face_values(iterator.Key().x)->Component(axis)(iterator.Key().y);}}}
    return A.Solve_Linear_System(rhs);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> int Vertices(const POLYGON<TV>& e)
{return e.X.Size();}
template<class T> int Vertices(const SEGMENT_2D<T>& e)
{return 2;}
template<class T> int Vertices(const POINT_SIMPLEX_1D<T>& e)
{return 1;}
template<class T_GRID> VORONOI_INTERPOLATION_CHIMERA<T_GRID>::
VORONOI_INTERPOLATION_CHIMERA(T_LAPLACE_GRID& laplace_grid_input,const ARRAY<T_FACE_ARRAYS_SCALAR*>& face_values_input,const ARRAY<T>& voronoi_face_values_input)
    :laplace_grid(laplace_grid_input),face_values(face_values_input),voronoi_face_values(voronoi_face_values_input)
{
    int n_global_grids=laplace_grid.n_global_grids;
    int n_voronoi_faces=laplace_grid.voronoi_faces.Size();
    
    voronoi_face_centroids.Resize(n_voronoi_faces);
    voronoi_face_centroid_values.Resize(n_voronoi_faces);
    voronoi_face_vertex_values.Resize(n_voronoi_faces);

    laplace_grid.voronoi.vectors.Remove_All();

    T threshold_factor=(T)1.001;

    //compute voronoi face velocities
    for(int voronoi_face_index=1;voronoi_face_index<=n_voronoi_faces;voronoi_face_index++){
        int n_vertices=Vertices(laplace_grid.voronoi_faces_geometry(voronoi_face_index));
        voronoi_face_vertex_values(voronoi_face_index).Resize(n_vertices);
        //voronoi_face_centroids(voronoi_face_index)=laplace_grid.voronoi_faces_geometry(voronoi_face_index).Center();
        
        TV centroid;
        TV centroid_value;

        TV location_1,location_2,location,normal;T distance;
        laplace_grid.Voronoi_Face(voronoi_face_index,location_1,location_2,location,normal,distance);
        
        for(int vertex_index=1;vertex_index<=n_vertices;vertex_index++){
            TV& x=laplace_grid.voronoi_faces_geometry(voronoi_face_index).X(vertex_index);
            T search_radius=threshold_factor*(x-location_1).Magnitude();
            voronoi_face_vertex_values(voronoi_face_index)(vertex_index)=Compute_Least_Squares_Vector(x,search_radius,face_values,voronoi_face_values);
            laplace_grid.voronoi.vectors.Append(VECTOR<TV,2>(x,voronoi_face_vertex_values(voronoi_face_index)(vertex_index)));///////////DEBUG
            centroid_value+=voronoi_face_vertex_values(voronoi_face_index)(vertex_index);
            centroid+=x;}

        centroid_value/=(T)n_vertices;
        centroid_value+=(voronoi_face_values(voronoi_face_index)-TV::Dot_Product(normal,centroid_value))*normal;
        centroid/=(T)n_vertices;
        voronoi_face_centroid_values(voronoi_face_index)=centroid_value;
        voronoi_face_centroids(voronoi_face_index)=centroid;
        laplace_grid.voronoi.vectors.Append(VECTOR<TV,2>(centroid,centroid_value));/////////////DEBUG
    }
    
    //compute cartesian face velocities
    cell_centroids.Resize(n_global_grids);
    cell_centroid_values.Resize(n_global_grids);
    face_centroid_values.Resize(n_global_grids);
    node_values.Resize(n_global_grids);
    for(int grid_index=1;grid_index<=n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index)){
            T_GRID& grid=laplace_grid.Grid(grid_index);
            FRAME<TV> frame=laplace_grid.Frame(grid_index);
            int n_boundary_cells=laplace_grid.boundary_cell_indices(grid_index).Size();
            cell_centroid_values(grid_index).Resize(n_boundary_cells);
            cell_centroids(grid_index).Resize(n_boundary_cells);
            
            for(int boundary_cell_index=1;boundary_cell_index<=n_boundary_cells;boundary_cell_index++){
                TV_INT cell_index=laplace_grid.boundary_cell_indices(grid_index)(boundary_cell_index);
                TV cell_location=frame*grid.X(cell_index);

                TV centroid;
                TV centroid_value;

                T face_weight=0;
                
                for(int axis=1;axis<=TV::dimension;axis++)
                    for(int side=1;side<=2;side++){
                        D_FACE_INDEX face_index(axis,cell_index+(side-1)*TV_INT::Axis_Vector(axis));
                        if(laplace_grid.Chimera_Face(grid_index,face_index)){
                            TV face_value;
                            TV normal=frame.r.Rotate(TV::Axis_Vector(axis));
                            for(int i=1;i<=T_GRID::number_of_nodes_per_face;i++){
                                TV_INT node_index=grid.Face_Node_Index(face_index.axis,face_index.index,i);
                                if(!node_values(grid_index).Contains(node_index)){
                                    TV node_location=frame*grid.Node(node_index);
                                    T search_radius=threshold_factor*(node_location-cell_location).Magnitude();
                                    TV node_value=Compute_Least_Squares_Vector(node_location,search_radius,face_values,voronoi_face_values);
                                    node_values(grid_index).Set(node_index,node_value);
                                    laplace_grid.voronoi.vectors.Append(VECTOR<TV,2>(node_location,node_value));///////////DEBUG
                                }
                                face_value+=node_values(grid_index).Get(node_index);}
                            
                            face_value/=T_GRID::number_of_nodes_per_face;
                            face_value+=(face_values(grid_index)->operator()(face_index)-TV::Dot_Product(normal,face_value))*normal;
                            face_centroid_values(grid_index).Set(face_index,face_value);

                            TV face_centroid=frame*grid.Face(face_index.axis,face_index.index);
                            T size=grid.Face_Size(axis);
                            centroid+=size*face_centroid;
                            centroid_value+=size*face_value;
                            face_weight+=size;

                            laplace_grid.voronoi.vectors.Append(VECTOR<TV,2>(face_centroid,face_value));///////////DEBUG
                        }}

                GRID_CELL_INDEX grid_cell_index(grid_index,cell_index);
                if(laplace_grid.cell_indices_to_incident_voronoi_face_indices.Contains(grid_cell_index)){
                    const ARRAY<int>& incident=laplace_grid.cell_indices_to_incident_voronoi_face_indices.Get(grid_cell_index);
                    for(int i=1;i<=incident.Size();i++){
                        T size=laplace_grid.voronoi_faces(incident(i)).y;
                        centroid+=size*voronoi_face_centroids(incident(i));
                        centroid_value+=size*voronoi_face_centroid_values(incident(i));
                        face_weight+=size;
                    }}

                centroid/=face_weight;
                centroid_value/=face_weight;
                cell_centroids(grid_index)(boundary_cell_index)=centroid;
                cell_centroid_values(grid_index)(boundary_cell_index)=centroid_value;
                
                /*D_MATRIX NNt;
                TV rhs;
                for(int axis=1;axis<=TV::dimension;axis++)
                    for(int side=1;side<=2;side++){
                        D_FACE_INDEX face_index(axis,cell_index+(side-1)*TV_INT::Axis_Vector(axis));
                        if(laplace_grid.Chimera_Face(grid_index,face_index)){
                            TV normal=laplace_grid.Frame(grid_index).r.Rotate(TV::Axis_Vector(axis));
                            NNt+=D_MATRIX::Outer_Product(normal,normal);
                            rhs+=normal*face_values(grid_index)->Component(face_index.axis)(face_index.index);}}
                GRID_CELL_INDEX grid_cell_index(grid_index,cell_index);
                if(laplace_grid.cell_indices_to_incident_voronoi_face_indices.Contains(grid_cell_index)){
                    const ARRAY<int>& incident=laplace_grid.cell_indices_to_incident_voronoi_face_indices.Get(grid_cell_index);
                    for(int i=1;i<=incident.Size();i++){
                        TV location_1,location_2,location,normal;T distance;
                        laplace_grid.Voronoi_Face(incident(i),location_1,location_2,location,normal,distance);
                        NNt+=D_MATRIX::Outer_Product(normal,normal);
                        rhs+=normal*voronoi_face_values(incident(i));}}
                
                        cell_center_values(grid_index)(boundary_cell_index)=NNt.Solve_Linear_System(rhs);*/
                
                laplace_grid.voronoi.vectors.Append(VECTOR<TV,2>(centroid,centroid_value));///////////DEBUG
            }}
}
//#####################################################################
// Function Nearest_Cell
//#####################################################################
template<class T_GRID> PAIR<int,typename T_GRID::VECTOR_INT> VORONOI_INTERPOLATION_CHIMERA<T_GRID>::Nearest_Cell(const TV& x) const
{
    T tolerance=1;//.001;
    GRID_CELL_INDEX nearest_cell;
    T nearest_distance=FLT_MAX;
    for(int grid_index=1;grid_index<=laplace_grid.n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index)){
            T_GRID& grid=laplace_grid.Grid(grid_index);
            TV x_local=laplace_grid.Frame(grid_index).Inverse_Times(x);
            TV x_clamped=grid.Clamp(x_local);
            T search_radius=tolerance*grid.DX().Magnitude();
            for(CELL_ITERATOR iterator(grid,grid.Clamp_To_Cell(RANGE<TV>(x_clamped-search_radius,x_clamped+search_radius),0));iterator.Valid();iterator.Next())
                if(laplace_grid.Chimera_Cell(grid_index,iterator.Cell_Index())){
                    T distance=(x_local-iterator.Location()).Magnitude();
                    if(distance<nearest_distance){
                        nearest_cell=GRID_CELL_INDEX(grid_index,iterator.Cell_Index());
                        nearest_distance=distance;}}}

    if(laplace_grid.Boundary_Cell(nearest_cell.x,nearest_cell.y))
        return nearest_cell;
    else
        return GRID_CELL_INDEX(0,TV_INT());
}
//#####################################################################
// Function From_Voronoi_Faces
//#####################################################################
#define WRAP(I,N) ((I-1)%N+1)
template<class T_GRID> typename T_GRID::VECTOR_T VORONOI_INTERPOLATION_CHIMERA<T_GRID>::From_Voronoi_Faces(const TV& x) const
{
    T thickness_factor=(T)1e-4;

    GRID_CELL_INDEX nearest_cell=Nearest_Cell(x);
    
    //if(x(2)>11.95 && x(2)<12.05)
    //    LOG::cout << "interpolating " << x << " " << nearest_cell.x << " " << nearest_cell.y << std::endl;

    if(nearest_cell.x){
        T_GRID& grid=laplace_grid.Grid(nearest_cell.x);
        FRAME<TV> frame=laplace_grid.Frame(nearest_cell.x);
        int boundary_cell_index=laplace_grid.boundary_cell_indices_to_linear_index(nearest_cell.x).Get(nearest_cell.y);
        
        TV cell_centroid=cell_centroids(nearest_cell.x)(boundary_cell_index);//frame*grid.X(nearest_cell.y);
        VECTOR<TV,TV::dimension+1> simplex_vertices;
        if(laplace_grid.cell_indices_to_incident_voronoi_face_indices.Contains(nearest_cell)){
            ARRAY<int> incident_indices=laplace_grid.cell_indices_to_incident_voronoi_face_indices.Get(nearest_cell);
            for(int i=1;i<=incident_indices.Size();i++){
                int voronoi_face_index=incident_indices(i);
                int n_vertices=Vertices(laplace_grid.voronoi_faces_geometry(voronoi_face_index));
                for(int vertex_index=1;vertex_index<=n_vertices;vertex_index++){
                    simplex_vertices(1)=cell_centroid;
                    simplex_vertices(2)=voronoi_face_centroids(incident_indices(i));
                    for(int j=0;j<(TV::dimension-1);j++)
                        simplex_vertices(j+3)=laplace_grid.voronoi_faces_geometry(voronoi_face_index).X(WRAP(vertex_index+j,n_vertices));
                    T_SIMPLEX simplex(simplex_vertices);
                    bool swapped=false;
                    if(simplex.Signed_Size()<0){
                        exchange(simplex_vertices(1),simplex_vertices(2));
                        simplex=T_SIMPLEX(simplex_vertices);
                        swapped=true;}
                    if(!simplex.Outside(x,thickness_factor*simplex.Minimum_Edge_Length())){
                        VECTOR<T,TV::dimension+1> bary_coordinates=simplex.Barycentric_Coordinates(x);
                        TV vector;
                        if(swapped) vector=bary_coordinates(2)*cell_centroid_values(nearest_cell.x)(boundary_cell_index)+bary_coordinates(1)*voronoi_face_centroid_values(incident_indices(i));
                        else vector=bary_coordinates(1)*cell_centroid_values(nearest_cell.x)(boundary_cell_index)+bary_coordinates(2)*voronoi_face_centroid_values(incident_indices(i));

                        //if(x(2)>11.95 && x(2)<12.05)
                        //    LOG::cout << "interpolating from voronoi face " << x << " " << simplex_vertices << " " << bary_coordinates << std::endl;

                        for(int j=0;j<(TV::dimension-1);j++)
                            vector+=bary_coordinates(j+3)*voronoi_face_vertex_values(voronoi_face_index)(WRAP(vertex_index+j,n_vertices));
                        return vector;}}}}

        for(int axis=1;axis<=TV::dimension;axis++)
            for(int side=1;side<=2;side++){
                D_FACE_INDEX face_index(axis,nearest_cell.y+(side-1)*TV_INT::Axis_Vector(axis));
                if(laplace_grid.Chimera_Face(nearest_cell.x,face_index)){
                    for(int vertex_index=1;vertex_index<=T_GRID::number_of_nodes_per_face;vertex_index++){
                        simplex_vertices(1)=cell_centroid;
                        simplex_vertices(2)=frame*grid.Face(face_index.axis,face_index.index);
                        for(int j=0;j<(TV::dimension-1);j++)
                            simplex_vertices(j+3)=frame*grid.Node(grid.Face_Node_Index(face_index.axis,face_index.index,WRAP(vertex_index+j,T_GRID::number_of_nodes_per_face)));
                        T_SIMPLEX simplex(simplex_vertices);
                        bool swapped=false;
                        if(simplex.Signed_Size()<0){
                            exchange(simplex_vertices(1),simplex_vertices(2));
                            simplex=T_SIMPLEX(simplex_vertices);
                        swapped=true;}
                        if(!simplex.Outside(x,thickness_factor*simplex.Minimum_Edge_Length())){
                            //if(x(2)>11.95 && x(2)<12.05)
                            //    LOG::cout << "interpolating from cartesian face " << x << " " << simplex_vertices(2) << std::endl;
                            VECTOR<T,TV::dimension+1> bary_coordinates=simplex.Barycentric_Coordinates(x);
                            TV vector;
                            if(swapped) vector=bary_coordinates(2)*cell_centroid_values(nearest_cell.x)(boundary_cell_index)+bary_coordinates(1)*face_centroid_values(nearest_cell.x).Get(face_index);
                            else vector=bary_coordinates(1)*cell_centroid_values(nearest_cell.x)(boundary_cell_index)+bary_coordinates(2)*face_centroid_values(nearest_cell.x).Get(face_index);
                            for(int j=0;j<(TV::dimension-1);j++)
                                vector+=bary_coordinates(j+3)*node_values(nearest_cell.x).Get(grid.Face_Node_Index(face_index.axis,face_index.index,WRAP(vertex_index+j,T_GRID::number_of_nodes_per_face)));
                                return vector;}}}}}
    
    return std::numeric_limits<T>::quiet_NaN()*TV::All_Ones_Vector();
}
//#####################################################################
// Function From_Cell_Faces
//#####################################################################
template<class T_GRID> typename T_GRID::VECTOR_T VORONOI_INTERPOLATION_CHIMERA<T_GRID>::From_Voronoi_Cell_Faces(const TV& x) const
{
    typedef MATRIX<T,2*TV::dimension,2*TV::dimension> D2_MATRIX;
    typedef VECTOR<T,2*TV::dimension> TV2;

    GRID_CELL_INDEX nearest_cell=Nearest_Cell(x);
    if(nearest_cell.x){
        T_GRID& grid=laplace_grid.Grid(nearest_cell.x);
        FRAME<TV> frame=laplace_grid.Frame(nearest_cell.x);
        //int boundary_cell_index=laplace_grid.boundary_cell_indices_to_linear_index(nearest_cell.x).Get(nearest_cell.y);
        
        TV2 rhs;
        D2_MATRIX A;
        
        if(laplace_grid.cell_indices_to_incident_voronoi_face_indices.Contains(nearest_cell)){
            ARRAY<int> incident_indices=laplace_grid.cell_indices_to_incident_voronoi_face_indices.Get(nearest_cell);
            for(int i=1;i<=incident_indices.Size();i++){
                int voronoi_face_index=incident_indices(i);
                TV location_1,location_2,location,normal;T distance;
                laplace_grid.Voronoi_Face(voronoi_face_index,location_1,location_2,location,normal,distance);
                TV offset=location-x;
                TV2 vector;
                for(int axis=1;axis<=TV::dimension;axis++){
                    vector(axis)=normal(axis);
                    vector(axis+TV::dimension)=normal(axis)*offset(axis);}
                A+=D2_MATRIX::Outer_Product(vector,vector);
                rhs+=vector*voronoi_face_values(voronoi_face_index);}}
        
        for(int axis=1;axis<=TV::dimension;axis++)
            for(int side=1;side<=2;side++){
                D_FACE_INDEX face_index(axis,nearest_cell.y+(side-1)*TV_INT::Axis_Vector(axis));
                if(laplace_grid.Chimera_Face(nearest_cell.x,face_index)){
                    TV location=frame*grid.Face(face_index.axis,face_index.index);
                    TV normal=frame.r.Rotate(TV::Axis_Vector(face_index.axis));
                    TV offset=location-x;
                    TV2 vector;
                    for(int axis=1;axis<=TV::dimension;axis++){
                        vector(axis)=normal(axis);
                        vector(axis+TV::dimension)=normal(axis)*offset(axis);}
                    A+=D2_MATRIX::Outer_Product(vector,vector);
                    rhs+=vector*face_values(nearest_cell.x)->operator()(face_index);}}

        /*TV2 divergence_constraint;
        for(int axis=1;axis<=TV::dimension;axis++)
            divergence_constraint(axis+TV::dimension)=1;
            A+=D2_MATRIX::Outer_Product(divergence_constraint,(T)100*divergence_constraint);*/
        
        TV2 solution=A.Solve_Linear_System(rhs);
        TV velocity;

        //LOG::cout << "interpolated " << x << " " << solution << std::endl;

        for(int axis=1;axis<=TV::dimension;axis++)
            velocity(axis)=solution(axis);
        return velocity;}

    return std::numeric_limits<T>::quiet_NaN()*TV::All_Ones_Vector();
}
//#####################################################################
template class VORONOI_INTERPOLATION_CHIMERA<GRID<VECTOR<float,1> > >;
template class VORONOI_INTERPOLATION_CHIMERA<GRID<VECTOR<float,2> > >;
template class VORONOI_INTERPOLATION_CHIMERA<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VORONOI_INTERPOLATION_CHIMERA<GRID<VECTOR<double,1> > >;
template class VORONOI_INTERPOLATION_CHIMERA<GRID<VECTOR<double,2> > >;
template class VORONOI_INTERPOLATION_CHIMERA<GRID<VECTOR<double,3> > >;
#endif
