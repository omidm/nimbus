//#####################################################################
// Copyright 2011, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace CUT_CELL_COMPUTATIONS
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/HEAPIFY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_SEGMENT_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Collisions_And_Grids/OBJECTS_IN_CELL.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/CUT_CELL.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
namespace PhysBAM{
namespace CUT_CELL_COMPUTATIONS{
//#####################################################################
// Helper Functions
//#####################################################################
namespace {
template<class T,int d>
bool Is_Occluded_Cell(const VECTOR<T,d> centroid,const VECTOR<T,d> cell_center,OBJECTS_IN_CELL<GRID<VECTOR<T,d> >,COLLISION_GEOMETRY_ID>& objects_in_cell,const COLLISION_GEOMETRY_COLLECTION<VECTOR<T,d> >& collision_geometry_collection,const VECTOR<int,d>& cell_index,const VECTOR<int,d>& neighbor_index)
{ARRAY<COLLISION_GEOMETRY_ID> collision_objects;objects_in_cell.Get_Objects_For_Cells(cell_index,neighbor_index,collision_geometry_collection.bodies.m,collision_objects);return collision_geometry_collection.Intersection_Between_Points(centroid,cell_center,&collision_objects);}
template<class T,int d>
bool Is_Occluded_Cell_Center(const VECTOR<T,d> centroid,const VECTOR<T,d> cell_center,OBJECTS_IN_CELL<GRID<VECTOR<T,d> >,COLLISION_GEOMETRY_ID>& objects_in_cell,const COLLISION_GEOMETRY_COLLECTION<VECTOR<T,d> >& collision_geometry_collection,const VECTOR<int,d>& cell_index)
{ARRAY<COLLISION_GEOMETRY_ID> collision_objects;objects_in_cell.Get_Objects_For_Cell(cell_index,collision_objects);return collision_geometry_collection.Intersection_Between_Points(centroid,cell_center,&collision_objects);}
};
//#####################################################################
// Function Compute_Cut_Geometries 1-D
//#####################################################################
template<class T>
void Compute_Cut_Geometries(const GRID<VECTOR<T,1> >& grid,const int num_ghost_cells,typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<T,1> > >::GRID_BASED_COLLISION_GEOMETRY& collision_bodies_affecting_fluid,ARRAY<CUT_CELLS<T,1>*,VECTOR<int,1> >& cut_cells){
    typedef VECTOR<T,1> TV;
    typedef VECTOR<int,1> TV_INT;
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid,num_ghost_cells);iterator.Valid();iterator.Next()){
        TV_INT index=iterator.Cell_Index();COLLISION_GEOMETRY_ID body_id;
        ARRAY<COLLISION_GEOMETRY_ID> collision_objects;collision_bodies_affecting_fluid.objects_in_cell.Get_Objects_For_Cell(index,collision_objects);
        if(!collision_objects.Size()) cut_cells(index)=0;
        else{
            cut_cells(index) = new CUT_CELLS<T,1>();
            RANGE<TV> cell_volume=iterator.Bounding_Box();
            RAY<TV> l_to_r_ray(cell_volume.min_corner,cell_volume.max_corner-cell_volume.min_corner,true);l_to_r_ray.t_max=l_to_r_ray.direction.Normalize();l_to_r_ray.semi_infinite=false;
            RAY<TV> r_to_l_ray(cell_volume.max_corner,cell_volume.min_corner-cell_volume.max_corner,true);r_to_l_ray.t_max=r_to_l_ray.direction.Normalize();r_to_l_ray.semi_infinite=false;
            if(collision_bodies_affecting_fluid.Intersection_With_Any_Simplicial_Object(l_to_r_ray,body_id,&collision_objects)){
                int poly_index=cut_cells(index)->geometry.Append(POLYGON<TV>(RANGE<TV>(cell_volume.min_corner, l_to_r_ray.Point(l_to_r_ray.t_max))));
                TV centroid=ARRAYS_COMPUTATIONS::Average(cut_cells(index)->geometry(poly_index).X);

                cut_cells(index)->visibility.Resize(poly_index);
                if(!Is_Occluded_Cell_Center<T,1>(centroid,grid.Center(index),collision_bodies_affecting_fluid.objects_in_cell,collision_bodies_affecting_fluid.collision_geometry_collection,index)){
                    cut_cells(index)->dominant_element=poly_index;
                    cut_cells(index)->visibility(poly_index).Append(index);}
                if(!cut_cells(index)->dominant_element)
                    for(int node=1;!cut_cells(index)->dominant_element && node<=cut_cells(index)->geometry(poly_index).X.Size();++node)
                        if(!Is_Occluded_Cell_Center<T,1>(cut_cells(index)->geometry(poly_index).X(node),grid.Center(index),collision_bodies_affecting_fluid.objects_in_cell,collision_bodies_affecting_fluid.collision_geometry_collection,index)){
                            cut_cells(index)->dominant_element=poly_index;
                            cut_cells(index)->visibility(poly_index).Append(index);}

                TV_INT neighbor_cell_index=index-TV_INT::Axis_Vector(1);
                if(!Is_Occluded_Cell<T,1>(centroid,grid.Center(neighbor_cell_index),collision_bodies_affecting_fluid.objects_in_cell,collision_bodies_affecting_fluid.collision_geometry_collection,index,neighbor_cell_index))
                    cut_cells(index)->visibility(poly_index).Append(neighbor_cell_index);}
            if(collision_bodies_affecting_fluid.Intersection_With_Any_Simplicial_Object(r_to_l_ray,body_id,&collision_objects)){
                int poly_index=cut_cells(index)->geometry.Append(POLYGON<TV>(RANGE<TV>(r_to_l_ray.Point(r_to_l_ray.t_max),cell_volume.max_corner)));
                TV centroid=ARRAYS_COMPUTATIONS::Average(cut_cells(index)->geometry(poly_index).X);

                cut_cells(index)->visibility.Resize(poly_index);
                if(!Is_Occluded_Cell_Center<T,1>(centroid,grid.Center(index),collision_bodies_affecting_fluid.objects_in_cell,collision_bodies_affecting_fluid.collision_geometry_collection,index)){
                    cut_cells(index)->dominant_element=poly_index;
                    cut_cells(index)->visibility(poly_index).Append(index);}
                if(!cut_cells(index)->dominant_element)
                    for(int node=1;!cut_cells(index)->dominant_element && node<=cut_cells(index)->geometry(poly_index).X.Size();++node)
                        if(!Is_Occluded_Cell_Center<T,1>(cut_cells(index)->geometry(poly_index).X(node),grid.Center(index),collision_bodies_affecting_fluid.objects_in_cell,collision_bodies_affecting_fluid.collision_geometry_collection,index)){
                            cut_cells(index)->dominant_element=poly_index;
                            cut_cells(index)->visibility(poly_index).Append(index);}

                TV_INT neighbor_cell_index=index+TV_INT::Axis_Vector(1);
                if(!Is_Occluded_Cell<T,1>(centroid,grid.Center(neighbor_cell_index),collision_bodies_affecting_fluid.objects_in_cell,collision_bodies_affecting_fluid.collision_geometry_collection,index,neighbor_cell_index))
                    cut_cells(index)->visibility(poly_index).Append(neighbor_cell_index);}
            if(!cut_cells(index)->geometry.Size()){delete cut_cells(index);cut_cells(index)=0;}}}

#if 0
    // for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,1);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
    //     if(!cut_cells(cell_index)) continue;
    //     std::stringstream ss;
    //     ss<<"Cut Cells in "<<cell_index<<"; "<<cut_cells(cell_index)->dominant_element<<std::endl;
    //     for(int i=1;i<=cut_cells(cell_index)->geometry.Size();++i){
    //         ss<<"Polygon "<<i<<":";for(int j=1;j<=cut_cells(cell_index)->geometry(i).X.Size();++j) ss<<" "<<cut_cells(cell_index)->geometry(i).X(j);ss<<std::endl;
    //         ss<<"\tCan See:";      for(int j=1;j<=cut_cells(cell_index)->visibility(i).Size();++j) ss<<" "<<cut_cells(cell_index)->visibility(i)(j);ss<<std::endl;}
    // LOG::filecout(ss.str());}
#endif
}
//#####################################################################
// Function Compute_Cut_Geometries 2-D
//#####################################################################
namespace{
template<class T>
struct CLIP_ENTRY
{
    int i;
    T a,b;
};

template<class T>
void Compute_Candidate_Nodes(const GRID<VECTOR<T,1> >& grid,const ARRAY<CUT_CELLS<T,1>*,VECTOR<int,1> >& cut_cells,const VECTOR<int,1>& cell_index,const int poly_index,ARRAY<VECTOR<T,1> >& candidate_nodes,typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<T,1> > >::GRID_BASED_COLLISION_GEOMETRY& collision_bodies_affecting_fluid)
{if(cut_cells(cell_index)) candidate_nodes.Append(ARRAYS_COMPUTATIONS::Average(cut_cells(cell_index)->geometry(poly_index).X)); else candidate_nodes.Append(grid.Center(cell_index));}

template<class T>
void Compute_Candidate_Nodes(const GRID<VECTOR<T,2> >& grid,const ARRAY<CUT_CELLS<T,2>*,VECTOR<int,2> >& cut_cells,const VECTOR<int,2>& cell_index,const int poly_index,ARRAY<VECTOR<T,2> >& candidate_nodes,typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<T,2> > >::GRID_BASED_COLLISION_GEOMETRY& collision_bodies_affecting_fluid)
{
    COLLISION_GEOMETRY_ID dummy_index;int dummy_aggregate_id;
    if(cut_cells(cell_index)){
        candidate_nodes.Append(ARRAYS_COMPUTATIONS::Average(cut_cells(cell_index)->geometry(poly_index).X));
        for(int i=1;i<=cut_cells(cell_index)->geometry(poly_index).X.Size();++i) candidate_nodes.Append(cut_cells(cell_index)->geometry(poly_index).X(i));
        for(int i=candidate_nodes.Size();i>=1;--i)
            if(collision_bodies_affecting_fluid.Inside_Any_Simplex_Of_Any_Body(candidate_nodes(i),dummy_index,dummy_aggregate_id)) candidate_nodes.Remove_Index(i);
        for(int i=candidate_nodes.Size();i>=1;--i) if(!cut_cells(cell_index)->geometry(poly_index).Inside_Polygon(candidate_nodes(i))) candidate_nodes.Remove_Index(i);

        if(!candidate_nodes.Size()){
            // std::stringstream ss;ss<<"WARNING: Unable to find a usual candidate node for cell "<<cell_index<<" poly "<<poly_index<<"; trying edge-intermediate nodes"<<std::endl;LOG::filecout(ss.str());
            for(int i=1;i<=cut_cells(cell_index)->geometry(poly_index).X.Size();++i){
                VECTOR<T,2> average_point=(T).5*(cut_cells(cell_index)->geometry(poly_index).X(i)+cut_cells(cell_index)->geometry(poly_index).X((i%cut_cells(cell_index)->geometry(poly_index).X.Size())+1));
                if(!collision_bodies_affecting_fluid.Inside_Any_Simplex_Of_Any_Body(average_point,dummy_index,dummy_aggregate_id)) candidate_nodes.Append(average_point);}}

    }
    else candidate_nodes.Append(grid.Center(cell_index));
}

template<class T>
void Compute_Candidate_Nodes(const GRID<VECTOR<T,3> >& grid,const ARRAY<CUT_CELLS<T,3>*,VECTOR<int,3> >& cut_cells,const VECTOR<int,3>& cell_index,const int poly_index,ARRAY<VECTOR<T,3> >& candidate_nodes,typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<T,3> > >::GRID_BASED_COLLISION_GEOMETRY& collision_bodies_affecting_fluid)
{PHYSBAM_NOT_IMPLEMENTED("Cut cells not supported in 3-D");}

template<class T>
void Compute_Divided_Polygons(ARRAY<POLYGON<VECTOR<T,2> > >& polygons,const int active_polygon,ARRAY<VECTOR<T,2> >& vertices_to_add,const T thickness_over_two)
{
    ARRAY<VECTOR<T,2> > old_polygon_positions(polygons(active_polygon).X);polygons(active_polygon).X.Resize(0);
    int cut_poly_start=0,cut_poly_end=0;
    for(int i=1;i<=old_polygon_positions.Size() && !(cut_poly_start && cut_poly_end);++i){
        SEGMENT_2D<T> poly_edge(old_polygon_positions(i),old_polygon_positions((i%old_polygon_positions.Size())+1));
        if(poly_edge.Inside(vertices_to_add(1),thickness_over_two)) cut_poly_start=i;
        if(poly_edge.Inside(vertices_to_add(vertices_to_add.Size()),thickness_over_two)) cut_poly_end=i;}

    if(cut_poly_start == cut_poly_end && vertices_to_add.Size()>2 && (old_polygon_positions(cut_poly_start)-vertices_to_add(1)).Magnitude_Squared() > (old_polygon_positions(cut_poly_start)-vertices_to_add(vertices_to_add.Size())).Magnitude_Squared())
        ARRAYS_COMPUTATIONS::Reverse_In_Place(vertices_to_add);
    const int new_poly_index=polygons.Append(POLYGON<VECTOR<T,2> >());
    if(cut_poly_start < cut_poly_end || (cut_poly_start==cut_poly_end && vertices_to_add.Size()>2)){
        for(int i=1;                     i<=cut_poly_start;              ++i) polygons(active_polygon).X.Append(old_polygon_positions(i));
        for(int i=1;                     i<=vertices_to_add.Size();      ++i) polygons(active_polygon).X.Append(vertices_to_add(i));
        for(int i=cut_poly_end+1;        i<=old_polygon_positions.Size();++i) polygons(active_polygon).X.Append(old_polygon_positions(i));
        for(int i=cut_poly_start+1;      i<=cut_poly_end;                ++i) polygons(new_poly_index).X.Append(old_polygon_positions(i));
        for(int i=vertices_to_add.Size();i>=1;                           --i) polygons(new_poly_index).X.Append(vertices_to_add(i));}
    else if(cut_poly_end < cut_poly_start){
        for(int i=1;                     i<=cut_poly_end;                ++i) polygons(active_polygon).X.Append(old_polygon_positions(i));
        for(int i=vertices_to_add.Size();i>=1;                           --i) polygons(active_polygon).X.Append(vertices_to_add(i));
        for(int i=cut_poly_start+1;      i<=old_polygon_positions.Size();++i) polygons(active_polygon).X.Append(old_polygon_positions(i));
        for(int i=cut_poly_end+1;        i<=cut_poly_start;              ++i) polygons(new_poly_index).X.Append(old_polygon_positions(i));
        for(int i=1;                     i<=vertices_to_add.Size();      ++i) polygons(new_poly_index).X.Append(vertices_to_add(i));}
    else{
        polygons(active_polygon).X=old_polygon_positions;
        polygons.Remove_End();
        return;}
}

template<class T>
bool Inside_Any_Simplex(const VECTOR<T,2>& min_corner,const VECTOR<T,2>& max_corner,const VECTOR<T,2>& X,const T thickness_over_two,int& side)
{
    SEGMENT_2D<T> segment(min_corner,min_corner);
    segment.x2.x=max_corner.x;segment.x2.y=min_corner.y;if(segment.Inside(X,thickness_over_two)){side=1;return true;}
    segment.x1=max_corner;                              if(segment.Inside(X,thickness_over_two)){side=2;return true;}
    segment.x2.x=min_corner.x;segment.x2.y=max_corner.y;if(segment.Inside(X,thickness_over_two)){side=3;return true;}
    segment.x1=min_corner;                              if(segment.Inside(X,thickness_over_two)){side=4;return true;}
    return false;
}

template<class T>
bool Robust_Inside_Polygon(const POLYGON<VECTOR<T,2> >& poly,const VECTOR<T,2>& X,const T thickness_over_two)
{
    for(int i=1;i<=poly.X.Size();++i){
        SEGMENT_2D<T> poly_edge(poly.X(i),poly.X((i%poly.X.Size())+1));
        if(poly_edge.Inside(X,thickness_over_two)) return true;}
    return poly.Inside_Polygon(X);
}
};

template<class T>
void Compute_Cut_Geometries(const GRID<VECTOR<T,2> >& grid,const int num_ghost_cells,typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<T,2> > >::GRID_BASED_COLLISION_GEOMETRY& collision_bodies_affecting_fluid,ARRAY<CUT_CELLS<T,2>*,VECTOR<int,2> >& cut_cells)
{
    LOG::SCOPE scope("Compute_Cut_Geometries");
    std::stringstream ss;
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;
    ARRAY<TV_INT> neighbor_offsets(8);
    neighbor_offsets(1)=TV_INT( 1,-1); neighbor_offsets(2)=TV_INT( 1, 0); neighbor_offsets(3)=TV_INT( 1, 1);
    neighbor_offsets(4)=TV_INT( 0,-1);                                    neighbor_offsets(5)=TV_INT( 0, 1);
    neighbor_offsets(6)=TV_INT(-1,-1); neighbor_offsets(7)=TV_INT(-1, 0); neighbor_offsets(8)=TV_INT(-1, 1);

    const T collision_thickness_over_two=(T).5*collision_bodies_affecting_fluid.collision_thickness;
    const T sqr_collision_thickness_over_two=collision_thickness_over_two*collision_thickness_over_two;
    for(COLLISION_GEOMETRY_ID body(1);body<=collision_bodies_affecting_fluid.collision_geometry_collection.bodies.Size();++body){
        if(!collision_bodies_affecting_fluid.Is_Active(body)) continue;
        COLLISION_GEOMETRY<TV>& collision_geometry=collision_bodies_affecting_fluid.collision_geometry_collection(body);
    
        HASHTABLE<TV_INT,ARRAY<CLIP_ENTRY<T> > > cell_simplices;
        for(int i=1;i<=collision_geometry.Number_Of_Simplices();i++){SEGMENT_2D<T> simplex=collision_geometry.World_Space_Simplex(i);
            RANGE<TV> box(simplex.Bounding_Box());TV_INT min_index=grid.Clamp_To_Cell(box.Minimum_Corner()-collision_thickness_over_two,num_ghost_cells),max_index=grid.Clamp_To_Cell(box.Maximum_Corner()+collision_thickness_over_two,num_ghost_cells);
            for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,RANGE<TV_INT>(min_index,max_index));iterator.Valid();iterator.Next()){T start,end;
                if(simplex.Clip_To_Box(iterator.Bounding_Box(),start,end,collision_thickness_over_two)){
                    CLIP_ENTRY<T> entry;entry.i=i;entry.a=start;entry.b=end;cell_simplices.Get_Or_Insert(iterator.Cell_Index()).Append(entry);}}}

        for(typename HASHTABLE<TV_INT,ARRAY<CLIP_ENTRY<T> > >::ITERATOR iterator(cell_simplices);iterator.Valid();iterator.Next()){const TV_INT& cell_index=iterator.Key();
            ARRAY<CLIP_ENTRY<T> >& segments(iterator.Data());
            if(!segments.Size()) continue;

            bool was_fixed=false;
            do{was_fixed=false;
                for(int i=1;i<=segments.Size() && !was_fixed;++i){
                    SEGMENT_2D<T> segment1=collision_geometry.World_Space_Simplex(segments(i).i);
                    for(int j=1;j<=segments.Size() && !was_fixed;++j){
                        SEGMENT_2D<T> segment2=collision_geometry.World_Space_Simplex(segments(j).i);
                        if(segment1.x2 == segment2.x1 && j != i+1){was_fixed=true;
                            if(i > j){exchange(segments(i),segments(j));exchange(segments(j+1),segments(i));}
                            else if(i < j) exchange(segments(i+1),segments(j));
                            else ss<<"STRANGENESS MYAR!"<<std::endl;}}}
            }while(was_fixed);

            if(!cut_cells(cell_index)){cut_cells(cell_index)=new CUT_CELLS<T,2>();cut_cells(cell_index)->geometry.Append(POLYGON<TV>(RANGE<TV>(grid.Node(cell_index),grid.Node(cell_index)+grid.dX)));}

            ARRAY<PAIR<TV,int> > notable_points;
            for(int i=1;i<=segments.Size();++i){SEGMENT_2D<T> segment(collision_geometry.World_Space_Simplex(segments(i).i));
                if(i==1 || segments(i).a==0) notable_points.Append(PAIR<TV,int>(segment.Point_From_Barycentric_Coordinates(segments(i).a),0));
                notable_points.Append(PAIR<TV,int>(segment.Point_From_Barycentric_Coordinates(segments(i).b),0));}
            for(int i=1;i<=notable_points.Size();++i) Inside_Any_Simplex<T>(grid.Node(cell_index),grid.Node(cell_index)+grid.dX,notable_points(i).x,collision_thickness_over_two,notable_points(i).y);

            while(notable_points.Size() && !notable_points(1).y) notable_points.Remove_Index(1);
            while(notable_points.Size() && !notable_points(notable_points.Size()).y) notable_points.Remove_Index(notable_points.Size());

            for(int i=notable_points.Size();i > 1;--i)
                if((notable_points(i).x-notable_points(i-1).x).Magnitude_Squared() < (T)1e-4*sqr_collision_thickness_over_two){
                    if(notable_points(i-1).y) notable_points.Remove_Index(i);
                    else notable_points.Remove_Index(i-1);}
                else if(notable_points(i-1).y && notable_points(i).y && notable_points(i-1).y == notable_points(i).y) notable_points.Remove_Index(i-1);

            ARRAY<POLYGON<TV> >& polygons(cut_cells(cell_index)->geometry);
            for(int index=1;index<notable_points.Size();++index){assert(notable_points(index).y);
                int next_index=index+1;while(!notable_points(next_index).y) ++next_index;
                for(int geom=polygons.Size();geom>=1;--geom)
                    if(Robust_Inside_Polygon(polygons(geom),notable_points(index).x,collision_thickness_over_two) && Robust_Inside_Polygon(polygons(geom),notable_points(next_index).x,collision_thickness_over_two)){
                        ARRAY<TV> vertices_to_add(next_index-index+1);
                        for(int i=0;i<=next_index-index;++i) vertices_to_add(i+1) = notable_points(index+i).x;
#if 0
                        ss<<"Calling Compute_Divided_Polygons on cell "<<cell_index<<std::endl;
                        for(int i=1;i<=vertices_to_add.Size();++i) ss<<"\t"<<vertices_to_add(i)<<std::endl;
#endif
                        ARRAY<TV> saved_old_polygon(cut_cells(cell_index)->geometry(geom).X);
                        Compute_Divided_Polygons(cut_cells(cell_index)->geometry,geom,vertices_to_add,collision_thickness_over_two);
                        ARRAY<TV> test_nodes;Compute_Candidate_Nodes(grid,cut_cells,cell_index,cut_cells(cell_index)->geometry.Size(),test_nodes,collision_bodies_affecting_fluid);
                        if(!test_nodes.Size()){
                            ss<<"WARNING: A cut cell in cell "<<cell_index<<" was computed which has no visible neighbors; discarding!"<<std::endl;
                            CUT_CELLS<T,2>::Print_Debug_Information(cut_cells(cell_index)->geometry);cut_cells(cell_index)->geometry(geom).X=saved_old_polygon;cut_cells(cell_index)->geometry.Remove_End();}
                        index=next_index;
                        if(index < notable_points.Size() && !notable_points(index+1).y) --index;
                        continue;}}}}

    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(!cut_cells(cell_index)) continue;
        for(int i=cut_cells(cell_index)->geometry.Size();i>=1;--i) if(cut_cells(cell_index)->geometry(i).Area() < (T)1e-4*grid.Cell_Size()) cut_cells(cell_index)->geometry.Remove_Index(i);

        cut_cells(cell_index)->visibility.Resize(cut_cells(cell_index)->geometry.Size());
        cut_cells(cell_index)->visibility_nodes.Resize(cut_cells(cell_index)->geometry.Size());
        for(int poly_index=1;poly_index<=cut_cells(cell_index)->geometry.Size();++poly_index){
            Compute_Candidate_Nodes(grid,cut_cells,cell_index,poly_index,cut_cells(cell_index)->visibility_nodes(poly_index),collision_bodies_affecting_fluid);
            ARRAY<TV>& test_nodes(cut_cells(cell_index)->visibility_nodes(poly_index));

            bool found_visible_candidate=false;
            for(int i=1;i<=test_nodes.Size() && !found_visible_candidate;++i){
                if(!Is_Occluded_Cell_Center<T,2>(test_nodes(i),grid.Center(cell_index),collision_bodies_affecting_fluid.objects_in_cell,
                                                 collision_bodies_affecting_fluid.collision_geometry_collection,cell_index)){
                    found_visible_candidate=true;
                    cut_cells(cell_index)->dominant_element=poly_index;
                    cut_cells(cell_index)->visibility(poly_index).Append(cell_index);}}}}

    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(!cut_cells(cell_index)) continue;
        for(int poly_index=1;poly_index<=cut_cells(cell_index)->geometry.Size();++poly_index){
            ARRAY<TV>& test_nodes(cut_cells(cell_index)->visibility_nodes(poly_index));
            for(int n=1;n<=neighbor_offsets.Size();++n){TV_INT neighbor_cell_index=cell_index+neighbor_offsets(n);
                if(!cut_cells.Valid_Index(neighbor_cell_index)) continue;
                bool found_visible_candidate=false;
                for(int i=1;i<=test_nodes.Size() && !found_visible_candidate;++i){
                    if(cut_cells(neighbor_cell_index) && !cut_cells(neighbor_cell_index)->dominant_element) continue;
                    if(!Is_Occluded_Cell<T,2>(test_nodes(i),grid.Center(neighbor_cell_index),collision_bodies_affecting_fluid.objects_in_cell,
                                              collision_bodies_affecting_fluid.collision_geometry_collection,cell_index,neighbor_cell_index)){
                        cut_cells(cell_index)->visibility(poly_index).Append(neighbor_cell_index);
                        found_visible_candidate=true;}}}}}

    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(cut_cells(cell_index)) for(int poly=cut_cells(cell_index)->geometry.Size();poly>=1;--poly)
            if(!cut_cells(cell_index)->visibility(poly).Size()){
                ss<<"WARNING: Trimming polygon "<<poly<<" from cut cells in "<<cell_index<<"; no visible neighbors!"<<std::endl;
                cut_cells(cell_index)->geometry.Remove_Index(poly);cut_cells(cell_index)->visibility.Remove_Index(poly);cut_cells(cell_index)->visibility_nodes.Remove_Index(poly);
                if(cut_cells(cell_index)->dominant_element > poly) cut_cells(cell_index)->dominant_element--;
                else if(cut_cells(cell_index)->dominant_element == poly) cut_cells(cell_index)->dominant_element=0;}}

#if 0 // Debug output
    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(!cut_cells(cell_index)) continue;
        ss<<"Cut Cells in "<<cell_index<<"; "<<cut_cells(cell_index)->dominant_element<<std::endl;
        for(int i=1;i<=cut_cells(cell_index)->geometry.Size();++i){
            ss<<"Polygon "<<i<<", "<<cut_cells(cell_index)->geometry(i).Area()<<"\t:";for(int j=1;j<=cut_cells(cell_index)->geometry(i).X.Size();++j) ss<<" "<<cut_cells(cell_index)->geometry(i).X(j);ss<<std::endl;
            ss<<"\tCan See:";      for(int j=1;j<=cut_cells(cell_index)->visibility(i).Size();++j) ss<<" "<<cut_cells(cell_index)->visibility(i)(j);ss<<std::endl;
            if(cut_cells(cell_index)->geometry.Size()==2)
                for(int j=1;j<=cut_cells(cell_index)->geometry.Size();++j) if(i!=j){
                    for(int a=1;a<=cut_cells(cell_index)->visibility(i).Size();++a) for(int b=1;b<=cut_cells(cell_index)->visibility(j).Size();++b)
                        if(cut_cells(cell_index)->visibility(i)(a) == cut_cells(cell_index)->visibility(j)(b)){
                            ss<<"ERROR: Two cut cells in cell "<<cell_index<<" can see the same cell ("<<cut_cells(cell_index)->visibility(j)(b)<<"); polygons "<<i<<" and "<<j<<std::endl;
                            CUT_CELLS<T,TV::dimension>::Print_Debug_Information(cut_cells(cell_index)->geometry);}
                }}}
#endif
    LOG::filecout(ss.str());
}
//#####################################################################
// Function Compute_Cut_Geometries 3-D
//#####################################################################
template<class T>
void Compute_Cut_Geometries(const GRID<VECTOR<T,3> >& grid,const int num_ghost_cells,typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<T,3> > >::GRID_BASED_COLLISION_GEOMETRY& collision_bodies_affecting_fluid,ARRAY<CUT_CELLS<T,3>*,VECTOR<int,3> >& cut_cells)
{PHYSBAM_NOT_IMPLEMENTED();}
//#####################################################################
// Function Compute_Temporal_Neighbors
//#####################################################################
template<class T,int d>
void Compute_Temporal_Neighbors(const GRID<VECTOR<T,d> >& grid,const int num_ghost_cells,const T dt,const ARRAY<VECTOR<int,d> >& neighbor_offsets,
                                typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<T,d> > >::GRID_BASED_COLLISION_GEOMETRY& collision_bodies_affecting_fluid,
                                const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >& cut_cells_n,const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >& cut_cells_np1,const ARRAY<bool,VECTOR<int,d> >& psi_np1,
                                HASHTABLE<PAIR<VECTOR<int,d>,int>,ARRAY<PAIR<VECTOR<int,d>,int> > >& visibility_graph)
{
    LOG::SCOPE scope("Compute_Temporal_Neighbors");
    std::stringstream ss;
    typedef VECTOR<T,d> TV;typedef VECTOR<int,d> TV_INT;
    visibility_graph.Clean_Memory();

    for(typename GRID<TV>::CELL_ITERATOR iter_n(grid);iter_n.Valid();iter_n.Next()){const TV_INT& cell_index=iter_n.Cell_Index();
        bool borders_cut_cell=(cut_cells_n(cell_index)!=0 || cut_cells_np1(cell_index)!=0);
        for(int i=1;i<=neighbor_offsets.Size() && !borders_cut_cell;++i){const TV_INT neighbor_index=cell_index+neighbor_offsets(i);
            if(!cut_cells_n.Valid_Index(neighbor_index)) continue;
            borders_cut_cell |= (cut_cells_n(neighbor_index)!=0 || cut_cells_np1(neighbor_index)!=0);}
        if(!borders_cut_cell) continue;

        // TODO(jontg): Precompute candidate nodes only once per cell
        if(cut_cells_n(cell_index)){
            for(int poly_index_n=1;poly_index_n<=cut_cells_n(cell_index)->geometry.Size();++poly_index_n){
                visibility_graph.Get_Or_Insert(PAIR<TV_INT,int>(cell_index,poly_index_n));
                ARRAY<TV>& test_nodes_n(cut_cells_n(cell_index)->visibility_nodes(poly_index_n));
                for(int n=1;n<=neighbor_offsets.Size();++n){const TV_INT neighbor_index=cell_index+neighbor_offsets(n);
                    if(!cut_cells_n.Valid_Index(neighbor_index)) continue;
                    if(cut_cells_np1(neighbor_index)){
                        for(int poly_index_np1=1;poly_index_np1<=cut_cells_np1(neighbor_index)->geometry.Size();++poly_index_np1){
                            ARRAY<TV>& test_nodes_np1(cut_cells_np1(neighbor_index)->visibility_nodes(poly_index_np1));
                            bool found_connection=false;
                            for(int i=1;i<=test_nodes_n.Size() && !found_connection;++i) for(int j=1;j<=test_nodes_np1.Size() && !found_connection;++j)
                                if(!collision_bodies_affecting_fluid.Any_Simplex_Crossover(test_nodes_n(i),test_nodes_np1(j),dt)) found_connection=true;//ss<<"Found a connection between "<<cell_index<<", "<<poly_index_n<<" and "<<neighbor_index<<", "<<poly_index_np1<<std::endl;
                            if(found_connection){visibility_graph.Get_Or_Insert(PAIR<TV_INT,int>(cell_index,poly_index_n)).Append(PAIR<TV_INT,int>(neighbor_index,poly_index_np1));
                            }}}
                    else{
                        if(!psi_np1(neighbor_index)) continue;
                        ARRAY<TV> test_nodes_np1;test_nodes_np1.Append(grid.Center(neighbor_index));
                        bool found_connection=false;
                        for(int i=1;i<=test_nodes_n.Size() && !found_connection;++i) for(int j=1;j<=test_nodes_np1.Size() && !found_connection;++j)
                            if(!collision_bodies_affecting_fluid.Any_Simplex_Crossover(test_nodes_n(i),test_nodes_np1(j),dt)) found_connection=true;
                        if(found_connection) visibility_graph.Get_Or_Insert(PAIR<TV_INT,int>(cell_index,poly_index_n)).Append(PAIR<TV_INT,int>(neighbor_index,0));}}}}
        else{
            visibility_graph.Get_Or_Insert(PAIR<TV_INT,int>(cell_index,0));
            ARRAY<TV> test_nodes_n;test_nodes_n.Append(grid.Center(cell_index));
            for(int n=1;n<=neighbor_offsets.Size();++n){const TV_INT neighbor_index=cell_index+neighbor_offsets(n);
                if(!cut_cells_n.Valid_Index(neighbor_index)) continue;
                if(cut_cells_np1(neighbor_index)){
                    for(int poly_index_np1=1;poly_index_np1<=cut_cells_np1(neighbor_index)->geometry.Size();++poly_index_np1){
                        ARRAY<TV>& test_nodes_np1(cut_cells_np1(neighbor_index)->visibility_nodes(poly_index_np1));
                        bool found_connection=false;
                        for(int i=1;i<=test_nodes_n.Size() && !found_connection;++i) for(int j=1;j<=test_nodes_np1.Size() && !found_connection;++j)
                            if(!collision_bodies_affecting_fluid.Any_Simplex_Crossover(test_nodes_n(i),test_nodes_np1(j),dt)) found_connection=true;
                        if(found_connection) visibility_graph.Get_Or_Insert(PAIR<TV_INT,int>(cell_index,0)).Append(PAIR<TV_INT,int>(neighbor_index,poly_index_np1));}}
                else{
                    if(!psi_np1(neighbor_index)) continue;
                    ARRAY<TV> test_nodes_np1;test_nodes_np1.Append(grid.Center(neighbor_index));
                    bool found_connection=false;
                    for(int i=1;i<=test_nodes_n.Size() && !found_connection;++i) for(int j=1;j<=test_nodes_np1.Size() && !found_connection;++j)
                        if(!collision_bodies_affecting_fluid.Any_Simplex_Crossover(test_nodes_n(i),test_nodes_np1(j),dt)) found_connection=true;
                    if(found_connection) visibility_graph.Get_Or_Insert(PAIR<TV_INT,int>(cell_index,0)).Append(PAIR<TV_INT,int>(neighbor_index,0));}}}}

#if 0
    for(typename HASHTABLE<PAIR<VECTOR<int,d>,int>,ARRAY<PAIR<VECTOR<int,d>,int> > >::ITERATOR iter(visibility_graph);iter.Valid();iter.Next()){
        ARRAY<PAIR<TV_INT,int> > visible_np1_neighbors;
        ss<<"Cell "<<iter.Key().x<<", Poly "<<iter.Key().y<<std::endl;
        for(int i=1;i<=iter.Data().Size();++i)
            if(visible_np1_neighbors.Contains(iter.Data()(i))) ss<<"ERROR: Two distinct polygons can see the same Np1 neighbor!"<<std::endl;
            else visible_np1_neighbors.Append(iter.Data()(i));
        for(int i=1;i<=visible_np1_neighbors.Size();++i) ss<<"\tCell "<<visible_np1_neighbors(i).x<<", Poly "<<visible_np1_neighbors(i).y<<std::endl;}
#endif
    LOG::filecout(ss.str());
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template void Compute_Cut_Geometries(const GRID<VECTOR<T,d> >&,const int,GRID_BASED_COLLISION_GEOMETRY_UNIFORM<GRID<VECTOR<T,d> > >&,ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >&); \
    template void Compute_Temporal_Neighbors(const GRID<VECTOR<T,d> >&,const int,const T,const ARRAY<VECTOR<int,d> >&, \
                                             COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<VECTOR<T,d> > >::GRID_BASED_COLLISION_GEOMETRY&, \
                                             const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >&,const ARRAY<CUT_CELLS<T,d>*,VECTOR<int,d> >&, const ARRAY<bool,VECTOR<int,d> >&,\
                                             HASHTABLE<PAIR<VECTOR<int,d>,int>,ARRAY<PAIR<VECTOR<int,d>,int> > >&);

INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif
}
}
