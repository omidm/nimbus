//#####################################################################
// Copyright 2006-2008, Kevin Der, Geoffrey Irving, Andrew Selle, Eftychios Sifakis, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUTTING_TETRAHEDRALIZED_VOLUME
//##################################################################### 
#ifndef __CUTTING_TETRAHEDRALIZED_VOLUME__
#define __CUTTING_TETRAHEDRALIZED_VOLUME__

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Geometry/Topology/POLYGON_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/ROBUST_SIMPLEX_INTERACTIONS.h>
#include <PhysBAM_Dynamics/Fracture/CUTTING_PARTICLES.h>
#include <PhysBAM_Dynamics/Fracture/CUTTING_POLYGON.h>
#include <PhysBAM_Dynamics/Fracture/INTERSECTION_REGISTRY.h>
namespace PhysBAM{

template<class T> class STACK;
template<class T> class TETRAHEDRALIZED_VOLUME;
template<class T> class TRIANGULATED_SURFACE;
template<class T> class GEOMETRY_PARTICLES;

template<class T>
class CUTTING_TETRAHEDRALIZED_VOLUME
{
public:
    typedef VECTOR<T,3> TV;

    // CURRENT STATE
    TETRAHEDRALIZED_VOLUME<T>* current_tetrahedralized_volume;
    
    // VALID ONLY FOR THE CURRENT CUT ITERATION
    const TRIANGULATED_SURFACE<T>* cutting_triangulated_surface;
    TETRAHEDRALIZED_VOLUME<T>* next_tetrahedralized_volume;
    TRIANGLE_MESH final_duplicated_boundary_mesh;
    ARRAY<ARRAY<int> > final_parents_per_new_particle;
    ARRAY<ARRAY<T> > final_parent_weights_per_new_particle;
    ARRAY<int> old_particle_per_new_collapsed_particle;
    ARRAY<ARRAY<int> > new_collapsed_tet_particle_per_old_tet_particle;
    ARRAY<int> previous_particle_index_per_new_particle_index;

private:
    // replaced at the end of current iteration
    CUTTING_SIMPLICES<T,3>* cutting_simplices;
    INTERSECTION_REGISTRY<T,3>* intersection_registry;
    POLYGON_MESH polygon_mesh;
    ARRAY<ARRAY<int> > simplices_per_current_tet;
    ARRAY<CUTTING_POLYGON> current_cutting_polygons;
    ARRAY<ARRAY<ARRAY<int> > > regions_per_tet;
    HASHTABLE<VECTOR<int,3>,VECTOR<int,3> > sorted_to_original_boundary_nodes;
    CUTTING_PARTICLES cutting_particles;

    // overwritten at end of current iteration and used in the next
    ARRAY<ARRAY<int> > cutting_polygons_per_cutting_simplex;

    // must be refreshed when first used in the current iteration
    int intersection_counter;
    ARRAY<int> last_old_simplex_index_per_current_tet;
    ARRAY<ARRAY<int> > intersecting_simplices_per_simplex;
    ARRAY<ARRAY<int> > polygons_per_element;
    HASHTABLE<PAIR<int,int>,int> hash_all_new_uncollapsed_particles; // <dup tet, orig particle> -> dup particle
    ARRAY<ARRAY<int> > uncollapsed_new_particles_per_all_current_particle_ids;
    ARRAY<ARRAY<int> > new_parents_per_new_particle;
    ARRAY<ARRAY<T> > new_parent_weights_per_new_particle;
    ARRAY<int> current_particle_id_per_uncollapsed_new_particle;
    ARRAY<ARRAY<int> > new_tets_per_current;
    ARRAY<int> current_particle_id_per_collapsed_new_particle;
    int num_new_tets,num_new_particles;
    UNION_FIND<> union_vertices;
    ARRAY<int> new_particle_indices;    
    HASHTABLE<int,int> dup_tet_before_to_after_collapse;
    HASHTABLE<int,int> dup_tet_after_to_before_collapse;

    // persistent
    T intersection_thickness;
    bool verbose;
    T abs_tol_for_level1_barycentric_coordinates;
    T abs_tol_for_level2_barycentric_coordinates;

public:

    CUTTING_TETRAHEDRALIZED_VOLUME(const bool verbose_input=false);    
    ~CUTTING_TETRAHEDRALIZED_VOLUME();

    const CUTTING_SIMPLICES<T,3>& Cutting_Simplices() const
    {return *cutting_simplices;}

    const INTERSECTION_REGISTRY<T,3>& Intersection_Registry() const
    {return *intersection_registry;}

    const CUTTING_PARTICLES& Cutting_Particles() const
    {return cutting_particles;}

    const POLYGON_MESH& Polygon_Mesh() const
    {return polygon_mesh;}

    const ARRAY< CUTTING_POLYGON >& Cutting_Polygons() const
    {return current_cutting_polygons;}

    void Set_Intersection_Thickness(const T intersection_thickness_input=0)  // Used only for hierarchy intersections; exact intersections performed at the simplex level
    {intersection_thickness=intersection_thickness_input;}

//#####################################################################
    void Initialize_Original_Embedding(const TETRAHEDRALIZED_VOLUME<T>& original_tetrahedralized_volume_input);
    void Cut_Material(const TRIANGULATED_SURFACE<T>& cutting_triangulated_surface_input,const int first_new=1);
    void Local_Planar_Coordinates_For_Intersection(int intersection_index, int cutting_simplex_index, VECTOR<T,2>& planar_coordinates) const;
    bool World_Coordinates_For_Intersection(int intersection_index, TV& world_coordinates) const; // returns true if on cut surface, false if tet node
    bool Is_Cutting_Particle(int intersection_index) const; // returns true if on cut surface, false if tet node
    void Write_Consistent_State(const STREAM_TYPE stream_type,const std::string& output_directory,const int suffix) const;
    void Read_Consistent_State(const STREAM_TYPE stream_type,const std::string& input_directory,const int suffix);
private:
    void Find_Triangle_Tet_Intersections(const int first_new);
    void Find_Triangle_Triangle_Intersections();
    void Fill_Intersection_Registry_With_Cuts();
    void Perform_Smart_Simplex_Intersection(const VECTOR<int,3>& simplices);
    void Split_Existing_Polygon_Edges();
    void Split_Existing_Polygons();
    void Add_Polygons_On_New_Simplices();
    void Find_Material_Regions();
    void Divide_Polygon_Particles_With_New_Segments(ARRAY<VECTOR<int,2> >& all_segments,const ARRAY<VECTOR<int,2> >& possible_segments_to_add,const ARRAY<ARRAY<int> >& polygon_particles,
        const int cutting_simplex,const bool flipped,ARRAY<ARRAY<ARRAY<int > > >& final_polygon_element_particles) const;
    void Get_Unoriented_Segments_On_Two_Simplices(const int simplex_1,const int simplex_2,ARRAY<VECTOR<int,2> >& new_unoriented_segments_on_simplex) const;
    void Get_Oriented_Segments_On_Edge_Of_Two_Simplices(const int simplex_1,const int simplex_2,const int shared_node_1,const int shared_node_2,
        ARRAY<VECTOR<int,2> >& new_oriented_segments_on_simplex) const;
    bool Potential_Segment_Should_Be_Added_To_Polygon(const ARRAY<ARRAY<int > >& particles_for_polygon,const bool flipped,const int simplex,const VECTOR<int,2>& nodes) const;
    bool Point_Is_Inside_Unoriented_Polygon(const ARRAY<int>& polygon_particles,const int simplex_owner,const int point) const;
    void Two_Dimensional_Region_Finding_On_Cutting_Simplex(const int cutting_simplex,const bool flipped,const ARRAY<VECTOR<int,2> >& segments,
        ARRAY<ARRAY<VECTOR<int,2> > >& unconnected_polygonal_regions,const bool survive_incomplete_regions=false) const;
    void Inside_Outside_Determination_For_Unconnected_Polygonal_Regions(const int simplex,const bool flipped,const ARRAY<ARRAY<VECTOR<int,2> > >& unconnected_polygonal_regions,
        ARRAY<ARRAY<ARRAY<int > > >& final_polygon_element_particles) const;
    void Determine_Duplicate_Tets_And_Duplicate_Particles();
    void Compute_Parents_And_Weights_For_Nodes(const int orig_intersection,const VECTOR<int,4>& vertices_for_tet,const VECTOR<int,4>& duplicate_particles_for_nodes,
        ARRAY<int>& new_parents,ARRAY<T>& parent_weights) const;
    void Duplicate_And_Merge_Elements();
    VECTOR<T,3> Compute_World_Space_Position_Of_Uncollapsed_Particle(const int uncollapsed_particle) const;
    void Create_Boundary_Surface();

    void Build_New_Cutting_Simplices_And_Intersection_Registry();
    void Add_Cuts_To_New_Cutting_Simplices(CUTTING_SIMPLICES<T,3>& new_cutting_simplices,HASHTABLE<int,int>& old_cut_to_new_cut) const;
    void Add_Boundaries_To_New_Cutting_Simplices(CUTTING_SIMPLICES<T,3>& new_cutting_simplices,HASHTABLE<VECTOR<int,3>,int>& volume_triangle_mesh_element_to_new_boundary,
        HASHTABLE<VECTOR<int,3>,VECTOR<int,3> >& new_sorted_to_original_boundary_nodes) const;
    void Add_New_Child_Simplices_And_Create_New_Polygon_Mesh(POLYGON_MESH& new_polygon_mesh,CUTTING_SIMPLICES<T,3>& new_cutting_simplices,ARRAY<ARRAY<int> >& simplices_per_new_tet,
        ARRAY<CUTTING_POLYGON>& new_cutting_polygons,ARRAY<ARRAY<ARRAY<int> > >& new_regions_per_tet,HASHTABLE<int,int>& old_cut_to_new_cut,
        HASHTABLE<VECTOR<int,3>,int>& volume_triangle_mesh_element_to_new_boundary,HASHTABLE<PAIR<int,int>,int>& old_child_simplex_to_new_child_simplex,
        HASHTABLE<int,int>& new_child_simplex_to_old_child_simplex,HASHTABLE<VECTOR<int,3>,VECTOR<int,3> >& new_sorted_to_original_boundary_nodes,
        ARRAY<ARRAY<int> >& new_cutting_polygons_per_cutting_simplex) const;
    void Build_New_Intersection_Registry(POLYGON_MESH& new_polygon_mesh,ARRAY<CUTTING_POLYGON>& new_cutting_polygons,HASHTABLE<int,int>& old_cut_to_new_cut,
        HASHTABLE<VECTOR<int,3>,int>& volume_triangle_mesh_element_to_new_boundary,HASHTABLE<PAIR<int,int>,int>& old_child_simplex_to_new_child_simplex,
        CUTTING_SIMPLICES<T,3>& new_cutting_simplices,INTERSECTION_REGISTRY<T,3>& new_intersection_registry,const HASHTABLE<int,int>& new_child_simplex_to_old_child_simplex,
        CUTTING_PARTICLES& new_cutting_particles) const;
    
    // HELPERS
    int Cutting_Polygon_To_Element(const int index) const;
    template<class T_ARRAY> int Add_Non_Tet_Node_Intersection_To_Registry(const T_ARRAY& simplices,const typename REBIND<T_ARRAY,VECTOR<T,2> >::TYPE& weights);
    template<class T_ARRAY> void Add_Non_Tet_Node_Intersection_To_Registry(const T_ARRAY& simplices,const typename REBIND<T_ARRAY,VECTOR<T,2> >::TYPE& weights,const int index);
    void Get_Particles_On_Simplices(const VECTOR<int,2>& simplices,ARRAY<int>& particles) const;
    VECTOR<T,2> Get_Node_Weights_From_Triangle(const int vertex,const VECTOR<int,3>& triangle_vertices) const;
    VECTOR<TV,3> Get_Face_Weights_From_Tet(const VECTOR<int,3>& face_vertices,const VECTOR<int,4>& tet_vertices) const;
    VECTOR<T,2> Get_Weight_From_Segment_On_Face(T bary,const int node_1,const int node_2,VECTOR<int,3> triangle) const;
    VECTOR<TV,3> Get_Triangle_Weights_From_Tet(const VECTOR<TV,3>& triangle_vertex_positions,const VECTOR<int,4>& tet_vertices) const;
    void Get_Simplex_Weights_For_Edge_Triangle_Intersection(const VECTOR<int,3>& simplices,const int triangle_array_index,const VECTOR<int,2>& shared_edges, VECTOR<VECTOR<T,2>,3>& all_weights) const;
    bool Intersects_And_Get_Simplex_Weights_For_Edge_Triangle_Intersection(const VECTOR<int,3>& simplices,const int triangle_array_index,const VECTOR<int,2>& shared_edges, VECTOR<VECTOR<T,2>,3>& all_weights) const;
    bool Intersects_And_Get_Simplex_Weights_For_Edge_Triangle_Intersection_Helper(const VECTOR<int,3>& simplices,const int triangle_array_index,const VECTOR<int,2>& shared_edges, VECTOR<VECTOR<T,2>,3>& all_weights, bool check_if_intersects) const;
    bool Segment_Reversed_In_Triangle(const VECTOR<int,3>& triangle_nodes,const VECTOR<int,2>& segment_nodes) const;
    bool Face_Reversed_In_Tetrahedron(const VECTOR<int,4>& tetrahedron_nodes,const VECTOR<int,3>& triangle_nodes,const GEOMETRY_PARTICLES<TV>& particles) const;
    void Get_Polygon_Edges(const int polygon_element_index,ARRAY<VECTOR<int,2> >& polygonal_segments) const;
    VECTOR<T,3> Get_Particle_Position_In_Local_Tet_Space(const int particle,const int child_simplex) const;
    VECTOR<T,3> Get_Child_Cutting_Simplex_Local_Normal(const int simplex) const;
    void Draw_Polygon(const int simplex,const bool flipped,const ARRAY<ARRAY<VECTOR<int,2> > >& unconnected_polygonal_regions) const;
//#####################################################################
};
}
#endif
