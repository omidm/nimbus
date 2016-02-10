//#####################################################################
// Copyright 2007, Kevin Der, Geoffrey Irving, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __CUTTING_GEOMETRY__
#define __CUTTING_GEOMETRY__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/STACK.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Fracture/CUTTING_SIMPLEX.h>
#include <PhysBAM_Geometry/Topology/POLYGON_MESH.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/ROBUST_SIMPLEX_INTERACTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Dynamics/Fracture/CUTTING_PARTICLES.h>
#include <PhysBAM_Dynamics/Fracture/CUTTING_POLYGON.h>
#include <PhysBAM_Dynamics/Fracture/CUTTING_SIMPLICES.h>
#include <PhysBAM_Dynamics/Fracture/INTERSECTION_REGISTRY.h>
namespace PhysBAM{

template<class TV,int d_input>
class CUTTING_GEOMETRY
{
public:
    enum WORKAROUND {d_ambient=TV::dimension,d_embed=d_input,d_cut=d_input-1};

    typedef TV VECTOR_T;
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d_embed>::OBJECT T_EMBEDDING_OBJECT;
    typedef typename BASIC_SIMPLEX_POLICY<TV,d_embed>::SIMPLEX T_EMBEDDING_SIMPLEX;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d_cut>::OBJECT T_CUT_OBJECT;
    typedef typename BASIC_SIMPLEX_POLICY<TV,d_cut>::SIMPLEX T_CUT_SIMPLEX;
    typedef typename MESH_POLICY<d_embed>::MESH T_EMBEDDING_MESH;
    typedef typename MESH_POLICY<d_cut>::MESH T_CUT_MESH;typedef T_CUT_MESH T_FACE_MESH;
    typedef CUTTING_SIMPLICES<T,d_embed> T_CUTTING_SIMPLICES;
    typedef CUTTING_SIMPLEX<T,d_embed> T_CUTTING_SIMPLEX;
    typedef INTERSECTION_REGISTRY<T,d_embed> T_INTERSECTION_REGISTRY;
    typedef VECTOR<T,d_embed> T_EMBEDDING_WEIGHTS; // an embedded point's weights (i.e. 2 for triangle )
    typedef VECTOR<T,d_cut> T_CUT_WEIGHTS;
    typedef VECTOR<T_EMBEDDING_WEIGHTS,d_cut+1> T_CUTTING_SIMPLEX_WEIGHTS;

    STATIC_ASSERT(d_input==2 || d_input==3);
    struct TWO_DIMENSIONAL{};struct THREE_DIMENSIONAL{};
    typedef typename IF<INTS_EQUAL<d_input,2>::value,TWO_DIMENSIONAL,THREE_DIMENSIONAL>::TYPE DIMENSION;

    // CURRENT STATE
    T_EMBEDDING_OBJECT* current_embedding;

    // VALID ONLY FOR THE CURRENT CUT ITERATION
    T_CUT_OBJECT cutting; // TODO: make into an instance owned by this class
    T_EMBEDDING_OBJECT* next_embedding;
    T_CUT_MESH final_duplicated_boundary_mesh;
    ARRAY<ARRAY<int> > final_parents_per_new_particle;
    ARRAY<ARRAY<T> > final_parent_weights_per_new_particle;
    ARRAY<int> old_particle_per_new_collapsed_particle;
    ARRAY<ARRAY<int> > new_collapsed_embedding_simplex_particle_per_old_embedding_simplex_particle;

public: // would be private but helper functions need to access
    // replaced at the end of current iteration
    T_CUTTING_SIMPLICES* cutting_simplices; // TODO: try to make an instance
    T_INTERSECTION_REGISTRY* intersection_registry; // TODO: try to make an instance
    POLYGON_MESH polygon_mesh;
    ARRAY<ARRAY<int> > simplices_per_current_embedding_simplex;
    ARRAY<CUTTING_POLYGON> current_cutting_polygons;
    ARRAY<ARRAY<ARRAY<int> > > regions_per_embedding_simplex;
    HASHTABLE<VECTOR<int,d_cut+1>,VECTOR<int,d_cut+1> > sorted_to_original_boundary_nodes;
    CUTTING_PARTICLES cutting_particles;

    // overwritten at end of current iteration and used in the next
    ARRAY<ARRAY<int> > cutting_polygons_per_cutting_simplex;

    // must be refreshed when first used in the current iteration
    ARRAY<STACK<int> > simplex_stack_per_current_embedding_simplex;
    ARRAY<ARRAY<int> > polygons_per_element;
#if 0
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
#endif

    // persistent
    T intersection_thickness;
    bool verbose;
    int first_new_cut_element;
public:

    CUTTING_GEOMETRY(const bool verbose_input=false);    
    virtual ~CUTTING_GEOMETRY();

    void Set_Intersection_Thickness(const T intersection_thickness_input=0)  // Used only for hierarchy intersections; exact intersections performed at the simplex level
    {intersection_thickness=intersection_thickness_input;}

    template<int d2> void
    Get_Particles_On_Simplices(const VECTOR<int,d2>& simplices,ARRAY<int>& particles) const
    {bool is_all_cutting_simplices=false;
    for(int i=1;i<=simplices.m;i++) if(cutting_simplices->simplices(simplices(i)).type!=T_CUTTING_SIMPLEX::LOCAL_CUT_FACE){is_all_cutting_simplices=false;break;}
    intersection_registry->Intersection_List(simplices,particles);
    // also grab particles on cuts if all cutting simplices
    if(is_all_cutting_simplices){VECTOR<int,d2> converted_simplices;ARRAY<int> particles_on_cuts;
        for(int i=1;i<=simplices.m;i++){
            converted_simplices(i)=cutting_simplices->simplices(simplices(i)).parent;assert((cutting_simplices->simplices(converted_simplices(i)).type==T_CUTTING_SIMPLEX::GLOBAL_CUT_FACE));}
        int embedding_simplex=cutting_simplices->simplices(simplices(1)).element_owner;const VECTOR<int,d_embed+1>& tet_nodes=current_embedding->mesh.elements(embedding_simplex);
        intersection_registry->Intersection_List_For_Cuts(converted_simplices,tet_nodes,particles_on_cuts); // only keep particle if simplices contained in our tet
        particles.Append_Elements(particles_on_cuts);}}

    template<class T_ARRAY> int Register_Cut_Intersection(const T_ARRAY& simplices,const typename REBIND<T_ARRAY,T_CUT_WEIGHTS>::TYPE& weights,const int index_input)
    {assert(simplices.Size()==weights.Size());
    int index=index_input;
    if(!index){
        index=intersection_registry->Number_Of_Intersections()+1; // TODO: check if this is the right index
        cutting_particles.Add_Intersection_Id(index);}        
    intersection_registry->Register_Intersection(simplices,weights,index);
    if(verbose){std::stringstream ss;ss<<"Intersection: "<<index<<std::endl;
        for(int i=1;i<=simplices.m;i++) ss<<"simplex "<<simplices(i)<<" has weights "<<weights(i)<<std::endl;ss<<std::endl;
    LOG::filecout(ss.str());}
    return index;}

//#####################################################################
    void Initialize_Original_Embedding(const T_EMBEDDING_OBJECT& original_embedding);
    void Cut_Material(T_EMBEDDING_OBJECT& next_embedding_input);
protected:
    virtual void Split_Existing_Polygon_Edges()=0;
    virtual void Split_Existing_Polygons()=0;
    virtual void Intersect_Simplex_With_Old_Simplices_In_Embedding(const int embedding_simplex,const int new_simplex)=0;
    void Fill_Intersection_Registry_With_Cuts();
    T_CUTTING_SIMPLEX_WEIGHTS Get_Face_Weights_From_Embedding_Simplex(const VECTOR<int,d_cut+1>& face_nodes,const VECTOR<int,d_embed+1>& embedding_simplex_nodes) const;
public:
    T_CUT_WEIGHTS Get_Node_Weights_From_Face(const int node,const VECTOR<int,d_cut+1>& face_nodes) const;
    //T_CUTTING_SIMPLEX_WEIGHTS Get_Cutting_Simplex_Weights_From_Embedding_Simplex(const VECTOR<TV,d_cut+1>& cutting_simplex_X,const VECTOR<TV,d_embed+1>& embedding_simplex_X) const;
 protected:
    void Get_Polygon_Edges(const int polygon_element_index,ARRAY<VECTOR<int,2> >& polygonal_segments) const;
//#####################################################################
};
}
#endif
