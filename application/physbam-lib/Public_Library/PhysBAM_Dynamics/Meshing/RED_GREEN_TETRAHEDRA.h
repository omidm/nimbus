//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RED_GREEN_TETRAHEDRA
//##################################################################### 
#ifndef __RED_GREEN_TETRAHEDRA__
#define __RED_GREEN_TETRAHEDRA__

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
namespace PhysBAM{

template<class T>
class RED_GREEN_TETRAHEDRA:public NONCOPYABLE
{
public:
    static const int number_of_red_children=8;
    TETRAHEDRALIZED_VOLUME<T>& object;
    ARRAY<VECTOR<int,2> > leaf_levels_and_indices; // (level,index) of each (leaf) tet in object
    ARRAY<ARRAY<int>*> leaf_number; // for each tet on each level it's leaf number (0 if it's not a leaf)
    ARRAY<TETRAHEDRON_MESH*> meshes;
    ARRAY<ARRAY<int>*> parent;
    ARRAY<ARRAY<VECTOR<int,number_of_red_children> >*> children; // length 8 but padded with zeros at the end for irregular refinement
    SEGMENT_MESH segment_mesh; // all of the segments in all the levels of the mesh
    ARRAY<int> segment_midpoints; // for each segment, its midpoint or 0
    ARRAY<int>* segment_index_from_midpoint_index; // inverse map of segment_midpoints
    ARRAY<VECTOR<T,3> >* rest_position;
    HASHTABLE<int,bool> refinement_target_segments;
private:
    ARRAY<VECTOR<int,2> > stack; // level and index of tetrahedra to be processed
    ARRAY<ARRAY<int>*> index_in_stack; // 0 if not in stack, or -1 if can guarantee it never needs to be in stack

public:
    RED_GREEN_TETRAHEDRA(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume_input);
    RED_GREEN_TETRAHEDRA(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume_input,ARRAY<VECTOR<T,3> >& rest_position_input);
    ~RED_GREEN_TETRAHEDRA();


    bool Leaf(const int level,const int tet) const
    {return !(*children(level))(tet)(1);}

    bool Regularly_Refined(const int level,const int tet) const
    {return (*children(level))(tet)(8) != 0;}

    bool Red(const int level,const int tet) const
    {return level == 1 || Regularly_Refined(level-1,(*parent(level))(tet));}

//#####################################################################
    void Initialize_Segment_Index_From_Midpoint_Index(); // TODO: check that this is correct
    void Print() const;
    void Initialize();
    void Clean_Memory();
    void Refine_Simplex_List(const ARRAY<int>& tetrahedron_list);
    void Subdivide_Segment_List(const ARRAY<VECTOR<int,2> >& segment_list);
    bool Find_Edge_And_Test_If_Green(const int node1,const int node2,int* segment_output=0,int* level_output=0,int* tet_output=0);
    void Coarsen_Green_Refinements(TETRAHEDRON_MESH& final_mesh,ARRAY<int>& t_junctions,ARRAY<VECTOR<int,2> >& t_junction_parents) const;
    void Coarsen_Complete_Refinements_Of_Subset(TETRAHEDRON_MESH& final_mesh,ARRAY<bool>& keep_tet_flag,ARRAY<int>& t_junctions,ARRAY<VECTOR<int,2> >& t_junction_parents,
        bool allow_red_coarsening=false,ARRAY<bool>* node_is_uncoarsenable=0) const;
    void Remove_Simplex_List(const ARRAY<int>& tetrahedron_list,ARRAY<HASHTABLE<int,int> >* level_simplex_maps=0); // TODO: fill this in
private:
    void Resolve_Stack();
    void Refine_If_Necessary(int const level,const int tet);
    void Regularly_Refine_Tet(const int level,const int tet);
    void Get_Existing_Subindices(const int level,const int tet,ARRAY<int>& midpoints,ARRAY<int>& subedges);
    void Ensure_Level_Exists(const int level);
    void Delete_Children(const int level,const int tet,ARRAY<int>& deleted_tet_indices,ARRAY<int>& deleted_edge_indices);
    int Add_Midpoint(const int segment,const int level,const int tet);
    void Add_Segment(ARRAY<int>& free_edge_indices,const int node1,const int node2);
    void Add_Tetrahedron(ARRAY<int>& free_tet_indices,const int level,const int i,const int j,const int k,const int l,const int parent_index);
    void Rebuild_Object();
    void Unmark_Parents(ARRAY<ARRAY<bool> >& tetrahedron_kept,int level,int tet) const;
    void Unmark_Children(ARRAY<ARRAY<bool> >& tetrahedron_kept,int level,int tet) const;
//#####################################################################
};
}
#endif

