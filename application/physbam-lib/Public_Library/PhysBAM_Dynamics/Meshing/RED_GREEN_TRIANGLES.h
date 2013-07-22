//#####################################################################
// Copyright 2006-2008, Avi Robinson-Mosher, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RED_GREEN_TRIANGLES
//##################################################################### 
#ifndef __RED_GREEN_TRIANGLES__
#define __RED_GREEN_TRIANGLES__

#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/OPERATION_HASH.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{

template<class TV>
class RED_GREEN_TRIANGLES:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::TRIANGULATED_OBJECT T_TRIANGULATED_OBJECT;
public:
    static const int number_of_red_children=4;
    T_TRIANGULATED_OBJECT& object;
    ARRAY<VECTOR<int,2> > leaf_levels_and_indices; // (level,index) of each (leaf) triangle in triangulated_object
    ARRAY<ARRAY<int>*> leaf_number; // for each triangle on each level it's leaf number (0 if it's not a leaf)
    ARRAY<TRIANGLE_MESH*> meshes;
    ARRAY<ARRAY<int>*> parent;
    ARRAY<ARRAY<VECTOR<int,number_of_red_children> >*> children; // length 4 but padded with zeros at the end for irregular refinement
    SEGMENT_MESH segment_mesh; // all of the segments in all the levels of the mesh
    ARRAY<int> segment_midpoints; // for each segment, its midpoint or 0
    ARRAY<TV>* rest_position;
    ARRAY<int>* segment_index_from_midpoint_index; // inverse map of segment_midpoints
    HASHTABLE<VECTOR<int,2>,int>* free_segment_midpoints;
    ARRAY<ARRAY<VECTOR<int,3> > > element_edges; // edge indices are indices into segment_mesh
private:
    ARRAY<VECTOR<int,2> > stack; // level and index of triangles to be processed
    ARRAY<ARRAY<int>*> index_in_stack; // 0 if not in stack, or -1 if can guarantee it never needs to be in stack

public:
    RED_GREEN_TRIANGLES(T_TRIANGULATED_OBJECT& triangulated_object_input)
        :object(triangulated_object_input),rest_position(0),segment_index_from_midpoint_index(0),free_segment_midpoints(0)
    {Initialize();}

    RED_GREEN_TRIANGLES(T_TRIANGULATED_OBJECT& triangulated_object_input,ARRAY<TV>& rest_position_input)
        :object(triangulated_object_input),rest_position(&rest_position_input)
    {Initialize();}

    ~RED_GREEN_TRIANGLES()
    {Clean_Memory();delete free_segment_midpoints;}

    bool Leaf(const int level,const int tri) const
    {return !(*children(level))(tri)(1);}

    bool Regularly_Refined(const int level,const int tri) const
    {return (*children(level))(tri)(4) != 0;}

    bool Red(const int level,const int tri) const
    {return level == 1 || Regularly_Refined(level-1,(*parent(level))(tri));}

    void Initialize_Segment_Index_From_Midpoint_Index()
    {if(segment_index_from_midpoint_index) delete segment_index_from_midpoint_index;
    segment_index_from_midpoint_index=new ARRAY<int>(object.particles.array_collection->Size());
    for(int s=1;s<=segment_midpoints.m;s++) if(segment_midpoints(s)) (*segment_index_from_midpoint_index)(segment_midpoints(s))=s;}

    void Add_Free_Segment_Midpoint(const VECTOR<int,2>& endpoints,const int midpoint)
    {if(!free_segment_midpoints) free_segment_midpoints=new HASHTABLE<VECTOR<int,2>,int>;
    free_segment_midpoints->Insert(endpoints.Sorted(),midpoint);}
    
//#####################################################################
    void Initialize();
    void Clean_Memory();
    template<class T_ARRAY> void Refine_Simplex_List(const T_ARRAY& triangle_list);
    void Unrefined_Parents(const int node,ARRAY<int>& parents,ARRAY<T>& weights) const;
    void Remove_Simplex_List(const ARRAY<int>& triangle_list,ARRAY<HASHTABLE<int,int> >* level_simplex_maps=0);
private:
    void Resolve_Stack();
    void Refine_If_Necessary(int const level,const int tri);
    void Regularly_Refine_Triangle(const int level,const int tri);
    void Get_Existing_Subindices(const int level,const int tri,ARRAY<int>& midpoints,ARRAY<int>& subedges);
    void Ensure_Level_Exists(const int level);
    void Delete_Children(const int level,const int tri,ARRAY<int>& deleted_tri_indices,ARRAY<int>& deleted_edge_indices);
    void Add_Midpoints_If_Not_Already_There(const int level,const int tri);
    int Add_Midpoint(const int segment,const int level,const int tri);
    void Add_Segment(ARRAY<int>& free_edge_indices,const int node1,const int node2);
    void Add_Triangle(ARRAY<int>& free_triangle_indices,const int level,const int i,const int j,const int k,const int parent_index);
    void Rebuild_Object();
    void Print() const;
//#####################################################################
};
}
#endif

