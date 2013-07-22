//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Craig Schroeder, Tamar Shinar, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class JOINT_MESH
//#####################################################################
#ifndef __JOINT_MESH__
#define __JOINT_MESH__

#include <PhysBAM_Tools/Data_Structures/DYNAMIC_LIST.h>
#include <PhysBAM_Tools/Data_Structures/UNDIRECTED_GRAPH.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_ID.h>
namespace PhysBAM{

template<class TV> class JOINT;

template<class TV>
class JOINT_MESH
{
    typedef typename TV::SCALAR T;
    DYNAMIC_LIST<JOINT<TV>,JOINT_ID> dynamic_list;
public:
    const ARRAY<JOINT<TV>*>& joints; // reference to array for traditional naming
    UNDIRECTED_GRAPH<int,JOINT_ID> undirected_graph;

    JOINT_MESH();
    ~JOINT_MESH();

    void Add_Articulation(const int parent_id,const int child_id,JOINT<TV>* joint)
    {undirected_graph.Add_Edge(parent_id,child_id,Add_Joint(joint));}
    
    void Add_Unowned_Articulation(const int parent_id,const int child_id,JOINT<TV>* joint)
    {dynamic_list.Add_Element(joint);undirected_graph.Add_Edge(parent_id,child_id,joint->id_number);}

    void Remove_Articulation(const JOINT_ID joint_id)
    {undirected_graph.Remove_Edge(joint_id);Remove_Joint(joint_id);}
    
    void Remove_Unowned_Articulation(const JOINT_ID joint_id)
    {undirected_graph.Remove_Edge(joint_id);Remove_Joint_Without_Destroying(joint_id);}

    void Deactivate_Articulation(const JOINT_ID joint_id)
    {undirected_graph.Remove_Edge(joint_id);Remove_Joint_Without_Destroying(joint_id);}
    
    void Reactivate_Articulation(const int parent_id,const int child_id,JOINT<TV>* joint)
    {dynamic_list.Reactivate_Element(joint,joint->id_number);undirected_graph.Add_Edge(parent_id,child_id,joint->id_number);}
    
    void Swap_Articulation(const int parent_id,const int child_id,JOINT<TV>* joint,const JOINT_ID joint_id)
    {undirected_graph.Remove_Edge(joint_id);dynamic_list.Swap_Elements(joint,joint->id_number,joint_id);undirected_graph.Add_Edge(parent_id,child_id,joint->id_number);}
    
    int Parent_Id(JOINT_ID joint_id) const
    {return undirected_graph.Edges(joint_id).x;}

    int Child_Id(JOINT_ID joint_id) const
    {return undirected_graph.Edges(joint_id).y;}

    JOINT<TV>* operator()(const JOINT_ID id)
    {return dynamic_list.Element(id);}

    const JOINT<TV>* operator()(const JOINT_ID id) const
    {return dynamic_list.Element(id);}

    JOINT_ID Size() const
    {return dynamic_list.Size();}

    bool Is_Active(JOINT_ID &joint_id) const
    {return (joint_id<=undirected_graph.Last_Edge() && undirected_graph.Edges(joint_id).x!=0 && undirected_graph.Edges(joint_id).y!=0);}
    
    JOINT_ID Add_Previously_Used_Joint(JOINT<TV>* joint)
    {return dynamic_list.Reactivate_Element(joint,joint->id_number);}

    void Remove_Joint_Without_Destroying(const JOINT_ID id)
    {dynamic_list.Deactivate_Element(id,false);}

    void Remove_All()
    {undirected_graph.Reset();dynamic_list.Remove_All();}

    void Clean_Memory()
    {dynamic_list.Clean_Memory();}

    int Joint_Index_From_Id(JOINT_ID id)
    {return dynamic_list.Element_Index(id);}

    void Remove_Joint(JOINT_ID joint_id)
    {dynamic_list.Remove_Element(joint_id);}

    const ARRAY<JOINT_ID>& Active_Elements() const
    {return dynamic_list.Active_Elements();}

//#####################################################################
    JOINT_ID Add_Joint(JOINT<TV>* new_joint_input);
    void Read(TYPED_ISTREAM& input,const std::string& output_directory,const int frame);
    void Write(TYPED_OSTREAM& output,const std::string& output_directory,const int frame) const;
//#####################################################################
};
}
#endif
