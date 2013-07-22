//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Neil Molino, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VIRTUAL_NODES
//##################################################################### 
#ifndef __VIRTUAL_NODES__
#define __VIRTUAL_NODES__

#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/VIRTUAL_NODE.h>
namespace PhysBAM{

class VIRTUAL_NODES
{
public:
    ARRAY<VIRTUAL_NODE> virtual_nodes;
    // TODO: make replicas a hash table
    ARRAY<ARRAY<int> > replicas; // indexes into virtual_nodes

    VIRTUAL_NODES()
    {}

    VIRTUAL_NODE& operator()(const int i)
    {return virtual_nodes(i);}

    int Add_Virtual_Node(const int corresponding_real_node)
    {return virtual_nodes.Append(VIRTUAL_NODE(corresponding_real_node));}

    void Initialize_Replicas()
    {int max_real_node=0;for(int i=1;i<=virtual_nodes.m;i++) max_real_node=max(max_real_node,virtual_nodes(i).corresponding_real_node);
    replicas.Remove_All();replicas.Resize(max_real_node);
    for(int i=1;i<=virtual_nodes.m;i++) replicas(virtual_nodes(i).corresponding_real_node).Append(i);}

    int Donor_Node(const int recipient,const int real_node) const
    {if(real_node<=replicas.m) for(int k=1;k<=replicas(real_node).m;k++){
        const VIRTUAL_NODE& virtual_node=virtual_nodes(replicas(real_node)(k));
        if(virtual_node.Is_Received_By(recipient)) return virtual_node.index;}
    return real_node;}

//##################################################################### 
};
}
#endif
