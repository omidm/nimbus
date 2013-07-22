//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Neil Molino, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VIRTUAL_NODE
//##################################################################### 
#ifndef __VIRTUAL_NODE__
#define __VIRTUAL_NODE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
namespace PhysBAM{

class VIRTUAL_NODE
{
public:
    int index;
    int corresponding_real_node;
    ARRAY<int> recipients;

    VIRTUAL_NODE()
    {}

    VIRTUAL_NODE(const int corresponding_real_node_input)
        :index(0),corresponding_real_node(corresponding_real_node_input)
    {}

    bool Is_Received_By(const int node) const
    {return recipients.Contains(node);}
    
//##################################################################### 
};
}
#endif
