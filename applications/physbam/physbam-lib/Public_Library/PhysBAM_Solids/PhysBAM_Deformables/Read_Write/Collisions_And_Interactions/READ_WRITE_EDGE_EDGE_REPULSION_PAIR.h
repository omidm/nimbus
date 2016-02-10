//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_EDGE_EDGE_REPULSION_PAIR
//#####################################################################
#ifndef __READ_WRITE_EDGE_EDGE_REPULSION_PAIR__
#define __READ_WRITE_EDGE_EDGE_REPULSION_PAIR__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<EDGE_EDGE_REPULSION_PAIR<VECTOR<T,2> >,RW>
{
public:
    static void Read(std::istream& input,EDGE_EDGE_REPULSION_PAIR<VECTOR<T,2> >& object)
    {Read_Binary<RW>(input,object.nodes,object.distance,object.normal);}

    static void Write(std::ostream& output,const EDGE_EDGE_REPULSION_PAIR<VECTOR<T,2> >& object)
    {Write_Binary<RW>(output,object.nodes,object.distance,object.normal);}
};
template<class RW,class T>
class Read_Write<EDGE_EDGE_REPULSION_PAIR<VECTOR<T,3> >,RW>
{
public:
    static void Read(std::istream& input,EDGE_EDGE_REPULSION_PAIR<VECTOR<T,3> >& object)
    {Read_Binary<RW>(input,object.nodes,object.distance,object.weights,object.normal);}

    static void Write(std::ostream& output,const EDGE_EDGE_REPULSION_PAIR<VECTOR<T,3> >& object)
    {Write_Binary<RW>(output,object.nodes,object.distance,object.weights,object.normal);}
};
}
#endif
