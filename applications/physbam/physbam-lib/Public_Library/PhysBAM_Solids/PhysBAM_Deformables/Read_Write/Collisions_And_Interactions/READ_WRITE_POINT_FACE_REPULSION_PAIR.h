//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_POINT_FACE_REPULSION_PAIR
//#####################################################################
#ifndef __READ_WRITE_POINT_FACE_REPULSION_PAIR__
#define __READ_WRITE_POINT_FACE_REPULSION_PAIR__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<POINT_FACE_REPULSION_PAIR<TV>,RW>
{
public:
    static void Read(std::istream& input,POINT_FACE_REPULSION_PAIR<TV>& object)
    {Read_Binary<RW>(input,object.nodes,object.distance,object.weights,object.normal);}

    static void Write(std::ostream& output,const POINT_FACE_REPULSION_PAIR<TV>& object)
    {Write_Binary<RW>(output,object.nodes,object.distance,object.weights,object.normal);}
};
}
#endif
