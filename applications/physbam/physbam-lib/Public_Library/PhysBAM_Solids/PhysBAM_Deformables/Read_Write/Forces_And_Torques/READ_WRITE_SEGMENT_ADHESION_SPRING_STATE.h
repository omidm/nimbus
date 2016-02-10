//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_SEGMENT_ADHESION_SPRING_STATE
//#####################################################################
#ifndef __READ_WRITE_SEGMENT_ADHESION_SPRING_STATE__
#define __READ_WRITE_SEGMENT_ADHESION_SPRING_STATE__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_ADHESION.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<SEGMENT_ADHESION_SPRING_STATE<TV>,RW>
{
public:
    static void Read(std::istream& input,SEGMENT_ADHESION_SPRING_STATE<TV>& object)
    {Read_Binary<RW>(input,object.nodes,object.weights,object.normal,object.distance,object.damping);}

    static void Write(std::ostream& output,const SEGMENT_ADHESION_SPRING_STATE<TV>& object)
    {Write_Binary<RW>(output,object.nodes,object.weights,object.normal,object.distance,object.damping);}
};
}
#endif
