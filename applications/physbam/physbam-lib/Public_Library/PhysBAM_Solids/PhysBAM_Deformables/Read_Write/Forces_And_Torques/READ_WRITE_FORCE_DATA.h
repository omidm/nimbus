//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_FORCE_DATA
//#####################################################################
#ifndef __READ_WRITE_FORCE_DATA__
#define __READ_WRITE_FORCE_DATA__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<FORCE_DATA<TV>,RW>
{
public:
    static void Read(std::istream& input,FORCE_DATA<TV>& object)
    {Read_Binary<RW>(input,object.name,object.state,object.first_action_point,object.second_action_point);}

    static void Write(std::ostream& output,const FORCE_DATA<TV>& object)
    {Write_Binary<RW>(output,object.name,object.state,object.first_action_point,object.second_action_point);}
};
}
#endif
