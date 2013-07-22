//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_JOINT
//#####################################################################
#ifndef __READ_WRITE_JOINT__
#define __READ_WRITE_JOINT__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<JOINT<TV>,RW>
{
public:
    static void Read(std::istream& input,JOINT<TV>& object)
    {Read_Binary<RW>(input,object.id_number,object.frame_pj,object.frame_jp,object.frame_cj,object.frame_jc,object.J,object.J_inverse,object.name);}

    static void Write(std::ostream& output,const JOINT<TV>& object)
    {Write_Binary<RW>(output,object.id_number,object.frame_pj,object.frame_jp,object.frame_cj,object.frame_jc,object.J,object.J_inverse,object.name);}
};
}
#endif
