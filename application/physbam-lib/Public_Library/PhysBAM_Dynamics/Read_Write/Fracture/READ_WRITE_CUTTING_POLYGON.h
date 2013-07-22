//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_CUTTING_POLYGON
//#####################################################################
#ifndef __READ_WRITE_CUTTING_POLYGON__
#define __READ_WRITE_CUTTING_POLYGON__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Dynamics/Fracture/CUTTING_POLYGON.h>
namespace PhysBAM{

template<class RW>
class Read_Write<CUTTING_POLYGON,RW>
{
public:
    static void Read(std::istream& input,CUTTING_POLYGON& object)
    {Read_Binary<RW>(input,object.polygon_index,object.simplex_owner,object.flipped,object.polygon_type);}

    static void Write(std::ostream& output,const CUTTING_POLYGON& object)
    {Write_Binary<RW>(output,object.polygon_index,object.simplex_owner,object.flipped,object.polygon_type);}
};
}
#endif
