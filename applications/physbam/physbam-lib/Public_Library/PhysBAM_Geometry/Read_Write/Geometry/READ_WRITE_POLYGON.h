//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_POLYGON
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_POLYGON__
#define __READ_WRITE_POLYGON__

#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/POLYGON.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<POLYGON<TV>,RW>
{
public:
    static void Read(std::istream& input,POLYGON<TV>& object)
    {Read_Binary<RW>(input,object.X);}

    static void Write(std::ostream& output,const POLYGON<TV>& object)
    {Write_Binary<RW>(output,object.X);}
};
}
#endif
#endif
