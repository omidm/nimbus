//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_SEGMENT_2D
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_SEGMENT_2D__
#define __READ_WRITE_SEGMENT_2D__

#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<SEGMENT_2D<T>,RW>
{
public:
    static void Read(std::istream& input,SEGMENT_2D<T>& object)
    {Read_Binary<RW>(input,object.x1,object.x2);}

    static void Write(std::ostream& output,const SEGMENT_2D<T>& object)
    {Write_Binary<RW>(output,object.x1,object.x2);}
};
}
#endif
#endif
