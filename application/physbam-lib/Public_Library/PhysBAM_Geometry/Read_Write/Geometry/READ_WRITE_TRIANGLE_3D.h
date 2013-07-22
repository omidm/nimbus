//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_TRIANGLE_3D
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_TRIANGLE_3D__
#define __READ_WRITE_TRIANGLE_3D__

#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<TRIANGLE_3D<T>,RW>
{
public:
    static void Read(std::istream& input,TRIANGLE_3D<T>& object)
    {Read_Binary<RW>(input,object.x1,object.x2,object.x3);}

    static void Write(std::ostream& output,const TRIANGLE_3D<T>& object)
    {Write_Binary<RW>(output,object.x1,object.x2,object.x3);}
};
}
#endif
#endif
