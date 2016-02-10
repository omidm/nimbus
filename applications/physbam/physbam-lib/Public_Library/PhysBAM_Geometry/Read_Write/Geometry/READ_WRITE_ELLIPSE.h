//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_ELLIPSE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_ELLIPSE__
#define __READ_WRITE_ELLIPSE__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Basic_Geometry/ELLIPSE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<ELLIPSE<T>,RW>
{
public:
    static void Read(std::istream& input,ELLIPSE<T>& object)
    {Read_Binary<RW>(input,object.center,object.radii(1),object.radii(2));}

    static void Write(std::ostream& output,const ELLIPSE<T>& object)
    {Write_Binary<RW>(output,object.center,object.radii(1),object.radii(2));}
};
}
#endif
#endif
