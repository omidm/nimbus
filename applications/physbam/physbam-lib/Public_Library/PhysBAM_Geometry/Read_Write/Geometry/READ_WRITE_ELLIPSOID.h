//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_ELLIPSOID
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_ELLIPSOID__
#define __READ_WRITE_ELLIPSOID__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Basic_Geometry/ELLIPSOID.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<ELLIPSOID<T>,RW>
{
public:
    static void Read(std::istream& input,ELLIPSOID<T>& object)
    {Read_Binary<RW>(input,object.center,object.radii(1),object.radii(2),object.radii(3));}

    static void Write(std::ostream& output,const ELLIPSOID<T>& object)
    {Write_Binary<RW>(output,object.center,object.radii(1),object.radii(2),object.radii(3));}
};
}
#endif
#endif
