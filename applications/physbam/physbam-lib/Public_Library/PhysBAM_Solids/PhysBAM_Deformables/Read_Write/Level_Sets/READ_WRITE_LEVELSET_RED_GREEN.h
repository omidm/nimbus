#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_LEVELSET_RED_GREEN
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_LEVELSET_RED_GREEN__
#define __READ_WRITE_LEVELSET_RED_GREEN__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Red_Green/READ_WRITE_RED_GREEN_GRID_2D.h>
#include <PhysBAM_Geometry/Read_Write/Red_Green/READ_WRITE_RED_GREEN_GRID_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Level_Sets/LEVELSET_RED_GREEN.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<LEVELSET_RED_GREEN<TV>,RW>
{
public:
    static void Read(std::istream& input,LEVELSET_RED_GREEN<TV>& object)
    {Read_Binary<RW>(input,object.grid,object.phi);}

    static void Write(std::ostream& output,const LEVELSET_RED_GREEN<TV>& object)
    {Write_Binary<RW>(output,object.grid,object.phi);}
};
}
#endif
#endif
#endif
