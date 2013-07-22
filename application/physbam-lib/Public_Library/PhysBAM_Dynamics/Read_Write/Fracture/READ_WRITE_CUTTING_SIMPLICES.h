//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_CUTTING_SIMPLICES
//#####################################################################
#ifndef __READ_WRITE_CUTTING_SIMPLICES__
#define __READ_WRITE_CUTTING_SIMPLICES__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Dynamics/Fracture/CUTTING_SIMPLICES.h>
namespace PhysBAM{

template<class RW,class T,int d>
class Read_Write<CUTTING_SIMPLICES<T,d>,RW>
{
public:
    static void Read(std::istream& input,CUTTING_SIMPLICES<T,d>& object)
    {Read_Binary<RW>(input,object.simplices,object.index_for_last_old_cutting_simplex);}

    static void Write(std::ostream& output,const CUTTING_SIMPLICES<T,d>& object)
    {Write_Binary<RW>(output,object.simplices,object.index_for_last_old_cutting_simplex);}
};
}
#endif
