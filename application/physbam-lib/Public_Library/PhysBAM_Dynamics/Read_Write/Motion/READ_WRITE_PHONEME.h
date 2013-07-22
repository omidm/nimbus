//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_PHONEME
//#####################################################################
#ifndef __READ_WRITE_PHONEME__
#define __READ_WRITE_PHONEME__

#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_ND.h>
#include <PhysBAM_Dynamics/Motion/PHONEME.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<PHONEME<T>,RW>
{
public:
    static void Read(std::istream& input,PHONEME<T>& object)
    {Read_Binary<RW>(input,object.name,object.previous_name,object.next_name,object.frame_length,object.time_length,object.controls);}

    static void Write(std::ostream& output,const PHONEME<T>& object)
    {Write_Binary<RW>(output,object.name,object.previous_name,object.next_name,object.frame_length,object.time_length,object.controls);}
};
}

#endif
