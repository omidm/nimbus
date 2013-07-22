//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_MOTION_SEQUENCE
//#####################################################################
#ifndef __READ_WRITE_MOTION_SEQUENCE__
#define __READ_WRITE_MOTION_SEQUENCE__

#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_FRAME.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Dynamics/Motion/MOTION_SEQUENCE.h>
namespace PhysBAM{

template<class RW,class T,class T2>
class Read_Write<MOTION_SEQUENCE<T,T2>,RW>
{
public:
    static void Read(std::istream& input,MOTION_SEQUENCE<T,T2>& object)
    {Read_Binary<RW>(input,object.time_grid,object.trajectories,object.valid,object.names);object.Update_Name_Lookup();}

    static void Write(std::ostream& output,const MOTION_SEQUENCE<T,T2>& object)
    {Write_Binary<RW>(output,object.time_grid,object.trajectories,object.valid,object.names);}
};
}

#endif
