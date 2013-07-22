//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_INTERSECTION_REGISTRY
//#####################################################################
#ifndef __READ_WRITE_INTERSECTION_REGISTRY__
#define __READ_WRITE_INTERSECTION_REGISTRY__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Dynamics/Fracture/INTERSECTION_REGISTRY.h>
namespace PhysBAM{

template<class RW,class T,int d>
class Read_Write<INTERSECTION_REGISTRY<T,d>,RW>
{
public:
    static void Read(std::istream& input,INTERSECTION_REGISTRY<T,d>& object)
    {Read_Binary<RW>(input,object.cutting_simplices,object.intersections_on_simplex,object.simplices_on_intersection,object.simplex_weights_on_intersection,object.index_for_last_old_intersection);}

    static void Write(std::ostream& output,const INTERSECTION_REGISTRY<T,d>& object)
    {Write_Binary<RW>(output,object.cutting_simplices,object.intersections_on_simplex,object.simplices_on_intersection,object.simplex_weights_on_intersection,object.index_for_last_old_intersection);}
};
}
#endif
