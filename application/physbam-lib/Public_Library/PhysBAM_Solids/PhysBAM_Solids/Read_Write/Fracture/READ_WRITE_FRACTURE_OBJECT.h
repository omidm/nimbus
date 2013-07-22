//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_FRACTURE_OBJECT
//#####################################################################
#ifndef __READ_WRITE_FRACTURE_OBJECT__
#define __READ_WRITE_FRACTURE_OBJECT__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fracture/FRACTURE_OBJECT.h>
namespace PhysBAM{

template<class RW,class T,int d>
class Read_Write<FRACTURE_OBJECT<T,d>,RW>
{
public:
    static void Read(std::istream& input,FRACTURE_OBJECT<T,d>& object)
    {Read_Binary<RW>(input,object.number_of_fracture_initiations,object.corresponding_node_in_reference,object.corresponding_simplex_in_reference);}

    static void Write(std::ostream& output,const FRACTURE_OBJECT<T,d>& object)
    {Write_Binary<RW>(output,object.number_of_fracture_initiations,object.corresponding_node_in_reference,object.corresponding_simplex_in_reference);}
};
}
#endif
