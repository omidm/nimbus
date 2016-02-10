//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_EMBEDDING
//#####################################################################
#ifndef __READ_WRITE_EMBEDDING__
#define __READ_WRITE_EMBEDDING__

#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_STRUCTURE.h>
#include <PhysBAM_Geometry/Read_Write/Topology/READ_WRITE_SIMPLEX_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDING.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<EMBEDDING<TV>,RW>:public Read_Write<STRUCTURE<TV>,RW>
{
public:
    static void Read_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
    {EMBEDDING<TV>& object=dynamic_cast<EMBEDDING<TV>&>(structure_object);
    object.material_surface.Clean_Memory();Read_Binary<RW>(input,object.material_surface_mesh);}

    static void Read_Structure_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
    {Read_Helper(input,structure_object);}

    static void Write_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
    {const EMBEDDING<TV>& object=dynamic_cast<const EMBEDDING<TV>&>(structure_object);
    Write_Binary<RW>(output,object.material_surface_mesh);}

    static void Write_Structure_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
    {Write_Helper(output,structure_object);}
};
}
#endif
