//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_EMBEDDED_MATERIAL_SURFACE
//#####################################################################
#ifndef __READ_WRITE_EMBEDDED_MATERIAL_SURFACE__
#define __READ_WRITE_EMBEDDED_MATERIAL_SURFACE__

#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_EMBEDDED_TRIANGULATED_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Topology/READ_WRITE_SIMPLEX_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDED_TRIANGULATED_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_MATERIAL_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Read_Write/Fracture/READ_WRITE_EMBEDDING.h>
namespace PhysBAM{

template<class RW,class TV,int d>
class Read_Write<EMBEDDED_MATERIAL_SURFACE<TV,d>,RW>:public Read_Write<STRUCTURE<TV>,RW>
{
public:
    static void Read_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
    {EMBEDDED_MATERIAL_SURFACE<TV,d>& object=dynamic_cast<EMBEDDED_MATERIAL_SURFACE<TV,d>&>(structure_object);
    object.material_surface.Clean_Memory();Read_Binary<RW>(input,object.material_surface_mesh,object.previously_perturbed);
    Read_Write<typename EMBEDDED_MATERIAL_SURFACE<TV,d>::T_EMBEDDED_OBJECT,RW>::Read_Structure(input,object.embedded_object);}

    static void Read_Structure_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
    {Read_Helper(input,structure_object);}

    static void Write_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
    {const EMBEDDED_MATERIAL_SURFACE<TV,d>& object=dynamic_cast<const EMBEDDED_MATERIAL_SURFACE<TV,d>&>(structure_object);
    Write_Binary<RW>(output,object.material_surface_mesh,object.previously_perturbed);
    Read_Write<typename EMBEDDED_MATERIAL_SURFACE<TV,d>::T_EMBEDDED_OBJECT,RW>::Write_Structure(output,object.embedded_object);}

    static void Write_Structure_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
    {Write_Helper(output,structure_object);}
};
}
#endif
