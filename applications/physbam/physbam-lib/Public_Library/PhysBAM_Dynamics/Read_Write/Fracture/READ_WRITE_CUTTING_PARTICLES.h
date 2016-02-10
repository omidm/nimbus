//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_CUTTING_PARTICLES
//#####################################################################
#ifndef __READ_WRITE_CUTTING_PARTICLES__
#define __READ_WRITE_CUTTING_PARTICLES__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Dynamics/Fracture/CUTTING_PARTICLES.h>
namespace PhysBAM{

template<class RW>
class Read_Write<CUTTING_PARTICLES,RW>
{
public:
    static void Read(std::istream& input,CUTTING_PARTICLES& object)
    {Read_Binary<RW>(input,object.tet_node_indices,object.intersection_indices,object.particle_ids_types,object.intersection_to_particle_id,object.tet_node_to_particle_id);}

    static void Write(std::ostream& output,const CUTTING_PARTICLES& object)
    {Write_Binary<RW>(output,object.tet_node_indices,object.intersection_indices,object.particle_ids_types,object.intersection_to_particle_id,object.tet_node_to_particle_id);}
};
}
#endif
