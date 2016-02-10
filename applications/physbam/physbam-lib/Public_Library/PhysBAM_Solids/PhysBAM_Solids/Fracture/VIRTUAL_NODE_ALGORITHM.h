//#####################################################################
// Copyright 2003-2006, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Neil Molino, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VIRTUAL_NODE_ALGORITHM
//#####################################################################
#ifndef __VIRTUAL_NODE_ALGORITHM__
#define __VIRTUAL_NODE_ALGORITHM__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/EMBEDDING_POLICY.h>
namespace PhysBAM{
namespace VIRTUAL_NODE_ALGORITHM{

template<class TV,int d> void
Rebuild_Embedded_Object(EMBEDDED_OBJECT<TV,d>& embedded_object,ARRAY<int>& map_to_old_particles,ARRAY<int>& map_to_old_embedded_particles,
    ARRAY<int>& map_to_old_simplices,const bool verbose=true);

};
}
#endif
