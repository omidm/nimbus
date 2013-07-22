//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE
//#####################################################################
#ifndef __READ_WRITE_EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE__
#define __READ_WRITE_EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE__

#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Read_Write/Fracture/READ_WRITE_EMBEDDED_MATERIAL_SURFACE.h>
namespace PhysBAM{

template<class RW,class T>
class Read_Write<EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>,RW>:public Read_Write<EMBEDDED_MATERIAL_SURFACE<VECTOR<T,3>,3>,RW>
{
};
}
#endif
