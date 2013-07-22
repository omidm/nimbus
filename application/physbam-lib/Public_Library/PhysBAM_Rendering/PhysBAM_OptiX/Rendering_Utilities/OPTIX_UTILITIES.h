//#####################################################################
// Copyright 2011, Valeria Nikolaenko, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// OPTIX_UTILITIES class
//#####################################################################
#ifdef USE_OPTIX
#ifndef __OPTIX_UTILITIES__
#define __OPTIX_UTILITIES__
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <optix_world.h>
#include <optixu/optixpp_namespace.h>

namespace PhysBAM{
using namespace optix;

class OPTIX_UTILITIES 
{
public:
    template<class T> static float3 Get_Float3(VECTOR<T,3> vector) 
    {
        float3 float3_vector={(float)vector.x,(float)vector.y,(float)vector.z};
        return float3_vector;
    }
};
}
#endif /* _OPTIX_UTILITIES_H_ */
#endif
