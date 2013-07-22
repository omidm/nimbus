//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_TO_SOLID_INTERPOLATION
//#####################################################################
#ifndef __FLUID_TO_SOLID_INTERPOLATION__
#define __FLUID_TO_SOLID_INTERPOLATION__
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/COUPLING_CONSTRAINT_ID.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION_BASE.h>

namespace PhysBAM{

template<class TV>
class FLUID_TO_SOLID_INTERPOLATION:public FLUID_TO_SOLID_INTERPOLATION_BASE<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    using FLUID_TO_SOLID_INTERPOLATION_BASE<TV>::index_map;

    struct ENTRY
    {
        T w;
        int i;
    };

public:
    T max_dist;
    ARRAY<VECTOR<ARRAY<ENTRY>,TV::m> > entries;
    ARRAY<int> coupled_particles;
    const PARTICLES<TV>& particles;

    FLUID_TO_SOLID_INTERPOLATION(const COLLISION_AWARE_INDEX_MAP<TV>& map,const PARTICLES<TV>& particles_input);
    virtual ~FLUID_TO_SOLID_INTERPOLATION();

//#####################################################################
    void Compute(const int ghost_cells) PHYSBAM_OVERRIDE;
    void Compute_Weights(const TV& X,int axis,ARRAY<ENTRY>& array);
    void Times_Add(const VECTOR_ND<T>& fluid_velocity,GENERALIZED_VELOCITY<TV>& solid_velocity) const PHYSBAM_OVERRIDE;
    void Transpose_Times_Add(const GENERALIZED_VELOCITY<TV>& solid_force,VECTOR_ND<T>& fluid_force) const PHYSBAM_OVERRIDE;
    void Print_Each_Matrix(int n,int fluid_faces,GENERALIZED_VELOCITY<TV>& G) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
