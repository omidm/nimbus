//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_TO_SOLID_INTERPOLATION_PHI
//#####################################################################
#ifndef __FLUID_TO_SOLID_INTERPOLATION_PHI__
#define __FLUID_TO_SOLID_INTERPOLATION_PHI__
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/FLUID_TO_SOLID_INTERPOLATION_CUT.h>

namespace PhysBAM{

template<class TV>
class FLUID_TO_SOLID_INTERPOLATION_PHI:public FLUID_TO_SOLID_INTERPOLATION_CUT<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef FLUID_TO_SOLID_INTERPOLATION_CUT<TV> BASE;
    using BASE::index_map;using BASE::cut_cells;using BASE::fluid_mass;using BASE::gradient;using BASE::curve;using BASE::psi_N;using BASE::outside_fluid;using BASE::Remove_Degeneracy;

public:
    const ARRAY<T,TV_INT>& phi;
    int cut_order;

    FLUID_TO_SOLID_INTERPOLATION_PHI(const COLLISION_AWARE_INDEX_MAP<TV>& map,const ARRAY<T,TV_INT>& phi,SEGMENTED_CURVE_2D<T>& curve,T density_input);
    virtual ~FLUID_TO_SOLID_INTERPOLATION_PHI();

//#####################################################################
    void Setup_Mesh();
    void Setup_Before_Compute(ARRAY<bool,TV_INT>& outside_fluid_input,const ARRAY<bool,FACE_INDEX<TV::m> >& psi_N_input);
    T Linear_Average(const FACE_INDEX<TV::m>& face);
    T Node_Average(const TV_INT& cell);
    T Compute_Cut3(const FACE_INDEX<TV::m>& face,T phi0,T phi1);
    T Compute_Cut4(const FACE_INDEX<TV::m>& face,T phi0,T phi1);
//#####################################################################
};
}
#endif
