//#####################################################################
// Copyright 2004-2006, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MELTING_PARAMETERS
//#####################################################################
#ifndef __MELTING_PARAMETERS__
#define __MELTING_PARAMETERS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

template<class TV,int d> class LEVELSET_SIMPLICIAL_OBJECT;

template<class TV,int d>
class MELTING_PARAMETERS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    enum BODY_TYPE {DEFORMABLE,RIGID};
    ARRAY<BODY_TYPE> body_type;
    ARRAY<int> body_index;
    int maximum_depth;
    bool refine_near_interface;
    T interface_refinement_multiplier;
    bool refine_for_high_deformation;
    T successive_refinement_multiplier;
    T expansive_refinement_threshold,compressive_refinement_threshold;
    bool write_overlay_levelset;
    bool use_constant_melting_speed;
    T constant_melting_speed;
    ARRAY<LEVELSET_SIMPLICIAL_OBJECT<TV,d>*> levelsets;
    ARRAY<FRAME<TV> > rigid_body_grid_frames;
    ARRAY<ARRAY<T>*> temperature;
    ARRAY<ARRAY<T>*> reaction;
    T rigid_body_coupling_density;
    bool reinitialize_after_melting_step;

    MELTING_PARAMETERS()
        :maximum_depth(1),refine_near_interface(false),interface_refinement_multiplier(3),refine_for_high_deformation(false),successive_refinement_multiplier(2),
        expansive_refinement_threshold((T).2),compressive_refinement_threshold((T).2),write_overlay_levelset(true),use_constant_melting_speed(false),constant_melting_speed(0),
        rigid_body_coupling_density(1),reinitialize_after_melting_step(true)
    {}

    virtual ~MELTING_PARAMETERS()
    {}

//#####################################################################
};
}
#endif
