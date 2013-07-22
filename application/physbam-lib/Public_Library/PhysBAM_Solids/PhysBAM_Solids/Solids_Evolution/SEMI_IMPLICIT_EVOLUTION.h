//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEMI_IMPLICIT_EVOLUTION
//#####################################################################
#ifndef __SEMI_IMPLICIT_EVOLUTION__
#define __SEMI_IMPLICIT_EVOLUTION__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
namespace PhysBAM{
#if 0
template<class T_input>
class SEMI_IMPLICIT_EVOLUTION:public SOLIDS_EVOLUTION<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;
protected:
    typedef SOLIDS_EVOLUTION<TV> BASE;
    using BASE::X_save;
public:
    using BASE::deformable_object;using BASE::Euler_Step_Position;using BASE::Save_Position;using BASE::solids_parameters;

    SEMI_IMPLICIT_EVOLUTION(SOLIDS_PARAMETERS<TV>& solids_parameters)
        :SOLIDS_EVOLUTION<TV>(solids_parameters)
    {}

    bool Use_CFL() const PHYSBAM_OVERRIDE
    {return true;}

//#####################################################################
    void Advance_One_Time_Step(const T dt,const T time) PHYSBAM_OVERRIDE;
//#####################################################################
};
#endif
}
#endif
