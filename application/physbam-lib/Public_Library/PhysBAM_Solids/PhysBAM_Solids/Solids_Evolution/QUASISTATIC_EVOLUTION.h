//#####################################################################
// Copyright 2006-2007, Ron Fedkiw, Geoffrey Irving, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUASISTATIC_EVOLUTION
//#####################################################################
#ifndef __QUASISTATIC_EVOLUTION__
#define __QUASISTATIC_EVOLUTION__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
namespace PhysBAM{

template<class TV>
class QUASISTATIC_EVOLUTION:public SOLIDS_EVOLUTION<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef SOLIDS_EVOLUTION<TV> BASE;
    using BASE::solid_body_collection;using BASE::solids_parameters;
private:
    using BASE::F_full;using BASE::R_full;using BASE::S_full;using BASE::B_full;
    using BASE::rigid_F_full;using BASE::rigid_R_full;using BASE::rigid_S_full;using BASE::rigid_B_full;
    ARRAY<TV> dX_full;
public:
    bool balance_external_forces_only;

    QUASISTATIC_EVOLUTION(SOLIDS_PARAMETERS<TV>& solids_parameters_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input);

    bool Use_CFL() const PHYSBAM_OVERRIDE
    {return false;}

    void Initialize_Rigid_Bodies(const T frame_rate, const bool restart) PHYSBAM_OVERRIDE
    {PHYSBAM_NOT_IMPLEMENTED();}

//#####################################################################
    void One_Newton_Step_Toward_Steady_State(const T time,ARRAY<TV>& dX_full);
    void Advance_One_Time_Step_Position(const T dt,const T time,const bool solids) PHYSBAM_OVERRIDE;
    void Advance_One_Time_Step_Velocity(const T dt,const T time,const bool solids) PHYSBAM_OVERRIDE {}
//#####################################################################
};
}
#endif
