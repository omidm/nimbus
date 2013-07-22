//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BACKWARD_EULER_EVOLUTION
//#####################################################################
#ifndef __BACKWARD_EULER_EVOLUTION__
#define __BACKWARD_EULER_EVOLUTION__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
namespace PhysBAM{

template<class TV>
class BACKWARD_EULER_EVOLUTION:public SOLIDS_EVOLUTION<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef SOLIDS_EVOLUTION<TV> BASE;
    using BASE::solid_body_collection;using BASE::solids_parameters;using BASE::Euler_Step_Position;
private:
    using BASE::F_full;using BASE::R_full;using BASE::S_full;using BASE::B_full;using BASE::rigid_F_full;using BASE::rigid_R_full;using BASE::rigid_B_full;using BASE::rigid_S_full;
    mutable ARRAY<TV> dV_full;
public:

    BACKWARD_EULER_EVOLUTION(SOLIDS_PARAMETERS<TV>& solids_parameters_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input);

    bool Use_CFL() const PHYSBAM_OVERRIDE
    {return true;}

    void Initialize_Rigid_Bodies(const T frame_rate, const bool restart) PHYSBAM_OVERRIDE
    {PHYSBAM_NOT_IMPLEMENTED();}

//#####################################################################
    void Advance_One_Time_Step_Position(const T dt,const T time,const bool solids) PHYSBAM_OVERRIDE {PHYSBAM_NOT_IMPLEMENTED();}
    void Advance_One_Time_Step_Velocity(const T dt,const T time,const bool solids) PHYSBAM_OVERRIDE;
private:
    void One_Newton_Step_Backward_Euler(const T dt,const T time,ARRAY_VIEW<const TV> V_save,ARRAY<TV>& dV_full);
//#####################################################################
};
}
#endif
