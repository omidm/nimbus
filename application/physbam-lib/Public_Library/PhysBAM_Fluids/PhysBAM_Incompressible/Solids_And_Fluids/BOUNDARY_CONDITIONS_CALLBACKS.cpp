//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Solids_And_Fluids/BOUNDARY_CONDITIONS_CALLBACKS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BOUNDARY_CONDITIONS_CALLBACKS<TV>::
BOUNDARY_CONDITIONS_CALLBACKS()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BOUNDARY_CONDITIONS_CALLBACKS<TV>::
~BOUNDARY_CONDITIONS_CALLBACKS()
{
}
//#####################################################################
// Function Mark_Outside
//#####################################################################
template<class TV> void BOUNDARY_CONDITIONS_CALLBACKS<TV>::
Mark_Outside(ARRAY<bool,FACE_INDEX<d> >& outside)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Mark_Outside
//#####################################################################
template<class TV> void BOUNDARY_CONDITIONS_CALLBACKS<TV>::
Mark_Outside(ARRAY<bool,TV_INT>& outside)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Get_Boundary_Along_Ray
//#####################################################################
template<class TV> typename BOUNDARY_CONDITIONS_CALLBACKS<TV>::RAY_TYPE BOUNDARY_CONDITIONS_CALLBACKS<TV>::
Get_Boundary_Along_Ray(const FACE_INDEX<d>& f1,const FACE_INDEX<d>& f2,T& theta,T& value)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return unused;
}
//#####################################################################
// Function Get_Boundary_Along_Ray
//#####################################################################
template<class TV> typename BOUNDARY_CONDITIONS_CALLBACKS<TV>::RAY_TYPE BOUNDARY_CONDITIONS_CALLBACKS<TV>::
Get_Boundary_Along_Ray(const TV_INT& c1,const TV_INT& c2,T& theta,T& value)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return unused;
}
template class BOUNDARY_CONDITIONS_CALLBACKS<VECTOR<float,1> >;
template class BOUNDARY_CONDITIONS_CALLBACKS<VECTOR<float,2> >;
template class BOUNDARY_CONDITIONS_CALLBACKS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BOUNDARY_CONDITIONS_CALLBACKS<VECTOR<double,1> >;
template class BOUNDARY_CONDITIONS_CALLBACKS<VECTOR<double,2> >;
template class BOUNDARY_CONDITIONS_CALLBACKS<VECTOR<double,3> >;
#endif
