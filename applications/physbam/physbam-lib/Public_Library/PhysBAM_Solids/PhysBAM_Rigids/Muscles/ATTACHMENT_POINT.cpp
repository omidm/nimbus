//#####################################################################
// Copyright 2005-2007, Kevin Der, Ron Fedkiw, Eran Guendelman, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ATTACHMENT_POINT
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/ATTACHMENT_POINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> ATTACHMENT_POINT<TV>::
ATTACHMENT_POINT(RIGID_BODY<TV>& rigid_body_input,TV object_space_position_input)
    :rigid_body(rigid_body_input),object_space_position(object_space_position_input)
{
}
//#####################################################################
// Function Embedded_Position
//#####################################################################
template<class TV> TV ATTACHMENT_POINT<TV>::
Embedded_Position() const
{
    return rigid_body.World_Space_Point(object_space_position);
}
//#####################################################################
// Function Embedded_Velocity
//#####################################################################
template<class TV> TV ATTACHMENT_POINT<TV>::
Embedded_Velocity() const
{
    return rigid_body.Pointwise_Object_Velocity(Embedded_Position());
}
//#####################################################################
// Function Apply_Impulse
//#####################################################################
template<class TV> void ATTACHMENT_POINT<TV>::
Apply_Impulse(const TV& impulse)
{
    rigid_body.Apply_Impulse_To_Body(Embedded_Position(),impulse);
}
template class ATTACHMENT_POINT<VECTOR<float,1> >;
template class ATTACHMENT_POINT<VECTOR<float,2> >;
template class ATTACHMENT_POINT<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ATTACHMENT_POINT<VECTOR<double,1> >;
template class ATTACHMENT_POINT<VECTOR<double,2> >;
template class ATTACHMENT_POINT<VECTOR<double,3> >;
#endif

