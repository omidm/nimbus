//#####################################################################
// Copyright 2009, Elliot English
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_RIGID_BODY_BINDING
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Bindings/RIGID_BODY_BINDING.h>
namespace PhysBAM{

void Initialize_Rigid_Body_Binding()
{
#define HELPER(T,d) \
    BINDING_REGISTRY<VECTOR<T,d> >::Register<RIGID_BODY_BINDING<VECTOR<T,d> > >();

    HELPER(float,1)
    HELPER(float,2)
    HELPER(float,3)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    HELPER(double,1)
    HELPER(double,2)
    HELPER(double,3)
#endif
}
}
