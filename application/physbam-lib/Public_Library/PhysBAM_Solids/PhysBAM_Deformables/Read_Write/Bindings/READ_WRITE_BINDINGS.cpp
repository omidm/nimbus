//#####################################################################
// Copyright 2009, Elliot English
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_BINDINGS
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING_DYNAMIC.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/PARTICLE_BINDING.h>
namespace PhysBAM{

void Initialize_Bindings()
{
#define HELPER(T,d) \
    BINDING_REGISTRY<VECTOR<T,d> >::Register<PARTICLE_BINDING<VECTOR<T,d> > >(); \
    BINDING_REGISTRY<VECTOR<T,d> >::Register<LINEAR_BINDING<VECTOR<T,d>,2> >(); \
    BINDING_REGISTRY<VECTOR<T,d> >::Register<LINEAR_BINDING<VECTOR<T,d>,3> >(); \
    BINDING_REGISTRY<VECTOR<T,d> >::Register<LINEAR_BINDING<VECTOR<T,d>,4> >(); \
    BINDING_REGISTRY<VECTOR<T,d> >::Register<LINEAR_BINDING_DYNAMIC<VECTOR<T,d> > >();

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
