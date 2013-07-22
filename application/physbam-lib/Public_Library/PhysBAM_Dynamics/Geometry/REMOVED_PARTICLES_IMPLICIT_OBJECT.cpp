//#####################################################################
// Copyright 2006, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REMOVED_PARTICLES_IMPLICIT_OBJECT
//#####################################################################
#include <PhysBAM_Dynamics/Geometry/REMOVED_PARTICLES_IMPLICIT_OBJECT.h>
using namespace PhysBAM;
template REMOVED_PARTICLES_IMPLICIT_OBJECT<VECTOR<float,3> >::REMOVED_PARTICLES_IMPLICIT_OBJECT(REMOVED_PARTICLES_PROCESSING<float>* particle_processing_input,LEVELSET_IMPLICIT_OBJECT<VECTOR<float,3> >* levelset_implicit);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template REMOVED_PARTICLES_IMPLICIT_OBJECT<VECTOR<double,3> >::REMOVED_PARTICLES_IMPLICIT_OBJECT(REMOVED_PARTICLES_PROCESSING<double>* particle_processing_input,LEVELSET_IMPLICIT_OBJECT<VECTOR<double,3> >* levelset_implicit);
#endif
