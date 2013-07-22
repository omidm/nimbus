//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSTITUTIVE_MODELS_FORWARD 
//#####################################################################
#ifndef __CONSTITUTIVE_MODELS_FORWARD__
#define __CONSTITUTIVE_MODELS_FORWARD__

#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class TV,int d> class STRAIN_MEASURE;
template<class T> class STRAIN_MEASURE_HEXAHEDRONS;
template<class T,int d> class CONSTITUTIVE_MODEL;
template<class T,int d> class ISOTROPIC_CONSTITUTIVE_MODEL;
template<class T,int d> class ANISOTROPIC_CONSTITUTIVE_MODEL;
template<class T,int d> class PLASTICITY_MODEL;
template<class T,int d> class DIAGONALIZED_STRESS_DERIVATIVE;
template<class T,int d> class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE;

}
#endif
