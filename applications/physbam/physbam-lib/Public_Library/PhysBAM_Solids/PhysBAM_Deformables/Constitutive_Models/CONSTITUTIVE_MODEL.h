//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Igor Neverov, Andrew Selle, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSTITUTIVE_MODEL
//##################################################################### 
#ifndef __CONSTITUTIVE_MODEL__
#define __CONSTITUTIVE_MODEL__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/CONSTITUTIVE_MODELS_FORWARD.h>
namespace PhysBAM{

template<class T,int d> class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE;
template<class T,int d>
class CONSTITUTIVE_MODEL:public NONCOPYABLE
{
    typedef VECTOR<T,d> TV;
public:
    bool enforce_definiteness;
    T constant_lambda,constant_mu; // Lame coefficients (used by almost all derived models)
    T constant_alpha,constant_beta; // isotropic damping parameters (used by all current derived models)
    ARRAY<T> lambda,mu; // spatially varying Lame coefficients
    ARRAY<T> alpha,beta; // spatially varying damping parameters

private:
    CONSTITUTIVE_MODEL();

    // all constitutive models should derive from one of these
    friend class ISOTROPIC_CONSTITUTIVE_MODEL<T,d>;
    friend class ANISOTROPIC_CONSTITUTIVE_MODEL<T,d>;
public:

    virtual ~CONSTITUTIVE_MODEL();

    virtual T Maximum_Elastic_Stiffness(const int simplex) const; // for elastic CFL computation
    virtual T Maximum_Damping_Stiffness(const int simplex) const; // for damping CFL computation

//#####################################################################
    virtual MATRIX<T,d> P_From_Strain_Rate(const DIAGONAL_MATRIX<T,d>& F,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const=0;
    virtual void Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,d>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d>& dPi_dF,const int simplex) const;
    virtual int P_From_Strain_Rate_Forces_Size() const;
    virtual void P_From_Strain_Rate_First_Half(const DIAGONAL_MATRIX<T,d>& F,ARRAY_VIEW<T> aggregate,const MATRIX<T,d>& F_dot,const T scale,const int simplex) const;
    virtual MATRIX<T,d> P_From_Strain_Rate_Second_Half(const DIAGONAL_MATRIX<T,d>& F,const ARRAY_VIEW<const T> aggregate,const T scale,const int simplex) const;
//#####################################################################
};
}
#endif
