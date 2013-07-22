//#####################################################################
// Copyright 2003-2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLASTICITY_CONTROL_3D
//##################################################################### 
#ifndef __PLASTICITY_CONTROL_3D__
#define __PLASTICITY_CONTROL_3D__

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/PLASTICITY_MODEL.h>
namespace PhysBAM{

template<class T>
class PLASTICITY_CONTROL_3D:public PLASTICITY_MODEL<T,3>
{
    typedef VECTOR<T,3> TV;
public:
    typedef PLASTICITY_MODEL<T,3> BASE;
    using BASE::Fp_inverse;

    ARRAY<SYMMETRIC_MATRIX<T,3> > log_Fp,log_Fp_goal;
    T yield_ratio;
private:
    T sqr_log_yield_ratio;

public:
    PLASTICITY_CONTROL_3D(const STRAIN_MEASURE<TV,3>& strain_measure,const ARRAY<VECTOR<T,3> >& X_goal,const T yield_ratio_input)
        :PLASTICITY_MODEL<T,3>(strain_measure.Dm_inverse.m),yield_ratio(yield_ratio_input)
    {
        assert(yield_ratio>0);sqr_log_yield_ratio=sqr(log(yield_ratio));
        log_Fp.Exact_Resize(strain_measure.Dm_inverse.m);
        log_Fp_goal.Exact_Resize(strain_measure.Dm_inverse.m);
        for(int t=1;t<=strain_measure.Dm_inverse.m;t++){
            MATRIX<T,3> F=strain_measure.Ds(X_goal.array,t)*strain_measure.Dm_inverse(t);
            if(F.Determinant()<=0) PHYSBAM_FATAL_ERROR("Cannot control towards inverting elements");
            MATRIX<T,3> U,V;DIAGONAL_MATRIX<T,3> Fp_goal_hat;F.Fast_Singular_Value_Decomposition(U,Fp_goal_hat,V);
            log_Fp_goal(t)=SYMMETRIC_MATRIX<T,3>::Conjugate(V,log(Fp_goal_hat));}
    }
    
    void Print_Extreme_Deformations()
    {T tension=-1e10,compression=1e10;
    for(int t=1;t<=log_Fp_goal.m;t++){
        DIAGONAL_MATRIX<T,3> eigenvalues=log_Fp_goal(t).Fast_Eigenvalues();
        tension=max(tension,eigenvalues.x11);compression=min(compression,eigenvalues.x22);}
    {std::stringstream ss;ss<<"Maximum goal expansion = "<<exp(tension)<<", compression = "<<exp(compression)<<std::endl;LOG::filecout(ss.str());}}
    
    bool Project_Fe(const DIAGONAL_MATRIX<T,3>& Fe_trial,DIAGONAL_MATRIX<T,3>& Fe_project) const PHYSBAM_OVERRIDE
    {DIAGONAL_MATRIX<T,3> Fe_log=log(Fe_trial.Max((T)1e-5));T sqr_norm=Fe_log.Frobenius_Norm_Squared();
    if(sqr_norm<=sqr_log_yield_ratio)return false;Fe_log*=sqrt(sqr_log_yield_ratio/sqr_norm);
    Fe_project=exp(Fe_log);return true;}
    
    void Project_Fp(const int tetrahedron,const MATRIX<T,3>& Fp_trial) PHYSBAM_OVERRIDE
    {MATRIX<T,3> U,V;DIAGONAL_MATRIX<T,3> Fp_trial_hat;Fp_trial.Fast_Singular_Value_Decomposition(U,Fp_trial_hat,V);
    SYMMETRIC_MATRIX<T,3> log_Fp_trial=SYMMETRIC_MATRIX<T,3>::Conjugate(V,log(Fp_trial_hat.Max((T)1e-5)));
    SYMMETRIC_MATRIX<T,3> log_Fp_trial_minus_log_Fp=log_Fp_trial-log_Fp(tetrahedron);
    T goal_dot_trial=SYMMETRIC_MATRIX<T,3>::Inner_Product(log_Fp_goal(tetrahedron)-log_Fp(tetrahedron),log_Fp_trial_minus_log_Fp);if(goal_dot_trial<=0)return;
    T t=goal_dot_trial/log_Fp_trial_minus_log_Fp.Frobenius_Norm_Squared();
    if(t<1)log_Fp_trial=t*log_Fp_trial_minus_log_Fp+log_Fp(tetrahedron);
    log_Fp(tetrahedron)=log_Fp_trial;Fp_inverse(tetrahedron)=exp(-log_Fp_trial);}
    
    void Read_State(TYPED_ISTREAM& input) PHYSBAM_OVERRIDE
    {PLASTICITY_MODEL<T,3>::Read_State(input);
    for(int t=1;t<=Fp_inverse.m;t++)log_Fp(t)=-log(Fp_inverse(t));}
    
//#####################################################################
};
}
#endif
