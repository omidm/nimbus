//#####################################################################
// Copyright 2005, Eftychios Sifakis
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMPLEX_MUSCLE_3D
//##################################################################### 
#ifndef __COMPLEX_MUSCLE_3D__
#define __COMPLEX_MUSCLE_3D__

#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/TRANSVERSE_ISOTROPY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
namespace PhysBAM{

template<class T>
class COMPLEX_MUSCLE_3D:public TRANSVERSE_ISOTROPY_3D<T>
{
    typedef VECTOR<T,3> TV;
public:
    typedef TRANSVERSE_ISOTROPY_3D<T> BASE;
    using BASE::constant_mu;using BASE::constant_lambda;using BASE::constant_alpha;using BASE::constant_beta;using BASE::failure_threshold;using BASE::Hessian_Index;

    T G_1_muscle,G_2_muscle,K_muscle,P_1,P_2,P_3,P_4,lambda_star_muscle,sigma_max,lambda_ofl;
    T G_1_tendon,G_2_tendon,K_tendon,L_1,L_2,L_3,L_4,lambda_star_tendon;
    ARRAY<bool>& is_tendon;
    T& activation;

    COMPLEX_MUSCLE_3D(ARRAY<bool>& is_tendon_input,T& activation_input,const T youngs_modulus=3e6,const T poissons_ratio=.475,const T Rayleigh_coefficient=.05,
        const T failure_threshold_input=.25)
        :TRANSVERSE_ISOTROPY_3D<T>(failure_threshold_input),is_tendon(is_tendon_input),activation(activation_input)
    {
        // Damping parameters
        assert(poissons_ratio>-1&&poissons_ratio<.5);
        constant_lambda=youngs_modulus*poissons_ratio/((1+poissons_ratio)*(1-2*poissons_ratio));
        constant_mu=youngs_modulus/(2*(1+poissons_ratio));
        constant_alpha=Rayleigh_coefficient*constant_lambda;
        constant_beta=Rayleigh_coefficient*constant_mu; 
        // Muscle parameters
        G_1_muscle=(T)5e2;
        G_2_muscle=(T)5e2;
        K_muscle=(T)1e7;
        P_1=(T)0.05;
        P_2=(T)6.6;
        lambda_star_muscle=(T)1.4;
        sigma_max=(T)3e5;
        lambda_ofl=(T)1.4;
        P_3=P_1*P_2*exp(P_2*(lambda_star_muscle/lambda_ofl-(T)1));
        P_4=P_1*(P_3/(P_1*P_2)-(T)1)-P_3*lambda_star_muscle/lambda_ofl;
        // Aponeurosis and external fascia parameters
        G_1_tendon=(T)5e4;
        G_2_tendon=(T)5e4;
        K_tendon=(T)1e8;
        L_1=(T)2.7e6;
        L_2=(T)46.4;
        lambda_star_tendon=(T)1.03;
        L_3=L_1*L_2*exp(L_2*(lambda_star_tendon-(T)1));
        L_4=L_3/L_2-L_1-L_3*lambda_star_tendon;
    }

    inline T Sqr_Acosh_Robust_Derivative(const T x) const
    {if(x<(T)1.001)return (T)2; else return (T)2*log(x+sqrt(sqr(x)-(T)1))/sqrt(sqr(x)-(T)1);}
    
    inline T Sqr_Acosh_Robust_Second_Derivative(const T x) const
    {if(x<(T)1.001)return (T)-two_thirds; else return ((T)2-x*Sqr_Acosh_Robust_Derivative(x))/(sqr(x)-(T)1);}

    T Muscle_Tension(const T stretch) const
    {T strain=stretch/lambda_ofl-(T)1,strain_abs=abs(strain),active_tension=0,passive_tension=0;
    if(stretch>lambda_star_muscle)passive_tension=P_3*stretch/lambda_ofl+P_4;else if(stretch>lambda_ofl)passive_tension=P_1*(exp(P_2*strain)-(T)1);
    if(strain_abs<(T).4)active_tension=((T)1-(T)4*sqr(strain));else if(strain_abs<(T).6)active_tension=(T)9*sqr(strain_abs-(T).6);
    return sigma_max*(activation*active_tension+passive_tension);}
        
    T Tendon_Tension(const T stretch) const
    {if(stretch>lambda_star_tendon) return L_3+L_4/stretch;
    if(stretch>(T)1) return L_1*(exp(L_2*(stretch-(T)1))-(T)1)/stretch;
    return 0;}
    
    T Muscle_Tension_Derivative(const T stretch) const
    {T strain=stretch/lambda_ofl-(T)1,strain_abs=abs(strain),active_tension_derivative=0,passive_tension_derivative=0;
    if(stretch>lambda_star_muscle)passive_tension_derivative=P_3/lambda_ofl ;else if(stretch>lambda_ofl)passive_tension_derivative=P_1*P_2*(exp(P_2*strain)-1)/lambda_ofl;
    if(strain_abs<(T).4)active_tension_derivative=-(T)8*strain;else if(strain_abs<(T).6)active_tension_derivative=(T)18*(strain-sign(strain)*(T).6);
    return sigma_max*(activation*active_tension_derivative+passive_tension_derivative);}

    T Tendon_Tension_Derivative(const T stretch) const
    {if(stretch>lambda_star_tendon) return -L_4/sqr(stretch);
    if(stretch>(T)1) return L_1*((L_2*stretch-(T)1)*exp(L_2*(stretch-(T)1))+(T)1)/sqr(stretch);
    return 0;}
    
    void Energy_Gradient(VECTOR_ND<T>& energy_gradient,const VECTOR_ND<T>& invariants,const int tetrahedron_index=0) const PHYSBAM_OVERRIDE
    {T G_1,G_2,K,dW3_dlambda;
    T stretch=pow(invariants(3),-(T)one_sixth)*sqrt(invariants(4));
    if(is_tendon(tetrahedron_index)){G_1=G_1_tendon;G_2=G_2_tendon;K=K_tendon;dW3_dlambda=Tendon_Tension(stretch);}
    else{G_1=G_1_muscle;G_2=G_2_muscle;K=K_muscle;dW3_dlambda=Muscle_Tension(stretch);}
    energy_gradient=VECTOR_ND<T>(5);
    // Along-fiber shear
    energy_gradient(4)-=(T)2*G_1*invariants(5)/cube(invariants(4));
    energy_gradient(5)+=G_1/sqr(invariants(4));
    // Cross-fiber shear
    T omega=(T).5*(invariants(1)*invariants(4)-invariants(5))/sqrt(invariants(3)*invariants(4));
    T dW2_domega=G_2*Sqr_Acosh_Robust_Derivative(omega);
    T omega_1=(T).5*sqrt(invariants(4)/invariants(3));
    T omega_3=(T)-.5*omega/invariants(3);
    T omega_4=(T).25*(invariants(1)*invariants(4)+invariants(5))/sqrt(invariants(3)*cube(invariants(4)));
    T omega_5=(T)-.5/sqrt(invariants(3)*invariants(4));
    energy_gradient(1)+=dW2_domega*omega_1;energy_gradient(3)+=dW2_domega*omega_3;energy_gradient(4)+=dW2_domega*omega_4;energy_gradient(5)+=dW2_domega*omega_5;
    // Along-fiber stretch
    T lambda_3=(T)-one_sixth*stretch/invariants(3);
    T lambda_4=(T).5*stretch/invariants(4);
    energy_gradient(3)+=dW3_dlambda*lambda_3;
    energy_gradient(4)+=dW3_dlambda*lambda_4;
    // Volumetric component
    energy_gradient(3)+=(T).25*K*log(invariants(3))/invariants(3);}

    void Energy_Hessian(VECTOR_ND<T>& energy_hessian,const VECTOR_ND<T>& invariants,const int tetrahedron_index=0) const PHYSBAM_OVERRIDE
    {T G_1,G_2,K,dW3_dlambda,ddW3_ddlambda;
    T stretch=pow(invariants(3),-(T)one_sixth)*sqrt(invariants(4));
    if(is_tendon(tetrahedron_index)){G_1=G_1_tendon;G_2=G_2_tendon;K=K_tendon;dW3_dlambda=Tendon_Tension(stretch);ddW3_ddlambda=Tendon_Tension_Derivative(stretch);}
    else{G_1=G_1_muscle;G_2=G_2_muscle;K=K_muscle;dW3_dlambda=Muscle_Tension(stretch);ddW3_ddlambda=Muscle_Tension_Derivative(stretch);}
    energy_hessian=VECTOR_ND<T>(15);
    // Along-fiber shear
    energy_hessian(Hessian_Index(4,4))+=(T)6*G_1*invariants(5)/pow(invariants(3),4);
    energy_hessian(Hessian_Index(5,4))-=(T)2*G_1/cube(invariants(4));
    // Cross-fiber shear
    T omega=(T).5*(invariants(1)*invariants(4)-invariants(5))/sqrt(invariants(3)*invariants(4));
    T dW2_domega=G_2*Sqr_Acosh_Robust_Derivative(omega);
    T ddW2_ddomega=G_2*Sqr_Acosh_Robust_Second_Derivative(omega);
    T omega_1=(T).5*sqrt(invariants(4)/invariants(3));
    T omega_3=(T)-.5*omega/invariants(3);
    T omega_4=(T).25*invariants(1)/sqrt(invariants(3)*invariants(4));
    T omega_5=(T)-.5/sqrt(invariants(3)*invariants(4));
    T omega_31=(T)-.25*sqrt(invariants(4)/cube(invariants(3)));
    T omega_41=(T).25/sqrt(invariants(3)*invariants(4));
    T omega_33=(T).75*omega/sqr(invariants(3));
    T omega_43=(T)-.125*(invariants(1)*invariants(4)+invariants(5))/sqrt(cube(invariants(3)*invariants(4)));
    T omega_53=(T).25/sqrt(cube(invariants(3))*invariants(4));
    T omega_44=(T)-.125*(invariants(1)*invariants(4)+(T)3*invariants(5))/sqrt(invariants(3)*pow(invariants(4),5));
    T omega_54=(T).25/sqrt(invariants(3)*cube(invariants(4)));
    energy_hessian(Hessian_Index(1,1))+=ddW2_ddomega*omega_1*omega_1;
    energy_hessian(Hessian_Index(3,1))+=ddW2_ddomega*omega_3*omega_1+dW2_domega*omega_31;
    energy_hessian(Hessian_Index(4,1))+=ddW2_ddomega*omega_4*omega_1+dW2_domega*omega_41;
    energy_hessian(Hessian_Index(5,1))+=ddW2_ddomega*omega_5*omega_1;
    energy_hessian(Hessian_Index(3,3))+=ddW2_ddomega*omega_3*omega_3+dW2_domega*omega_33;
    energy_hessian(Hessian_Index(4,3))+=ddW2_ddomega*omega_4*omega_3+dW2_domega*omega_43;
    energy_hessian(Hessian_Index(5,3))+=ddW2_ddomega*omega_5*omega_3+dW2_domega*omega_53;
    energy_hessian(Hessian_Index(4,4))+=ddW2_ddomega*omega_4*omega_4+dW2_domega*omega_44;
    energy_hessian(Hessian_Index(5,4))+=ddW2_ddomega*omega_5*omega_4+dW2_domega*omega_54;
    energy_hessian(Hessian_Index(5,5))+=ddW2_ddomega*omega_5*omega_5;
    // Along-fiber stretch
    T lambda_3=(T)-one_sixth*stretch/invariants(3);
    T lambda_4=(T).5*stretch/invariants(4);
    T lambda_33=((T)7/(T)36)*stretch*sqr(invariants(3));
    T lambda_43=(T)-one_twelfth*stretch/(invariants(3)*invariants(4));
    T lambda_44=(T)-.25*stretch/sqr(invariants(4));
    energy_hessian(Hessian_Index(3,3))+=ddW3_ddlambda*lambda_3*lambda_3+dW3_dlambda*lambda_33;
    energy_hessian(Hessian_Index(4,3))+=ddW3_ddlambda*lambda_4*lambda_3+dW3_dlambda*lambda_43;
    energy_hessian(Hessian_Index(4,4))+=ddW3_ddlambda*lambda_4*lambda_4+dW3_dlambda*lambda_44;
    // Volumetric component
    energy_hessian(Hessian_Index(3,3))+=(T).25*K*((T)1-log(invariants(3)))/sqr(invariants(3));}

//#####################################################################
};

template<class T> FINITE_VOLUME<VECTOR<T,3>,3>*
Create_Complex_Muscle(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,const ARRAY<VECTOR<T,3> >& fiber_field,ARRAY<bool>& is_tendon,
    T& activation,const T youngs_modulus=(T)3e6,const T poissons_ratio=(T).475,const T Rayleigh_coefficient=(T).05,const T failure_threshold=(T).25,const bool verbose=true)
{
    COMPLEX_MUSCLE_3D<T>* constitutive_model=new COMPLEX_MUSCLE_3D<T>(is_tendon,activation,youngs_modulus,poissons_ratio,Rayleigh_coefficient,failure_threshold);
    FINITE_VOLUME<VECTOR<T,3>,3>* fvm=Create_Finite_Volume(tetrahedralized_volume,*constitutive_model);
    constitutive_model->Initialize_Fiber_Field_From_Current_State(fvm->strain_measure,fiber_field);
    return fvm;
}

template<class T> FINITE_VOLUME<VECTOR<T,3>,3>*
Create_Quasistatic_Complex_Muscle(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,const ARRAY<VECTOR<T,3> >& fiber_field,
    ARRAY<bool>& is_tendon,T& activation,const T failure_threshold=(T).25,const bool verbose=true,const bool precompute_stiffness_matrix=true)
{   
    COMPLEX_MUSCLE_3D<T>* constitutive_model=new COMPLEX_MUSCLE_3D<T>(is_tendon,activation,0,0,0,failure_threshold);
    FINITE_VOLUME<VECTOR<T,3>,3>* fvm=Create_Quasistatic_Finite_Volume(tetrahedralized_volume,*constitutive_model,precompute_stiffness_matrix);
    constitutive_model->Initialize_Fiber_Field_From_Current_State(fvm->strain_measure,fiber_field);
    return fvm;
}

}
#endif
