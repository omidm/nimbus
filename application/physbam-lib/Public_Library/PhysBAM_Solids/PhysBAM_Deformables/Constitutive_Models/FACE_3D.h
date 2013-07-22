//#####################################################################
// Copyright 2004-2007, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_3D
//##################################################################### 
#ifndef __FACE_3D__
#define __FACE_3D__

#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ANISOTROPIC_CONSTITUTIVE_MODEL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/STRAIN_MEASURE_HEXAHEDRONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME_HEXAHEDRONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/INCOMPRESSIBLE_FINITE_VOLUME.h>
namespace PhysBAM{

template<class T>
class FACE_3D:public ANISOTROPIC_CONSTITUTIVE_MODEL<T,3>
{
    typedef VECTOR<T,3> TV;
public:
    typedef ANISOTROPIC_CONSTITUTIVE_MODEL<T,3> BASE;
    using BASE::enforce_definiteness;using BASE::constant_lambda;using BASE::constant_mu;using BASE::constant_alpha;using BASE::constant_beta;
    using BASE::use_isotropic_component_of_stress_derivative_only;

    T constant_mu_10,constant_mu_01,constant_kappa;
    ARRAY<T> *tet_mu_10,*tet_mu_01,*tet_kappa;
    T fiber_p1,fiber_p2,fiber_p3,fiber_p4,fiber_cutoff;
    ARRAY<T> fiber_max_stress;
    T failure_threshold;
    ARRAY<T>& muscle_activations;
    ARRAY<ARRAY<int> > tet_muscles;
    ARRAY<ARRAY<VECTOR<T,3> > > tet_fibers;
    ARRAY<ARRAY<T> > tet_densities;
    ARRAY<ARRAY<T> > tension,tension_derivative,active_tension_unit_activation;
    int* single_activation_used_for_force_derivative;

    FACE_3D(const T mu_10_input,const T mu_01_input,const T kappa_input,const T Rayleigh_coefficient,ARRAY<T>& muscle_activations_input,
        int* single_activation_used_for_force_derivative_input,const ARRAY<T>* peak_isometric_stress,const T failure_threshold_input)
        :constant_mu_10(mu_10_input),constant_mu_01(mu_01_input),constant_kappa(kappa_input),tet_mu_10(0),tet_mu_01(0),tet_kappa(0),failure_threshold(failure_threshold_input),
        muscle_activations(muscle_activations_input),single_activation_used_for_force_derivative(single_activation_used_for_force_derivative_input)
    {
        use_isotropic_component_of_stress_derivative_only=true;
        constant_mu=2*(constant_mu_10+constant_mu_01);
        constant_lambda=4*constant_mu_10-four_thirds*constant_mu_01+constant_kappa;
        constant_alpha=Rayleigh_coefficient*constant_lambda;
        constant_beta=Rayleigh_coefficient*constant_mu;
        fiber_p1=(T).05;
        fiber_p2=(T)6.6;
        fiber_cutoff=(T)1.4;
        T cutoff_scaled=fiber_p2*(fiber_cutoff-1);
        if(peak_isometric_stress) {assert(peak_isometric_stress->m==muscle_activations.m);fiber_max_stress=*peak_isometric_stress;} 
        else {fiber_max_stress.Resize(muscle_activations.m);ARRAYS_COMPUTATIONS::Fill(fiber_max_stress,(T)3e5);}
        fiber_p3=fiber_p1*fiber_p2*(exp(cutoff_scaled)-1);
        fiber_p4=fiber_p1*(exp(cutoff_scaled)*(1-fiber_p2*fiber_cutoff)+fiber_p2-1);
    }

    void Initialize_Inhomogeneous_Material_Properties(const int number_of_tetrahedrons)
    {tet_mu_10=new ARRAY<T>(number_of_tetrahedrons);ARRAYS_COMPUTATIONS::Fill(*tet_mu_10,constant_mu_10);
    tet_mu_01=new ARRAY<T>(number_of_tetrahedrons);ARRAYS_COMPUTATIONS::Fill(*tet_mu_01,constant_mu_01);
    tet_kappa=new ARRAY<T>(number_of_tetrahedrons);ARRAYS_COMPUTATIONS::Fill(*tet_kappa,constant_kappa);}

    void Set_Material_Properties_Of_Subset(const ARRAY<int>& tetrahedron_list,const T mu_10,const T mu_01,const T kappa)
    {assert(tet_mu_10&&tet_mu_01&&tet_kappa);
    for(int i=1;i<=tetrahedron_list.m;i++) {(*tet_mu_10)(tetrahedron_list(i))=mu_10;(*tet_mu_01)(tetrahedron_list(i))=mu_01;(*tet_kappa)(tetrahedron_list(i))=kappa;}}

    void Initialize_Fiber_Data(const STRAIN_MEASURE<TV,3>& strain_measure,const ARRAY<ARRAY<int> >& muscle_tets,const ARRAY<ARRAY<VECTOR<T,3> > >& muscle_fibers,
        const ARRAY<ARRAY<T> >& muscle_densities)
    {int n=strain_measure.mesh.elements.m;
    tet_muscles.Resize(n);tet_fibers.Resize(n);tet_densities.Resize(n);tension.Resize(n);tension_derivative.Resize(n);active_tension_unit_activation.Resize(n);
    for(int m=1;m<=muscle_tets.m;m++) for(int t=1;t<=muscle_tets(m).m;t++){
        tet_muscles(muscle_tets(m)(t)).Append(m);tet_densities(muscle_tets(m)(t)).Append(muscle_densities(m)(t));
        tet_fibers(muscle_tets(m)(t)).Append(strain_measure.F(muscle_tets(m)(t)).Transpose_Times(muscle_fibers(m)(t)));
        tension(muscle_tets(m)(t)).Append(0);active_tension_unit_activation(muscle_tets(m)(t)).Append(0);tension_derivative(muscle_tets(m)(t)).Append(0);}}

    void Initialize_Fiber_Data(const STRAIN_MEASURE_HEXAHEDRONS<T>& strain_measure,const ARRAY<ARRAY<int> >& muscle_tets,const ARRAY<ARRAY<VECTOR<T,3> > >& muscle_fibers,
        const ARRAY<ARRAY<T> >& muscle_densities)
    {int n=8*strain_measure.hexahedron_mesh.hexahedrons.m;
    tet_muscles.Resize(n);tet_fibers.Resize(n);tet_densities.Resize(n);tension.Resize(n);tension_derivative.Resize(n);active_tension_unit_activation.Resize(n);
    for(int m=1;m<=muscle_tets.m;m++) for(int t=1;t<=muscle_tets(m).m;t++){
        tet_muscles(muscle_tets(m)(t)).Append(m);tet_densities(muscle_tets(m)(t)).Append(muscle_densities(m)(t));
        tet_fibers(muscle_tets(m)(t)).Append(strain_measure.F((muscle_tets(m)(t)-1)%8+1,(muscle_tets(m)(t)-1)/8+1).Transpose_Times(muscle_fibers(m)(t)));
        tension(muscle_tets(m)(t)).Append(0);active_tension_unit_activation(muscle_tets(m)(t)).Append(0);tension_derivative(muscle_tets(m)(t)).Append(0);}}

    T Tension(const int tet_index,const int tet_muscle_index,const T stretch) const
    {T strain=stretch-1,strain_abs=abs(strain),activation=muscle_activations(tet_muscles(tet_index)(tet_muscle_index)),density=tet_densities(tet_index)(tet_muscle_index);
    T active_tension=0,passive_tension=0,scale=(T)25/(T)6;
    if(stretch>fiber_cutoff)passive_tension=fiber_p3*stretch+fiber_p4;else if(stretch>1)passive_tension=fiber_p1*(exp(fiber_p2*strain)-fiber_p2*strain-1);
    if(strain_abs<.4)active_tension=activation*density*(1-scale*sqr(strain));else if(strain_abs<.6)active_tension=2*scale*activation*density*sqr(strain_abs-(T).6);
    return fiber_max_stress(tet_muscles(tet_index)(tet_muscle_index))*(active_tension+passive_tension);}
        
    T Active_Tension_Unit_Activation(const int tet_index,const int tet_muscle_index,const T stretch) const
    {T strain=stretch-1,strain_abs=abs(strain),density=tet_densities(tet_index)(tet_muscle_index),active_tension=0,scale=(T)25/(T)6;
    if(strain_abs<.4)active_tension=density*(1-scale*sqr(strain));else if(strain_abs<.6)active_tension=2*scale*density*sqr(strain_abs-(T).6);
    return fiber_max_stress(tet_muscles(tet_index)(tet_muscle_index))*active_tension;}
        
    T Tension_Derivative(const int tet_index,const int tet_muscle_index,const T stretch) const
    {T strain=stretch-1,strain_abs=abs(strain),activation=muscle_activations(tet_muscles(tet_index)(tet_muscle_index)),density=tet_densities(tet_index)(tet_muscle_index);
    T active_tension_derivative=0,passive_tension_derivative=0,scale=(T)25/(T)6;
    if(stretch>fiber_cutoff)passive_tension_derivative=fiber_p3;else if(stretch>1)passive_tension_derivative=fiber_p1*fiber_p2*(exp(fiber_p2*strain)-1);
    if(strain_abs<.4)active_tension_derivative=-2*scale*activation*density*strain;else if(strain_abs<.6)active_tension_derivative=4*scale*activation*density*(strain-sign(strain)*(T).6);
    return fiber_max_stress(tet_muscles(tet_index)(tet_muscle_index))*(active_tension_derivative+passive_tension_derivative);}

    MATRIX<T,3> P_From_Strain(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& V,const T scale,const int tetrahedron_index) const PHYSBAM_OVERRIDE
    {T mu_10=(tet_mu_10)?(*tet_mu_10)(tetrahedron_index):constant_mu_10,mu_01=(tet_mu_01)?(*tet_mu_01)(tetrahedron_index):constant_mu_01,kappa=(tet_kappa)?(*tet_kappa)(tetrahedron_index):constant_kappa;
    if(single_activation_used_for_force_derivative&&(*single_activation_used_for_force_derivative))return P_From_Strain_Unit_Activation(F,V,scale,tetrahedron_index,*single_activation_used_for_force_derivative);
    DIAGONAL_MATRIX<T,3> F_threshold=F.Max(failure_threshold),C=F_threshold*F_threshold,F_cube=C*F_threshold,F_inverse=F_threshold.Inverse(),isotropic_part;
    MATRIX<T,3> anisotropic_part;T I_C=C.Trace(),II_C=(C*C).Trace(),J=F_threshold.Determinant(),Jcc=pow(J,-(T)two_thirds);
    isotropic_part=(2*Jcc*(mu_10+Jcc*mu_01*I_C))*F_threshold-(2*Jcc*Jcc*mu_01)*F_cube+((kappa*log(J)-(T)two_thirds*Jcc*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))))*F_inverse;
    for(int m=1;m<=tet_muscles(tetrahedron_index).m;m++){
        VECTOR<T,3> fiber=V.Transpose_Times(tet_fibers(tetrahedron_index)(m)),F_fiber=F_threshold*fiber;
        anisotropic_part+=(tension(tetrahedron_index)(m)/F_fiber.Magnitude())*MATRIX<T,3>::Outer_Product(F_fiber,fiber);}
    return scale*(anisotropic_part+isotropic_part);}
    
    MATRIX<T,3> P_From_Strain_Unit_Activation(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& V,const T scale,const int tetrahedron_index,const int muscle_id) const
    {assert(1<=muscle_id&&muscle_id<=muscle_activations.m);DIAGONAL_MATRIX<T,3> F_threshold=F.Max(failure_threshold);
    int tet_muscle_index=tet_muscles(tetrahedron_index).Find(muscle_id);
    if(!tet_muscle_index)return MATRIX<T,3>();
    VECTOR<T,3> fiber=V.Transpose_Times(tet_fibers(tetrahedron_index)(tet_muscle_index)),F_fiber=F_threshold*fiber;
    return scale*(active_tension_unit_activation(tetrahedron_index)(tet_muscle_index)/F_fiber.Magnitude())*MATRIX<T,3>::Outer_Product(F_fiber,fiber);}

    MATRIX<T,3> P_From_Strain_Rate(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& F_dot,const T scale,const int tetrahedron) const PHYSBAM_OVERRIDE
    {SYMMETRIC_MATRIX<T,3> strain_rate=F_dot.Symmetric_Part(); 
    return 2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace();}

    void Isotropic_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& F,DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,const int tetrahedron_index=0) const PHYSBAM_OVERRIDE
    {T mu_10=(tet_mu_10)?(*tet_mu_10)(tetrahedron_index):constant_mu_10,mu_01=(tet_mu_01)?(*tet_mu_01)(tetrahedron_index):constant_mu_01,kappa=(tet_kappa)?(*tet_kappa)(tetrahedron_index):constant_kappa;
    DIAGONAL_MATRIX<T,3> F_threshold=F.Max(failure_threshold),C=F_threshold*F_threshold,F_cube=C*F_threshold,F_inverse=F_threshold.Inverse();
    T I_C=C.Trace(),II_C=(C*C).Trace(),J=F_threshold.Determinant(),Jcc=pow(J,-(T)two_thirds);
    SYMMETRIC_MATRIX<T,3> alpha;
    alpha.x11=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x11-C.x11));alpha.x21=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x22-C.x11));alpha.x31=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x33-C.x11));
    alpha.x22=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x22-C.x22));alpha.x32=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x33-C.x22));alpha.x33=2*Jcc*(mu_10+Jcc*mu_01*(I_C-C.x33-C.x33));
    SYMMETRIC_MATRIX<T,3> beta;
    beta.x11=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x11*F_inverse.x11-2*Jcc*Jcc*mu_01*F_threshold.x11*F_threshold.x11;
    beta.x21=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x22*F_inverse.x11-2*Jcc*Jcc*mu_01*F_threshold.x22*F_threshold.x11;
    beta.x31=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x33*F_inverse.x11-2*Jcc*Jcc*mu_01*F_threshold.x33*F_threshold.x11;
    beta.x22=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x22*F_inverse.x22-2*Jcc*Jcc*mu_01*F_threshold.x22*F_threshold.x22;
    beta.x32=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x33*F_inverse.x22-2*Jcc*Jcc*mu_01*F_threshold.x33*F_threshold.x22;
    beta.x33=(Jcc*(T)two_thirds*(mu_10*I_C+Jcc*mu_01*(I_C*I_C-II_C))-kappa*log(J))*F_inverse.x33*F_inverse.x33-2*Jcc*Jcc*mu_01*F_threshold.x33*F_threshold.x33;
    SYMMETRIC_MATRIX<T,3> eta;
    eta.x11=4*Jcc*Jcc*mu_01;
    eta.x31=-Jcc*(T)one_third*(4*mu_10+8*Jcc*mu_01*I_C);
    eta.x32=8*(T)one_third*Jcc*Jcc*mu_01;
    eta.x33=Jcc*(T)one_ninth*(4*mu_10*I_C+Jcc*8*mu_01*(I_C*I_C-II_C))+kappa;
    MATRIX<T,3> F_base(F_threshold.x11,F_threshold.x22,F_threshold.x33,F_cube.x11,F_cube.x22,F_cube.x33,F_inverse.x11,F_inverse.x22,F_inverse.x33);
    SYMMETRIC_MATRIX<T,3> gamma=SYMMETRIC_MATRIX<T,3>::Conjugate(F_base,eta);
    dPi_dF.x1111=alpha.x11+beta.x11+gamma.x11;dPi_dF.x2222=alpha.x22+beta.x22+gamma.x22;dPi_dF.x3333=alpha.x33+beta.x33+gamma.x33;
    dPi_dF.x2211=gamma.x21;dPi_dF.x3311=gamma.x31;dPi_dF.x3322=gamma.x32;
    dPi_dF.x2121=alpha.x21;dPi_dF.x3131=alpha.x31;dPi_dF.x3232=alpha.x32;
    dPi_dF.x2112=beta.x21;dPi_dF.x3113=beta.x31;dPi_dF.x3223=beta.x31;
    if(enforce_definiteness) dPi_dF.Enforce_Definiteness();}

    void Stress_Derivative(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& V,DIAGONALIZED_STRESS_DERIVATIVE<T,3>& dP_dF,const int simplex) const PHYSBAM_OVERRIDE
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}

    MATRIX<T,3> dP_From_dF(const MATRIX<T,3>& dF,const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& V,const DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>& dPi_dF,const T scale,
        const int tetrahedron_index) const
    {MATRIX<T,3> dP=scale*dPi_dF.Differential(dF);
    DIAGONAL_MATRIX<T,3> F_threshold=F.Max(failure_threshold);
    for(int m=1;m<=tet_muscles(tetrahedron_index).m;m++){
        VECTOR<T,3> fiber=V.Transpose_Times(tet_fibers(tetrahedron_index)(m)),F_fiber=F_threshold*fiber,dF_fiber=dF*fiber;
        T stretch_squared=F_fiber.Magnitude_Squared(),stretch=sqrt(stretch_squared);
        T c1=tension(tetrahedron_index)(m)/stretch,c2=(tension_derivative(tetrahedron_index)(m)-c1)/stretch_squared;
        VECTOR<T,3> dPw=c1*dF_fiber+c2*VECTOR<T,3>::Dot_Product(dF_fiber,F_fiber)*F_fiber;
        dP+=scale*MATRIX<T,3>::Outer_Product(dPw,fiber);}
    return dP;}

    void Update_State_Dependent_Auxiliary_Variables(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& V,const int tetrahedron_index) PHYSBAM_OVERRIDE
    {DIAGONAL_MATRIX<T,3> F_threshold=F.Max(failure_threshold);
    for(int m=1;m<=tet_muscles(tetrahedron_index).m;m++){
        VECTOR<T,3> fiber=V.Transpose_Times(tet_fibers(tetrahedron_index)(m)),F_fiber=F_threshold*fiber;T stretch=F_fiber.Magnitude();
        tension(tetrahedron_index)(m)=Tension(tetrahedron_index,m,stretch);
        active_tension_unit_activation(tetrahedron_index)(m)=Active_Tension_Unit_Activation(tetrahedron_index,m,stretch);
        if(enforce_definiteness) tension_derivative(tetrahedron_index)(m)=max((T)0,Tension_Derivative(tetrahedron_index,m,stretch));
        else tension_derivative(tetrahedron_index)(m)=Tension_Derivative(tetrahedron_index,m,stretch);}}

//#####################################################################
};

template<class T> FINITE_VOLUME<VECTOR<T,3>,3>*
Create_Face(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,const ARRAY<ARRAY<int> >& muscle_tets,
    const ARRAY<ARRAY<VECTOR<T,3> > >& muscle_fibers,const ARRAY<ARRAY<T> >& muscle_densities,ARRAY<T>& muscle_activations,const ARRAY<T>* peak_isometric_stress=0,
    const T mu_10=(T)6e4,const T mu_01=(T)2e4,const T kappa=(T)6e4,const T Rayleigh_coefficient=(T).05,const T failure_threshold=(T).25,const bool limit_time_step_by_strain_rate=true,
    const T max_strain_per_time_step=(T).1,const bool use_rest_state_for_strain_rate=true,const bool precompute_stiffness_matrix=true,const bool verbose=true)
{
    FACE_3D<T>* constitutive_model=new FACE_3D<T>(mu_10,mu_01,kappa,Rayleigh_coefficient,muscle_activations,0,peak_isometric_stress,failure_threshold);
    FINITE_VOLUME<VECTOR<T,3>,3>* fvm=Create_Finite_Volume(tetrahedralized_volume,constitutive_model,limit_time_step_by_strain_rate,max_strain_per_time_step,
        use_rest_state_for_strain_rate,precompute_stiffness_matrix,verbose);
    constitutive_model->Initialize_Fiber_Data(fvm->strain_measure,muscle_tets,muscle_fibers,muscle_densities);
    return fvm;
}  

template<class T> INCOMPRESSIBLE_FINITE_VOLUME<VECTOR<T,3>,3>*
Create_Incompressible_Face(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,const ARRAY<ARRAY<int> >& muscle_tets,
    const ARRAY<ARRAY<VECTOR<T,3> > >& muscle_fibers,const ARRAY<ARRAY<T> >& muscle_densities,ARRAY<T>& muscle_activations,const ARRAY<T>* peak_isometric_stress=0,
    const T mu_10=(T)6e4,const T mu_01=(T)2e4,const T kappa=(T)6e4,const T Rayleigh_coefficient=(T).05,const T failure_threshold=(T).25,const bool limit_time_step_by_strain_rate=true,
    const T max_strain_per_time_step=(T).1,const bool use_rest_state_for_strain_rate=true,const bool precompute_stiffness_matrix=true,const bool verbose=true)
{
    {std::stringstream ss;ss<<"Creating Incompressible_Face"<<std::endl;LOG::filecout(ss.str());}
    FACE_3D<T>* constitutive_model=new FACE_3D<T>(mu_10,mu_01,kappa,Rayleigh_coefficient,muscle_activations,0,peak_isometric_stress,failure_threshold);
    INCOMPRESSIBLE_FINITE_VOLUME<VECTOR<T,3>,3>* fvm=Create_Incompressible_Finite_Volume(tetrahedralized_volume,constitutive_model,verbose);
    constitutive_model->Initialize_Fiber_Data(fvm->strain_measure,muscle_tets,muscle_fibers,muscle_densities);
    return fvm;
} 

template<class T> FINITE_VOLUME<VECTOR<T,3>,3>*
Create_Quasistatic_Face(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,const ARRAY<ARRAY<int> >& muscle_tets,
    const ARRAY<ARRAY<VECTOR<T,3> > >& muscle_fibers,const ARRAY<ARRAY<T> >& muscle_densities,ARRAY<T>& muscle_activations,
    int* single_activation_used_for_force_derivative=0,const ARRAY<T>* peak_isometric_stress=0,CONSTITUTIVE_MODEL<T,3> **face_constitutive_model=0,const T mu_10=(T)6e4,
    const T mu_01=(T)2e4,const T kappa=(T)6e4,const T failure_threshold=(T).25,const bool precompute_stiffness_matrix=true,const bool verbose=false)
{   
    FACE_3D<T>* constitutive_model=new FACE_3D<T>(mu_10,mu_01,kappa,0,muscle_activations,single_activation_used_for_force_derivative,peak_isometric_stress,
        failure_threshold);
    FINITE_VOLUME<VECTOR<T,3>,3>* fvm=Create_Quasistatic_Finite_Volume(tetrahedralized_volume,constitutive_model,precompute_stiffness_matrix,verbose);
    constitutive_model->Initialize_Fiber_Data(fvm->strain_measure,muscle_tets,muscle_fibers,muscle_densities);
    if(face_constitutive_model) *face_constitutive_model=constitutive_model;
    return fvm;
}

template<class T> FINITE_VOLUME_HEXAHEDRONS<T>*
Create_Quasistatic_Face(HEXAHEDRALIZED_VOLUME<T>& hexahedralized_volume,const ARRAY<ARRAY<int> >& muscle_tets,
    const ARRAY<ARRAY<VECTOR<T,3> > >& muscle_fibers,const ARRAY<ARRAY<T> >& muscle_densities,ARRAY<T>& muscle_activations,
    CONSTITUTIVE_MODEL<T,3>*& face_constitutive_model,int* single_activation_used_for_force_derivative=0,const ARRAY<T>* peak_isometric_stress=0,const T mu_10=(T)6e4,
    const T mu_01=(T)2e4,const T kappa=(T)6e4,const T failure_threshold=(T).25,const bool verbose=false)
{
    FACE_3D<T>* constitutive_model=new FACE_3D<T>(mu_10,mu_01,kappa,0,muscle_activations,single_activation_used_for_force_derivative,peak_isometric_stress,
        failure_threshold);
    FINITE_VOLUME_HEXAHEDRONS<T>* fvm=Create_Quasistatic_Finite_Volume(hexahedralized_volume,constitutive_model);
    constitutive_model->Initialize_Fiber_Data(fvm->strain_measure,muscle_tets,muscle_fibers,muscle_densities);
    face_constitutive_model=constitutive_model;
    return fvm;
}

}
#endif
