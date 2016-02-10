//#####################################################################
// Copyright 2002-2004, Silvia Salinas-Blemker, Geoffrey Irving, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MUSCLE_3D
//#####################################################################
#ifndef __MUSCLE_3D__
#define __MUSCLE_3D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_NXN.h>
#include <PhysBAM_Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ANISOTROPIC_CONSTITUTIVE_MODEL.h>
namespace PhysBAM{

template<class T>
class MUSCLE_3D:public ANISOTROPIC_CONSTITUTIVE_MODEL<T,3>
{
    typedef VECTOR<T,3> TV;
public:
    typedef ANISOTROPIC_CONSTITUTIVE_MODEL<T,3> BASE;
    using BASE::constant_alpha;using BASE::constant_beta;

    ARRAY<T> inp_activation; // between 0 and 1
    T inp_matrix_c1,inp_matrix_c2,inp_bulk; // in Pa
    T inp_pass_c1,inp_pass_c2; // in Pa
    T inp_pass_starlam; // unitless
    T peak_isometric_stress; // in Pa
    ARRAY<TV> fiber_field;
    ARRAY<int> material_numbers;
    ARRAY<int> motor_units;
    T tendon_matrix_c1;
    T tendon_matrix_c2;
    T tendon_pass_c1;
    T tendon_pass_c2;
    T tendon_pass_starlam;
    T tendon_bulk;
    T muscle_pass_c3,muscle_pass_c4;
    T tendon_pass_c3,tendon_pass_c4;
    T failure_threshold;
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume;

    MUSCLE_3D(T alpha_input,T beta_input,T inp_activation_input,T inp_matrix_c1_input,T inp_matrix_c2_input,T inp_bulk_input,T inp_pass_c1_input,
        T inp_pass_c2_input,T inp_pass_starlam_input,T peak_isometric_stress_input,T failure_threshold_input,
        int number_of_motor_units,TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume_input)
        :inp_matrix_c1(inp_matrix_c1_input),inp_matrix_c2(inp_matrix_c2_input),
        inp_bulk(inp_bulk_input),inp_pass_c1(inp_pass_c1_input),inp_pass_c2(inp_pass_c2_input),inp_pass_starlam(inp_pass_starlam_input),
        peak_isometric_stress(peak_isometric_stress_input),tendon_matrix_c1(T(25000)),tendon_matrix_c2(T(0)),tendon_pass_c1(T(0)),
        tendon_pass_c2(T(0)),tendon_pass_starlam(T(1)),tendon_bulk(T(50000)),failure_threshold(failure_threshold_input),tetrahedralized_volume(tetrahedralized_volume_input)
    {
        constant_alpha=alpha_input;constant_beta=beta_input;
        inp_activation.Resize(number_of_motor_units);ARRAYS_COMPUTATIONS::Fill(inp_activation,inp_activation_input);
        if(number_of_motor_units==1)
        {motor_units.Resize(tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m);ARRAYS_COMPUTATIONS::Fill(motor_units,1);}
        muscle_pass_c3=inp_pass_c1*inp_pass_c2*exp(inp_pass_c2*(inp_pass_starlam-1));
        muscle_pass_c4=inp_pass_c1*(exp(inp_pass_c2*(inp_pass_starlam-1))-1)-muscle_pass_c3*inp_pass_starlam;
        tendon_pass_c3=tendon_pass_c1*tendon_pass_c2*exp(tendon_pass_c2*(tendon_pass_starlam-1));
        tendon_pass_c4=tendon_pass_c1*(exp(tendon_pass_c2*(tendon_pass_starlam-1))-1)-tendon_pass_c3*tendon_pass_starlam;
        material_numbers.Resize(tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m);ARRAYS_COMPUTATIONS::Fill(material_numbers,0);
        fiber_field.Resize(tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m);ARRAYS_COMPUTATIONS::Fill(fiber_field,TV(1,0,0));
    }

    static MUSCLE_3D<T>* Create(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,const T failure_threshold)
    {return new MUSCLE_3D<T>(2500,1000,0,50000,0,50000,0,0,1,60000,failure_threshold,1,tetrahedralized_volume);}

    void Initialize_Fibers_Tendons_And_Motor_Units(ARRAY<TV>& fibers,ARRAY<bool>& tendon,ARRAY<int>& motor_units_input)
    {fiber_field=fibers.m;
    material_numbers.Resize(fibers.m);for(int t=1;t<=tendon.m;t++) if(tendon(t)) material_numbers(t)=1;else material_numbers(t)=0;
    if(motor_units_input.m) motor_units=motor_units_input;
    else{motor_units.Resize(fibers.m);ARRAYS_COMPUTATIONS::Fill(motor_units,1);}}

    void Correct_Fiber_Field_For_QR(STRAIN_MEASURE<TV,3>& strain_measure)
    {for(int t=1;t<=fiber_field.m;t++)fiber_field(t)=strain_measure.F(t).Transposed()*fiber_field(t);}

    T Fiber_Tension(const int tetrahedron,const T stretch) const
    {T pass_c1=inp_pass_c1,pass_c2=inp_pass_c2,pass_c3=muscle_pass_c3,pass_c4=muscle_pass_c4,pass_starlam=inp_pass_starlam,activation=inp_activation(motor_units(tetrahedron));
    if(material_numbers(tetrahedron)){pass_c1=tendon_pass_c1;pass_c2=tendon_pass_c2;pass_c3=tendon_pass_c3;pass_c4=tendon_pass_c4;pass_starlam=tendon_pass_starlam;activation=0;}
    T passive_fiber_tension=0,active_fiber_tension=0;
    T strain=stretch-1,stretch_inv=1/stretch,strain_abs=abs(strain);
    // Passive Stretch:
    if(stretch>1){
        if(stretch<pass_starlam)passive_fiber_tension=pass_c1*stretch_inv*(exp(pass_c4*strain)-1);
        else passive_fiber_tension=pass_c3+pass_c4*stretch_inv;}
    if(strain_abs>.6)active_fiber_tension=0;
    else if(strain_abs>.4)active_fiber_tension=stretch_inv*peak_isometric_stress*activation*9*sqr(strain_abs-(T).6);
    else active_fiber_tension=stretch_inv*peak_isometric_stress*activation*(1-4*sqr(strain_abs));
    return (T).5*stretch_inv*passive_fiber_tension+active_fiber_tension;}

    MATRIX<T,3> P_From_Strain(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& V,const T scale,const int tetrahedron) const PHYSBAM_OVERRIDE
    {DIAGONAL_MATRIX<T,3> F_threshold=F.Max(failure_threshold);
    T matrix_c1=inp_matrix_c1,matrix_c2=inp_matrix_c2,bulk=inp_bulk;
    if(material_numbers(tetrahedron)){matrix_c1=tendon_matrix_c1;matrix_c2=tendon_matrix_c2;bulk=tendon_bulk;}
    TV fiber_m=V.Transpose_Times(fiber_field(tetrahedron)),fiber_s=F_threshold*fiber_m;
    DIAGONAL_MATRIX<T,3> F_inv=F_threshold.Inverse(),C=F_threshold*F_threshold,F_cube=C*F_threshold;
    T J=F_threshold.Determinant(),Jc=pow(J,-(T)one_third),Jcc=Jc*Jc,scaler=4*Jcc;
    T C_trace=C.Trace(),CC_trace=(C*C).Trace(),I1=Jcc*C_trace;
    T stretch_squared=fiber_s.Magnitude_Squared(),fiber_stretch=Jc*sqrt(stretch_squared);
    T W1_scaled=scaler*matrix_c1,W2_scaled=scaler*matrix_c2,W12_scaled=W1_scaled+I1*W2_scaled,tension_scaled=scaler*Fiber_Tension(tetrahedron,fiber_stretch);
    T pressure_J=bulk*log(J);
    T false_pressure=(T)one_third*(W12_scaled*C_trace-Jcc*W2_scaled*CC_trace+tension_scaled*stretch_squared);
    MATRIX<T,3> fiber_tensor=MATRIX<T,3>::Outer_Product(fiber_s,fiber_m);
    DIAGONAL_MATRIX<T,3> isotropic_part=scale*(W12_scaled*F_threshold-Jcc*W2_scaled*F_cube+(pressure_J-false_pressure)*F_inv);
    return isotropic_part+scale*tension_scaled*fiber_tensor;}

    MATRIX<T,3> P_From_Strain_Rate(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& F_dot,const T scale,const int tetrahedron) const PHYSBAM_OVERRIDE
    {SYMMETRIC_MATRIX<T,3> strain_rate=F_dot.Symmetric_Part();
    return 2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace();}

    SYMMETRIC_MATRIX<T,3> Sigma_From_Strain(const int tetrahedron,const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& V,const T scale) const
    {return (T(1)/F.Determinant())*MATRIX<T,3>::Right_Multiply_With_Symmetric_Result(P_From_Strain(F,V,scale,tetrahedron),F);}

    T Elasticity_From_Strain(const int tetrahedron,const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& V) const
    {T matrix_c1=inp_matrix_c1,matrix_c2=inp_matrix_c2,bulk=inp_bulk;
    if(material_numbers(tetrahedron)){matrix_c1=tendon_matrix_c1;matrix_c2=tendon_matrix_c2;bulk=tendon_bulk;} // tendon

    SYMMETRIC_MATRIX<T,3> sigma=Sigma_From_Strain(tetrahedron,F,V,1);
    T J=F.Determinant(),Jc=pow(J,-(T)one_third),Jc_sqr=Jc*Jc;

    DIAGONAL_MATRIX<T,3> C=F*F,Ct=C*Jc_sqr,Cts=Ct*Ct;

    T I1=Ct.Trace();
    T Cts_trace=(Ct*Ct).Trace();
    T I2=(T).5*(I1*I1-Cts_trace);

    T w1=matrix_c1; T w2=matrix_c2;
    T  pr=bulk*log(J)/J;

    T ch4=-w2;
    T cc4=T(4)/J; T cc3=T(-1./3.)*cc4; T cc2=T(-2./3.);

    T cc1=T(-2./9.)*( w1*I1+T(2)*w2*I2);
    T cc5= T(4./9.)*( I1*(I1*w2)+Cts_trace*(-w2));

    T caci2=I1*w2;

    T aci1=Ct.x11*caci2-Cts.x11*w2;
    T aci2=Ct.x22*caci2-Cts.x22*w2;
    T aci3=Ct.x33*caci2-Cts.x33*w2;

    aci1=cc2*sigma.x11+cc3*aci1;
    aci2=cc2*sigma.x22+cc3*aci2;
    aci3=cc2*sigma.x33+cc3*aci3;

    SYMMETRIC_MATRIX_NXN<T> elasticity(6);
    elasticity(1,1)=-pr-2*cc1+2*aci1+cc5;
    elasticity(2,2)=-pr-2*cc1+2*aci2+cc5;
    elasticity(3,3)=-pr-2*cc1+2*aci3+cc5;
    cc1=-(T)1.5*cc1;ch4=(T).5*ch4;
    elasticity(4,4)=-pr+cc1+cc4*(ch4*Ct.x11*Ct.x22);
    elasticity(5,5)=-pr+cc1+cc4*(ch4*Ct.x33*Ct.x22);
    elasticity(6,6)=-pr+cc1+cc4*(ch4*Ct.x33*Ct.x11);

    return elasticity.Trace();}

    T CFL_Elastic(const T density,const T altitude,const int tetrahedron) const
    {T max_elastic_trace=0;
    for(int j=1;j<tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m;j++){
        max_elastic_trace= max(abs(Elasticity_From_Strain(j,DIAGONAL_MATRIX<T,3>::Identity_Matrix(),DIAGONAL_MATRIX<T,3>::Identity_Matrix())),max_elastic_trace);}
    T elastic_sound_speed=sqrt(max_elastic_trace/density);
    return altitude/elastic_sound_speed;}

//#####################################################################
};
}
#endif


