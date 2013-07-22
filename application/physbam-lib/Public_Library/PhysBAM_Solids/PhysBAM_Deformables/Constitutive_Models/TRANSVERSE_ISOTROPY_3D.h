//#####################################################################
// Copyright 2005, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRANSVERSE_ISOTROPY_3D
//##################################################################### 
#ifndef __TRANSVERSE_ISOTROPY_3D__
#define __TRANSVERSE_ISOTROPY_3D__

#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ANISOTROPIC_CONSTITUTIVE_MODEL.h>
namespace PhysBAM{

template<class T>
class TRANSVERSE_ISOTROPY_3D:public ANISOTROPIC_CONSTITUTIVE_MODEL<T,3>
{
    typedef VECTOR<T,3> TV;
public:
    typedef ANISOTROPIC_CONSTITUTIVE_MODEL<T,3> BASE;
    using BASE::enforce_definiteness;using BASE::constant_alpha;using BASE::constant_beta;using BASE::use_isotropic_component_of_stress_derivative_only;
    
    T failure_threshold;
    ARRAY<VECTOR<T,3> > fiber_field;

    TRANSVERSE_ISOTROPY_3D(const T failure_threshold_input)
        :failure_threshold(failure_threshold_input)
    {
        use_isotropic_component_of_stress_derivative_only=false;
    }

    void Initialize_Fiber_Field_From_Current_State(const STRAIN_MEASURE<TV,3>& strain_measure,const ARRAY<VECTOR<T,3> >& fiber_field_input)
    {int n=strain_measure.tetrahedron_mesh.tetrahedrons.m;fiber_field.Resize(n);for(int t=1;t<=n;t++) fiber_field(t)=strain_measure.F(t).Transpose_Times(fiber_field_input(t));}

    inline int Hessian_Index(const int m,const int n) const
    {assert(m>=n);return m+(n-1)*(10-n)/2;}

    void Invariants(VECTOR_ND<T>& invariants,const DIAGONAL_MATRIX<T,3>& C,const VECTOR<T,3> V_fiber) const
    {invariants(1)=C.Trace();invariants(2)=(C*C).Trace();invariants(3)=C.Determinant();
    invariants(4)=VECTOR<T,3>::Dot_Product(V_fiber,C*V_fiber);invariants(5)=(C*V_fiber).Magnitude_Squared();}

    SYMMETRIC_MATRIX<T,3> S(const DIAGONAL_MATRIX<T,3>& C,const VECTOR<T,3>& V_fiber,const VECTOR_ND<T>& invariants,const VECTOR_ND<T> energy_gradient) const
    {SYMMETRIC_MATRIX<T,3> result;
    if(energy_gradient(1)) result+=(T)2*energy_gradient(1);
    if(energy_gradient(2)) result+=(T)4*energy_gradient(2)*C;
    if(energy_gradient(3)) result+=(T)2*energy_gradient(3)*invariants(3)*C.Inverse();
    if(energy_gradient(4)) result+=(T)2*energy_gradient(4)*SYMMETRIC_MATRIX<T,3>::Outer_Product(V_fiber);
    if(energy_gradient(5)) result+=(T)2*energy_gradient(5)*(C*SYMMETRIC_MATRIX<T,3>::Outer_Product(V_fiber)).Twice_Symmetric_Part();
    return result;}

    MATRIX<T,3> P_From_Strain(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& V,const T scale,const int tetrahedron) const PHYSBAM_OVERRIDE
    {DIAGONAL_MATRIX<T,3> F_threshold=F.Max(failure_threshold),C=F_threshold*F_threshold;VECTOR<T,3> V_fiber=V.Transpose_Times(fiber_field(tetrahedron));
    VECTOR_ND<T> invariants(5),energy_gradient(5);SYMMETRIC_MATRIX<T,3> result;Invariants(invariants,C,V_fiber);Energy_Gradient(energy_gradient,invariants,tetrahedron);
    return scale*F_threshold*S(C,V_fiber,invariants,energy_gradient);}

    MATRIX<T,3> P_From_Strain_Rate(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& F_dot,const T scale,const int tetrahedron) const PHYSBAM_OVERRIDE
    {SYMMETRIC_MATRIX<T,3> strain_rate=F_dot.Symmetric_Part(); // Use linear damping by default
    return 2*scale*constant_beta*strain_rate+scale*constant_alpha*strain_rate.Trace();}

    void Stress_Derivative(const DIAGONAL_MATRIX<T,3>& F,const MATRIX<T,3>& V,DIAGONALIZED_STRESS_DERIVATIVE<T,3>& dP_dF,const int tetrahedron=0) const PHYSBAM_OVERRIDE
    {DIAGONAL_MATRIX<T,3> F_threshold=F.Max(failure_threshold),C=F_threshold*F_threshold,C_inverse=C.Inverse();VECTOR<T,3> V_fiber=V.Transpose_Times(fiber_field(tetrahedron));
    SYMMETRIC_MATRIX<T,3> C_inverse_outer=SYMMETRIC_MATRIX<T,3>::Outer_Product(VECTOR<T,3>(C_inverse.x11,C_inverse.x22,C_inverse.x33));
    SYMMETRIC_MATRIX<T,3> V_fiber_outer=SYMMETRIC_MATRIX<T,3>::Outer_Product(V_fiber);
    VECTOR_ND<T> invariants(5),energy_gradient(5),energy_hessian(15);Invariants(invariants,C,V_fiber);
    Energy_Gradient(energy_gradient,invariants,tetrahedron);Energy_Hessian(energy_hessian,invariants,tetrahedron);
    ARRAY<VECTOR<T,3> > J_d(5),J_s(5);
    J_d(1)=VECTOR<T,3>((T)1,(T)1,(T)1);
    J_d(2)=(T)2*VECTOR<T,3>(C.x11,C.x22,C.x33);
    J_d(3)=invariants(3)*VECTOR<T,3>(C_inverse.x11,C_inverse.x22,C_inverse.x33);
    J_d(4)=VECTOR<T,3>(V_fiber_outer.x11,V_fiber_outer.x22,V_fiber_outer.x33);
    J_d(5)=(T)2*VECTOR<T,3>(C.x11*V_fiber_outer.x11,C.x22*V_fiber_outer.x22,C.x33*V_fiber_outer.x33);
    J_s(4)=(T)root_two*VECTOR<T,3>(V_fiber_outer.x21,V_fiber_outer.x31,V_fiber_outer.x32);
    J_s(5)=(T)root_two*VECTOR<T,3>((C.x11+C.x22)*V_fiber_outer.x21,(C.x11+C.x33)*V_fiber_outer.x31,(C.x22+C.x33)*V_fiber_outer.x32);
    dP_dF.dSdC_d=dP_dF.dSdC_s=SYMMETRIC_MATRIX<T,3>();dP_dF.dSdC_ds=MATRIX<T,3>();dP_dF.F=F_threshold;dP_dF.S=S(C,V_fiber,invariants,energy_gradient);
    for(int n=1;n<=5;n++) for(int m=n;m<=5;m++){
        int hessian_index=Hessian_Index(m,n);if(!energy_hessian(hessian_index))continue;
        if(m==n){
            dP_dF.dSdC_d+=(T)2*energy_hessian(hessian_index)*SYMMETRIC_MATRIX<T,3>::Outer_Product(J_d(m));
            dP_dF.dSdC_s+=(T)2*energy_hessian(hessian_index)*SYMMETRIC_MATRIX<T,3>::Outer_Product(J_s(m));
            dP_dF.dSdC_ds+=(T)2*energy_hessian(hessian_index)*MATRIX<T,3>::Outer_Product(J_d(m),J_s(m));}
        else{
            dP_dF.dSdC_d+=(T)2*energy_hessian(hessian_index)*MATRIX<T,3>::Outer_Product(J_d(m),J_d(n)).Twice_Symmetric_Part();
            dP_dF.dSdC_s+=(T)2*energy_hessian(hessian_index)*MATRIX<T,3>::Outer_Product(J_s(m),J_s(n)).Twice_Symmetric_Part();
            dP_dF.dSdC_ds+=(T)2*energy_hessian(hessian_index)*MATRIX<T,3>::Outer_Product(J_d(m),J_s(n)).Twice_Symmetric_Part();}}
    if(energy_gradient(2)){
        dP_dF.dSdC_d+=(T)4*energy_gradient(2);
        dP_dF.dSdC_s+=(T)4*energy_gradient(2);}
    if(energy_gradient(3)){
        dP_dF.dSdC_d+=(T)2*energy_gradient(3)*invariants(3)*SYMMETRIC_MATRIX<T,3>(0,C_inverse_outer.x21,C_inverse_outer.x31,0,C_inverse_outer.x32,0);
        dP_dF.dSdC_s-=(T)2*energy_gradient(3)*invariants(3)*DIAGONAL_MATRIX<T,3>(C_inverse_outer.x21,C_inverse_outer.x31,C_inverse_outer.x32);}
    if(energy_gradient(5)){
        dP_dF.dSdC_d+=(T)4*energy_gradient(5)*DIAGONAL_MATRIX<T,3>(V_fiber_outer.x11,V_fiber_outer.x22,V_fiber_outer.x33);
        dP_dF.dSdC_s+=(T)2*energy_gradient(5)*SYMMETRIC_MATRIX<T,3>(V_fiber_outer.x11+V_fiber_outer.x22,V_fiber_outer.x32,V_fiber_outer.x31,
                                                                      V_fiber_outer.x11+V_fiber_outer.x33,V_fiber_outer.x21,V_fiber_outer.x22+V_fiber_outer.x33);
        dP_dF.dSdC_ds+=(T)2*(T)root_two*energy_gradient(5)*MATRIX<T,3>(V_fiber_outer.x21,V_fiber_outer.x21,0,V_fiber_outer.x31,0,V_fiber_outer.x31,0,V_fiber_outer.x32,V_fiber_outer.x32);}
    if(enforce_definiteness)dP_dF.Enforce_Definiteness();}

//#####################################################################
    virtual void Energy_Gradient(VECTOR_ND<T>& energy_gradient,const VECTOR_ND<T>& invariants,const int tetrahedron) const=0;
    virtual void Energy_Hessian(VECTOR_ND<T>& energy_hessian,const VECTOR_ND<T>& invariants,const int tetrahedron) const=0;
//#####################################################################
};
}
#endif
