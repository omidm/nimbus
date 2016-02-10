//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BW_SHEAR_FORCES
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BW_SHEAR_FORCES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <cfloat>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BW_SHEAR_FORCES<TV>::
BW_SHEAR_FORCES(PARTICLES<TV>& particles,TRIANGLE_MESH& triangle_mesh_input,const T stiffness_coefficient_input,const T damping_coefficient_input)
    :BW_MATERIAL_SPACE_FORCES<TV,1>(particles,triangle_mesh_input,stiffness_coefficient_input,damping_coefficient_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BW_SHEAR_FORCES<TV>::
~BW_SHEAR_FORCES()
{
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void BW_SHEAR_FORCES<TV>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    for(TRIANGLE_ITERATOR iterator(force_simplices);iterator.Valid();iterator.Next()){int s=iterator.Data();
        typename BASE::STATE& state=states(s);
        typename BASE::MATERIAL_FORCE_STATE& material_force_state=material_force_states(s);
        Compute_UV_Deformation(s);

        state.C=VECTOR<T,1>(material_force_state.rest_state_triangle_area*TV::Dot_Product(material_force_state.w_u*material_force_state.w_u_magnitude,
                material_force_state.w_v*material_force_state.w_v_magnitude));
        state.C_dot=VECTOR<T,1>();
        for(int i=1;i<=3;i++){
            state.dC_dx(i)=MATRIX<T,3,1>(material_force_state.rest_state_triangle_area*(material_force_state.dwu_dx(i)*material_force_state.w_u*material_force_state.w_u_magnitude+
                    material_force_state.w_v*material_force_state.w_v_magnitude*material_force_state.dwv_dx(i)));
            state.C_dot+=state.dC_dx(i).Transposed()*particles.V(state.nodes[i]);}
        for(int i=1;i<=3;i++) for(int j=1;j<=3;j++){
            MATRIX<T,3> dC_dxi_dxj;
            T diagonal_term=material_force_state.rest_state_triangle_area*(material_force_state.dwu_dx(i)*material_force_state.dwv_dx(j)+
                material_force_state.dwu_dx(j)*material_force_state.dwv_dx(i));
            for(int s=1;s<=3;s++) for(int t=1;t<=3;t++) if(s==t) dC_dxi_dxj(s,t)=diagonal_term; else dC_dxi_dxj(s,t)=0;
            state.dC_dxi_dxj_times_C(i,j)=dC_dxi_dxj*state.C(1);
            state.dC_dxi_dxj_times_C_dot(i,j)=dC_dxi_dxj*state.C_dot(1);}}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR BW_SHEAR_FORCES<TV>::
Potential_Energy(int s,const T time) const
{
    return (T).5*states(s).stiffness_coefficient*sqr(states(s).C(1));
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR BW_SHEAR_FORCES<TV>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    const_cast<BW_SHEAR_FORCES<TV>* >(this)->Update_Position_Based_State(time,true);
    for(TRIANGLE_ITERATOR iterator(force_simplices);iterator.Valid();iterator.Next()){int s=iterator.Data();
        potential_energy+=Potential_Energy(s,time);}
    return potential_energy;
}
//#####################################################################
// Function Create_BW_Shear_Force
//#####################################################################
template<class TV> BW_SHEAR_FORCES<TV>* PhysBAM::
Create_BW_Shear_Force(PARTICLES<TV>& particles,TRIANGLE_MESH& triangle_mesh,const typename TV::SCALAR stiffness_coefficient_input,const typename TV::SCALAR damping_coefficient_input)
{
    BW_SHEAR_FORCES<TV>* sf=new BW_SHEAR_FORCES<TV>(particles,triangle_mesh,stiffness_coefficient_input,damping_coefficient_input);
    return sf;
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template class BW_SHEAR_FORCES<VECTOR<T,3> >; \
    template BW_SHEAR_FORCES<VECTOR<T,3> >* PhysBAM::Create_BW_Shear_Force<VECTOR<T,3> >(PARTICLES<VECTOR<T,3> >&,TRIANGLE_MESH&,const T,const T);

INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
