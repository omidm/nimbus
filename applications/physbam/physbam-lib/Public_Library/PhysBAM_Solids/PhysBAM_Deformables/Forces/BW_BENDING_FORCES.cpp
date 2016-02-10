//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BW_BENDING_FORCES
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/cyclic_shift.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Functions.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BW_BENDING_FORCES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <cfloat>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BW_BENDING_FORCES<TV>::
BW_BENDING_FORCES(PARTICLES<TV>& particles,TRIANGLE_MESH& triangle_mesh_input,const T stiffness_coefficient_input,const T damping_coefficient_input)
    :BW_FORCES<TV,1,4>(particles,triangle_mesh_input,stiffness_coefficient_input,damping_coefficient_input),assume_constant_normal_length(true)
{
    // constant matrices
    VECTOR<VECTOR<T,3>,3> pos=VECTOR<VECTOR<T,3>,3>(VECTOR<T,3>(1,0,0),VECTOR<T,3>(0,1,0),VECTOR<T,3>(0,0,1));
    VECTOR<VECTOR<T,3>,3> neg=VECTOR<VECTOR<T,3>,3>(VECTOR<T,3>(-1,0,0),VECTOR<T,3>(0,-1,0),VECTOR<T,3>(0,0,-1));
    dq_a_i_j_t(1,2)=neg;dq_a_i_j_t(1,3)=pos;dq_a_i_j_t(2,3)=neg;dq_a_i_j_t(2,1)=pos;dq_a_i_j_t(3,1)=neg;dq_a_i_j_t(3,2)=pos;
    dq_b_i_j_t(2,4)=neg;dq_b_i_j_t(2,3)=pos;dq_b_i_j_t(3,2)=neg;dq_b_i_j_t(3,4)=pos;dq_b_i_j_t(4,3)=neg;dq_b_i_j_t(4,2)=pos;


    bool adjacent_triangles_defined=(triangle_mesh.adjacent_elements!=0);if(!adjacent_triangles_defined) triangle_mesh.Initialize_Adjacent_Elements();
    int number_quadruples=0;
    for(int t=1;t<=triangle_mesh.elements.m;t++) for(int a=1;a<=(*triangle_mesh.adjacent_elements)(t).m;a++) if((*triangle_mesh.adjacent_elements)(t)(a)>t) number_quadruples++;
    constraint_particles.Resize(number_quadruples,false);
    int index=0; // reset number
    for(int t=1;t<=triangle_mesh.elements.m;t++){
        int t1,t2,t3;triangle_mesh.elements(t).Get(t1,t2,t3);
        for(int a=1;a<=(*triangle_mesh.adjacent_elements)(t).m;a++){
            int s=(*triangle_mesh.adjacent_elements)(t)(a);
            if(s>t){
                int s1,s2,s3;triangle_mesh.elements(s).Get(s1,s2,s3);
                if(t1==s1 || t1==s2 || t1==s3){cyclic_shift(t1,t2,t3);if(t1==s1 || t1==s2 || t1==s3) cyclic_shift(t1,t2,t3);}
                constraint_particles(++index).Set(t2,t3,t1,triangle_mesh.Other_Node(t2,t3,s));}}}
    if(!adjacent_triangles_defined){delete triangle_mesh.adjacent_elements;triangle_mesh.adjacent_elements=0;}

    states.Resize(number_quadruples,false,false);
    bending_states.Resize(number_quadruples,false,false);
    ARRAY<bool> particle_is_simulated(particles.array_collection->Size());ARRAYS_COMPUTATIONS::Fill(particle_is_simulated,true);
    Update_Mpi(particle_is_simulated,0);
    for(CONSTRAINT_ITERATOR iterator(force_simplices);iterator.Valid();iterator.Next()){int q=iterator.Data();
        int x2,x1,x0,x3;constraint_particles(q).Get(x2,x1,x0,x3);
        typename BASE::STATE& state=states(q);
        state.nodes=VECTOR<int,4>(x0,x1,x2,x3);
        state.stiffness_coefficient=stiffness_coefficient_input;
        state.damping_coefficient=damping_coefficient_input;}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BW_BENDING_FORCES<TV>::
~BW_BENDING_FORCES()
{
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void BW_BENDING_FORCES<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    for(int q=1;q<=constraint_particles.m;q++)
        for(int i=1;i<=3;i++) for(int j=i+1;j<=4;j++) dependency_mesh.Add_Element_If_Not_Already_There(VECTOR<int,2>(constraint_particles(q)[i],constraint_particles(q)[j]));
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void BW_BENDING_FORCES<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    force_simplices.Update(constraint_particles,particle_is_simulated);
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void BW_BENDING_FORCES<TV>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    for(CONSTRAINT_ITERATOR iterator(force_simplices);iterator.Valid();iterator.Next()){int q=iterator.Data();
        typename BASE::STATE& state=states(q);
        BENDING_STATE& bending_state=bending_states(q);
        TV x31=particles.X(state.nodes(3))-particles.X(state.nodes(1)),x21=particles.X(state.nodes(2))-particles.X(state.nodes(1));
        TV x24=particles.X(state.nodes(2))-particles.X(state.nodes(4)),x34=particles.X(state.nodes(3))-particles.X(state.nodes(4));
        bending_state.n_a=TV::Cross_Product(x31,x21);
        bending_state.n_b=TV::Cross_Product(x24,x34);
        bending_state.e=particles.X(state.nodes(2))-particles.X(state.nodes(3));
        bending_state.n_a_magnitude=bending_state.n_a.Normalize();
        bending_state.n_b_magnitude=bending_state.n_b.Normalize();
        bending_state.e_magnitude=bending_state.e.Normalize();
        
        T cos_theta=TV::Dot_Product(bending_state.n_a,bending_state.n_b);
        T sin_theta=TV::Dot_Product(TV::Cross_Product(bending_state.n_a,bending_state.n_b),bending_state.e);
        state.C=VECTOR<T,1>(atan2(sin_theta,cos_theta));
        state.C_dot=VECTOR<T,1>();

        VECTOR<VECTOR<T,3>,4> q_a=VECTOR<VECTOR<T,3>,4>(particles.X(state.nodes(3))-particles.X(state.nodes(2)),
            particles.X(state.nodes(1))-particles.X(state.nodes(3)),particles.X(state.nodes(2))-particles.X(state.nodes(1)),VECTOR<T,3>());
        VECTOR<VECTOR<T,3>,4> q_b=VECTOR<VECTOR<T,3>,4>(VECTOR<T,3>(),particles.X(state.nodes(3))-particles.X(state.nodes(4)),
            particles.X(state.nodes(4))-particles.X(state.nodes(2)),particles.X(state.nodes(2))-particles.X(state.nodes(3)));
        VECTOR<T,4> q_e=VECTOR<T,4>(0,1,-1,0);

        VECTOR<VECTOR<T,3>,4> dsin_theta_dx;
        VECTOR<VECTOR<T,3>,4> dcos_theta_dx;
        VECTOR<VECTOR<VECTOR<T,3>,3>,4> dn_a_dx;
        VECTOR<VECTOR<VECTOR<T,3>,3>,4> dn_b_dx;
        VECTOR<VECTOR<VECTOR<T,3>,3>,4> de_dx;
        MATRIX<MATRIX<T,3>,4> d2sin_theta_dx_dy;
        MATRIX<MATRIX<T,3>,4> d2cos_theta_dx_dy;
        for(int i=1;i<=4;i++){
            MATRIX<T,3> S_q_a=MATRIX<T,3>::Cross_Product_Matrix(q_a(i));
            MATRIX<T,3> S_q_b=MATRIX<T,3>::Cross_Product_Matrix(q_b(i));
            for(int j=1;j<=3;j++){
                dn_a_dx(i)(j)=-S_q_a.Column(j)/bending_state.n_a_magnitude;
                dn_b_dx(i)(j)=-S_q_b.Column(j)/bending_state.n_b_magnitude;
                de_dx(i)(j)=(q_e(i)*TV::Axis_Vector(j))/bending_state.e_magnitude;
                if(!assume_constant_normal_length){
                    dn_a_dx(i)(j) -= (bending_state.n_a/bending_state.n_a_magnitude)*TV::Dot_Product(bending_state.n_a,-S_q_a.Column(j));
                    dn_b_dx(i)(j) -= (bending_state.n_b/bending_state.n_b_magnitude)*TV::Dot_Product(bending_state.n_b,-S_q_b.Column(j));
                    de_dx(i)(j) -= (bending_state.e/bending_state.e_magnitude)*TV::Dot_Product(bending_state.e,(q_e(i)*TV::Axis_Vector(j)));}
                dcos_theta_dx(i)(j)=TV::Dot_Product(dn_a_dx(i)(j),bending_state.n_b)+TV::Dot_Product(bending_state.n_a,dn_b_dx(i)(j));                
                dsin_theta_dx(i)(j)=TV::Dot_Product(TV::Cross_Product(dn_a_dx(i)(j),bending_state.n_b)+TV::Cross_Product(bending_state.n_a,dn_b_dx(i)(j)),bending_state.e)+
                    TV::Dot_Product(TV::Cross_Product(bending_state.n_a,bending_state.n_b),de_dx(i)(j));
                state.dC_dx(i)(j,1)=cos_theta*dsin_theta_dx(i)(j)-sin_theta*dcos_theta_dx(i)(j);}
            state.C_dot+=state.dC_dx(i).Transposed()*particles.V(state.nodes[i]);}
        for(int i=1;i<=4;i++){
            for(int j=1;j<=4;j++){
                for(int s=1;s<=3;s++) for(int t=1;t<=3;t++){
                    TV dn_a_dxs_dyt=-MATRIX<T,3>::Cross_Product_Matrix(dq_a_i_j_t(i,j)(t)).Column(s)/bending_state.n_a_magnitude;
                    TV dn_b_dxs_dyt=-MATRIX<T,3>::Cross_Product_Matrix(dq_b_i_j_t(i,j)(t)).Column(s)/bending_state.n_b_magnitude;
                    TV de_dx_dy=TV();
                    if(!assume_constant_normal_length){
                        MATRIX<T,3> S_q_a=MATRIX<T,3>::Cross_Product_Matrix(q_a(i));
                        MATRIX<T,3> S_q_b=MATRIX<T,3>::Cross_Product_Matrix(q_b(i));
                        MATRIX<T,3> S_q_a_j=MATRIX<T,3>::Cross_Product_Matrix(q_a(j));
                        MATRIX<T,3> S_q_b_j=MATRIX<T,3>::Cross_Product_Matrix(q_b(j));
                        dn_a_dxs_dyt=-(MATRIX<T,3>::Outer_Product(dn_a_dx(j)(t),bending_state.n_a)+
                            MATRIX<T,3>::Outer_Product(bending_state.n_a,dn_a_dx(j)(t)))*(-S_q_a.Column(s)/bending_state.n_a_magnitude)+
                            ((MATRIX<T,3>::Identity_Matrix()-MATRIX<T,3>::Outer_Product(bending_state.n_a,bending_state.n_a))/sqr(bending_state.n_a_magnitude))*
                            (bending_state.n_a_magnitude*(-MATRIX<T,3>::Cross_Product_Matrix(dq_a_i_j_t(i,j)(t)).Column(s))-
                                (TV::Dot_Product(-S_q_a_j.Column(t),bending_state.n_a)*-S_q_a.Column(s)));
                        dn_b_dxs_dyt = -(MATRIX<T,3>::Outer_Product(dn_b_dx(j)(t),bending_state.n_b)+
                            MATRIX<T,3>::Outer_Product(bending_state.n_b,dn_b_dx(j)(t)))*(-S_q_b.Column(s)/bending_state.n_b_magnitude) +
                            ((MATRIX<T,3>::Identity_Matrix()-MATRIX<T,3>::Outer_Product(bending_state.n_b,bending_state.n_b))/sqr(bending_state.n_b_magnitude))*
                            (bending_state.n_b_magnitude*(-MATRIX<T,3>::Cross_Product_Matrix(dq_b_i_j_t(i,j)(t)).Column(s))-
                                (TV::Dot_Product(-S_q_b_j.Column(t),bending_state.n_b)*-S_q_b.Column(s)));
                        de_dx_dy = -(MATRIX<T,3>::Outer_Product((q_e(j)*TV::Axis_Vector(t))/bending_state.e_magnitude,bending_state.e)+
                            MATRIX<T,3>::Outer_Product(bending_state.e,(q_e(j)*TV::Axis_Vector(t))/bending_state.e_magnitude))*
                            (-MATRIX<T,3>::Cross_Product_Matrix(q_e(i)*TV::All_Ones_Vector()).Column(s)/bending_state.e_magnitude)-
                            ((MATRIX<T,3>::Identity_Matrix()-MATRIX<T,3>::Outer_Product(bending_state.e,bending_state.e))/sqr(bending_state.e_magnitude))*
                            q_e(i)*q_e(j)*bending_state.e(t)*TV::Axis_Vector(s);}
                    T first_term=TV::Dot_Product(TV::Cross_Product(dn_a_dxs_dyt,bending_state.n_b)+TV::Cross_Product(dn_a_dx(i)(s),dn_b_dx(j)(t))+TV::Cross_Product(dn_a_dx(j)(t),dn_b_dx(i)(s))+
                        TV::Cross_Product(bending_state.n_a,dn_b_dxs_dyt),bending_state.e);
                    d2sin_theta_dx_dy(i,j)(s,t)=first_term+TV::Dot_Product(TV::Cross_Product(dn_a_dx(i)(s),bending_state.n_b)+TV::Cross_Product(bending_state.n_a,dn_b_dx(i)(s)),de_dx(j)(t))+
                        TV::Dot_Product(TV::Cross_Product(dn_a_dx(j)(t),bending_state.n_b)+TV::Cross_Product(bending_state.n_a,dn_b_dx(j)(t)),de_dx(i)(s))+
                        TV::Dot_Product(TV::Cross_Product(bending_state.n_a,bending_state.n_b),de_dx_dy);
                    d2cos_theta_dx_dy(i,j)(t,s)=TV::Dot_Product(dn_a_dxs_dyt,bending_state.n_b)+TV::Dot_Product(dn_b_dx(j)(t),dn_a_dx(i)(s))+
                        TV::Dot_Product(dn_a_dx(j)(t),dn_b_dx(i)(s))+TV::Dot_Product(bending_state.n_a,dn_b_dxs_dyt);}}}

        for(int i=1;i<=4;i++){
            for(int j=1;j<=4;j++){
                MATRIX<T,3> dC_dxi_dxj;
                for(int s=1;s<=3;s++) for(int t=1;t<=3;t++)
                    dC_dxi_dxj(s,t)=dcos_theta_dx(j)(t)*dsin_theta_dx(i)(s)+cos_theta*d2sin_theta_dx_dy(i,j)(s,t)-dsin_theta_dx(j)(t)*dcos_theta_dx(i)(s)-sin_theta*d2cos_theta_dx_dy(i,j)(s,t);
                state.dC_dxi_dxj_times_C(i,j)=dC_dxi_dxj*state.C(1);
                state.dC_dxi_dxj_times_C_dot(i,j)=dC_dxi_dxj*state.C_dot(1);}}}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR BW_BENDING_FORCES<TV>::
Potential_Energy(int s,const T time) const
{
    return (T).5*states(s).stiffness_coefficient*states(s).C.Magnitude_Squared();
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR BW_BENDING_FORCES<TV>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    const_cast<BW_BENDING_FORCES<TV>* >(this)->Update_Position_Based_State(time,true);
    for(CONSTRAINT_ITERATOR iterator(force_simplices);iterator.Valid();iterator.Next()){int s=iterator.Data();
        potential_energy+=Potential_Energy(s,time);}
        return potential_energy;
}
//#####################################################################
// Function Create_BW_Bending_Force
//#####################################################################
template<class TV> BW_BENDING_FORCES<TV>* PhysBAM::
Create_BW_Bending_Force(PARTICLES<TV>& particles,TRIANGLE_MESH& triangle_mesh,const typename TV::SCALAR stiffness_coefficient_input,const typename TV::SCALAR damping_coefficient_input)
{
    BW_BENDING_FORCES<TV>* sf=new BW_BENDING_FORCES<TV>(particles,triangle_mesh,stiffness_coefficient_input,damping_coefficient_input);
    return sf;
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template class BW_BENDING_FORCES<VECTOR<T,3> >; \
    template BW_BENDING_FORCES<VECTOR<T,3> >* PhysBAM::Create_BW_Bending_Force<VECTOR<T,3> >(PARTICLES<VECTOR<T,3> >&,TRIANGLE_MESH&,const T,const T);

INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
