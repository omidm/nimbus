//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BW_FORCES
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BW_FORCES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <cfloat>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d,int m> BW_FORCES<TV,d,m>::
BW_FORCES(PARTICLES<TV>& particles,TRIANGLE_MESH& triangle_mesh_input,const T stiffness_coefficient_input,const T damping_coefficient_input)
    :DEFORMABLES_FORCES<TV>(particles),triangle_mesh(triangle_mesh_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d,int m> BW_FORCES<TV,d,m>::
~BW_FORCES()
{
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV,int d,int m> void BW_FORCES<TV,d,m>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(SIMPLEX_ITERATOR iterator(force_simplices);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        // f_not
        for(int i=1;i<=m;i++) F(state.nodes[i])+=-state.stiffness_coefficient*state.dC_dx(i)*state.C;

        // f_not (damping)
        for(int i=1;i<=m;i++) F(state.nodes[i])+=-state.damping_coefficient*state.dC_dx(i)*state.C_dot;}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV,int d,int m> void BW_FORCES<TV,d,m>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(SIMPLEX_ITERATOR iterator(force_simplices);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        // df/dv (damping)
        for(int i=1;i<=m;i++) for(int j=1;j<=m;j++){
            MATRIX<T,3> K=-state.damping_coefficient*(state.dC_dx(i)*state.dC_dx(j).Transposed());
            F(state.nodes[i])+=(K*V(state.nodes[j]));}}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV,int d,int m> void BW_FORCES<TV,d,m>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(SIMPLEX_ITERATOR iterator(force_simplices);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        // df_dx * v_0
        for(int i=1;i<=m;i++) for(int j=1;j<=m;j++){
            MATRIX<T,3> K=-state.stiffness_coefficient*(state.dC_dx(i)*state.dC_dx(j).Transposed()+state.dC_dxi_dxj_times_C(i,j));
            F(state.nodes[i])+=K*V(state.nodes[j]);}

        // df_dx * v_0 (damping)
        for(int i=1;i<=m;i++) for(int j=1;j<=m;j++){
            MATRIX<T,3> K=-state.damping_coefficient*(state.dC_dxi_dxj_times_C_dot(i,j));
            F(state.nodes[i])+=K*V(state.nodes[j]);}}
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV,int d,int m> void BW_FORCES<TV,d,m>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV,int d,int m> typename TV::SCALAR BW_FORCES<TV,d,m>::
CFL_Strain_Rate() const
{
    return 0;
}
//#####################################################################
#define INSTANTIATION_HELPER(T,d,m) \
    template class BW_FORCES<VECTOR<T,3>,d,m>;

INSTANTIATION_HELPER(float,1,3)
INSTANTIATION_HELPER(float,2,3)
INSTANTIATION_HELPER(float,1,4)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1,3)
INSTANTIATION_HELPER(double,2,3)
INSTANTIATION_HELPER(double,1,4)
#endif
