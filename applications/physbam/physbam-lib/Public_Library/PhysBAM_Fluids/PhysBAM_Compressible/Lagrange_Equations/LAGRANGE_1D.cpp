//#####################################################################
// Copyright 2002, Ronald Fedkiw
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/LAGRANGE_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
// providing time is optional, mass is never updated
template<class T> void LAGRANGE_1D<T>::
Euler_Step(const T dt,const T time)
{  
    int i;
    int m=grid.m;

    // find nodal masses
    ARRAY<T,VECTOR<int,1> > M_node(1,m);
    for(i=1;i<=m-1;i++){T mass_over_2=mass(i)/2;M_node(i)+=mass_over_2;M_node(i+1)+=mass_over_2;}
    
    // get grid lengths
    ARRAY<T,VECTOR<int,1> > L(1,m-1);grid.Get_Lengths(L);
    
    // forces F, and heating rates H = -F dot V 
    ARRAY<T,VECTOR<int,1> > F(1,m),H(1,m-1);

    // pressue
    for(i=1;i<=m-1;i++){
        T pressure=eos.p(mass(i)/L(i),energy(i));
        F(i)+=-pressure;H(i)+=pressure*velocity(i);
        F(i+1)+=pressure;H(i)+=-pressure*velocity(i+1);}

    // boundary conditions - external forces 
    for(i=1;i<=m;i++) F(i)+=external_force(i);

    // artificial viscosity  
    ARRAY<T,VECTOR<int,1> > Q(1,m-1);
    artificial_viscosity->Get_Artificial_Viscosity(eos,grid,mass,velocity,energy,Q);
    for(i=1;i<=m-1;i++){
        F(i)+=-Q(i);H(i)+=Q(i)*velocity(i);
        F(i+1)+=Q(i);H(i)+=-Q(i)*velocity(i+1);}

    // material strength
    if(material_strength)
        for(i=1;i<=m-1;i++){
            T S=-stiffness*(L(i)/L_not(i)-1);
            F(i)+=-S;H(i)+=S*velocity(i);
            F(i+1)+=S;H(i)+=-S*velocity(i+1);}

    // update grid
    grid.Euler_Step(velocity,dt);

    // update velocity
    for(i=1;i<=m;i++) velocity(i)+=dt*F(i)/M_node(i);

    // boundary conditions - external velocities 
    for(i=1;i<=m;i++) if(fixed_velocity(i)) velocity(i)=external_velocity(i);
    
    // update interal energy 
    for(i=1;i<=m-1;i++) energy(i)+=dt*H(i)/mass(i);
}
//#####################################################################
// Function CFL
//#####################################################################
// includes the effects of artificial viscosity
template<class T> T LAGRANGE_1D<T>::
CFL()
{
    int i;
    int m=grid.m;

    ARRAY<T,VECTOR<int,1> > L(1,m-1);grid.Get_Lengths(L);
    
    ARRAY<T,VECTOR<int,1> > Q(1,m-1);
    artificial_viscosity->Get_Artificial_Viscosity(eos,grid,mass,velocity,energy,Q);

    ARRAY<T,VECTOR<int,1> > S(1,m-1);
    if(material_strength) for(i=1;i<=m-1;i++) S(i)=-stiffness*(L(i)/L_not(i)-1);

    ARRAY<T,VECTOR<int,1> > timestep(1,m-1);
    for(i=1;i<=m-1;i++){
        T density=mass(i)/L(i);         
        T P=eos.p(density,energy(i))+Q(i)+abs(S(i));
        T E=eos.e_From_p_And_rho(P,density);
        timestep(i)=eos.c(density,E)/L(i);}

    T dt_convect=timestep.Max();
    return 1/dt_convect;
}
//#####################################################################
template class LAGRANGE_1D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LAGRANGE_1D<double>;
#endif
