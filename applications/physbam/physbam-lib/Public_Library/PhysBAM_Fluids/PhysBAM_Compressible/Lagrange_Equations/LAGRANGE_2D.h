//#####################################################################
// Copyright 2002, Ronald Fedkiw
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAGRANGE_2D  
//##################################################################### 
//
// Input an EOS class, that is a class that inherits the virtual base class EOS.
// Input a GRID_LAGRANGE_2D class.
// Input mass as (1,m-1,1,n-1) as the mass of each cell.
// Input velcoity as (1,m,1,n) as the velocity of each node.
// Input energy as (1,m-1,1,n-1) as the internal energy of each cell.
//
// Use external_force and external_velocity to set boundary conditions.  
//
//#####################################################################
#ifndef __LAGRANGE_2D__
#define __LAGRANGE_2D__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/ARTIFICIAL_VISCOSITY_2D.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/ARTIFICIAL_VISCOSITY_VNR_2D.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/GRID_LAGRANGE_2D.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/LAGRANGE.h>
namespace PhysBAM{

template<class T>
class LAGRANGE_2D:public LAGRANGE<T>
{
public:
    using LAGRANGE<T>::eos;

    GRID_LAGRANGE_2D<T>& grid; // lagrangian grid
    ARRAY<T,VECTOR<int,2> >& mass;    // mass of each cell
    ARRAY<T,VECTOR<int,2> >& u;          // x-velocity of each node
    ARRAY<T,VECTOR<int,2> >& v;          // y-velocity of each node
    ARRAY<T,VECTOR<int,2> >& energy; // internal energy of each cell
    ARRAY<T,VECTOR<int,2> > external_force_x; // boundary condition
    ARRAY<T,VECTOR<int,2> > external_force_y; // boundary condition
    ARRAY<int,VECTOR<int,2> > fixed_velocity;       // (1) use external velocity, (0) compute velocity 
    ARRAY<T,VECTOR<int,2> > external_u;          // boundary condition
    ARRAY<T,VECTOR<int,2> > external_v;          // boundary condition
    ARTIFICIAL_VISCOSITY_2D<T>* artificial_viscosity;
    ARTIFICIAL_VISCOSITY_VNR_2D<T> artificial_viscosity_default; 
    int material_strength;    // (0) no, (1) material strength on
    T stiffness;               // for material strength
    ARRAY<T,VECTOR<int,2> > L1_not,L2_not; // length at which there is no restoring force

public:
    LAGRANGE_2D(EOS<T>& eos_input,GRID_LAGRANGE_2D<T>& grid_input,ARRAY<T,VECTOR<int,2> >& mass_input,ARRAY<T,VECTOR<int,2> >& u_input,ARRAY<T,VECTOR<int,2> >& v_input,ARRAY<T,VECTOR<int,2> >& energy_input)
        :LAGRANGE<T>(eos_input),grid(grid_input),mass(mass_input),u(u_input),v(v_input),energy(energy_input),external_force_x(1,grid.m,1,grid.n),external_force_y(1,grid.m,1,grid.n),
        fixed_velocity(1,grid.m,1,grid.n),external_u(1,grid.m,1,grid.n),external_v(1,grid.m,1,grid.n),L1_not(1,grid.m-1,1,grid.n),L2_not(1,grid.m,1,grid.n-1)
    {
        material_strength=0;
        artificial_viscosity=&artificial_viscosity_default;
    }
    
    void Set_Custom_Artificial_Viscosity(ARTIFICIAL_VISCOSITY_2D<T>& artificial_viscosity_input)
    {artificial_viscosity=&artificial_viscosity_input;}

    void Set_Material_Strength(const T stiffness_input)
    {material_strength=1;stiffness=stiffness_input;grid.Get_Lengths(L1_not,L2_not);}
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    T CFL();
//#####################################################################
};
}    
#endif
