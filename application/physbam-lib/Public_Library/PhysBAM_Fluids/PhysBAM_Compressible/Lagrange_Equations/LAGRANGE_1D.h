//#####################################################################
// Copyright 2002, Ronald Fedkiw
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LAGRANGE_1D  
//##################################################################### 
//
// Input an EOS class, that is a class that inherits the virtual base class EOS.
// Input a GRID_LAGRANGE_1D class.
// Input mass as (1,m-1) as the mass of each cell.
// Input velcoity as (1,m) as the velocity of each node.
// Input energy as (1,m-1) as the internal energy of each cell.
//
// Use external_force and external_velocity to set boundary conditions. 
//
//#####################################################################
#ifndef __LAGRANGE_1D__
#define __LAGRANGE_1D__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/ARTIFICIAL_VISCOSITY_1D.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/ARTIFICIAL_VISCOSITY_VNR_1D.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/GRID_LAGRANGE_1D.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/LAGRANGE.h>
namespace PhysBAM{

template<class T>
class LAGRANGE_1D:public LAGRANGE<T>
{
public:
    using LAGRANGE<T>::eos;

    GRID_LAGRANGE_1D<T>& grid; // lagrangian grid
    ARRAY<T,VECTOR<int,1> >& mass;     // mass of each cell
    ARRAY<T,VECTOR<int,1> >& velocity; // velocity of each node
    ARRAY<T,VECTOR<int,1> >& energy;   // internal energy of each cell
    ARRAY<T,VECTOR<int,1> > external_force;     // boundary condition
    ARRAY<int,VECTOR<int,1> > fixed_velocity;        // (1) use external velocity, (0) compute velocity 
    ARRAY<T,VECTOR<int,1> > external_velocity; // boundary condition
    ARTIFICIAL_VISCOSITY_1D<T>* artificial_viscosity;
    ARTIFICIAL_VISCOSITY_VNR_1D<T> artificial_viscosity_default;
    int material_strength;    // (0) no, (1) material strength on
    T stiffness;            // for material strength
    ARRAY<T,VECTOR<int,1> > L_not; // length at which there is no restoring force

public:
    LAGRANGE_1D(EOS<T>& eos_input,GRID_LAGRANGE_1D<T>& grid_input,ARRAY<T,VECTOR<int,1> >& mass_input,ARRAY<T,VECTOR<int,1> >& velocity_input,ARRAY<T,VECTOR<int,1> >& energy_input)
        :LAGRANGE<T>(eos_input),grid(grid_input),mass(mass_input),velocity(velocity_input),energy(energy_input),external_force(1,grid.m),
        fixed_velocity(1,grid.m),external_velocity(1,grid.m),L_not(1,grid.m-1)
    {
        material_strength=0;
        artificial_viscosity=&artificial_viscosity_default;
    }

    void Set_Custom_Artificial_Viscosity(ARTIFICIAL_VISCOSITY_1D<T>& artificial_viscosity_input)
    {artificial_viscosity=&artificial_viscosity_input;}

    void Set_Material_Strength(const T stiffness_input)
    {material_strength=1;stiffness=stiffness_input;grid.Get_Lengths(L_not);}

//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    T CFL();
//#####################################################################
}; 
}  
#endif
