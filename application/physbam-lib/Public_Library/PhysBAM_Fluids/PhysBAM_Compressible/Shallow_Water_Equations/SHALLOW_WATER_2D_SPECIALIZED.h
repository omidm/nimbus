//#####################################################################
// Copyright 2002-2004, Ron Fedkiw, Eran Guendelman, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_2D_SPECIALIZED  
//##################################################################### 
//
// Input U as 3 by (1,m) by (1,n) for h, h*u, and h*v.
//
//#####################################################################
#ifndef __SHALLOW_WATER_2D_SPECIALIZED__
#define __SHALLOW_WATER_2D_SPECIALIZED__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Shallow_Water_Equations/SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Shallow_Water_Equations/SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_G.h>
namespace PhysBAM{

template<class T>
class SHALLOW_WATER_2D_SPECIALIZED
{
    typedef VECTOR<T,2> TV;typedef VECTOR<T,3> TV_DIMENSION;
    enum {d=3};
public:
    BOUNDARY_UNIFORM<GRID<TV>,TV_DIMENSION>* boundary;
    CONSERVATION<GRID<TV>,2>* conservation;
private:
    BOUNDARY_UNIFORM<GRID<TV>,TV_DIMENSION> boundary_default;
    CONSERVATION_ENO_LLF<GRID<TV>,2> conservation_default;
    CONSERVATION_ENO_LLF<GRID<TV>,3> internal_conservation;
public:
    T gravity;
    T min_height;
    T epsilon_u,epsilon_v;
    GRID<TV>& grid;
    ARRAY<TV_DIMENSION,VECTOR<int,2> >& U; // h, u, and v
    ARRAY<T,VECTOR<int,2> >* ground_ghost;
    ARRAY<T,VECTOR<int,2> > eta_ghost;
    ARRAY<VECTOR<T,6> ,VECTOR<int,2> > postprocessed_flux_x,postprocessed_flux_y; // for debugging
    SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F<T> eigensystem_F;
    SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_G<T> eigensystem_G;

    SHALLOW_WATER_2D_SPECIALIZED(GRID<TV>& grid_input,ARRAY<TV_DIMENSION,VECTOR<int,2> >& U_input,const T gravity_input=9.8,const T min_height_input=1e-3)
        :boundary(&boundary_default),conservation(&conservation_default),
        gravity(gravity_input),min_height(min_height_input),grid(grid_input),U(U_input),ground_ghost(0),eta_ghost(-2,grid.counts.x+3,-2,grid.counts.y+3),
        postprocessed_flux_x(0,grid.counts.x,1,grid.counts.y),postprocessed_flux_y(1,grid.counts.x,0,grid.counts.y),
        eigensystem_F(gravity_input,eta_ghost,min_height_input),eigensystem_G(gravity_input,eta_ghost,min_height_input)
    {
        Set_Viscosity_Coefficients(0,0);
    }

    void Set_Custom_Boundary(BOUNDARY_UNIFORM<GRID<TV>,TV_DIMENSION>& boundary_input)
    {boundary=&boundary_input;}

    void Set_Custom_Conservation(CONSERVATION<GRID<TV>,2>& conservation_input)
    {conservation=&conservation_input;}
    
    void Set_Ground_Ghost(ARRAY<T,VECTOR<int,2> >& ground_ghost_input)
    {ground_ghost=&ground_ghost_input;}

    void Set_Viscosity_Coefficients(const T epsilon_u_input,const T epsilon_v_input)
    {epsilon_u=epsilon_u_input;epsilon_v=epsilon_v_input;}
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    T CFL();
//#####################################################################
};   
}
#endif
