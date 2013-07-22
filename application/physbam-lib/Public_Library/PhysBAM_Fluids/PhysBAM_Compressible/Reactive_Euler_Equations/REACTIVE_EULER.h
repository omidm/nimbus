//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REACTIVE_EULER.h  
//##################################################################### 
//
// Inherited by REACTIVE_EULER_1D, REACTIVE_EULER_2D, and REACTIVE_EULER_3D.
// Initializes the equation of state information.
//
//#####################################################################
#ifndef __REACTIVE_EULER__
#define __REACTIVE_EULER__    

#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/REACTIVE_EOS.h>
namespace PhysBAM{

template<class T_GRID>
class REACTIVE_EULER
{
    typedef typename T_GRID::SCALAR T;typedef VECTOR<T,T_GRID::dimension+3> TV_DIMENSION;
public:
    REACTIVE_EOS<T>& eos; // needed for equation of state functions
    BOUNDARY_UNIFORM<T_GRID,TV_DIMENSION>* boundary;
    CONSERVATION<T_GRID,T_GRID::dimension+3>* conservation;
protected:
    int cut_out_grid; // (1) cut out grid, (0) no cut out grid 
private:
    BOUNDARY_UNIFORM<T_GRID,TV_DIMENSION> boundary_default;
    CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+3> conservation_default;

protected:
    REACTIVE_EULER(REACTIVE_EOS<T>& eos_input)
        :eos(eos_input)
    {
        cut_out_grid=0;
        boundary=&boundary_default;
        conservation=&conservation_default;
    }

public:
    void Set_Custom_Boundary(BOUNDARY_UNIFORM<T_GRID,TV_DIMENSION>& boundary_input)
    {boundary=&boundary_input;}
    
    void Set_Custom_Conservation(CONSERVATION<T_GRID,T_GRID::dimension+3>& conservation_input)
    {conservation=&conservation_input;}

//#####################################################################
    T e(const T rho,const T rho_u,const T E);
    T e(const T rho,const T rho_u,const T rho_v,const T E);
    T e(const T rho,const T rho_u,const T rho_v,const T rho_w,const T E);
//#####################################################################
};   
}
#endif
