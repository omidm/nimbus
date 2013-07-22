//#####################################################################
// Copyright 2002-2003, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER
//##################################################################### 
//
// Inherited by SHALLOW_WATER_1D, SHALLOW_WATER_2D, and SHALLOW_WATER_3D.
//
//#####################################################################
#ifndef __SHALLOW_WATER__
#define __SHALLOW_WATER__    

#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
namespace PhysBAM{

template<class T_GRID>
class SHALLOW_WATER
{
    typedef typename T_GRID::SCALAR T;typedef VECTOR<T,T_GRID::dimension+1> TV_DIMENSION;
public:
    BOUNDARY_UNIFORM<T_GRID,TV_DIMENSION>* boundary;
    CONSERVATION<T_GRID,T_GRID::dimension+1>* conservation;
private:
    BOUNDARY_UNIFORM<T_GRID,TV_DIMENSION> boundary_default;
    CONSERVATION_ENO_LLF<T_GRID,T_GRID::dimension+1> conservation_default;

protected:
    SHALLOW_WATER()
    {
        boundary=&boundary_default;
        conservation=&conservation_default;
    }
public:

    void Set_Custom_Boundary(BOUNDARY_UNIFORM<T_GRID,TV_DIMENSION>& boundary_input)
    {boundary=&boundary_input;}
    
    void Set_Custom_Conservation(CONSERVATION<T,T_GRID::dimension+1>& conservation_input)
    {conservation=&conservation_input;}

//#####################################################################
};
}    
#endif

