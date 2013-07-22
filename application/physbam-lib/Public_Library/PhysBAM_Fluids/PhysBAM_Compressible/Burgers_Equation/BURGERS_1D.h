//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BURGERS_1D  
//##################################################################### 
//
// Input U as 1 by (1,m).
//
//#####################################################################
#ifndef __BURGERS_1D__
#define __BURGERS_1D__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Burgers_Equation/BURGERS_1D_EIGENSYSTEM_F.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <cfloat>
namespace PhysBAM{

template<class T>
class BURGERS_1D
{
    typedef VECTOR<T,1> TV;typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
public:
    BOUNDARY_UNIFORM<GRID<TV>,TV>* boundary;
    CONSERVATION<GRID<TV>,1>* conservation;
protected:
    GRID<TV>& grid;
    ARRAY<TV,VECTOR<int,1> >& U;
    BOUNDARY_UNIFORM<GRID<TV>,TV> boundary_default;
    CONSERVATION_ENO_LLF<GRID<TV>,1> conservation_default;
    BURGERS_1D_EIGENSYSTEM_F<T> eigensystem_F;
private:
     T max_time_step;

public:
    BURGERS_1D(GRID<TV>& grid_input,ARRAY<TV,VECTOR<int,1> >& U_input)
        :grid(grid_input),U(U_input)
    {
        boundary=&boundary_default;
        conservation=&conservation_default;
        Set_Max_Time_Step();
    }
    
    void Set_Custom_Boundary(BOUNDARY_UNIFORM<GRID<TV>,TV>& boundary_input)
    {boundary=&boundary_input;}
    
    void Set_Custom_Conservation(CONSERVATION<GRID<TV>,1>& conservation_input)
    {conservation=&conservation_input;}
    
    void Set_Max_Time_Step(const T max_time_step_input=FLT_MAX)
    {max_time_step=max_time_step_input;}
    
//#####################################################################
    void Euler_Step(const T dt,const T time=0);
    T CFL();
//#####################################################################
};   
}
#endif
