//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Jon Gretarsson, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER  
//##################################################################### 
//
// Inherited by EULER_1D, EULER_2D, and EULER_3D.
// Initializes the equation of state information.
//
//#####################################################################
#ifndef __EULER__
#define __EULER__    

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/BOUNDARY_OBJECT_EULER.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/BOUNDARY_OBJECT_SOLID_VELOCITY.h>
namespace PhysBAM{

template<class T_GRID>
class EULER
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::VECTOR_T TV;
    typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;typedef typename T_ARRAYS_BASE::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_BASE;
    enum {d=T_GRID::dimension+2};
public:
    EOS<T>* eos;
    BOUNDARY_UNIFORM<T_GRID,TV_DIMENSION>* boundary;
    CONSERVATION<T_GRID,d>* conservation;
    T cfl_number;
    VECTOR<bool,2*T_GRID::dimension> open_boundaries;
    bool use_solid_velocity_in_ghost_cells;

protected:
    bool cut_out_grid; // (1) cut out grid, (0) no cut out grid 
    T max_time_step;
    bool use_force;
    T gravity;
    TV downward_direction;
private:
    BOUNDARY_UNIFORM<T_GRID,TV_DIMENSION> boundary_default;
    CONSERVATION_ENO_LLF<T_GRID,d> conservation_default;
    EOS_GAMMA<T> eos_default;

protected:
    EULER();
    virtual ~EULER();

public:
    static VECTOR<T,1> Get_Velocity(const ARRAYS_ND_BASE<VECTOR<VECTOR<T,3>,1> >& U,const VECTOR<int,1>& cell)
    {return VECTOR<T,1>(U(cell)(2))/U(cell)(1);}

    static VECTOR<T,2> Get_Velocity(const ARRAYS_ND_BASE<VECTOR<VECTOR<T,4>,2> >& U,const VECTOR<int,2>& cell)
    {return VECTOR<T,2>(U(cell)(2),U(cell)(3))/U(cell)(1);}

    static VECTOR<T,3> Get_Velocity(const ARRAYS_ND_BASE<VECTOR<VECTOR<T,5>,3> >& U,const VECTOR<int,3>& cell)
    {return VECTOR<T,3>(U(cell)(2),U(cell)(3),U(cell)(4))/U(cell)(1);}

    static VECTOR<T,1> Get_Velocity(const VECTOR<T,3>& u)
    {return VECTOR<T,1>(u(2))/u(1);}

    static VECTOR<T,2> Get_Velocity(const VECTOR<T,4>& u)
    {return VECTOR<T,2>(u(2),u(3))/u(1);}

    static VECTOR<T,3> Get_Velocity(const VECTOR<T,5>& u)
    {return VECTOR<T,3>(u(2),u(3),u(4))/u(1);}

    static T Get_Velocity_Component(const T_ARRAYS_DIMENSION_BASE& U,const TV_INT& cell,const int axis)
    {assert(axis>=1 && axis<=T_GRID::dimension);return U(cell)(axis+1)/U(cell)(1);}

    static T Get_Velocity_Component(const TV_DIMENSION& U,const int axis)
    {assert(axis>=1 && axis<=T_GRID::dimension);return U(axis+1)/U(1);}

    static T Get_Density(const T_ARRAYS_DIMENSION_BASE& U,const TV_INT& cell)
    {return U(cell)(1);}

    static T Get_Total_Energy(const T_ARRAYS_DIMENSION_BASE& U,const TV_INT& cell)
    {return U(cell)(T_GRID::dimension+2);}

    static void Set_Euler_State_From_rho_velocity_And_internal_energy(T_ARRAYS_DIMENSION_BASE& U,const TV_INT& cell,const T rho,const TV& velocity,const T e)
    {U(cell)(1)=rho;
    for(int k=1;k<=T_GRID::dimension;k++) U(cell)(k+1)=rho*velocity[k];
    U(cell)(T_GRID::dimension+2)=rho*(e+velocity.Magnitude_Squared()*(T).5);}

    static VECTOR<T,T_GRID::dimension+2> Get_Euler_State_From_rho_velocity_And_internal_energy(const T rho,const TV& velocity,const T e)
    {VECTOR<T,T_GRID::dimension+2> U;
     U(1)=rho;
     for(int k=1;k<=T_GRID::dimension;k++) U(k+1)=rho*velocity[k];
     U(T_GRID::dimension+2)=rho*(e+(T).5*velocity.Magnitude_Squared());
     return U;}

    static T e(const ARRAYS_ND_BASE<VECTOR<VECTOR<T,3>,1> >& U,const VECTOR<int,1>& cell)
    {return e(U(cell)(1),U(cell)(2),U(cell)(3));}

    static T e(const ARRAYS_ND_BASE<VECTOR<VECTOR<T,4>,2> >& U,const VECTOR<int,2>& cell)
    {return e(U(cell)(1),U(cell)(2),U(cell)(3),U(cell)(4));}

    static T e(const ARRAYS_ND_BASE<VECTOR<VECTOR<T,5>,3> >& U,const VECTOR<int,3>& cell)
    {return e(U(cell)(1),U(cell)(2),U(cell)(3),U(cell)(4),U(cell)(5));}

    static T e(const VECTOR<T,3>& u)
    {return e(u(1),u(2),u(3));}

    static T e(const VECTOR<T,4>& u)
    {return e(u(1),u(2),u(3),u(4));}

    static T e(const VECTOR<T,5>& u)
    {return e(u(1),u(2),u(3),u(4),u(5));}

    static T e(const T rho,const T rho_u,const T E)
    {return E/rho-sqr(rho_u/rho)/2;}

    static T e(const T rho,const T rho_u,const T rho_v,const T E)
    {return E/rho-(sqr(rho_u/rho)+sqr(rho_v/rho))/2;}

    static T e(const T rho,const T rho_u,const T rho_v,const T rho_w,const T E)
    {return E/rho-(sqr(rho_u/rho)+sqr(rho_v/rho)+sqr(rho_w/rho))/2;}

    static T p(EOS<T>* eos_input,const TV_DIMENSION& u)
    {return eos_input->p(u(1),e(u));}

    static T entropy(EOS<T>* eos,const TV_DIMENSION& u)
    {T internal_energy=e(u);T rho=u(1);return eos->S(rho,internal_energy);}

    static T enthalpy(const EOS<T>& eos,const VECTOR<T,T_GRID::dimension+2>& u)
    {T internal_energy=e(u);T rho=u(1);T p=eos.p(rho,internal_energy);return internal_energy+p/rho;}

    T enthalpy(const VECTOR<T,T_GRID::dimension+2>& u) const
    {T internal_energy=e(u);T rho=u(1);T p=eos->p(rho,internal_energy);return internal_energy+p/rho;}

    void Set_Custom_Boundary(BOUNDARY_UNIFORM<T_GRID,TV_DIMENSION>& boundary_input)
    {boundary=&boundary_input;
    for(int axis=1;axis<=T_GRID::dimension;axis++){
        open_boundaries(2*axis-1)=boundary->Constant_Extrapolation(2*axis-1);
        open_boundaries(2*axis)=boundary->Constant_Extrapolation(2*axis);}}
    
    void Set_Custom_Conservation(CONSERVATION<T_GRID,d>& conservation_input)
    {conservation=&conservation_input;
        if(use_solid_velocity_in_ghost_cells) conservation->Set_Custom_Object_Boundary(*new BOUNDARY_OBJECT_SOLID_VELOCITY<T_GRID>);
        else conservation->Set_Custom_Object_Boundary(*new BOUNDARY_OBJECT_EULER<T_GRID>);}

    void Set_Custom_Equation_Of_State(EOS<T>& eos_input)
    {eos=&eos_input;}

    void Set_Max_Time_Step(const T max_time_step_input=1e8)
    {max_time_step=max_time_step_input;}

    void Set_Gravity(const T gravity_input=9.8)
    {gravity=gravity_input;downward_direction=T_GRID::dimension>1?-TV::Axis_Vector(2):TV();}

    void Set_Gravity(const T gravity_input,const TV& downward_direction_input)
    {gravity=gravity_input;downward_direction=downward_direction_input;}

    void Set_CFL_Number(const T cfl_number_input=.5)
    {cfl_number=cfl_number_input;}

//#####################################################################
    virtual void Log_Parameters() const;
//#####################################################################
};
//#####################################################################
}
#endif
