//#####################################################################
// Copyright 2002-2009, Doug Enright, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Michael Lentine, Frank Losasso, Duc Nguyen, Nick Rasmussen, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INCOMPRESSIBLE_FLUID_EVOLUTION  
//#####################################################################
#ifndef __INCOMPRESSIBLE_FLUID_EVOLUTION__
#define __INCOMPRESSIBLE_FLUID_EVOLUTION__

#include <PhysBAM_Tools/Advection/ADVECTION_FORWARD.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Utilities/Find_Type.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{
template<class T_GRID> class INCOMPRESSIBLE_FLUIDS_FORCES;
template<class T_GRID> class BOUNDARY_MAC_GRID_SOLID_WALL_SLIP;
template<class T_GRID,class T> class BOUNDARY_UNIFORM;

template<class T_GRID>
class INCOMPRESSIBLE_FLUID_EVOLUTION:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
public:
    
    T_GRID grid;
    ADVECTION<T_GRID,T>* advection;
    BOUNDARY_UNIFORM<T_GRID,T>* boundary;
    ARRAY<INCOMPRESSIBLE_FLUIDS_FORCES<T_GRID>*> fluids_forces;
    T max_time_step;
protected:               
    BOUNDARY_MAC_GRID_SOLID_WALL_SLIP<T_GRID>& boundary_default;
public:

    INCOMPRESSIBLE_FLUID_EVOLUTION(const T_GRID& grid_input);
    virtual ~INCOMPRESSIBLE_FLUID_EVOLUTION();

    void Set_Custom_Advection(ADVECTION<T_GRID,T>& advection_input)
    {advection=&advection_input;}

    void Set_Custom_Boundary(BOUNDARY_UNIFORM<T_GRID,T>& boundary_input)
    {boundary=&boundary_input;}

    template<class T_FORCE> T_FORCE
    Find_Force(const int index=1)
    {return Find_Type<T_FORCE>(fluids_forces,index);}

    template<class T_FORCE> const T_FORCE
    Find_Force(const int index=1) const
    {return Find_Type<T_FORCE>(fluids_forces,index);}

//#####################################################################
    void Advance_One_Time_Step_Convection(const T dt,const T time,const T_FACE_ARRAYS_SCALAR& advecting_face_velocities,T_FACE_ARRAYS_SCALAR& face_velocities_to_advect,const int number_of_ghost_cells);
    void Advance_One_Time_Step_Forces(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const int number_of_ghost_cells);
    void Advance_One_Time_Step_Implicit_Part(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Initialize_Grids(const T_GRID& grid_input);
    T CFL(T_FACE_ARRAYS_SCALAR& face_velocities) const;
    int Add_Force(INCOMPRESSIBLE_FLUIDS_FORCES<T_GRID>* force);
//#####################################################################
};
}
#endif

