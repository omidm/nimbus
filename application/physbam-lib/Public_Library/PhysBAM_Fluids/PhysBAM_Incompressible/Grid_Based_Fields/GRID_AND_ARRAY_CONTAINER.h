//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_AND_ARRAY_CONTAINER
//#####################################################################
#ifndef __GRID_AND_ARRAY_CONTAINER__    
#define __GRID_AND_ARRAY_CONTAINER__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Grids_Dyadic_Advection/ADVECTION_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_RLE_Advection/ADVECTION_POLICY_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_UTILITIES.h>
namespace PhysBAM{

template<class T_GRID> struct BOUNDARY_POLICY;
template<class T_GRID> struct GRID_ARRAYS_POLICY;

template<class T_GRID,class T2>
class GRID_AND_ARRAY_CONTAINER:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,T_GRID::dimension> TV_SIDES;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename REBIND<T_ARRAYS_SCALAR,T2>::TYPE T_ARRAYS_T2;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename REBIND<T_ARRAYS_SCALAR,TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename ADVECTION_POLICY<T_GRID>::ADVECTION_SEMI_LAGRANGIAN_SCALAR T_ADVECTION_SEMI_LAGRANGIAN_SCALAR;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY_SCALAR;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_REFLECTION T_BOUNDARY_REFLECTION;
public:
    T_GRID& grid;
    T_ARRAYS_T2 array;
    ADVECTION<T_GRID,T>* advection;
    ADVECTION<T_GRID,T>* advection_maccormack;
    T_BOUNDARY_SCALAR* boundary;
private:
    T_ADVECTION_SEMI_LAGRANGIAN_SCALAR& advection_default;
protected:
    T_BOUNDARY_SCALAR& boundary_default; 
    const T_FACE_ARRAYS_SCALAR* face_velocities;
    const T_ARRAYS_VECTOR* cell_velocities;
public:

    GRID_AND_ARRAY_CONTAINER(T_GRID& grid_input);
    virtual ~GRID_AND_ARRAY_CONTAINER();

    void Clean_Memory()
    {array.Clean_Memory();}
    
    virtual void Initialize_Array(const int ghost_cells=0,const bool initialize_new_elements=true,const bool copy_existing_elements=true)
    {array.Resize(grid.Cell_Indices(ghost_cells),initialize_new_elements,copy_existing_elements);}
  
    void Initialize_Domain_Boundary_Conditions(const TV_SIDES& domain_walls=TV_SIDES::Constant_Vector(TV_BOOL2::Constant_Vector(true)))
    {boundary->Set_Constant_Extrapolation(VECTOR_UTILITIES::Complement(domain_walls));}

    void Set_Custom_Boundary(T_BOUNDARY_SCALAR& boundary_input)
    {boundary=&boundary_input;}

    void Set_Custom_Advection(ADVECTION<T_GRID,T>& advection_input)
    {advection=&advection_input;}

    void Set_To_Constant_Value(const T2& value)
    {array.Fill(value);}

    void Set_Velocity(const T_FACE_ARRAYS_SCALAR* face_velocities_input)
    {face_velocities=face_velocities_input;}
    
    void Set_Velocity(const T_ARRAYS_VECTOR* cell_velocities_input)
    {cell_velocities=cell_velocities_input;}

    template<class T_ARRAYS_BOOL>
    void Use_Maccormack_Advection(const T_ARRAYS_BOOL& cell_mask)
    {advection_maccormack=new ADVECTION_MACCORMACK_UNIFORM<T_GRID,T,ADVECTION<T_GRID,T> >(*advection,0,&cell_mask,0);
    Set_Custom_Advection(*advection_maccormack);}

//#####################################################################
    virtual void Euler_Step(const T dt,const T time,const int number_of_ghost_cells);
//#####################################################################
};
}
#endif
