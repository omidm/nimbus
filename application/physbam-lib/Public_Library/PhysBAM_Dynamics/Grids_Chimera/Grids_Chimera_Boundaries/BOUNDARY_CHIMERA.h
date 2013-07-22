//#####################################################################
// Copyright 2011, Linhai Qiu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_CHIMERA
//#####################################################################
#ifndef __BOUNDARY_CHIMERA__
#define __BOUNDARY_CHIMERA__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Boundaries/BOUNDARY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GRID.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>

namespace PhysBAM{

template<class T_GRID,class T2> class CHIMERA_GRID;

template<class T_GRID,class T2=typename T_GRID::SCALAR>
class BOUNDARY_CHIMERA:public REBIND<typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR,T2>::TYPE
{
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::VECTOR_T TV;
    typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,T_GRID::dimension> TV_SIDES;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;typedef typename T_ARRAYS_BASE::template REBIND<T2>::TYPE T_ARRAYS_T2;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename REBIND<T_FACE_ARRAYS_SCALAR,T2>::TYPE T_FACE_ARRAYS_T2;typedef typename REBIND<T_FACE_ARRAYS_SCALAR,int>::TYPE T_FACE_ARRAYS_INT;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename REBIND<typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR,T2>::TYPE T_BOUNDARY_T2;
    typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;

public:
    CHIMERA_GRID<T_GRID,T2>* chimera_grid;
    T_BOUNDARY_T2& boundary;
public:
    BOUNDARY_CHIMERA(CHIMERA_GRID<T_GRID,T2>* chimera_grid_input,T_BOUNDARY_T2& boundary_input);
    virtual ~BOUNDARY_CHIMERA();
//#####################################################################
    void Set_Constant_Extrapolation(const TV_SIDES& constant_extrapolation_input=TV_SIDES::Constant_Vector(TV_BOOL2::Constant_Vector(true)));
    bool Constant_Extrapolation(const int side) const PHYSBAM_OVERRIDE;
    void Set_Fixed_Boundary(const bool use_fixed_boundary_input=true,const T2 fixed_boundary_value_input=T2());
    void Apply_Boundary_Condition(const T_GRID& grid,T_ARRAYS_T2& u,const T time) PHYSBAM_OVERRIDE;
    void Apply_Boundary_Condition_Face(const T_GRID& grid,T_FACE_ARRAYS_T2& u,const T time) PHYSBAM_OVERRIDE;
    void Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells_input) PHYSBAM_OVERRIDE;
    void Fill_Ghost_Cells_Face(const T_GRID& grid,const T_FACE_ARRAYS_T2& u,T_FACE_ARRAYS_T2& u_ghost,const T time,const int number_of_ghost_cells_input) PHYSBAM_OVERRIDE;
    static T Interpolation_From_Delaunay_Triangle_Vertices(const VECTOR<T,1> location,const ARRAY<VECTOR<T,1> > points,const ARRAY<T> values,bool& selected,const T tolerance=1e-6,const ARRAY<T>* const phi_values=0,const T interface_value=0);
    static T Interpolation_From_Delaunay_Triangle_Vertices(const VECTOR<T,2> location,const ARRAY<VECTOR<T,2> > points,const ARRAY<T> values,bool& selected,const T tolerance=1e-6,const ARRAY<T>* const phi_values=0,const T interface_value=0);
    static T Interpolation_From_Delaunay_Triangle_Vertices(const VECTOR<T,3> location,const ARRAY<VECTOR<T,3> > points,const ARRAY<T> values,bool& selected,const T tolerance=1e-6,const ARRAY<T>* const phi_values=0,const T interface_value=0);
//#####################################################################
};
}
#endif
