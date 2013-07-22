#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTION_RLE
//#####################################################################
#ifndef __PROJECTION_RLE__
#define __PROJECTION_RLE__

#include <PhysBAM_Tools/Grids_PDE_Linear/PROJECTION.h>
#include <PhysBAM_Geometry/Grids_RLE_PDE_Linear/LAPLACE_COLLIDABLE_RLE.h>
namespace PhysBAM{

template<class T_GRID>
class PROJECTION_RLE:public PROJECTION<typename T_GRID::SCALAR>
{
public:
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_Y_ITERATOR FACE_Y_ITERATOR;

    using PROJECTION<T>::use_non_zero_divergence;

    const T_GRID& grid;
    ARRAY<T> p;
    ARRAY<T>& V;
    ARRAY<T>& face_velocities; // for templatization purposes only
    LAPLACE_COLLIDABLE_RLE<T_GRID> laplace;
    ARRAY<T> divergence; // use this to set up a non-zero divergence

    PROJECTION_RLE(const T_GRID& grid_input,ARRAY<T>& V_input);
    virtual ~PROJECTION_RLE();

    virtual void Initialize_Grid()
    {laplace.Initialize_Grid();
    Use_Non_Zero_Divergence(use_non_zero_divergence); // call this since the grid changed
    p.Resize(grid.number_of_cells);}

    void Use_Non_Zero_Divergence(const bool use_non_zero_divergence_input=true)
    {use_non_zero_divergence=use_non_zero_divergence_input;
    if(use_non_zero_divergence) divergence.Resize(grid.number_of_cells);else divergence.Resize(0);}

//#####################################################################
    virtual void Make_Divergence_Free(const T time,const ARRAY<T>* phi_ghost);
    virtual void Enforce_Velocity_Compatibility();
protected:
    struct Compute_Horizontal_Divergence{template<class T_FACE> static void Apply(PROJECTION_RLE<T_GRID>& projection);};
    void Compute_Vertical_Divergence();
    struct Apply_Horizontal_Pressure_Gradients{template<class T_FACE> static void Apply(const PROJECTION_RLE<T_GRID>& projection,const ARRAY<T>* phi_ghost);};
    void Apply_Vertical_Pressure_Gradients(const ARRAY<T>* phi_ghost);
    struct Compute_Boundary_Areas{template<class T_FACE> static void Apply(const PROJECTION_RLE<T_GRID>& projection,ARRAY<T>& boundary_area);};
    struct Zero_Out_Compatibility_Errors{template<class T_FACE> static void Apply(PROJECTION_RLE<T_GRID>& projection,const ARRAY<T>& compatibility_fraction);};
    struct Zero_Out_Neumann_Regions{template<class T_FACE> static void Apply(PROJECTION_RLE<T_GRID>& projection);};
public:
    void Transfer_Pressure(const T_GRID& new_grid);
//#####################################################################
};
}
#endif
#endif
