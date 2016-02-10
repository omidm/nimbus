//#####################################################################
// Copyright 2010, Mridul Aanjaneya, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_CAVITATION_UNIFORM  
//#####################################################################
#ifndef __EULER_CAVITATION_UNIFORM__
#define __EULER_CAVITATION_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_POLICY_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID> class EULER_UNIFORM;
template<class T_GRID> class LAPLACE_COLLIDABLE_UNIFORM;

template<class TV>
class EULER_CAVITATION_UNIFORM
{
    typedef GRID<TV> T_GRID;
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_FACE_ARRAYS_DIMENSION_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef GRID<VECTOR<T,TV::dimension-1> > T_GRID_LOWER_DIM;typedef typename T_GRID_LOWER_DIM::CELL_ITERATOR CELL_ITERATOR_LOWER_DIM;
    typedef typename T_GRID_LOWER_DIM::VECTOR_INT TV_INT_LOWER_DIM;
    typedef typename GRID<VECTOR<T,1> >::CELL_ITERATOR CELL_ITERATOR_1D;
    
private:
    EULER_UNIFORM<T_GRID>& euler;
    T_ARRAYS_SCALAR clamped_momentum_divergence;
    T_ARRAYS_SCALAR clamped_internal_energy_divergence;
    T epsilon;
    bool clamp_density;        // true if for density, false otherwise

public:
    T_ARRAYS_SCALAR p_cavitation;
    LAPLACE_COLLIDABLE_UNIFORM<T_GRID>* elliptic_solver;

    EULER_CAVITATION_UNIFORM(EULER_UNIFORM<T_GRID>& euler_input, const bool clamp_density_input, const T epsilon_input);
    ~EULER_CAVITATION_UNIFORM();

    void Apply_Cavitation_Correction(const T dt,const T time, T_FACE_ARRAYS_SCALAR& face_velocities);
    bool Is_Density_Clamped();

private:
    void Compute_Face_Pressure_From_Cell_Pressures(const T_GRID& face_grid,T_FACE_ARRAYS_SCALAR& p_face,const T_ARRAYS_SCALAR& p_cell);
    void Compute_Pressure(const T dt,const T time);
    void Compute_Clamped_Momentum_Divergence(const T dt);
    void Compute_Clamped_Internal_Energy_Divergence(const T dt);
    void Apply_Pressure(const T dt,const T time, T_FACE_ARRAYS_SCALAR& face_velocities);
    void Apply_Pressure_To_Density(const T dt);
    void Apply_Pressure_To_Internal_Energy(const T dt);
    void Initialize_Grid();
    void Fill_Ghost_Pressures_Along_Neumann_Boundaries();
    void Compute_Right_Hand_Side(const T dt);
    void Log_Parameters() const;
//#####################################################################
};
}
#endif
