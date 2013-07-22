//#####################################################################
// Copyright 2007, Jon Gretarsson, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_PROJECTION_UNIFORM  
//#####################################################################
#ifndef __EULER_PROJECTION_UNIFORM__
#define __EULER_PROJECTION_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_HAMILTON_JACOBI_ENO.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_REFLECTION_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/BOUNDARY_OBJECT_REFLECTION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_PROJECTION.h>
namespace PhysBAM{

template<class T_GRID> class EULER_UNIFORM;
template<class T_GRID> class EULER_LAPLACE;
template<class T_GRID> class POISSON_COLLIDABLE_UNIFORM;
template<class TV> class INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS;

template<class T_GRID>
class EULER_PROJECTION_UNIFORM:public EULER_PROJECTION<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_INT TV_INT;typedef VECTOR<T,T_GRID::dimension+2> TV_DIMENSION;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_FACE_ARRAYS_DIMENSION_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef GRID<VECTOR<T,TV::dimension-1> > T_GRID_LOWER_DIM;typedef typename T_GRID_LOWER_DIM::CELL_ITERATOR CELL_ITERATOR_LOWER_DIM;
    typedef typename T_GRID_LOWER_DIM::VECTOR_INT TV_INT_LOWER_DIM;
    typedef typename GRID<VECTOR<T,1> >::CELL_ITERATOR CELL_ITERATOR_1D;
    
public:

    T_ARRAYS_SCALAR p; // p should always store the real pressure.
    T_ARRAYS_SCALAR p_save_for_projection;
    T_ARRAYS_SCALAR p_advected;
    T_ARRAYS_SCALAR density_scaling;
    T_ARRAYS_SCALAR one_over_rho_c_squared;
    bool is_pressure_scaled;
    T dt_scale_pressure;
    T_ARRAYS_SCALAR p_dirichlet;
    T_FACE_ARRAYS_SCALAR face_velocities,face_velocities_save;
    T_FACE_ARRAYS_DIMENSION_SCALAR *fluxes;
    EULER_UNIFORM<T_GRID>* euler;
    EULER_LAPLACE<POISSON_COLLIDABLE_UNIFORM<T_GRID> >* elliptic_solver;
    POISSON_COLLIDABLE_UNIFORM<T_GRID>* poisson;
    T_ARRAYS_SCALAR divergence; // use this to set up a non-zero divergence
    BOUNDARY_UNIFORM<T_GRID,T>* pressure_boundary;
    BOUNDARY_REFLECTION_UNIFORM<T_GRID,T> pressure_boundary_default;
    bool save_fluxes,use_exact_neumann_face_location,use_neumann_condition_for_outflow_boundaries;
#if 1
    ADVECTION_HAMILTON_JACOBI_ENO<GRID<VECTOR<T,1> >,T> pressure_advection_HJ;
#else
    ADVECTION_HAMILTON_JACOBI_ENO<T_GRID,T> pressure_advection_HJ;
#endif
    BOUNDARY_OBJECT_REFLECTION<T_GRID,T> pressure_object_boundary;
    int hj_eno_order;

private:
    INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS<TV>* incompressible_coupling_callbacks;
    bool transition_to_using_implicit_pressure;
    
public:

    EULER_PROJECTION_UNIFORM(EULER_UNIFORM<T_GRID>* euler_input);
    ~EULER_PROJECTION_UNIFORM();

    void Save_State(T_FACE_ARRAYS_SCALAR& face_velocities_s)
    {T_FACE_ARRAYS_SCALAR::Copy(face_velocities,face_velocities_s);}

    void Restore_State(T_FACE_ARRAYS_SCALAR& face_velocities_s)
    {T_FACE_ARRAYS_SCALAR::Copy(face_velocities_s,face_velocities);}

    void Scale_Pressure_By_Dt(const T dt)
    {assert(!is_pressure_scaled);p*=dt;dt_scale_pressure=dt;is_pressure_scaled=true;}

    void Unscale_Pressure_By_Dt(const T dt)
    {assert(is_pressure_scaled && (dt==dt_scale_pressure));p*=(1/dt);is_pressure_scaled=false;}

    void Set_Transition_To_Using_Implicit_Pressure(const bool transition_to_using_implicit_pressure_input)
    {transition_to_using_implicit_pressure=transition_to_using_implicit_pressure_input;}

    void Exchange_Pressures_For_Projection()
    {T_ARRAYS_SCALAR::Exchange_Arrays(p,p_save_for_projection);}

    void Set_Incompressible_Coupling_Callbacks(INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS<TV>* incompressible_coupling_callbacks_input)
    {incompressible_coupling_callbacks=incompressible_coupling_callbacks_input;}

    void Set_Use_Exact_Neumann_Face_Location(const bool use_exact_neumann_face_location_input)
    {use_exact_neumann_face_location=use_exact_neumann_face_location_input;
    euler->conservation->Set_Use_Exact_Neumann_Face_Location(use_exact_neumann_face_location);}

    void Set_Constant_Extrapolation_Pressure_Boundary()
    {pressure_boundary->Set_Constant_Extrapolation(euler->boundary->constant_extrapolation);}

    void Set_Custom_Pressure_Boundary(BOUNDARY_UNIFORM<T_GRID,T> &pressure_boundary_input)
    {pressure_boundary=&pressure_boundary_input;Set_Constant_Extrapolation_Pressure_Boundary();}

//#####################################################################
private:
    void Compute_Pressure(const T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    template<class FACE_LOOKUP> void Compute_Divergence(const FACE_LOOKUP &face_lookup);
public:
    virtual void Initialize_Grid();
    void Get_Pressure(T_ARRAYS_SCALAR& pressure) const;
    void Fill_Face_Weights_For_Projection(const T dt,const T time,T_FACE_ARRAYS_SCALAR& beta_face);
    void Get_Ghost_Density(const T dt,const T time,const int number_of_ghost_cells,T_ARRAYS_SCALAR& density_ghost) const;
    void Get_Ghost_Centered_Velocity(const T dt,const T time,const int number_of_ghost_cells,T_ARRAYS_VECTOR& centered_velocity_ghost) const;
    void Make_Boundary_Faces_Neumann(T_FACE_ARRAYS_BOOL& psi_N);
    void Project(const T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Get_Dirichlet_Boundary_Conditions(const T_ARRAYS_DIMENSION_SCALAR& U_dirichlet);
    void Set_Dirichlet_Boundary_Conditions(const T time);
    void Compute_Advected_Pressure(const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const T_FACE_ARRAYS_SCALAR face_velocities_for_solid_faces,const T dt);
    void Compute_Right_Hand_Side(const T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Compute_One_Over_rho_c_Squared();
    void Compute_Density_Weighted_Face_Velocities(const T dt,const T time,const T_FACE_ARRAYS_BOOL& psi_N);
    static void Compute_Density_Weighted_Face_Velocities(const T_GRID& face_grid,T_FACE_ARRAYS_SCALAR& face_velocities,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const T_ARRAYS_BOOL& psi,const T_FACE_ARRAYS_BOOL& psi_N);
    static void Compute_Face_Pressure_From_Cell_Pressures(const T_GRID& face_grid,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const T_ARRAYS_BOOL& psi,T_FACE_ARRAYS_SCALAR& p_face,const T_ARRAYS_SCALAR& p_cell);
    void Get_Ghost_Pressures(const T dt,const T time,const T_ARRAYS_BOOL& psi_D,const T_FACE_ARRAYS_BOOL& psi_N,
        const T_ARRAYS_SCALAR& pressure,T_ARRAYS_SCALAR& p_ghost);
    void Get_Pressure_At_Faces(const T dt,const T time,const T_ARRAYS_SCALAR& p_ghost,T_FACE_ARRAYS_SCALAR& p_face);
    void Apply_Pressure(const T_ARRAYS_SCALAR& p_ghost,const T_FACE_ARRAYS_SCALAR& p_face,const T_FACE_ARRAYS_SCALAR& face_velocities_star,
        const T_ARRAYS_BOOL& psi_D,const T_FACE_ARRAYS_BOOL& psi_N,const T dt,const T time);
    static void Apply_Pressure(const T_ARRAYS_SCALAR& p_ghost,const T_FACE_ARRAYS_SCALAR& p_face,const T_FACE_ARRAYS_SCALAR& face_velocities_star,
            const T_ARRAYS_BOOL& psi_D,const T_FACE_ARRAYS_BOOL& psi_N,const T dt,const T time,const T_ARRAYS_SCALAR& density_scaling,
            T_FACE_ARRAYS_DIMENSION_SCALAR *fluxes,EULER_UNIFORM<T_GRID>* euler);
    bool Consistent_Boundary_Conditions() const;
    void Log_Parameters() const;
//#####################################################################
};
}
#endif
