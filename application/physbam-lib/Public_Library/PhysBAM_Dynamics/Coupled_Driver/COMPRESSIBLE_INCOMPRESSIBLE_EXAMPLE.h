//#####################################################################
// Copyright 2012, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE
//#####################################################################
#ifndef __COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE__
#define __COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE__

#include <PhysBAM_Dynamics/Boundaries/BOUNDARY_PHI_WATER.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_CALLBACKS.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_EVOLUTION_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Boundaries/GEOMETRY_BOUNDARY_POLICY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/LINEAR_INTERPOLATION_COLLIDABLE_CELL_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_HAMILTON_JACOBI_ENO.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEMI_LAGRANGIAN_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/EXAMPLE.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM {

template<class T_GRID> class LEVELSET_MULTIPLE_UNIFORM;

template<class TV>
class COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE:public EXAMPLE<TV>, public LEVELSET_CALLBACKS<GRID<TV> >
{
    typedef EXAMPLE<TV> BASE;typedef typename TV::SCALAR T;typedef GRID<TV> T_GRID;
    typedef VECTOR<T,TV::dimension+2> TV_DIMENSION;typedef VECTOR<int,TV::dimension> TV_INT;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY_SCALAR;
    typedef typename REBIND_LENGTH<T_BOUNDARY_SCALAR,TV::dimension+2>::TYPE T_BOUNDARY_DIMENSION_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_FACE_ARRAYS_DIMENSION_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;
    typedef typename LEVELSET_POLICY<GRID<TV> >::PARTICLE_LEVELSET T_PARTICLE_LEVELSET;
    typedef typename GEOMETRY_BOUNDARY_POLICY<GRID<TV> >::BOUNDARY_PHI_WATER T_BOUNDARY_PHI_WATER;
    typedef EXTRAPOLATION_UNIFORM<T_GRID,T> T_EXTRAPOLATION_SCALAR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
    typedef typename T_LINEAR_INTERPOLATION_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_LINEAR_INTERPOLATION_DIMENSION;
    typedef typename LEVELSET_POLICY<T_GRID>::EXTRAPOLATION_SCALAR::template REBIND<TV_DIMENSION>::TYPE T_EXTRAPOLATION_DIMENSION_SCALAR;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<TV> >::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    enum {d=TV::dimension+2};

  public:
    using BASE::Write_Frame_Title;

    using BASE::output_directory;using BASE::first_frame;using BASE::restart;using BASE::stream_type;using BASE::parse_args;using BASE::write_substeps_level;using BASE::frame_title;using BASE::write_last_frame;
    
    T_GRID grid;
    MPI_UNIFORM_GRID<T_GRID> *mpi_grid;
    T_ARRAYS_DIMENSION_SCALAR U,U_ghost;
    T_FACE_ARRAYS_SCALAR face_velocities,compressible_face_velocities;
    T_ARRAYS_SCALAR pressure,one_over_rho_c_squared,compressible_pressure,entropy;
    T_ARRAYS_BOOL psi;
    T_FACE_ARRAYS_BOOL psi_N,valid_mask;
    VECTOR<VECTOR<bool,2>,TV::dimension> domain_walls;

    T cfl,rho_incompressible,viscosity,surface_tension;
    int scale,eno_order,rk_order;
    int output_number,current_frame;
    int number_of_cells_to_extrapolate,pls_particles_per_cell;
    bool write_debug_data,apply_viscosity,use_rf;
    EOS<T>* eos;
    T_BOUNDARY_SCALAR boundary_scalar;
    T_BOUNDARY_SCALAR *boundary,*phi_boundary;
    T_BOUNDARY_PHI_WATER phi_boundary_water;
    T_BOUNDARY_DIMENSION_SCALAR *compressible_boundary;
    T_BOUNDARY_DIMENSION_SCALAR *compressible_boundary_scalar;
    VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,TV::dimension> cfl_eigensystem;
    CONSERVATION<T_GRID,d>* conservation_law_solver;
    PARTICLE_LEVELSET_EVOLUTION_UNIFORM<T_GRID> particle_levelset_evolution;
    EULER_LAPLACE<POISSON_COLLIDABLE_UNIFORM<T_GRID> > projection;
    ADVECTION_HAMILTON_JACOBI_ENO<T_GRID,T> eno_advection;
    ADVECTION_SEMI_LAGRANGIAN_UNIFORM<T_GRID,T> advection_scalar;
    T_GRID_BASED_COLLISION_GEOMETRY collision_bodies_affecting_fluid;

    COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE(const STREAM_TYPE& stream_type);
    virtual ~COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE();

    void Get_Levelset_Velocity(const T_GRID& grid,T_LEVELSET& levelset,T_FACE_ARRAYS_SCALAR& V_levelset,const T time) const PHYSBAM_OVERRIDE
    {V_levelset=face_velocities;}
    void Get_Levelset_Velocity(const T_GRID& grid,LEVELSET_MULTIPLE_UNIFORM<T_GRID>& levelset_multiple,T_FACE_ARRAYS_SCALAR& V_levelset,const T time) const PHYSBAM_OVERRIDE {}
    void Set_Incompressible_Density(const T rho_incompressible_input)
    {rho_incompressible=rho_incompressible_input;}
    void Set_Viscosity(const T viscosity_input)
    {apply_viscosity=true;viscosity=viscosity_input;}
    void Set_Surface_Tension(const T surface_tension_input)
    {surface_tension=surface_tension_input;}
    void Set_PLS_Particles_Per_Cell(const int pls_particles_per_cell_input=0)
    {pls_particles_per_cell=pls_particles_per_cell_input;}
//#####################################################################
    virtual void Initialize_Grids();
    virtual void Initialize_Fluid_State() = 0;
    virtual void Initialize_Boundaries() = 0;
    virtual void Set_Boundary_Conditions(const T dt,const T time) = 0;
    virtual void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE;
//#####################################################################
    void Write_Substep(const std::string& title,const int substep,const int level=0);
    void Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time);
    void Set_Eigensystems();
    void Log_Parameters() const PHYSBAM_OVERRIDE;
    void Register_Options() PHYSBAM_OVERRIDE;
    void Parse_Options() PHYSBAM_OVERRIDE;
    void Parse_Late_Options() PHYSBAM_OVERRIDE;
    void Read_Output_Files(const int frame);
    void Compute_Divergence(const T_FACE_LOOKUP& face_lookup);
    void Compute_One_Over_Rho_C_Squared();
    void Compute_Right_Hand_Side(const T dt);
    void Fill_Face_Weights_For_Projection(const T dt,const T time);
    void Compute_Pressure(const T dt,const T time);
    void Project(const T dt,const T time);
    void Compute_Density_Weighted_Face_Velocities(const T dt,const T time);
    void Apply_Pressure(T_ARRAYS_SCALAR& p_ghost,T_FACE_ARRAYS_SCALAR& p_face,const T dt,const T time);
    void Compute_Face_Pressure_From_Cell_Pressures(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& p_face,const T_ARRAYS_SCALAR& p_cell,const T dt,const T time);
    template<class T_TYPE> void Extrapolate_State_One_Sided(T_TYPE& Z,T_ARRAYS_SCALAR& phi);
    void Extrapolate_Velocity_Across_Interface(T_FACE_ARRAYS_SCALAR& face_velocities_input,T_ARRAYS_SCALAR& phi_ghost,const T band_width=3);
    void Update_Psi_And_Extrapolate_State_Into_Incompressible_Flow(const T dt,const T time);
    void Limit_Dt(T& dt,const T time);
//#####################################################################
};
}
#endif // __COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE__
