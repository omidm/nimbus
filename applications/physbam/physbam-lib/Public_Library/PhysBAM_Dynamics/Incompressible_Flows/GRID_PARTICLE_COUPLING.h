//#####################################################################
// Copyright 2011, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID_PARTICLE_COUPLING
//#####################################################################
#ifndef __GRID_PARTICLE_COUPLING__
#define __GRID_PARTICLE_COUPLING__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SURFACE_TENSION_FORCE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Dynamics/Particles/FLUID_PARTICLES.h>
namespace PhysBAM{

template<class TV>
class GRID_PARTICLE_COUPLING:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    typedef FACE_INDEX<TV::m> T_FACE_INDEX;
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;

public:
    GRID<TV>& grid;
    FLUID_PARTICLES<TV>& particles;
    ARRAY<T,TV_INT> cell_weights,pressure;
    ARRAY<bool,TV_INT> psi_D,surface_cells;
    ARRAY<bool,T_FACE_INDEX> psi_N,surface_faces;
    ARRAY<T,T_FACE_INDEX> face_velocities,face_forces;
    ARRAY<T,TV_INT> phi;
    T_LEVELSET levelset;
    SURFACE_TENSION_FORCE<TV>* surface_tension_force;

    T number_particle_per_cell,pic_ratio,viscosity,density;
    bool output_matrices,use_constant_interpolation,output_face_velocities,fix_volume_error,use_log_based_volume_correction,use_IC_preconditioner,feed_back_for_momentum_conservation,not_feed_back_to_surface_faces,use_simplified_mass,apply_pressure,recompute_mass;
    int max_cg_iterations;

    SPARSE_MATRIX_FLAT_MXN<T> HT,G,K,KT,J,W,H,WH,C,WT,Gu;
    SPARSE_MATRIX_FLAT_NXN<T> A;
    VECTOR_ND<T> Mi,x,b,q,s,r,k,z,v,dv;
    ARRAY<int,TV_INT> cell_indices;
    ARRAY<int,T_FACE_INDEX> face_indices;
    ARRAY<TV_INT> index_cells;
    ARRAY<T_FACE_INDEX> index_faces;
    VECTOR_ND<T> particle_mass,face_mass;
    T record_time;

public:
    GRID_PARTICLE_COUPLING(GRID<TV>& grid_input,FLUID_PARTICLES<TV>& particles_input);
    ~GRID_PARTICLE_COUPLING();
    void Initialize_Grid();
    void Solve(const T& dt);
    void Update_Position_Based_State();
    void Output_Matrices()const;
    void Set_Surface_Tension_Force(SURFACE_TENSION_FORCE<TV>* input)
    {surface_tension_force=input;}
    void Add_Explicit_Forces(const T time,const T dt);
    void Add_External_Forces(const T dt);
    void Compute_Grid_To_Surface_Mesh_Interpolation_Matrix();

    int Flatten_Particle_Index(int p,int axis)
    {return (axis-1)*particles.Size()+p;}

    void Project_Out_Null_Space_Of_Nullspace();
    void Compute_Cell_Weights();
    void Compute_Levelset();
    void Compute_Mass();
    void Setup_Boundary_Conditions();
    void Build_Index_Mapping();
    void Compute_Particle_To_Grid_Interpolation_Matrix_Transpose_Linear();
    void Compute_Particle_To_Grid_Interpolation_Matrix_Transpose_Constant();
    void Compute_Grid_Gradient_Matrix();
    void Compute_Grid_Diffusion_Matrix(const T dt);
    void Compute_Full_Coupling_Matrix();
    void Solve_Coupling_System(const T dt);
    void Interpolate_From_Grid_To_Particles(const T dt);
    void Add_Volume_Error_To_RHS(VECTOR_ND<T>& b,const T dt);
    SPARSE_MATRIX_FLAT_NXN<T> Add_Diagonal_Terms_To_Matrix(const T diagonal_term,const int start_row,const int end_row,const SPARSE_MATRIX_FLAT_NXN<T>& matrix);
    SPARSE_MATRIX_FLAT_MXN<T> Add_Diagonal_Terms_To_Matrix(const T diagonal_term,const int start_row,const int end_row,const SPARSE_MATRIX_FLAT_MXN<T>& matrix);
    SPARSE_MATRIX_FLAT_MXN<T> Build_Matrix_From_Blocks(SPARSE_MATRIX_FLAT_MXN<T>& block11,SPARSE_MATRIX_FLAT_MXN<T>& block12,SPARSE_MATRIX_FLAT_MXN<T>& block21,SPARSE_MATRIX_FLAT_MXN<T>& block22);
    T Linear_Kernel(const TV& p1,const TV& p2,const T& one_over_r)
    {return clamp_min(TV::Constant_Vector(1)-abs(p1-p2)*one_over_r,(T)0).Product();}
//#####################################################################
};
}
#endif
