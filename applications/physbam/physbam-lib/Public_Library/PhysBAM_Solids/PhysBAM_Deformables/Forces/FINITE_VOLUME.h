//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Neil Molino, Igor Neverov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FINITE_VOLUME
//#####################################################################
#ifndef __FINITE_VOLUME__
#define __FINITE_VOLUME__

#include <PhysBAM_Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <PhysBAM_Tools/Math_Tools/FACTORIAL.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/CONSTITUTIVE_MODELS_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV,int d> class DIAGONALIZED_SEMI_IMPLICIT_ELEMENT;

template<class TV,int d>
class FINITE_VOLUME:public DEFORMABLES_FORCES<TV>,public FINITE_VOLUME_TAG
{
    typedef typename TV::SCALAR T;
    typedef MATRIX<T,TV::m,d> T_MATRIX;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d>::OBJECT T_OBJECT;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;using BASE::max_strain_per_time_step;using BASE::use_velocity_independent_forces;using BASE::use_velocity_dependent_forces;
    using BASE::cfl_number;using BASE::compute_half_forces;
    typedef typename FORCE_ELEMENTS::ITERATOR FORCE_ITERATOR;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

    STRAIN_MEASURE<TV,d>& strain_measure;
    CONSTITUTIVE_MODEL<T,d>& constitutive_model;
    PLASTICITY_MODEL<T,d>* plasticity_model;
    ARRAY<T> Be_scales;
    ARRAY<T> sqrt_Be_scales;
    ARRAY<T_MATRIX> U;
    ARRAY<MATRIX<T,d> > De_inverse_hat;
    ARRAY<DIAGONAL_MATRIX<T,d> > Fe_hat;
    ARRAY<DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,d> >* dPi_dFe;
    ARRAY<DIAGONALIZED_STRESS_DERIVATIVE<T,d> >* dP_dFe;
    ARRAY<T>* Be_scales_save;
    ARRAY<MATRIX<T,d> >* V;
    T twice_max_strain_per_time_step; // for asynchronous
    ARRAY<DIAGONALIZED_SEMI_IMPLICIT_ELEMENT<TV,d> >* semi_implicit_data;
    bool use_uniform_density;
    ARRAY<T>* density_list;
    T density;
protected:
    ISOTROPIC_CONSTITUTIVE_MODEL<T,d>* isotropic_model;
    ANISOTROPIC_CONSTITUTIVE_MODEL<T,d>* anisotropic_model;
    IMPLICIT_OBJECT<TV>* implicit_surface;
    ARRAY<SYMMETRIC_MATRIX<T,TV::m> >* node_stiffness;
    ARRAY<MATRIX<T,TV::m> >* edge_stiffness;
public:
    FORCE_ELEMENTS force_elements;
    FORCE_ELEMENTS* force_segments;
    FORCE_ELEMENTS force_particles;
protected:
    bool destroy_data;
public:
    int half_force_size;
    int total_half_force_size;

    FINITE_VOLUME(const bool use_uniform_density_input,STRAIN_MEASURE<TV,d>& strain_measure,CONSTITUTIVE_MODEL<T,d>& constitutive_model,PLASTICITY_MODEL<T,d>* plasticity_model=0);
    virtual ~FINITE_VOLUME();

    void Define_Inversion_Via_Implicit_Surface(IMPLICIT_OBJECT<TV>& implicit_surface_input)
    {if(TV::m==d) PHYSBAM_FATAL_ERROR();implicit_surface=&implicit_surface_input;}

//#####################################################################
    static FINITE_VOLUME* Create(T_OBJECT& object,CONSTITUTIVE_MODEL<T,d>* constitutive_model,PLASTICITY_MODEL<T,d>* plasticity_model,const bool verbose,const bool use_uniform_density);
    // initialization
    void Save_V();
    void Save_Stress_Derivative();
    void Use_Quasistatics();
    void Use_Stiffness_Matrix();
    void Update_Be_Scales();
    void Initialize_Be_Scales_Save();
    void Copy_Be_Scales_Save_Into_Be_Scales(const ARRAY<int>& map);
    // solids_forces functions
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    int Velocity_Dependent_Forces_Size() const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE;
    void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const PHYSBAM_OVERRIDE;
    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
    // asynchronous
    void Semi_Implicit_Impulse_Precomputation(const T time,const T max_dt,ARRAY<T>* time_plus_dt,const bool verbose);
    void Semi_Implicit_Recompute_Dt(const int element,T& time_plus_dt);
    void Add_Semi_Implicit_Impulse(const int element,const T dt,T* time_plus_dt);
    void Add_Semi_Implicit_Impulses(const T dt,ARRAY<T>* time_plus_dt);
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
    // unused
    void Store_Velocities() PHYSBAM_OVERRIDE {};
//#####################################################################
};

template<class T_OBJECT,class T,int d> FINITE_VOLUME<typename T_OBJECT::VECTOR_T,d>*
Create_Finite_Volume(T_OBJECT& object,CONSTITUTIVE_MODEL<T,d>* constitutive_model,const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,
    const bool use_rest_state_for_strain_rate=true,const bool precompute_stiffness_matrix=false,const bool verbose=true,PLASTICITY_MODEL<T,d>* plasticity_model=0,
    const bool use_uniform_density=false)
{
    typedef typename T_OBJECT::VECTOR_T TV;
    FINITE_VOLUME<TV,d>* fvm=FINITE_VOLUME<TV,d>::Create(object,constitutive_model,plasticity_model,verbose,use_uniform_density);
    fvm->Limit_Time_Step_By_Strain_Rate(limit_time_step_by_strain_rate,max_strain_per_time_step);
    fvm->Use_Rest_State_For_Strain_Rate(use_rest_state_for_strain_rate);
    if(precompute_stiffness_matrix) fvm->Use_Stiffness_Matrix();
    return fvm;
}

template<class T_OBJECT,class T,int d> FINITE_VOLUME<typename T_OBJECT::VECTOR_T,d>*
Create_Quasistatic_Finite_Volume(T_OBJECT& object,CONSTITUTIVE_MODEL<T,d>* constitutive_model,const bool precompute_stiffness_matrix=false,const bool verbose=true,
    const bool use_uniform_density=false)
{
    typedef typename T_OBJECT::VECTOR_T TV;
    FINITE_VOLUME<TV,d>* fvm=Create_Finite_Volume(object,constitutive_model,true,(T).1,true,precompute_stiffness_matrix,verbose,(PLASTICITY_MODEL<T,d>*)0,use_uniform_density);
    fvm->Use_Quasistatics();if(precompute_stiffness_matrix) fvm->Use_Stiffness_Matrix();return fvm;
}

template<class T_OBJECT,class T,int d> FINITE_VOLUME<typename T_OBJECT::VECTOR_T,d>*
Create_Backward_Euler_Finite_Volume(T_OBJECT& object,CONSTITUTIVE_MODEL<T,d>* constitutive_model,const bool limit_time_step_by_strain_rate=true,
    const T max_strain_per_time_step=.1,const bool use_rest_state_for_strain_rate=true,const bool precompute_stiffness_matrix=false,const bool verbose=true,
    const bool use_uniform_density=false)
{
    typedef typename T_OBJECT::VECTOR_T TV;
    FINITE_VOLUME<TV,d>* fvm=Create_Finite_Volume(object,constitutive_model,limit_time_step_by_strain_rate,max_strain_per_time_step,
        use_rest_state_for_strain_rate,precompute_stiffness_matrix,verbose,(PLASTICITY_MODEL<T,d>*)0,use_uniform_density);
    fvm->Use_Quasistatics();return fvm;
}

}
#endif
