//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_SPRINGS
//#####################################################################
#ifndef __LINEAR_SPRINGS__
#define __LINEAR_SPRINGS__

#include <PhysBAM_Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV> class LINEAR_SPRINGS_SYSTEM_VECTOR;

template<class TV>
class LINEAR_SPRINGS:public DEFORMABLES_FORCES<TV>,public SPRINGS_TAG
{
    typedef typename TV::SCALAR T;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;using BASE::Invalidate_CFL;using BASE::cfl_number;
    using BASE::max_strain_per_time_step;using BASE::use_rest_state_for_strain_rate;using BASE::use_implicit_velocity_independent_forces;using BASE::compute_half_forces;
    typedef typename FORCE_ELEMENTS::ITERATOR SEGMENT_ITERATOR;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;
    typedef MATRIX<T,TV::dimension> T_MATRIX;

public:
    SEGMENT_MESH& segment_mesh;
    ARRAY<T> youngs_modulus; // units of force (i.e. force per unit strain)
    T constant_youngs_modulus; // units of force (i.e. force per unit strain)
    ARRAY<T> restlength,visual_restlength; // visual restlength corresponds to length between particles; restlength may be larger than this to avoid zero/small restlength
    ARRAY<T> damping; // units of force*time (i.e. force per unit strain rate)
    T constant_damping; // units of force*time (i.e. force per unit strain rate)
    bool use_plasticity;
    ARRAY<T> plastic_yield_strain;
    ARRAY<T> plastic_hardening;
    ARRAY<T> plastic_visual_restlength;
    T plasticity_clamp_ratio;
    bool cache_strain;
    mutable ARRAY<VECTOR<T,2> > strains_of_segment; // VECTOR<T,2>(strain_rate, strain)

    ARRAY<T> energy_correction_forces;
    ARRAY<int> energy_correction_sign;
    ARRAY<T> extra_energy;
    ARRAY<T> potential_energy_save;
    ARRAY<T> force_correction;
    ARRAY<TV> previously_applied_forces;
    bool use_kinetic_energy_fix;
    bool relaxation_fraction;
    bool use_gauss_seidel_in_energy_correction;
    bool allow_kd_direction_flip;
    bool verbose;

    ARRAY<T> residual_PE;
    ARRAY<T> delta_PE;
    T total_delta_PE;
    ARRAY<T> force_estimates;
    ARRAY<ARRAY<int> > incident_nodes;
    T saved_constant_youngs_modulus;
    ARRAY<TV> saved_V;

protected:
    struct STATE{
        STATE()
            :coefficient(0),sqrt_coefficient(0)
        {}

        VECTOR<int,2> nodes; // copy of nodes to reduce associativity needed in cache
        T coefficient;
        T sqrt_coefficient;
        TV direction;
    };
    ARRAY<STATE> states;
    ARRAY<T> current_lengths;

public:
    FORCE_ELEMENTS force_segments;
    LINEAR_SPRINGS(PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh_input,const bool implicit);

    virtual ~LINEAR_SPRINGS();

    virtual void Set_Stiffness(const T youngs_modulus_input)
    {constant_youngs_modulus=youngs_modulus_input;youngs_modulus.Clean_Memory();Invalidate_CFL();}

    virtual void Set_Stiffness(ARRAY_VIEW<const T> youngs_modulus_input)
    {constant_youngs_modulus=0;youngs_modulus=youngs_modulus_input;Invalidate_CFL();}

    virtual void Set_Restlength(ARRAY_VIEW<const T> restlength_input)
    {restlength=restlength_input;visual_restlength=restlength;Invalidate_CFL();}

    virtual void Set_Damping(const T damping_input)
    {constant_damping=damping_input;damping.Clean_Memory();Invalidate_CFL();}

    virtual void Set_Damping(ARRAY_VIEW<const T> damping_input)
    {constant_damping=0;damping=damping_input;Invalidate_CFL();}

    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE
    {} // Add_Force_Differential always enforces definiteness

//#####################################################################
    template<class T_FIELD> void Enable_Plasticity(const T_FIELD& plastic_yield_strain_input,const T_FIELD& plastic_hardening_input,const T plasticity_clamp_ratio_input=4);
    virtual void Set_Restlength_From_Particles();
    virtual void Set_Restlength_From_Material_Coordinates(ARRAY_VIEW<TV> material_coordinates);
    virtual void Clamp_Restlength(const T clamped_restlength);
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    virtual void Set_Stiffness_Based_On_Reduced_Mass(const T scaling_coefficient); // assumes mass and restlength are already defined
    virtual void Set_Stiffness_Based_On_Reduced_Mass(ARRAY_VIEW<const T> scaling_coefficient); // assumes mass and restlength are already defined
    virtual void Set_Overdamping_Fraction(const T overdamping_fraction=1); // 1 is critically damped
    virtual void Set_Overdamping_Fraction(ARRAY_VIEW<const T> overdamping_fraction); // 1 is critically damped
    void Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction=1); // 1 is critically damped
    void Clamp_Restlength_With_Fraction_Of_Springs(const T fraction=.01);
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    int Velocity_Dependent_Forces_Size() const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const PHYSBAM_OVERRIDE;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
    T Average_Restlength() const;
    void Print_Restlength_Statistics() const;
    void Print_Deformation_Statistics() const;
    T Maximum_Compression_Or_Expansion_Fraction(int* index=0) const;
    T Potential_Energy(const int s,const T time) const;
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
    T Residual_Energy(const T time) const PHYSBAM_OVERRIDE;
    void Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name="") const PHYSBAM_OVERRIDE;

    void Prepare_Energy_Correction_Force();
    void Compute_Energy_Correction_Force(ARRAY_VIEW<const TV> velocity_save,const int max_particle_degree,const T time,const T dt);
    void Apply_Energy_Correction(const T time,const T dt);
    void Add_Connectivity(ARRAY<int>& particle_degree);
    TV Endpoint_Velocity(int s,int b) const;
    T Endpoint_Kinetic_Energy(int s,int b) const;
    T Endpoint_Kinetic_Energy(int s) const;
    typename TV::SCALAR Effective_Impulse_Factor(int s) const;
    T Compute_Spring_Work(int s,ARRAY_VIEW<const TV> velocity_save,const T time,const T dt) const;
    void Compute_Previously_Applied_Forces();
    void Setup_Set_Velocity_From_Positions(const T time,const bool is_position_update,const bool reset_alphas) PHYSBAM_OVERRIDE;
    int Get_Element_Count() PHYSBAM_OVERRIDE;
    FORCE_ELEMENTS* Get_Force_Elements() PHYSBAM_OVERRIDE;
    ARRAY<int>* Incident_Nodes(const int force_element) PHYSBAM_OVERRIDE;
    TV Get_Direction(const int force_element) PHYSBAM_OVERRIDE;
    T Get_Combined_One_Over_Mass(const int force_element) PHYSBAM_OVERRIDE;
    ARRAY<int>* Incident_Force_Elements(const int particle) PHYSBAM_OVERRIDE;
    void Choose_Solution(const bool use_orig_force,const int force_element,const T dt,const T alpha1,const T alpha2,ARRAY<T>& v_n_hats) PHYSBAM_OVERRIDE;
    TV Get_Force(const int force_element,const int particle,const bool use_original_force) PHYSBAM_OVERRIDE;
    T Get_Force(const int force_element) PHYSBAM_OVERRIDE;
    void Set_Force(const int force_element,const T force) PHYSBAM_OVERRIDE;
    void Get_Damping_Force(const int particle,TV& damping_force,const T dt,const bool use_coefficient) PHYSBAM_OVERRIDE;
    void Update_Residual_Energy(const int force_element,const T residual_energy,const T time) PHYSBAM_OVERRIDE;
    T Get_Residual_Energy(const int force_element) PHYSBAM_OVERRIDE;
    void Compute_Quadratic_Contribution_For_Force(T& A,T& a,T&c,const T dt,const int force_element,const T combined_one_over_mass,const bool ignore_PE_terms) PHYSBAM_OVERRIDE;
    void Compute_Quadratic_Contribution_For_Node(T& B,T& C,T&b,const T dt,const int node,const int force_element,const T combined_one_over_mass,const T v_n_correction,
        const bool ignore_PE_terms) PHYSBAM_OVERRIDE;
    void Compute_Quadratic_Contribution_For_Residual(T& B,T& C,T& a,T&b,T&c,const T dt,const T time,const int force_element,const bool ignore_PE_terms) PHYSBAM_OVERRIDE;
    void Store_Delta_PE(const T time) PHYSBAM_OVERRIDE;
    T Get_Total_Delta_PE() PHYSBAM_OVERRIDE;
    void Save_And_Reset_Elastic_Coefficient() PHYSBAM_OVERRIDE;
    void Restore_Elastic_Coefficient() PHYSBAM_OVERRIDE;
    void Store_Velocities();
    void Apply_Energy_Correction_Impulse(const int s,const T force,const T dt);
    void Save_Potential_Energy(const T time) PHYSBAM_OVERRIDE;
    void Compute_Energy_Error(ARRAY_VIEW<const TV> velocity_save,const T time,const T dt);
//#####################################################################
};

template<class TV> LINEAR_SPRINGS<TV>*
Create_Edge_Springs(PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh,const typename TV::SCALAR stiffness=2e3,
    const typename TV::SCALAR overdamping_fraction=1,const bool limit_time_step_by_strain_rate=true,const typename TV::SCALAR max_strain_per_time_step=.1,
    const bool use_rest_state_for_strain_rate=true,const typename TV::SCALAR restlength_enlargement_fraction=0,const bool verbose=true,const bool implicit=false);

template<class T_OBJECT> LINEAR_SPRINGS<typename T_OBJECT::VECTOR_T>*
Create_Edge_Springs(T_OBJECT& object,
    const typename T_OBJECT::SCALAR stiffness=2e3,const typename T_OBJECT::SCALAR overdamping_fraction=1,const bool limit_time_step_by_strain_rate=true,
    const typename T_OBJECT::SCALAR max_strain_per_time_step=.1,const bool use_rest_state_for_strain_rate=true,const typename T_OBJECT::SCALAR restlength_enlargement_fraction=0,
    const bool verbose=true,const bool implicit=false);

}
#endif
