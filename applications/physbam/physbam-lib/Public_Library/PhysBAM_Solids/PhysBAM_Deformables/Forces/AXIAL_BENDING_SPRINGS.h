//#####################################################################
// Copyright 2004-2008, Eran Guendelman, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class AXIAL_BENDING_SPRINGS
//#####################################################################
// A spring force between the edge connecting two triangles and the cross-edge connecting their two opposite vertices.
// The spring corresponds to the shortest line connecting the two edge lines (connections can be anywhere on the infinite lines
// supporting the edges, not restricted to being on the edges themselves.
//#####################################################################
#ifndef __AXIAL_BENDING_SPRINGS__
#define __AXIAL_BENDING_SPRINGS__

#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

class TRIANGLE_MESH;

template<class T_input>
class AXIAL_BENDING_SPRINGS:public DEFORMABLES_FORCES<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;using BASE::max_strain_per_time_step;using BASE::use_rest_state_for_strain_rate;using BASE::Invalidate_CFL;using BASE::cfl_number;
    typedef typename FORCE_ELEMENTS::ITERATOR SPRING_ITERATOR;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

    TRIANGLE_MESH& triangle_mesh;
    ARRAY<T> youngs_modulus; // units of force (i.e. force per unit strain)
    ARRAY<T> restlength,visual_restlength; // visual restlength corresponds to length between particles; restlength may be larger than this to avoid zero/small restlength
    ARRAY<T> damping; // units of force*time (i.e. force per unit strain rate)
    ARRAY<VECTOR<int,4> > spring_particles; // spring is shortest line between segment with particles (1,2) and segment with particles (3,4)
    ARRAY<T> attached_edge_restlength; // total rest length of edges of attached triangles which are not the center edge
    ARRAY<T> attached_edge_length; // total length of edges of attached triangles which are not the center edge
private:
    ARRAY<T> optimization_coefficient;
    ARRAY<T> optimization_current_length;
    ARRAY<TV> optimization_direction;
    ARRAY<VECTOR<T,4> > optimization_weights;
    FORCE_ELEMENTS force_springs;
public:
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
    ARRAY<ARRAY<int> > incident_elements;
    ARRAY<T> saved_youngs_modulus;

    AXIAL_BENDING_SPRINGS(PARTICLES<TV>& particles_input,TRIANGLE_MESH& triangle_mesh_input);

    virtual ~AXIAL_BENDING_SPRINGS();

    void Set_Stiffness(const T youngs_modulus_input)
    {ARRAYS_COMPUTATIONS::Fill(youngs_modulus,youngs_modulus_input);Invalidate_CFL();}

    void Set_Stiffness(ARRAY_VIEW<const T> youngs_modulus_input)
    {ARRAY<T>::Copy(youngs_modulus_input,youngs_modulus);Invalidate_CFL();}

    void Clamp_Restlength(const T clamped_restlength)
    {for(int i=1;i<=restlength.m;i++) restlength(i)=max(visual_restlength(i),clamped_restlength);}

    void Set_Damping(const T damping_input)
    {ARRAYS_COMPUTATIONS::Fill(damping,damping_input);Invalidate_CFL();}

    void Set_Damping(ARRAY_VIEW<const T> damping_input)
    {ARRAY<T>::Copy(damping_input,damping);Invalidate_CFL();}

//#####################################################################
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    void Initialize();
    void Set_Restlength_From_Particles();
    void Set_Overdamping_Fraction(const T overdamping_fraction=1); // 1 is critically damped
    void Set_Overdamping_Fraction(ARRAY_VIEW<const T> overdamping_fraction); // 1 is critically damped
    void Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction=1); // 1 is critically damped
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const;
    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
    void Set_Stiffness_Based_On_Reduced_Mass(const T scaling_coefficient); // assumes mass and restlength are already defined
    void Axial_Vector(const VECTOR<int,4>& nodes,T& axial_length,TV& axial_direction,VECTOR<T,2>& weights,T& attached_edge_length) const;

    T Potential_Energy(int s,const T time) const;
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
    void Prepare_Energy_Correction_Force();
    void Compute_Energy_Correction_Force(ARRAY_VIEW<const TV> velocity_save,const int max_particle_degree,const T time,const T dt);
    void Apply_Energy_Correction(const T time,const T dt);
    void Add_Connectivity(ARRAY<int>& particle_degree);
    T Endpoint_Mass(int s,int b) const;
    TV Endpoint_Velocity(int s,int b) const;
    TV Endpoint_Velocity(ARRAY_VIEW<const TV> velocity,int s,int b) const;
    T Endpoint_Kinetic_Energy(int s,int b) const;
    T Endpoint_Kinetic_Energy(ARRAY_VIEW<const TV> velocity,int s,int b) const;
    T Endpoint_Kinetic_Energy(int s) const;
    T Effective_Impulse_Factor(int s) const;
    T Compute_Spring_Work(int s,ARRAY_VIEW<const TV> velocity_save,const T time,const T dt) const;
    void Compute_Previously_Applied_Forces();
    void Apply_Energy_Correction_Impulse(const int s,const T force,const T dt);
    void Save_Potential_Energy(const T time);
    void Compute_Energy_Error(ARRAY_VIEW<const TV> velocity_save,const T time,const T dt);

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
//#####################################################################
};

template<class T> AXIAL_BENDING_SPRINGS<T>*
Create_Axial_Bending_Springs(PARTICLES<VECTOR<T,3> >& particles,TRIANGLE_MESH& triangle_mesh,const T clamped_restlength=(T).01,
    const T stiffness=2/(1+sqrt((T)2)),const T overdamping_fraction=2,const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,
    const bool use_rest_state_for_strain_rate=true,const bool verbose=true);

template<class T> AXIAL_BENDING_SPRINGS<T>*
Create_Axial_Bending_Springs(TRIANGULATED_SURFACE<T>& triangulated_surface,
    const T clamped_restlength=(T).01,const T stiffness=2/(1+sqrt((T)2)),const T overdamping_fraction=2,const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,
    const bool use_rest_state_for_strain_rate=true,const bool verbose=true);
}
#endif
