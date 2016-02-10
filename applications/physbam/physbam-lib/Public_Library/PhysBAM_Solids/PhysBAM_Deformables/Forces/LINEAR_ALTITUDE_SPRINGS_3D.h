//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Igor Neverov, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_ALTITUDE_SPRINGS_3D
//#####################################################################
#ifndef __LINEAR_ALTITUDE_SPRINGS_3D__
#define __LINEAR_ALTITUDE_SPRINGS_3D__

#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS.h>
namespace PhysBAM{

template<class T_input>
class LINEAR_ALTITUDE_SPRINGS_3D:public LINEAR_ALTITUDE_SPRINGS<VECTOR<T_input,3>,3>
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;typedef typename TV::SPIN T_SPIN;
    typedef LINEAR_ALTITUDE_SPRINGS<VECTOR<T,3>,3> BASE;typedef typename BASE::SPRING_STATE SPRING_STATE;typedef typename BASE::SPRING_PARAMETER SPRING_PARAMETER;
    using BASE::use_springs_compressed_beyond_threshold;using BASE::spring_compression_fraction_threshold;using BASE::print_number_used;
    using BASE::force_elements;
public:
    using BASE::particles;using BASE::Invalidate_CFL;
    using BASE::max_strain_per_time_step;using BASE::mesh;using BASE::use_rest_state_for_strain_rate;
    using BASE::spring_states;using BASE::spring_states_all_springs;
    using BASE::use_plasticity;using BASE::use_shortest_spring_only;
    using BASE::cache_strain;using BASE::strains_of_spring;using BASE::strains_of_spring_all_springs;
    using BASE::compute_half_forces;
    typedef typename BASE::ELEMENT_ITERATOR ELEMENT_ITERATOR;
    using BASE::parameters;

    ARRAY<T> potential_energy_save;
    ARRAY<T> delta_PE;
    T total_delta_PE;
    ARRAY<T> total_PE;
    ARRAY<T> residual_PE;
    ARRAY<T> force_estimates;
    ARRAY<ARRAY<int> > incident_nodes;
    ARRAY<VECTOR<T,4> > saved_youngs_modulus;

    LINEAR_ALTITUDE_SPRINGS_3D(PARTICLES<TV>& particles,TETRAHEDRON_MESH& tetrahedron_mesh);
    virtual ~LINEAR_ALTITUDE_SPRINGS_3D();

//#####################################################################
    void Set_Stiffness_Based_On_Reduced_Mass(const T scaling_coefficient); // assumes mass and restlength are already defined
    void Set_Restlength_From_Particles();
    void Set_Restlength_From_Material_Coordinates(ARRAY_VIEW<const TV> material_coordinates);
    void Set_Overdamping_Fraction(const T overdamping_fraction); // 1 is critically damped
    void Set_Overdamping_Fraction(const ARRAY<VECTOR<T,4> >& overdamping_fraction); // 1 is critically damped
    void Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction=1); // 1 is critically damped
    void Fill_Node_Indices(int i,int j,int k,int l,int isolated_node_number,int& node1,int& node2,int& node3,int& node4) const;
    bool Fill_Spring_State(int t,int isolated_node_number,int node1,int node2,int node3,int node4,SPRING_STATE& state);
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    int Velocity_Dependent_Forces_Size() const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    bool Compute_Strain_Rate_And_Strain(int t,int isolated_node_number,int node1,int node2,int node3,int node4,T& strain_rate,T& strain) const;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
    void Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name="") const PHYSBAM_OVERRIDE;

    T Potential_Energy(const int t,const T time) const;
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
    T Residual_Energy(const T time) const PHYSBAM_OVERRIDE;
    void Save_Potential_Energy(const T time) PHYSBAM_OVERRIDE;
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
template<class T> LINEAR_ALTITUDE_SPRINGS_3D<T>*
Create_Altitude_Springs(PARTICLES<VECTOR<T,3> >& particles,TETRAHEDRON_MESH& mesh,
    const T stiffness=200,const T overdamping_fraction=2,const bool use_compressed_by_threshold_only=true,const T fraction_compression=.1,const bool limit_time_step_by_strain_rate=true,
    const T max_strain_per_time_step=.1,const bool use_rest_state_for_strain_rate=true,const T restlength_enlargement_fraction=0,const bool verbose=true);
template<class T> LINEAR_ALTITUDE_SPRINGS_3D<T>*
Create_Altitude_Springs(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,
    const T stiffness=4/(1+sqrt((T)2)),const T overdamping_fraction=2,const bool use_compressed_by_threshold_only=true,const T fraction_compression=.1,
    const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,const bool use_rest_state_for_strain_rate=true,const T restlength_enlargement_fraction=0,
    const bool verbose=true);
}
#endif
