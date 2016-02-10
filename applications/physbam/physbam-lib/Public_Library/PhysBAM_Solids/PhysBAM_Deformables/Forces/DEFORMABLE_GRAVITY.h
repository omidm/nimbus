//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_GRAVITY
//#####################################################################
#ifndef __DEFORMABLE_GRAVITY__
#define __DEFORMABLE_GRAVITY__

#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/POINTWISE_DEFORMABLE_FORCE.h>
namespace PhysBAM{

template<class TV>
class DEFORMABLE_GRAVITY:public POINTWISE_DEFORMABLE_FORCE<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef POINTWISE_DEFORMABLE_FORCE<TV> BASE;
    using BASE::particles;using BASE::influenced_particles;using BASE::mpi_solids;using BASE::force_particles;using BASE::is_simulated;
    typedef typename FORCE_ELEMENTS::ITERATOR ELEMENT_ITERATOR;

    T gravity;
    TV downward_direction;

    ARRAY<T> residual_PE;
    ARRAY<T> force_estimates;
    ARRAY<ARRAY<int> > incident_nodes;
    ARRAY<ARRAY<int> > incident_elements;
    ARRAY<T> potential_energy_save;
    ARRAY<T> delta_PE;
    T total_delta_PE;
public:

    DEFORMABLE_GRAVITY(PARTICLES<TV>& particles_input,ARRAY<int>* influenced_particles_input,const T gravity_input=9.8,const TV& downward_direction_input=-TV::Axis_Vector(2-(TV::dimension==1)))
        :POINTWISE_DEFORMABLE_FORCE<TV>(particles_input,influenced_particles_input),gravity(gravity_input),downward_direction(downward_direction_input)
    {}

    DEFORMABLE_GRAVITY(PARTICLES<TV>& particles_input,const bool influence_all_particles_input,const T gravity_input=9.8,const TV& downward_direction_input=-TV::Axis_Vector(2-(TV::dimension==1)))
        :POINTWISE_DEFORMABLE_FORCE<TV>(particles_input,influence_all_particles_input),gravity(gravity_input),downward_direction(downward_direction_input)
    {}

    template<class T_MESH>
    DEFORMABLE_GRAVITY(PARTICLES<TV>& particles_input,const T_MESH& mesh,const T gravity_input=9.8,const TV& downward_direction_input=((TV::dimension==1)?((T)-1*TV::All_Ones_Vector()):TV(VECTOR<T,2>(0,-1))))
        :POINTWISE_DEFORMABLE_FORCE<TV>(particles_input,mesh),gravity(gravity_input),downward_direction(downward_direction_input)
    {
        mesh.elements.Flattened().Get_Unique(*influenced_particles);
    }

    virtual ~DEFORMABLE_GRAVITY()
    {}

    void Set_Gravity(const TV& gravity_vector)
    {downward_direction=gravity_vector;gravity=downward_direction.Normalize();}

    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE
    {}

    void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const PHYSBAM_OVERRIDE
    {}

    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE
    {}

    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE
    {}

    int Velocity_Dependent_Forces_Size() const PHYSBAM_OVERRIDE
    {return 0;}

    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const PHYSBAM_OVERRIDE
    {}

    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE
    {}

//#####################################################################
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    T Potential_Energy(int p,const T time) const;
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
    void Save_Potential_Energy(const T time) PHYSBAM_OVERRIDE;
    T Residual_Energy(const T time) const PHYSBAM_OVERRIDE;
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
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
