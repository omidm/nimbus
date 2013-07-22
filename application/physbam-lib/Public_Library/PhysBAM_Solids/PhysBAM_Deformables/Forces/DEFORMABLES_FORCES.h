//#####################################################################
// Copyright 2002-2008, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Igor Neverov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLES_FORCES
//#####################################################################
#ifndef __DEFORMABLES_FORCES__
#define __DEFORMABLES_FORCES__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLES_PARTICLES_FORWARD.h>
namespace PhysBAM{

template<class TV> class MPI_SOLIDS;
class SEGMENT_MESH;

template<class TV>
struct FORCE_DATA{
    typedef typename TV::SCALAR T;

    FORCE_DATA():state(0)
    {}

    std::string name;
    T state; // holds signed strech for springs
    TV first_action_point,second_action_point;
};

template<class TV>
class DEFORMABLES_FORCES:public NONCOPYABLE
{
public:
    typedef typename TV::SCALAR T;
    PARTICLES<TV>& particles;
protected:
    T cfl_number;
    bool allow_external_cfl_number;
public:
    struct FREQUENCY_DATA
    {
        FREQUENCY_DATA()
            :elastic_squared(0),damping(0)
        {}

        T elastic_squared,damping;
    };
    bool cfl_initialized;
    
    bool use_rest_state_for_strain_rate;
    bool limit_time_step_by_strain_rate;
    T max_strain_per_time_step; // for limiting the timestep in the CFL calculation
    bool use_velocity_independent_forces,use_velocity_dependent_forces,use_force_differential,use_implicit_velocity_independent_forces;
    bool use_position_based_state;
    int unique_id;
    bool compute_half_forces;

    DEFORMABLES_FORCES(PARTICLES<TV>& particles);
    virtual ~DEFORMABLES_FORCES();

    static int Get_Unique_Id()
    {static int next_unique_id=0;return ++next_unique_id;}

    void Set_CFL_Number(const T cfl_number_input)
    {if(allow_external_cfl_number) cfl_number=cfl_number_input;}

    bool CFL_Valid() const
    {return cfl_initialized;}

    void Invalidate_CFL()
    {cfl_initialized=false;}

    void Validate_CFL()
    {cfl_initialized=true;}

//#####################################################################
    virtual void Use_Rest_State_For_Strain_Rate(const bool use_rest_state_for_strain_rate_input=true);
    virtual void Limit_Time_Step_By_Strain_Rate(const bool limit_time_step_by_strain_rate_input=true,const T max_strain_per_time_step_input=.1);
    virtual void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const=0;
    virtual void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)=0;
    virtual void Update_Position_Based_State(const T time,const bool is_position_update){}
    virtual void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const=0;
    virtual void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const=0;
    virtual int Velocity_Dependent_Forces_Size() const;
    virtual void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const;
    virtual void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const;
    virtual void Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const;
    virtual void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const;
    virtual void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const;
    virtual void Enforce_Definiteness(const bool enforce_definiteness_input);
    virtual T CFL_Strain_Rate() const=0;
    virtual void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)=0;
    virtual T Potential_Energy(const T time) const;
    virtual T Residual_Energy(const T time) const;
    virtual void Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name="") const;
    virtual void Save_Potential_Energy(const T time);
    virtual void Compute_Energy_Error(ARRAY_VIEW<const TV> velocity_save,const T time,const T dt);
    virtual void Prepare_Energy_Correction_Force();
    virtual void Compute_Energy_Correction_Force(ARRAY_VIEW<const TV> velocity_save,const int max_particle_degree,const T time,const T dt);
    virtual void Apply_Energy_Correction(const T time,const T dt);
    virtual void Add_Connectivity(ARRAY<int>& particle_degree);
    virtual void Compute_Previously_Applied_Forces();
    virtual void Setup_Set_Velocity_From_Positions(const T time,const bool is_position_update,const bool reset_alphas);
    virtual int Get_Element_Count();
    virtual FORCE_ELEMENTS* Get_Force_Elements();
    virtual ARRAY<int>* Incident_Nodes(const int force_element);
    virtual TV Get_Direction(const int force_element);
    virtual T Get_Combined_One_Over_Mass(const int force_element);
    virtual ARRAY<int>* Incident_Force_Elements(const int particle);
    virtual void Set_Force(const int force_element,const T force);
    virtual void Choose_Solution(const bool use_orig_force,const int force_element,const T dt,const T alpha1,const T alpha2,ARRAY<T>& v_n_hats);
    virtual TV Get_Force(const int force_element,const int particle,const bool use_original_force);
    virtual T Get_Force(const int force_element);
    virtual void Get_Damping_Force(const int particle,TV& damping_force,const T dt,const bool use_coefficient);
    virtual void Update_Residual_Energy(const int force_element,const T residual_energy,const T time);
    virtual T Get_Residual_Energy(const int force_element);
    virtual void Compute_Quadratic_Contribution_For_Force(T& A,T& a,T&c,const T dt,const int force_element,const T combined_one_over_mass,const bool ignore_PE_terms);
    virtual void Compute_Quadratic_Contribution_For_Node(T& B,T& C,T&b,const T dt,const int node,const int force_element,const T combined_one_over_mass,const T v_n_correction,
        const bool ignore_PE_terms);
    virtual void Compute_Quadratic_Contribution_For_Residual(T& B,T& C,T& a,T&b,T&c,const T dt,const T time,const int force_element,const bool ignore_PE_terms);
    virtual void Store_Delta_PE(const T time);
    virtual T Get_Total_Delta_PE();
    virtual void Save_And_Reset_Elastic_Coefficient();
    virtual void Restore_Elastic_Coefficient();
    virtual void Store_Velocities();
//#####################################################################
};

// tag classes to allow iteration over certain types of forces
struct SPRINGS_TAG{};
struct FINITE_VOLUME_TAG{};

}
#include <PhysBAM_Solids/PhysBAM_Deformables/Read_Write/Forces_And_Torques/READ_WRITE_FORCE_DATA.h>
#endif
