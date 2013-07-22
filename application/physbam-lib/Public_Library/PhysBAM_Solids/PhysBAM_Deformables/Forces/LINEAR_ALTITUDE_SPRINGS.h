//#####################################################################
// Copyright 2004-2007, Ronald Fedkiw, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_ALTITUDE_SPRINGS
//#####################################################################
#ifndef __LINEAR_ALTITUDE_SPRINGS__
#define __LINEAR_ALTITUDE_SPRINGS__

#include <PhysBAM_Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SOLIDS_FORCES_POLICY.h>
namespace PhysBAM{

template<class TV,int d>
class LINEAR_ALTITUDE_SPRINGS:public DEFORMABLES_FORCES<TV>,public SPRINGS_TAG
{
    typedef typename TV::SCALAR T;
    typedef typename MESH_POLICY<d>::MESH T_MESH;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;using BASE::max_strain_per_time_step;using BASE::use_rest_state_for_strain_rate;using BASE::Invalidate_CFL;
protected:
    using BASE::cfl_number;
public:
    typedef typename FORCE_ELEMENTS::ITERATOR ELEMENT_ITERATOR;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

    T_MESH& mesh;

    struct SPRING_PARAMETER
    {
        SPRING_PARAMETER()
            :youngs_modulus(0),restlength(0),visual_restlength(0),damping(0)
        {}

        T youngs_modulus;
        T restlength,visual_restlength;
        T damping;
    };
    struct PLASTIC_PARAMETER
    {
        PLASTIC_PARAMETER()
            :yield_strain(0),hardening(0),visual_restlength(0)
        {}

        T yield_strain;
        T hardening;
        T visual_restlength;
    };

    ARRAY<VECTOR<SPRING_PARAMETER,d+1> > parameters;
    bool use_plasticity;
    ARRAY<VECTOR<PLASTIC_PARAMETER,d+1> > plastic_parameters;
    T plasticity_clamp_ratio;

    bool use_shortest_spring_only;

    bool cache_strain;
    mutable ARRAY<VECTOR<T,2> > strains_of_spring; // VECTOR<T,2>(strain_rate, strain)
    mutable ARRAY<VECTOR<VECTOR<T,2>,d+1> > strains_of_spring_all_springs; // VECTOR<T,2>(strain_rate, strain)

    FORCE_ELEMENTS force_elements;
protected:

    struct SPRING_STATE{
        SPRING_STATE()
            :node(0),coefficient(0),current_length(0),sqrt_coefficient(0)
        {}

        int node;
        T coefficient,current_length,sqrt_coefficient;
        TV direction;
        VECTOR<T,d> barycentric;
    };
    ARRAY<SPRING_STATE> spring_states;
    ARRAY<VECTOR<SPRING_STATE,4> > spring_states_all_springs;
    bool use_springs_compressed_beyond_threshold; // only use the springs compressed enough
    T spring_compression_fraction_threshold;
    bool print_number_used;
public:

    LINEAR_ALTITUDE_SPRINGS(PARTICLES<TV>& particles,T_MESH& mesh_input);

    virtual ~LINEAR_ALTITUDE_SPRINGS();

    void Print_Number_Used(const bool print_number_used_input=true)
    {print_number_used=print_number_used_input;}

    void Use_Springs_Compressed_Beyond_Threshold_Only(const bool use_springs_compressed_beyond_threshold_input=true,const T threshold_fraction=.25)
    {use_springs_compressed_beyond_threshold=use_springs_compressed_beyond_threshold_input;spring_compression_fraction_threshold=threshold_fraction;}

//#####################################################################
    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    void Compute_Plasticity(const int h,const int t,const T current_length);
    void Clamp_Restlength_With_Fraction_Of_Springs(const T fraction=.01);
    void Print_Restlength_Statistics() const;
    void Set_Stiffness(const T youngs_modulus_input);
    void Set_Stiffness(const ARRAY<VECTOR<T,d+1> >& youngs_modulus_input);
    void Set_Restlength(const ARRAY<VECTOR<T,d+1> >& restlength_input);
    virtual void Clamp_Restlength(const T clamped_restlength);
    void Set_Damping(const ARRAY<VECTOR<T,d+1> >& damping_input);
    void Set_Damping(const T damping_input);
    void Enable_Plasticity(const ARRAY<VECTOR<T,d+1> >& plastic_yield_strain_input,const ARRAY<VECTOR<T,d+1> >& plastic_hardening_input,const T plasticity_clamp_ratio_input=4);
    void Enable_Plasticity(const T plastic_yield_strain_input,const T plastic_hardening_input,const T plasticity_clamp_ratio_input=4);
//#####################################################################
};

template<class T,class TV,class T_MESH> typename SOLIDS_FORCES_POLICY<TV,T_MESH::dimension>::LINEAR_ALTITUDE_SPRINGS*
Create_Altitude_Springs_Base(PARTICLES<TV>& particles,T_MESH& mesh,const T stiffness,const T overdamping_fraction,
    const bool use_compressed_by_threshold_only,const T fraction_compression,const bool limit_time_step_by_strain_rate,const T max_strain_per_time_step,
    const bool use_rest_state_for_strain_rate,const T restlength_enlargement_fraction,const bool verbose);

}
#endif
