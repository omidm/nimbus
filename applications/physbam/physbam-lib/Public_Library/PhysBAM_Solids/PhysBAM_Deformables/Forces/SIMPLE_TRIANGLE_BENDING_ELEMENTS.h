//#####################################################################
// Copyright 2002-2007, Zhaosheng Bao, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SIMPLE_TRIANGLE_BENDING_ELEMENTS
//#####################################################################
#ifndef __SIMPLE_TRIANGLE_BENDING_ELEMENTS__
#define __SIMPLE_TRIANGLE_BENDING_ELEMENTS__

#include <PhysBAM_Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
namespace PhysBAM{

template<class TV> class BINDING_LIST;

template<class T_input>
class SIMPLE_TRIANGLE_BENDING_ELEMENTS:public LINEAR_SPRINGS<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef LINEAR_SPRINGS<TV> BASE;
    using BASE::particles;using BASE::max_strain_per_time_step;using BASE::use_rest_state_for_strain_rate;using BASE::Limit_Time_Step_By_Strain_Rate;
    typedef typename FORCE_ELEMENTS::ITERATOR QUADRUPLE_ITERATOR;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

    ARRAY<VECTOR<int,4> >& bending_quadruples; // for each quadruple, the bending axis is elements 2 and 3.
    ARRAY<T> bending_stiffness;
    ARRAY<T> damping;
private:
    ARRAY<VECTOR<int,4> > bending_quadruples_default;
    ARRAY<T> damping_coefficient;
    bool print_number_ignored;
    FORCE_ELEMENTS force_quadruples;
    ARRAY<VECTOR<int,2> > linear_bindings;
    SEGMENT_MESH spring_connectivity;
    BINDING_LIST<TV>& binding_list;
    ARRAY<T> fictitious_edge_length;
    ARRAY<T> base_youngs_modulus;
public:
    using BASE::force_segments;using BASE::states;using BASE::current_lengths;using BASE::restlength;using BASE::youngs_modulus;using BASE::visual_restlength;

    SIMPLE_TRIANGLE_BENDING_ELEMENTS(PARTICLES<TV>& particles,BINDING_LIST<TV>& binding_list_input,bool implicit);
    SIMPLE_TRIANGLE_BENDING_ELEMENTS(PARTICLES<TV>& particles,ARRAY<VECTOR<int,4> >& bending_quadruples_input,BINDING_LIST<TV>& binding_list_input,bool implicit);
    virtual ~SIMPLE_TRIANGLE_BENDING_ELEMENTS();

    void Print_Number_Ignored(const bool print_number_ignored_input=true)
    {print_number_ignored=print_number_ignored_input;}

    void Set_Stiffness(const T bending_stiffness_input)
    {ARRAYS_COMPUTATIONS::Fill(bending_stiffness,bending_stiffness_input);}

    void Set_Stiffness(ARRAY_VIEW<const T> bending_stiffness_input)
    {bending_stiffness=bending_stiffness_input;}

    void Set_Damping(const T damping_input)
    {ARRAYS_COMPUTATIONS::Fill(damping,damping_input);}

    void Set_Damping(ARRAY_VIEW<const T> damping_input)
    {damping=damping_input;}

//#####################################################################
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    void Set_Quadruples_From_Triangle_Mesh(TRIANGLE_MESH& mesh);
    void Set_Constants_From_Particles(const T material_stiffness,const T material_damping);
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
//#####################################################################
};

template<class T> SIMPLE_TRIANGLE_BENDING_ELEMENTS<T>*
Create_Simple_Bending_Elements(PARTICLES<VECTOR<T,3> >& particles,TRIANGLE_MESH& mesh,BINDING_LIST<VECTOR<T,3> >& binding_list,const T stiffness=1e-3,const T damping=1e-3,
    const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,const bool use_plasticity=false,const T plastic_yield=3,const T plastic_hardening=1,
    const T cutoff_fraction_of_minimum_area=0,const T cutoff_fraction_of_triangles=0,const bool verbose=true,const bool implicit=false);

template<class T> SIMPLE_TRIANGLE_BENDING_ELEMENTS<T>*
Create_Simple_Bending_Elements(TRIANGULATED_SURFACE<T>& triangulated_surface,BINDING_LIST<VECTOR<T,3> >& binding_list,const T stiffness=1e-3,
    const T damping=1e-3,const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,const bool use_plasticity=false,const T plastic_yield=3,const T plastic_hardening=1,
    const T cutoff_fraction_of_minimum_area=0,const T cutoff_fraction_of_triangles=0,const bool verbose=true,const bool implicit=false);

}
#endif
