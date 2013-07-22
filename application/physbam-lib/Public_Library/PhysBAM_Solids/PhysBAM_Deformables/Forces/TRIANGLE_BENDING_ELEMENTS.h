//#####################################################################
// Copyright 2002-2007, Zhaosheng Bao, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_BENDING_ELEMENTS
//#####################################################################
#ifndef __TRIANGLE_BENDING_ELEMENTS__
#define __TRIANGLE_BENDING_ELEMENTS__

#include <PhysBAM_Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
namespace PhysBAM{

template<class T_input>
class TRIANGLE_BENDING_ELEMENTS:public DEFORMABLES_FORCES<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;using BASE::max_strain_per_time_step;using BASE::use_rest_state_for_strain_rate;using BASE::Limit_Time_Step_By_Strain_Rate;
    typedef typename FORCE_ELEMENTS::ITERATOR QUADRUPLE_ITERATOR;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

    ARRAY<VECTOR<int,4> >& bending_quadruples; // for each quadruple, the bending axis is elements 2 and 3.
    ARRAY<T> bending_stiffness;
    ARRAY<T> sine_half_rest_angle;
    ARRAY<T> damping;
    ARRAY<T>* plastic_yield;
    ARRAY<T>* plastic_hardening;
    ARRAY<T>* sine_half_elastic_angle;
    T area_cutoff;
private:
    ARRAY<VECTOR<int,4> > bending_quadruples_default;
    ARRAY<T> elastic_s,damping_coefficient;
    ARRAY<VECTOR<VECTOR<T,3>,4> > force_directions; // n1,dj,dk,n2
    HASHTABLE<VECTOR<int,2>,int>* reference_bending_quadruples_hashtable;
    ARRAY<VECTOR<int,4> >* bending_quadruples_save;
    ARRAY<T>* reference_bending_stiffness;
    ARRAY<T>* reference_sine_half_rest_angle;
    ARRAY<T>* reference_damping;
    ARRAY<T>* plastic_yield_save;
    ARRAY<T>* reference_plastic_hardening;
    ARRAY<T>* sine_half_elastic_angle_save;
    bool print_number_ignored;
    FORCE_ELEMENTS force_quadruples;
public:

    TRIANGLE_BENDING_ELEMENTS(PARTICLES<TV>& particles)
        :DEFORMABLES_FORCES<TV>(particles),
        bending_quadruples(bending_quadruples_default),plastic_yield(0),plastic_hardening(0),sine_half_elastic_angle(0),area_cutoff(0),
        reference_bending_quadruples_hashtable(0),bending_quadruples_save(0),reference_sine_half_rest_angle(0),plastic_yield_save(0),sine_half_elastic_angle_save(0)
    {
        Print_Number_Ignored();
        Limit_Time_Step_By_Strain_Rate(false);
    }

    TRIANGLE_BENDING_ELEMENTS(PARTICLES<TV>& particles,ARRAY<VECTOR<int,4> >& bending_quadruples_input)
        :DEFORMABLES_FORCES<TV>(particles),bending_quadruples(bending_quadruples_input),bending_stiffness(bending_quadruples_input.m),
        sine_half_rest_angle(bending_quadruples_input.m),damping(bending_quadruples_input.m),plastic_yield(0),plastic_hardening(0),sine_half_elastic_angle(0),area_cutoff(0),
        reference_bending_quadruples_hashtable(0),bending_quadruples_save(0),reference_sine_half_rest_angle(0),plastic_yield_save(0),sine_half_elastic_angle_save(0)
    {
        Print_Number_Ignored();
        Limit_Time_Step_By_Strain_Rate(false);
    }

    ~TRIANGLE_BENDING_ELEMENTS()
    {delete plastic_yield;delete plastic_hardening;delete reference_bending_quadruples_hashtable;delete reference_sine_half_rest_angle;
    delete sine_half_elastic_angle_save;delete plastic_yield_save;delete bending_quadruples_save;}

    void Print_Number_Ignored(const bool print_number_ignored_input=true)
    {print_number_ignored=print_number_ignored_input;}

    void Set_Stiffness(const T bending_stiffness_input)
    {ARRAYS_COMPUTATIONS::Fill(bending_stiffness,bending_stiffness_input);}

    void Set_Stiffness(ARRAY_VIEW<const T> bending_stiffness_input)
    {bending_stiffness=bending_stiffness_input;}

    void Set_Sine_Half_Rest_Angle(const T sine_half_rest_angle_input)
    {ARRAYS_COMPUTATIONS::Fill(sine_half_rest_angle,sine_half_rest_angle_input);}

    void Set_Sine_Half_Rest_Angle(ARRAY_VIEW<const T> sine_half_rest_angle_input)
    {sine_half_rest_angle=sine_half_rest_angle_input;}

    void Set_Damping(const T damping_input)
    {ARRAYS_COMPUTATIONS::Fill(damping,damping_input);}

    void Set_Damping(ARRAY_VIEW<const T> damping_input)
    {damping=damping_input;}

    template<class T_FIELD>
    void Enable_Plasticity(const T_FIELD& plastic_yield_input,const T_FIELD& plastic_hardening_input)
    {if(!plastic_yield) plastic_yield=new ARRAY<T>;plastic_yield->Resize(bending_quadruples.m,false,false);
    if(!plastic_hardening) plastic_hardening=new ARRAY<T>;plastic_hardening->Resize(bending_quadruples.m,false,false);
    ARRAYS_COMPUTATIONS::Fill(*plastic_yield,plastic_yield_input);ARRAYS_COMPUTATIONS::Fill(*plastic_hardening,plastic_hardening_input);
    if(!sine_half_elastic_angle) sine_half_elastic_angle=new ARRAY<T>;*sine_half_elastic_angle=sine_half_rest_angle;}

    void Set_Area_Cutoff_From_Triangulated_Surface(TRIANGULATED_SURFACE<T>& triangulated_surface,const T cutoff_ratio=(T)1)
    {area_cutoff=cutoff_ratio*triangulated_surface.Minimum_Area();}

private:
    T Sine_Half_Angle_Between(const TV& n1,const TV& n2,const TV& e) const
    {T sine_half_psi=(T).5*(n1-n2).Magnitude();
    if(TV::Triple_Product(n1,n2,e)<0) sine_half_psi=-sine_half_psi;
    return sine_half_psi;}
public:

    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE
    {} // TODO: Implement Me?

//#####################################################################
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    void Set_Area_Cutoff_With_Fraction_Of_Triangles(TRIANGULATED_SURFACE<T>& triangulated_surface,const T fraction=.01);
    void Set_Quadruples_From_Triangle_Mesh(TRIANGLE_MESH& mesh);
    void Set_Quadruples_From_Reference_Triangle_Mesh(TRIANGLE_MESH& mesh,const ARRAY<int>& triangle_map_to_reference);
    void Set_Constants_From_Particles(const T material_stiffness,const T material_damping);
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
    T Compute_Discrete_Shell_Energy();
    void Initialize_Reference_Quantities(const int hash_multiple=2);
    void Copy_Back_Reference_Quantities(const ARRAY<int>& corresponding_node_in_reference_embedded_triangulated_surface);
    void Initialize_Save_Quantities();
    void Copy_Back_Save_Quantities(const ARRAY<int>& node_map_to_saved);
//#####################################################################
};

template<class T> TRIANGLE_BENDING_ELEMENTS<T>*
Create_Bending_Elements(PARTICLES<VECTOR<T,3> >& particles,TRIANGLE_MESH& mesh,const T stiffness=1e-3,const T damping=1e-3,
    const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,const bool use_plasticity=false,const T plastic_yield=3,const T plastic_hardening=1,
    const T cutoff_fraction_of_minimum_area=0,const T cutoff_fraction_of_triangles=0,const bool verbose=true)
{
    TRIANGLE_BENDING_ELEMENTS<T>* bend=new TRIANGLE_BENDING_ELEMENTS<T>(particles);
    bend->Set_Quadruples_From_Triangle_Mesh(mesh);
    bend->Set_Constants_From_Particles(stiffness,damping);
    bend->Limit_Time_Step_By_Strain_Rate(limit_time_step_by_strain_rate,max_strain_per_time_step);
    if(use_plasticity) bend->Enable_Plasticity(plastic_yield,plastic_hardening);
    if(cutoff_fraction_of_minimum_area){
        TRIANGULATED_SURFACE<T> triangulated_surface(mesh,particles);
        bend->Set_Area_Cutoff_From_Triangulated_Surface(triangulated_surface,cutoff_fraction_of_minimum_area);}
    if(cutoff_fraction_of_triangles){
        TRIANGULATED_SURFACE<T> triangulated_surface(mesh,particles);
        bend->Set_Area_Cutoff_With_Fraction_Of_Triangles(triangulated_surface,cutoff_fraction_of_triangles);}
    return bend;
}

template<class T> TRIANGLE_BENDING_ELEMENTS<T>*
Create_Bending_Elements(TRIANGULATED_SURFACE<T>& triangulated_surface,const T stiffness=1e-3,
    const T damping=1e-3,const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,const bool use_plasticity=false,const T plastic_yield=3,const T plastic_hardening=1,
    const T cutoff_fraction_of_minimum_area=0,const T cutoff_fraction_of_triangles=0,const bool verbose=true)
{
    return Create_Bending_Elements(dynamic_cast<PARTICLES<VECTOR<T,3> >&>(triangulated_surface.particles),triangulated_surface.mesh,stiffness,damping,limit_time_step_by_strain_rate,max_strain_per_time_step,use_plasticity,
        plastic_yield,plastic_hardening,cutoff_fraction_of_minimum_area,cutoff_fraction_of_triangles,verbose);
}

}
#endif
