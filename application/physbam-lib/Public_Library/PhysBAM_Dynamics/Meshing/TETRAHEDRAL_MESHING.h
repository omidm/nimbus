//#####################################################################
// Copyright 2003-2006, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Mike Rodgers, Tamar Shinar, Eftychios Sifakis, Joey Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TETRAHEDRAL_MESHING
//##################################################################### 
#ifndef __TETRAHEDRAL_MESHING__
#define __TETRAHEDRAL_MESHING__

#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
#include <PhysBAM_Dynamics/Meshing/LEVEL_SET_FORCES_AND_VELOCITIES.h>
namespace PhysBAM{

template<class T> class RED_GREEN_TETRAHEDRA;
template<class TV> class SOLID_BODY_COLLECTION;

    
class EXTRA_REFINEMENT_CRITERIA_HELPER
{
public:
    virtual ~EXTRA_REFINEMENT_CRITERIA_HELPER() {};
    virtual int operator()(int index)=0;
};

template<class T>
class TETRAHEDRAL_MESHING:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;
public:
    SOLIDS_PARAMETERS<TV>& solids_parameters; 
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    SOLIDS_EVOLUTION<TV>* solids_evolution;
    IMPLICIT_OBJECT<TV>* implicit_surface;
    LEVEL_SET_FORCES_AND_VELOCITIES<TV>* level_set_forces_and_velocities;
    const STREAM_TYPE stream_type;
    std::string output_directory;
    int frame;
    T curvature_subdivision_threshold,interpolation_error_subdivision_threshold;
    T maximum_boundary_edge_length;
    T density; // density of the material being meshed
    T boundary_mass_multiplicative_factor; // 1 does nothing
    T dynamic_ether_viscosity;
    bool use_finite_volume;
    T youngs_modulus,poissons_ratio,Rayleigh_coefficient;
    bool use_masses_and_springs;
    T edge_spring_stiffness,edge_spring_overdamping_fraction;
    T altitude_spring_stiffness,altitude_spring_overdamping_fraction;
    bool use_global_quality_criteria_for_early_exit;
    bool replace_green_refinement_with_embedded_t_junctions;
    bool allow_coarsening_to_non_graded_mesh;
    bool use_constant_mass;
    EXTRA_REFINEMENT_CRITERIA_HELPER *extra_refinement_criteria; // should return -1, 0, or 1 for never, maybe, or always refine
    bool symmetric_initial_grid;
private:
    ARRAY<int> map_from_nodes_to_boundary_list;
    ARRAY<ARRAY<int>*> layers; // each sublist is a list of nodes in one layer, with the first sublist being the boundary
    ARRAY<VECTOR<T,3> > boundary_mesh_normals;     
    T worst_boundary_quality,worst_interior_quality;
    ARRAY<ARRAY<int> >* dependent_nodes;
    TRIANGLE_MESH* boundary_mesh;
public:

    TETRAHEDRAL_MESHING(const STREAM_TYPE stream_type);
    ~TETRAHEDRAL_MESHING();

    void Set_Curvature_Subdivision_Threshold(const T curvature_subdivision_threshold_input=.7)
    {curvature_subdivision_threshold=curvature_subdivision_threshold_input;}

    void Set_Interpolation_Error_Subdivision_Threshold(const T interpolation_error_subdivision_threshold_input=.08)
    {interpolation_error_subdivision_threshold=interpolation_error_subdivision_threshold_input;}

    void Set_Maximum_Boundary_Edge_Length(const T maximum_boundary_edge_length_input=FLT_MAX)
    {maximum_boundary_edge_length=maximum_boundary_edge_length_input;}

    void Set_Density(const T density_input=1000)
    {density=density_input;}

    void Increase_Mass_On_Boundary(const T boundary_mass_multiplicative_factor_input=1) // 1 does nothing
    {boundary_mass_multiplicative_factor=boundary_mass_multiplicative_factor_input;}

    void Use_Dynamic_Ether_Viscosity(const T dynamic_ether_viscosity_input=0)
    {dynamic_ether_viscosity=dynamic_ether_viscosity_input;}

    void Use_Finite_Elements(const T youngs_modulus_input=500,const T poissons_ratio_input=.475,const T Rayleigh_coefficient_input=.1)
    {use_finite_volume=true;use_masses_and_springs=false;
    youngs_modulus=youngs_modulus_input;poissons_ratio=poissons_ratio_input;Rayleigh_coefficient=Rayleigh_coefficient_input;}

    void Use_Masses_And_Springs(const T edge_spring_stiffness_input=1e-4,const T edge_spring_overdamping_fraction_input=5,const T altitude_spring_stiffness_input=1e-4,
        const T altitude_spring_overdamping_fraction_input=5)
    {use_masses_and_springs=true;use_finite_volume=false;
    edge_spring_stiffness=edge_spring_stiffness_input;
    edge_spring_overdamping_fraction=edge_spring_overdamping_fraction_input;altitude_spring_stiffness=altitude_spring_stiffness_input;
    altitude_spring_overdamping_fraction=altitude_spring_overdamping_fraction_input;}

    void Use_Global_Quality_Criteria_For_Early_Exit(const bool use_global_quality_criteria_for_early_exit_input=true)
    {use_global_quality_criteria_for_early_exit=use_global_quality_criteria_for_early_exit_input;}

    void Replace_Green_Refinement_With_Embedded_T_Junctions(const bool replace_green_refinement_with_embedded_t_junctions_input=false,
        const bool allow_coarsening_to_non_graded_mesh_input=false)
    {replace_green_refinement_with_embedded_t_junctions=replace_green_refinement_with_embedded_t_junctions_input;
    allow_coarsening_to_non_graded_mesh=allow_coarsening_to_non_graded_mesh_input;}

//#####################################################################
    void Initialize(IMPLICIT_OBJECT<TV>* implicit_surface_input);
    void Snap_Nodes_To_Level_Set_Boundary(const int iterations=100);
    void Initialize_Optimization(const bool verbose=true);
    void Create_Final_Mesh_With_Optimization(const int number_of_initial_steps=4,const int number_of_final_steps=6,const bool verbose=true);
    void Optimization_Sweep(const T compression_fraction,const bool verbose=true);
private:
    void Optimize_Boundary_Layer(const T compression_fraction,const bool reverse=false);
    void Optimize_Interior_Layer(const int layer,const bool reverse=false);
    void Search_For_Best_Position(const int node,const ARRAY<VECTOR<T,3> >& direction,bool include_boundary_terms=false);
    T Quality_Of_Worst_Incident_Tetrahedron(const int node);
    T Quality_Of_Worst_Incident_Boundary_Triangle(const int node);
    T Quality_Of_Worst_Dependent_Tetrahedron(const int node);
    T Quality_Of_Worst_Dependent_Boundary_Triangle(const int node);
    void Compute_Boundary_Mesh_Normals();
    void Update_Dependent_Nodes(const int node);
public:
    void Initialize_Dynamics();
    void Create_Final_Mesh_With_Dynamics(const T time_step,const int number_of_force_steps,const int number_of_velocity_steps,const bool verbose=true);
private:
    void Advance_Dynamics(const T time,const T stopping_time,const bool verbose=true);
public:
    void Create_Initial_Mesh(const T bcc_lattice_cell_size=0,const bool use_adaptive_refinement=false,const int max_subdivision_levels=7,const bool discard_to_get_nice_topology=true,
        const bool verbose=true,const bool use_aggressive_tet_pruning_globally=false,const ARRAY<RANGE<TV> >* bounding_boxes_for_aggressive_pruning=0);
private:
    bool Tetrahedron_Refinement_Criteria(const int index) const;
    void Discard_To_Get_Nice_Topology(RED_GREEN_TETRAHEDRA<T>& redgreen,ARRAY<bool>& keep_tet_flag,const bool verbose=true);
    void Envelope_Interior_Nodes(ARRAY<bool>& keep_tet_flag);
public:
    void Write_Output_Files(const int frame);
//#####################################################################
};
}
#endif
