//#####################################################################
// Copyright 2004-2008, Eran Guendelman, Geoffrey Irving, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_COLLISIONS
//##################################################################### 
#ifndef __RIGID_BODY_COLLISIONS__
#define __RIGID_BODY_COLLISIONS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Collision_Detection/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_PARTICLE_INTERSECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_ID.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Parallel_Computation/MPI_RIGIDS.h>
#include <string>
namespace PhysBAM{

template<class TV> class RIGIDS_EVOLUTION;
template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES;
template<class TV> class RIGID_BODY_INTERSECTIONS;
template<class TV> struct RIGID_BODY_PARTICLE_INTERSECTION;
template<class TV> class RIGID_BODY_CLUSTER_BINDINGS;
template<class TV> class RIGID_BODY_CONTACT_GRAPH;
template<class TV> class RIGID_BODY;
template<class TV> class RIGID_TRIANGLE_COLLISIONS;
template<class TV> class RIGID_TRIANGLE_COLLISIONS_GEOMETRY;
template<class TV> class FRACTURE_PATTERN;
template<class TV> class RIGID_DEFORMABLE_COLLISIONS;
class RIGID_BODY_SKIP_COLLISION_CHECK;
class RIGID_BODY_COLLISION_MANAGER;
template<class TV> class ARTICULATED_RIGID_BODY;
template<class TV> class RIGID_BODY_STATE;
template<class TV> class RIGID_BODY_COLLISION_PARAMETERS;
template<class TV> class RIGIDS_COLLISION_CALLBACKS;

template<class TV>
class RIGID_BODY_COLLISIONS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;

    typedef bool (*UPDATE_ANALYTIC_CONTACT_PAIR_T)(RIGID_BODY_COLLISIONS<TV>&,RIGIDS_COLLISION_CALLBACKS<TV>&,const int,const int,IMPLICIT_OBJECT<TV>*,IMPLICIT_OBJECT<TV>*,const bool,const int,const T,const T,const T,const bool);
public:
    RIGID_BODY_COLLISION_PARAMETERS<TV>& parameters;
    RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks;
    bool verbose,prune_stacks_from_contact,prune_contact_using_velocity;
    RIGID_BODY_COLLISION_MANAGER* collision_manager;
    ARRAY<ARRAY<VECTOR<int,2> > > precomputed_contact_pairs_for_level;
    ARRAY<ARRAY<VECTOR<int,2> > > saved_contact_pairs_for_level; // updated contact pairs
    RIGID_BODY_SKIP_COLLISION_CHECK& skip_collision_check;

    int collision_pair_iterations;
    // TODO: Make ARRAYS so that this can be indexed safely.
    int contact_level_iterations,contact_pair_iterations;
    int shock_propagation_iterations,shock_propagation_level_iterations,shock_propagation_pair_iterations;
    int push_out_iterations,push_out_level_iterations,push_out_pair_iterations;
    bool use_freezing_with_push_out,use_gradual_push_out;
    T desired_separation_distance;
    bool rolling_friction;
    ARRAY<COLLISION_GEOMETRY_ID> object_indices; // rigid body particle indices
public:
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>& rigids_example_forces_and_velocities;
    COLLISION_GEOMETRY_SPATIAL_PARTITION<COLLISION_GEOMETRY<TV>,ARRAY<COLLISION_GEOMETRY<TV>*,COLLISION_GEOMETRY_ID>,COLLISION_GEOMETRY_ID>* spatial_partition;
    RIGID_BODY_INTERSECTIONS<TV>& intersections;
    RIGID_TRIANGLE_COLLISIONS<TV>* triangle_collisions;
    RIGID_TRIANGLE_COLLISIONS_GEOMETRY<TV>* triangle_collision_geometry;

    RIGID_BODY_CONTACT_GRAPH<TV>& contact_graph;
    ARRAY<FRAME<TV> > displacement;
    // TODO: make this per superfragment once rigid bodies can be in more than one superfragment
public:
    HASHTABLE<int,int> id_to_sid;
    HASHTABLE<TRIPLE<int,int,TV> > rigid_body_particle_intersections; // elements are body1,body2,body1_location
    bool store_collision_intersections_for_projection;
    bool use_static_body_masses;
    bool use_parent_normal;
    HASHTABLE<VECTOR<int,2> > pairs_processed_by_collisions;
    ARRAY<VECTOR<int,2> > contact_pairs_from_collisions;
    HASHTABLE<VECTOR<int,2> > pairs_processed_by_contact;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_body_cluster_bindings;
    ARRAY<ARRAY<int> > contact_stack;
    HASHTABLE<VECTOR<int,2>,TRIPLE<int,int,int> > pairs_scale;  
    HASHTABLE<VECTOR<std::string,2>,bool (*)(RIGID_BODY_COLLISIONS&,const int,const int,IMPLICIT_OBJECT<TV>*,IMPLICIT_OBJECT<TV>*,const T,const T,const bool)> analytic_collision_registry;
    HASHTABLE<VECTOR<std::string,2>,UPDATE_ANALYTIC_CONTACT_PAIR_T> analytic_contact_registry;
    FRACTURE_PATTERN<T>* fracture_pattern;
    VECTOR<ARRAY<int>,2> added_bodies;
    VECTOR<int,2> fractured_bodies;

    ARRAY<JOINT_ID> contact_joints;
    ARRAY<VECTOR<int,2> > saved_pairs;
    MPI_RIGIDS<TV>* mpi_rigids;
    ARRAY<TWIST<TV> > mpi_rigid_velocity_save;
    ARRAY<T_SPIN> mpi_rigid_angular_momentum_save;
    ARRAY<TV> mpi_rigid_X_save;
    ARRAY<ROTATION<TV> > mpi_rigid_rotation_save;

    OPERATION_HASH<COLLISION_GEOMETRY_ID> get_potential_collisions_already_added;

    RIGID_BODY_COLLISIONS(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,RIGID_BODY_COLLISION_PARAMETERS<TV>& parameters_input,
        RIGIDS_COLLISION_CALLBACKS<TV>& collision_callbacks_input,RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>& rigids_example_forces_and_velocities_input);
    virtual ~RIGID_BODY_COLLISIONS();

    void Set_Collision_Pair_Iterations(const int collision_pair_iterations_input=20)
    {collision_pair_iterations=collision_pair_iterations_input;}

    void Set_Contact_Level_Iterations(const int contact_level_iterations_input=2)
    {contact_level_iterations=contact_level_iterations_input;}

    void Set_Contact_Pair_Iterations(const int contact_pair_iterations_input=20)
    {contact_pair_iterations=contact_pair_iterations_input;}
    
    void Set_Shock_Propagation_Iterations(const int shock_propagation_iterations_input=1)
    {shock_propagation_iterations=shock_propagation_iterations_input;}

    void Set_Shock_Propagation_Level_Iterations(const int shock_propagation_level_iterations_input=2)
    {shock_propagation_level_iterations=shock_propagation_level_iterations_input;}

    void Set_Shock_Propagation_Pair_Iterations(const int shock_propagation_pair_iterations_input=10)
    {shock_propagation_pair_iterations=shock_propagation_pair_iterations_input;}
    
    void Set_Push_Out_Iterations(const int push_out_iterations_input=1)
    {push_out_iterations=push_out_iterations_input;}

    void Set_Push_Out_Level_Iterations(const int push_out_level_iterations_input=2)
    {push_out_level_iterations=push_out_level_iterations_input;}

    void Set_Push_Out_Pair_Iterations(const int push_out_pair_iterations_input=5)
    {push_out_pair_iterations=push_out_pair_iterations_input;}

    void Use_Gradual_Push_Out(const bool use_gradual_push_out_input=true)
    {use_gradual_push_out=use_gradual_push_out_input;}

    void Use_Freezing_With_Push_Out(const bool use_freezing_with_push_out_input=true)
    {use_freezing_with_push_out=use_freezing_with_push_out_input;}

    void Set_Desired_Separation_Distance(const T desired_separation_distance_input=0)
    {desired_separation_distance=desired_separation_distance_input;}

    void Use_Rolling_Friction(const bool rolling_friction_input=true)
    {rolling_friction=rolling_friction_input;}

//#####################################################################
    void Initialize_Data_Structures(const bool reset=true);
    static void Adjust_Bounding_Boxes(RIGID_BODY_COLLECTION<TV>& rigid_body_collection,const T false_thickness=1,const T extra_padding_distance=0);
    bool Get_First_Intersection_Point(const int id_1,const int id_2,T& smallest_value,int& smallest_index,
        TV& collision_location,TV& collision_normal,TV& collision_relative_velocity,const bool ignore_separating,const T desired_separation_distance,bool& ignored_separating,const T dt,const bool swap=false);
    bool Get_Deepest_Intersection_Point(const int id_1,const int id_2,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,T& smallest_value,int& smallest_index,
        TV& collision_location,TV& collision_normal,TV& collision_relative_velocity,const bool ignore_separating,const T desired_separation_distance,bool& ignored_separating)const;
    void Get_Contact_Pairs(const T dt,const T time,ARRAY<VECTOR<int,2> >& pairs);
    void Compute_Contact_Graph(const T dt,const T time,ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body);
    void Process_Push_Out_Legacy();
    void Process_Push_Out(const bool perform_collision_body_collisions,const T residual_push_out_depth);
    void Process_Push_Out_Projected_Gauss_Seidel();
    void Get_Rigid_Bodies_Intersecting_Rigid_Body(const int particle_index,ARRAY<int>& rigid_bodies,ARRAY<TV>& collision_locations,ARRAY<TV>& body_distances,
        ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const T residual_push_out_depth) const;
    bool Push_Out_From_Rigid_Body(RIGID_BODY<TV>& rigid_body,ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> >& particle_intersections,const T move_fraction,const T residual_push_out_depth);
    bool Check_For_Any_Interpenetration();
    ARRAY<VECTOR<int,2> > Find_All_Bounding_Box_Pairs(const T thickness);
    void Print_Interpenetration_Statistics();
    void Add_Elastic_Collisions(const T dt,const T time);
    void Process_Contact_Using_Graph(const T dt,const T time,ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body,const bool correct_contact_energy,const bool use_saved_pairs=false);
    void Shock_Propagation_Using_Graph(const T dt,const T time,ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body,const bool use_saved_pairs=false);
    bool Update_Collision_Pair(const int id_1,const int id_2,const T dt,const T time,const bool mpi_one_ghost);
    bool Update_Collision_Pair_Helper(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const T dt,const T time,const TV& collision_location,const TV& collision_normal,
        const TV& collision_relative_velocity,const bool mpi_one_ghost);
    bool Update_Levelset_Collision_Pair(const int id_1,const int id_2,const T dt,const T time,const bool mpi_one_ghost);
    bool Update_Surface_Collision_Pair(const int id_1,const int id_2,const T dt,const T time,const bool mpi_one_ghost);
    bool Update_Analytic_Multibody_Collision(RIGID_BODY<TV>& body1,RIGID_BODY<TV>& body2,const T dt,const T time,const bool mpi_one_ghost);
    bool Update_Analytic_Multibody_Collision(const int id_1,const int id_2,MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>& multibody,IMPLICIT_OBJECT<TV>& levelset,const T dt,const T time,const bool mpi_one_ghost);
    bool Either_Body_Collides_With_The_Other(const int rigid_body_id_1,const int rigid_body_id_2) const;
    bool Body_Collides_With_The_Other(const int rigid_body_id_1,const int rigid_body_id_2) const;
    void Apply_Stacking_Contact();
    int Get_Rigid_Body_Depth(int i,ARRAY<int,int>& depths);
    void Compute_Contact_Frequency();
    void Construct_Stacks();
    void Register_Analytic_Collisions();
    void Clean_Up_Fractured_Items_From_Lists(ARRAY<VECTOR<int,2> >& pairs,const int current_pair,const bool called_from_contactt);
    void Initialize_All_Contact_Projections(const bool enforce_rigid_rigid_contact_in_cg);
    void Remove_Contact_Joints();
    void Create_Contact_Joint(const RIGID_BODY<TV>& parent,const RIGID_BODY<TV>& child,const TV& location);
    void Get_Bounding_Box_Collision_Pairs_Of_Body(ARRAY<VECTOR<int,2> >& pairs,int id,const T thickness=0);
    void Get_Bounding_Box_Collision_Pairs(const T dt,const T time,ARRAY<VECTOR<int,2> >& pairs,const bool add_contact_pairs,const bool reinitialize_spatial_partition,const T thickness=0);
    void Euler_Step_Position(const int id,const T dt,const T time);

    void Apply_Prestabilization_To_Joint(const T dt,const T time,ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body,const JOINT_ID joint_id,const T epsilon_scale);
    void Set_Level_Temporarily_Static(const int level);
    void Clear_Temporarily_Static();
    friend class RIGID_DEFORMABLE_COLLISIONS<TV>;
//#####################################################################
};
}
#endif
