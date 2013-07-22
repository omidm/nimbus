//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARTICULATED_RIGID_BODY
//#####################################################################
#ifndef __ARTICULATED_RIGID_BODY__
#define __ARTICULATED_RIGID_BODY__

#include <PhysBAM_Tools/Data_Structures/DIRECTED_GRAPH.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_1D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_POLICY.h>
namespace PhysBAM{
template<class TV> class JOINT_MESH;
template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> class ARTICULATED_VECTOR;
template<class TV>
class ARTICULATED_RIGID_BODY_BASE:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;typedef typename TV::SPIN T_SPIN;
    typedef typename RIGID_BODY_POLICY<TV>::WORLD_SPACE_INERTIA_TENSOR T_WORLD_SPACE_INERTIA_TENSOR;
public:
    enum WORKAROUND {d=TV::m,s=T_SPIN::m};

    DIRECTED_GRAPH<int>* breadth_first_directed_graph;
    ARRAY<ARRAY<JOINT_ID> > process_list;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    JOINT_MESH<TV>& joint_mesh;
    MUSCLE_LIST<TV>* muscle_list;
    ARRAY<T> muscle_activations;
    ARRAY<T_SPIN> angular_momenta_save;
    ARRAY<TV> linear_velocities_save;
    T iterative_tolerance;
    int poststabilization_iterations,max_iterations,poststabilization_projection_iterations;
    bool use_epsilon_scale;
    int contact_level_iterations,shock_propagation_level_iterations,actuation_iterations;
    T line_search_interval_tolerance;
    int max_line_search_iterations;
    bool use_shock_propagation,do_final_pass;
    bool use_pd_actuators,use_muscle_actuators;
    bool enforce_nonnegative_activations;
    bool clamp_negative_activations;
    int activation_optimization_iterations;
    bool global_post_stabilization;
    bool constrain_pd_directions;
    bool verbose;
    bool use_angular_damping;
    T zero_row_tolerance;
    ARRAY<MATRIX_MXN<T>,JOINT_ID> v_to_lambda,lambda_to_delta_v; // components of the projections matrices for poststabilization
    bool use_krylov_poststab;
    bool use_krylov_prestab;
    bool use_poststab_in_cg;
    bool use_prestab;
    bool use_poststab;

    bool check_stale,is_stale;
    int last_read;
    ARRAY<int>* frame_list;

    ARTICULATED_RIGID_BODY_BASE(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input);
    virtual ~ARTICULATED_RIGID_BODY_BASE();

    void Set_Iterative_Tolerance(const T iterative_tolerance_input)
    {iterative_tolerance=iterative_tolerance_input;}

    void Set_Poststabilization_Iterations(const int poststabilization_iterations_input)
    {poststabilization_iterations=poststabilization_iterations_input;}

    void Set_Max_Iterations(const int max_iterations_input)
    {max_iterations=max_iterations_input;}

    void Use_Epsilon_Scale(const bool use_epsilon_scale_input)
    {use_epsilon_scale=use_epsilon_scale_input;}

    void Set_Contact_Level_Iterations(const int contact_level_iterations_input)
    {contact_level_iterations=contact_level_iterations_input;}

    void Set_Shock_Propagation_Level_Iterations(const int shock_propagation_level_iterations_input)
    {shock_propagation_level_iterations=shock_propagation_level_iterations_input;}

    void Set_Actuation_Iterations(const int actuation_iterations_input)
    {actuation_iterations=actuation_iterations_input;}

    void Set_Max_Line_Search_Iterations(const int max_line_search_iterations_input)
    {max_line_search_iterations=max_line_search_iterations_input;}

    void Set_Line_Search_Interval_Tolerance(const T line_search_interval_tolerance_input)
    {line_search_interval_tolerance=line_search_interval_tolerance_input;}

    void Set_Use_Shock_Propagation(const bool use_shock_propagation_input)
    {use_shock_propagation=use_shock_propagation_input;}

    void Set_Do_Final_Pass(const bool do_final_pass_input)
    {do_final_pass=do_final_pass_input;}

    void Use_PD_Actuators()
    {use_pd_actuators=true;use_muscle_actuators=false;}

    void Use_Muscle_Actuators()
    {use_muscle_actuators=true;use_pd_actuators=false;}

    void Use_No_Actuators() 
    {use_muscle_actuators=use_pd_actuators=false;}

    bool Has_Actuators() const
    {return use_muscle_actuators || use_pd_actuators;}

//#####################################################################
    void Remove_All();
    int Parent_Id(JOINT_ID joint_id) const;
    int Child_Id(JOINT_ID joint_id) const;
    RIGID_BODY<TV>* Parent(JOINT_ID joint_id);
    const RIGID_BODY<TV>* Parent(JOINT_ID joint_id) const;
    RIGID_BODY<TV>* Child(JOINT_ID joint_id);
    const RIGID_BODY<TV>* Child(JOINT_ID joint_id) const;
    TV Joint_Location(JOINT_ID joint_id);
    void Set_Consistent_Child_Frame(JOINT_ID joint_id);
    void Initialize_Breadth_First_Directed_Graph(const int root);
    void Update_With_Breadth_First_Directed_Graph(const int root,const int node_to_propagate_from=0);
    void Update_From_Edge_In_Breadth_First_Directed_Graph(const int root,const JOINT_ID edge);
    void Update_Joint_Frames();
    void Store_Velocities_And_Momenta();
    void Restore_Velocities_And_Momenta();
    virtual void Solve_Velocities_for_PD(const T time,const T dt,bool test_system,bool print_matrix)=0;
    virtual void Apply_Poststabilization_To_Joint(const JOINT_ID joint_id,const bool target_pd=false);
    virtual void Output_Articulation_Points(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const=0;
    virtual void Compute_Desired_PD_Velocity(const T dt,const T time){}
    template<class T_CONSTRAINT_FUNCTION> void Compute_Constraint_Correcting_Impulse(T_CONSTRAINT_FUNCTION& constraint_error_function,typename T_CONSTRAINT_FUNCTION::T_IMPULSE& j) const;
    void Apply_Prestabilization_To_Joint(const JOINT_ID joint_id,const T dt,const T epsilon_scale);
    void Delta_Relative_Twist(const JOINT_ID joint_id,const bool target_pd,TV& location,TWIST<TV>& delta_relative_twist);
    void Initialize_Poststabilization_Projection();
    void Poststabilization_Projection(ARRAY_VIEW<TWIST<TV> > twist,const bool symmetric=false);
    void Generate_Process_List_Using_Contact_Graph(const RIGID_BODY_CONTACT_GRAPH<TV>& contact_graph);
    void Apply_Poststabilization(bool test_system,bool print_matrix,const bool target_pd=false,const bool skip_global_post_stabilized_joints=false,const bool angular_damping_only=false);
    void Read(const STREAM_TYPE stream_type,const std::string& directory,const int frame);
    void Write(const STREAM_TYPE stream_type,const std::string& directory,const int frame);
    void Effective_Inertia_Inverse(MATRIX<T,d+s>& inertia_inverse,JOINT_ID joint_id) const;
    void Substitute_Joint_Parent_Body(JOINT_ID joint_id,int new_parent,const FRAME<TV>& frame);
    void Substitute_Joint_Child_Body(JOINT_ID joint_id,int new_child,const FRAME<TV>& frame);
    void Substitute_Joint_Parent_Body(JOINT_ID joint_id,int new_parent);
    void Substitute_Joint_Child_Body(JOINT_ID joint_id,int new_child);
    void Apply_Poststabilization_With_CG(T dt,bool correct_position,bool test_system,bool print_matrix);
    void Resize_System_Vector(ARTICULATED_VECTOR<TV>& v) const;
    TWIST<TV> Joint_Error(JOINT_ID joint_id) const;
private:
    void Poststabilization_Projection_Joint(const JOINT_ID joint_id,ARRAY_VIEW<TWIST<TV> > twist);
    void Update_Child_From_Parents(const int child_id,const ARRAY<int>& parents);
    void Update_Parent_Joint_States(const int child_id,const ARRAY<int>& parents);
//#####################################################################
};
}
//#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_1D.h>
#endif
