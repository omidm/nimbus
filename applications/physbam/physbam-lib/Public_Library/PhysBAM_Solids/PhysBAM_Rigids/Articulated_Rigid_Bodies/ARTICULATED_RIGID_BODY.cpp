//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Michael Lentine, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARTICULATED_RIGID_BODY
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Krylov_Solvers/SYMMQMR.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Octave/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_SYSTEM.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_CONTACT_GRAPH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/CONSTRAINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> ARTICULATED_RIGID_BODY_BASE<TV>::
ARTICULATED_RIGID_BODY_BASE(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
    :breadth_first_directed_graph(0),rigid_body_collection(rigid_body_collection_input),joint_mesh(*new JOINT_MESH<TV>),muscle_list(new MUSCLE_LIST<TV>(rigid_body_collection)),iterative_tolerance((T)1e-5),
    poststabilization_iterations(5),max_iterations(1000),poststabilization_projection_iterations(4),use_epsilon_scale(true),contact_level_iterations(5),shock_propagation_level_iterations(5),
    actuation_iterations(50),line_search_interval_tolerance((T)1e-6),max_line_search_iterations(100),use_shock_propagation(true),do_final_pass(true),use_pd_actuators(false),
    use_muscle_actuators(false),enforce_nonnegative_activations(true),clamp_negative_activations(true),activation_optimization_iterations(50),global_post_stabilization(false),
    constrain_pd_directions(false),verbose(false),use_angular_damping(false),zero_row_tolerance((T)1e-6),use_krylov_poststab(false),use_krylov_prestab(false),use_poststab_in_cg(true),
    use_prestab(true),use_poststab(true),check_stale(false),is_stale(true),last_read(-1)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> ARTICULATED_RIGID_BODY_BASE<TV>::
~ARTICULATED_RIGID_BODY_BASE()
{
    delete breadth_first_directed_graph;
    delete muscle_list;
    delete &joint_mesh;
}
//#####################################################################
// Function Set_Consistent_Child_Frame
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Set_Consistent_Child_Frame(JOINT_ID joint_id)
{
    JOINT<TV>& joint=*joint_mesh(joint_id);
    joint.Set_Joint_To_Child_Frame(Child(joint_id)->Frame().Inverse()*Parent(joint_id)->Frame()*joint.frame_pj);
}
//#####################################################################
// Function Initialize_Breadth_First_Directed_Graph
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Initialize_Breadth_First_Directed_Graph(const int root)
{
    delete breadth_first_directed_graph;breadth_first_directed_graph=new DIRECTED_GRAPH<int>(rigid_body_collection.rigid_body_particle.array_collection->Size());
    joint_mesh.undirected_graph.Breadth_First_Directed_Graph(root,*breadth_first_directed_graph);breadth_first_directed_graph->Generate_Levels();
}
//#####################################################################
// Function Update_With_Breadth_First_Directed_Graph
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Update_With_Breadth_First_Directed_Graph(const int root,const int node_to_propagate_from)
{
    if(!breadth_first_directed_graph || root!=breadth_first_directed_graph->Nodes_In_Level(1)(1)) Initialize_Breadth_First_Directed_Graph(root);
    if(!node_to_propagate_from) for(int k=1;k<=breadth_first_directed_graph->Number_Of_Levels();k++){
        int child=breadth_first_directed_graph->Nodes_In_Level(k)(1);
        Update_Child_From_Parents(child,breadth_first_directed_graph->Parents(child));}
    else{
        Update_Parent_Joint_States(node_to_propagate_from,breadth_first_directed_graph->Parents(node_to_propagate_from));
        Update_Child_From_Parents(node_to_propagate_from,breadth_first_directed_graph->Parents(node_to_propagate_from));
        ARRAY<bool> changed(breadth_first_directed_graph->Number_Of_Levels());
        int start_level=breadth_first_directed_graph->Level_Of_Node(node_to_propagate_from)+1;changed(start_level-1)=true;
        for(int k=start_level;k<=breadth_first_directed_graph->Number_Of_Levels();k++){
            int child=breadth_first_directed_graph->Nodes_In_Level(k)(1);ARRAY<int>& parents=breadth_first_directed_graph->Parents(child);
            bool parent_changed=false;for(int i=1;i<=parents.m;i++) if(changed(Value(parents(i)))){parent_changed=true;break;} // TODO: this looks broken!
            if(parent_changed) Update_Child_From_Parents(child,parents);}}
}
//#####################################################################
// Function Remove_All
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Remove_All()
{
    process_list.Remove_All();
    muscle_activations.Remove_All();
    angular_momenta_save.Remove_All();
    linear_velocities_save.Remove_All();
    muscle_list->Clean_Memory();
    joint_mesh.Remove_All();
    if(breadth_first_directed_graph) breadth_first_directed_graph->Reset();
}
//#####################################################################
// Function Parent_Id
//#####################################################################
template<class TV> int ARTICULATED_RIGID_BODY_BASE<TV>::
Parent_Id(JOINT_ID joint_id) const
{
    return joint_mesh.Parent_Id(joint_id);
}
//#####################################################################
// Function Child_Id
//#####################################################################
template<class TV> int ARTICULATED_RIGID_BODY_BASE<TV>::
Child_Id(JOINT_ID joint_id) const
{
    return joint_mesh.Child_Id(joint_id);
}
//#####################################################################
// Function Parent
//#####################################################################
template<class TV> RIGID_BODY<TV>* ARTICULATED_RIGID_BODY_BASE<TV>::
Parent(JOINT_ID joint_id)
{
    return &rigid_body_collection.Rigid_Body(Parent_Id(joint_id));
}
//#####################################################################
// Function Parent
//#####################################################################
template<class TV> const RIGID_BODY<TV>* ARTICULATED_RIGID_BODY_BASE<TV>::
Parent(JOINT_ID joint_id) const
{
    return &rigid_body_collection.Rigid_Body(Parent_Id(joint_id));
}
//#####################################################################
// Function Child
//#####################################################################
template<class TV> RIGID_BODY<TV>* ARTICULATED_RIGID_BODY_BASE<TV>::
Child(JOINT_ID joint_id)
{
    return &rigid_body_collection.Rigid_Body(Child_Id(joint_id));
}
//#####################################################################
// Function Child
//#####################################################################
template<class TV> const RIGID_BODY<TV>* ARTICULATED_RIGID_BODY_BASE<TV>::
Child(JOINT_ID joint_id) const
{
    return &rigid_body_collection.Rigid_Body(Child_Id(joint_id));
}
template<class TV> TV ARTICULATED_RIGID_BODY_BASE<TV>::
Joint_Location(JOINT_ID joint_id)
{
    return joint_mesh(joint_id)->Location(*Parent(joint_id),*Child(joint_id));
}
//#####################################################################
// Function Update_From_Edge_In_Breadth_First_Directed_Graph
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Update_From_Edge_In_Breadth_First_Directed_Graph(const int root,const JOINT_ID edge)
{
    if(!breadth_first_directed_graph || root!=breadth_first_directed_graph->Nodes_In_Level(1)(1)) Initialize_Breadth_First_Directed_Graph(root);
    int undirected_child_id=Child_Id(edge),undirected_parent_id=Parent_Id(edge);
    int directed_child_id=undirected_child_id;
    if(breadth_first_directed_graph->Level_Of_Node(undirected_child_id)<breadth_first_directed_graph->Level_Of_Node(undirected_parent_id)) directed_child_id=undirected_parent_id;
    Update_Child_From_Parents(directed_child_id,breadth_first_directed_graph->Parents(directed_child_id));
    Update_With_Breadth_First_Directed_Graph(root,directed_child_id);
}
//#####################################################################
// Function Update_Joint_Frames
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Update_Joint_Frames()
{
    for(int i=1;i<=joint_mesh.joints.m;i++){
        JOINT<TV>& joint=*joint_mesh.joints(i);
        joint.Set_Joint_Frame(joint.Compute_Current_Joint_Frame(*Parent(joint.id_number),*Child(joint.id_number)),false);}
}
//#####################################################################
// Function Store_Velocities_And_Momenta
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Store_Velocities_And_Momenta()
{
    ARRAY<int>& dynamic_rigid_body_particles=rigid_body_collection.dynamic_rigid_body_particles;
    linear_velocities_save.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size());angular_momenta_save.Resize(rigid_body_collection.rigid_body_particle.array_collection->Size());
    angular_momenta_save.Subset(dynamic_rigid_body_particles)=rigid_body_collection.rigid_body_particle.angular_momentum.Subset(dynamic_rigid_body_particles);
    for(int i=1;i<=dynamic_rigid_body_particles.m;i++){int p=dynamic_rigid_body_particles(i);
        linear_velocities_save(p)=rigid_body_collection.rigid_body_particle.V(p);}
}
//#####################################################################
// Function Restore_Velocities_And_Momenta
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Restore_Velocities_And_Momenta()
{
    ARRAY<int>& dynamic_rigid_body_particles=rigid_body_collection.dynamic_rigid_body_particles;
    rigid_body_collection.rigid_body_particle.angular_momentum.Subset(dynamic_rigid_body_particles)=angular_momenta_save.Subset(dynamic_rigid_body_particles);
    for(int i=1;i<=dynamic_rigid_body_particles.m;i++){int p=dynamic_rigid_body_particles(i);
        rigid_body_collection.rigid_body_particle.V(p)=linear_velocities_save(p);}
    rigid_body_collection.Update_Angular_Velocity(dynamic_rigid_body_particles);
}
//#####################################################################
// Function Update_Child_From_Parents
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Update_Child_From_Parents(const int child_id,const ARRAY<int>& parents)
{
    if(!parents.m) return;TV F_wc_sum_t;ROTATION<TV> current_hemisphere;
    ARRAY<ROTATION<TV> > rotations;
    for(int p=1;p<=parents.m;p++){
        JOINT_ID edge_id(0);
        for(int e=1;e<=joint_mesh.undirected_graph.Adjacent_Edges(child_id).m;e++){JOINT_ID edge=joint_mesh.undirected_graph.Adjacent_Edges(child_id)(e);
            if(Child_Id(edge)==parents(p) || Parent_Id(edge)==parents(p)){edge_id=edge;break;}}
        JOINT<TV>& joint=*joint_mesh(edge_id);
        FRAME<TV> F_wc=rigid_body_collection.Rigid_Body(parents(p)).Frame()*joint.F_pc(); // TODO: looks broken, should be line below?
        F_wc_sum_t+=F_wc.t;rotations.Append(F_wc.r);}
    rigid_body_collection.Rigid_Body(child_id).X()=F_wc_sum_t/(T)parents.m;
    rigid_body_collection.Rigid_Body(child_id).Rotation()=ROTATION<TV>::Average_Rotation(rotations);
}
//#####################################################################
// Function Update_Parent_Joint_States
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Update_Parent_Joint_States(const int child_id,const ARRAY<int>& parents)
{
    if(!parents.m) return;
    for(int p=1;p<=parents.m;p++) for(int e=1;e<=joint_mesh.undirected_graph.Adjacent_Edges(child_id).m;e++){JOINT_ID edge=joint_mesh.undirected_graph.Adjacent_Edges(child_id)(e);
        // J = (F_wj)^-1*F_wc*(F_jc)^-1 = F_jw*F_wc*F_cj
        if(Child_Id(edge)==parents(p)){
            JOINT<TV>& joint=*joint_mesh(edge);
            joint.Set_Joint_Frame(joint.Compute_Current_Joint_Frame(rigid_body_collection.Rigid_Body(child_id),rigid_body_collection.Rigid_Body(parents(p))));
            break;}
        else if(Parent_Id(edge)==parents(p)){
            JOINT<TV>& joint=*joint_mesh(edge);
            joint.Set_Joint_Frame(joint.Compute_Current_Joint_Frame(rigid_body_collection.Rigid_Body(parents(p)),rigid_body_collection.Rigid_Body(child_id)));
            break;}}
}
//####################################################################################
// Function Compute_Constraint_Correcting_Impulse
//####################################################################################
template<class TV> template<class T_CONSTRAINT_FUNCTION> void ARTICULATED_RIGID_BODY_BASE<TV>::
Compute_Constraint_Correcting_Impulse(T_CONSTRAINT_FUNCTION& constraint_error_function,typename T_CONSTRAINT_FUNCTION::T_IMPULSE& j) const
{
    typedef typename T_CONSTRAINT_FUNCTION::T_IMPULSE T_IMPULSE;
    typedef typename T_CONSTRAINT_FUNCTION::T_CONSTRAINT_ERROR T_CONSTRAINT_ERROR;

    T impulse_bound=0,iterative_tolerance_squared=sqr(iterative_tolerance),tau=(T).5*(sqrt((T)5)-1);

    j=T_IMPULSE();
    T_CONSTRAINT_ERROR f_of_j=constraint_error_function.F(j);
    T f_of_j_norm=constraint_error_function.Convergence_Norm_Squared(f_of_j);

    int iterations;
    for(iterations=1;iterations<=max_iterations;iterations++){
        if(f_of_j_norm<=iterative_tolerance_squared) break;

        T_IMPULSE delta_j=constraint_error_function.Jacobian(j).Solve_Linear_System(-f_of_j);
        if(iterations==1) impulse_bound=constraint_error_function.epsilon_scale*constraint_error_function.Magnitude(delta_j);

        // determine the step size z
        QUADRATIC<T> z_quadratic(constraint_error_function.Magnitude_Squared(delta_j),2*constraint_error_function.Inner_Product(j,delta_j),
            constraint_error_function.Magnitude_Squared(j)-sqr(impulse_bound));
        z_quadratic.Compute_Roots();
        if(z_quadratic.roots==0) break;
        T z_largest_root=z_quadratic.roots==1?z_quadratic.root1:z_quadratic.root2;
        if(z_largest_root<0) break;
        T z=min(z_largest_root,(T)1);

        // golden section search
        T_IMPULSE a=j,b=j+z*delta_j,j_1=a+(1-tau)*(b-a),j_2=a+tau*(b-a);
        T_CONSTRAINT_ERROR f_of_a=f_of_j,f_of_j_1=constraint_error_function.F(j_1),f_of_j_2=constraint_error_function.F(j_2),f_of_b=constraint_error_function.F(b);
        T f_of_a_norm=f_of_j_norm,f_of_j_1_norm=constraint_error_function.Convergence_Norm_Squared(f_of_j_1),
            f_of_j_2_norm=constraint_error_function.Convergence_Norm_Squared(f_of_j_2),f_of_b_norm=constraint_error_function.Convergence_Norm_Squared(f_of_b);

        int line_search_iterations=0;
        while(constraint_error_function.Magnitude(b-a)>line_search_interval_tolerance && line_search_iterations++<=max_line_search_iterations){
            // Throw out [a,j_1] if f(j_2) or f(b) is smallest
            if((f_of_j_2_norm<f_of_j_1_norm && f_of_j_2_norm<f_of_a_norm) || (f_of_b_norm<f_of_j_1_norm && f_of_b_norm<f_of_a_norm)){
                a=j_1;j_1=j_2;j_2=a+tau*(b-a);
                f_of_a=f_of_j_1;f_of_j_1=f_of_j_2;f_of_a_norm=f_of_j_1_norm;f_of_j_1_norm=f_of_j_2_norm;
                f_of_j_2=constraint_error_function.F(j_2);
                f_of_j_2_norm=constraint_error_function.Convergence_Norm_Squared(f_of_j_2);}
            else{
                b=j_2;j_2=j_1;j_1=a+(1-tau)*(b-a);
                f_of_b=f_of_j_2;f_of_j_2=f_of_j_1;f_of_b_norm=f_of_j_2_norm;f_of_j_2_norm=f_of_j_1_norm;
                f_of_j_1=constraint_error_function.F(j_1);
                f_of_j_1_norm=constraint_error_function.Convergence_Norm_Squared(f_of_j_1);}}

        // pick the min among a,j_1,j_2,b
        T_IMPULSE j_old=j;j=a;f_of_j=f_of_a;f_of_j_norm=f_of_a_norm;
        if(f_of_j_1_norm<f_of_j_norm){j=j_1;f_of_j=f_of_j_1;f_of_j_norm=f_of_j_1_norm;}
        if(f_of_j_2_norm<f_of_j_norm){j=j_2;f_of_j=f_of_j_2;f_of_j_norm=f_of_j_2_norm;}
        if(f_of_b_norm<f_of_j_norm){j=b;f_of_j=f_of_b;f_of_j_norm=f_of_b_norm;}
        if(j==j_old) break;}

    j/=constraint_error_function.dt;

    if(iterations>max_iterations) {std::stringstream ss;ss<<"PRESTABILIZATION: DID NOT CONVERGE (got "<<sqrt(f_of_j_norm)<<" vs. tolerance "<<iterative_tolerance<<")"<<std::endl;LOG::filecout(ss.str());}
    else if(f_of_j_norm>iterative_tolerance_squared && verbose)
        {std::stringstream ss;ss<<"PRESTABILIZATION: GOT IMPULSE "<<j<<" with tolerance "<<sqrt(f_of_j_norm)<<" vs. tolerance "<<iterative_tolerance<<std::endl;LOG::filecout(ss.str());}
}
//####################################################################################
// Function Apply_Prestabilization_To_Joint
//####################################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Apply_Prestabilization_To_Joint(const JOINT_ID joint_id,const T dt,const T epsilon_scale)
{
    if(!use_prestab || use_krylov_prestab) return;

    RIGID_BODY<TV> &parent=*Parent(joint_id),&child=*Child(joint_id);
    if(parent.Has_Infinite_Inertia() && child.Has_Infinite_Inertia()) return; // can't do anything

    TV location,jn;T_SPIN j_tau;JOINT<TV>& joint=*joint_mesh(joint_id);
    if(joint.Has_Prismatic_Constraint() && joint.Has_Angular_Constraint()){
        typename LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<TV>::T_IMPULSE j_combined;
        LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<TV> f_error(debug_cast<ARTICULATED_RIGID_BODY<TV>&>(*this),joint_id,dt,epsilon_scale,location);
        Compute_Constraint_Correcting_Impulse(f_error,j_combined);
        j_combined.Get_Subvector(1,jn);j_combined.Get_Subvector(TV::dimension+1,j_tau);
        RIGID_BODY<TV>::Apply_Impulse(parent,child,location,jn,j_tau);}
    else if(joint.Has_Prismatic_Constraint()){
        LINEAR_CONSTRAINT_FUNCTION<TV> f_error(debug_cast<ARTICULATED_RIGID_BODY<TV>&>(*this),joint_id,dt,epsilon_scale,location);
        Compute_Constraint_Correcting_Impulse(f_error,jn);
        RIGID_BODY<TV>::Apply_Impulse(parent,child,location,jn);}
    else if(joint.Has_Angular_Constraint()){
        ANGULAR_CONSTRAINT_FUNCTION<TV> f_error(debug_cast<ARTICULATED_RIGID_BODY<TV>&>(*this),joint_id,dt,epsilon_scale);
        if(TV::dimension==3) Compute_Constraint_Correcting_Impulse(f_error,j_tau);
        else if(TV::dimension==2) j_tau=-epsilon_scale*(f_error.Jacobian(T_SPIN()).Inverse()*f_error.F(T_SPIN())); // f(j_tau) is linear in 2D
        RIGID_BODY<TV>::Apply_Impulse(parent,child,TV(),TV(),j_tau);}
}
//####################################################################################
// Function Delta_Relative_Twist
//####################################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Delta_Relative_Twist(const JOINT_ID joint_id,const bool target_pd,TV& location,TWIST<TV>& delta_relative_twist)
{
    JOINT<TV>& joint=*joint_mesh(joint_id);
    RIGID_BODY<TV> &parent=*Parent(joint_id),&child=*Child(joint_id);
    location=joint.Location(parent,child);
    // TODO: need to fix translation to reflect correct translation although this should work for the basic examples
    parent.Update_Angular_Velocity();child.Update_Angular_Velocity();
    TWIST<TV> relative_twist=RIGID_BODY<TV>::Relative_Twist(parent,child,location),relative_twist_new=relative_twist;
    // IMPORTANT: joint function's relative velocity is child w.r.t. parent while rigid body uses parent w.r.t. child
    if(target_pd && joint.joint_function && joint.joint_function->active) relative_twist_new.angular=-joint.joint_function->desired_angular_velocity;
    joint.Constrain_Relative_Twist(parent.Frame(),relative_twist_new);
    delta_relative_twist=relative_twist_new-relative_twist;
}
//####################################################################################
// Function Apply_Poststabilization_To_Joint
//####################################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Apply_Poststabilization_To_Joint(const JOINT_ID joint_id,const bool target_pd)
{
    JOINT<TV>& joint=*joint_mesh(joint_id);
    if(!joint.Has_Prismatic_Constraint() && !joint.Has_Angular_Constraint()) return;

    RIGID_BODY<TV>& parent=*Parent(joint_id),&child=*Child(joint_id);
    if(parent.Has_Infinite_Inertia() && child.Has_Infinite_Inertia()) return;

    TV location;TWIST<TV> delta_relative_twist;
    Delta_Relative_Twist(joint_id,target_pd,location,delta_relative_twist);

    if(joint.Has_Prismatic_Constraint() && use_angular_damping && joint.angular_damping){
        if(joint.Has_Angular_Constraint() || (joint.joint_function && joint.joint_function->active)) PHYSBAM_NOT_IMPLEMENTED("NOT SUPPORTED FOR NOW... JUST DOING THIS FOR UNACTUATED NET");
        //TV delta_relative_angular_velocity=precomputed_desired_damped_angular_velocities(joint_index);
        delta_relative_twist.angular/=joint.angular_damping;
        MATRIX_MXN<T> angular_constraint_matrix(MATRIX<T,T_SPIN::dimension>::Identity_Matrix()),prismatic_constraint_matrix(MATRIX<T,T_SPIN::dimension>::Identity_Matrix());
        RIGID_BODY<TV>::Apply_Sticking_And_Angular_Sticking_Impulse(parent,child,location,delta_relative_twist,angular_constraint_matrix,prismatic_constraint_matrix);}
    else if(joint.Has_Prismatic_Constraint()){
        // get axes of constraint for rotation -- if PD is on, then all axes are assumed constrained
        MATRIX_MXN<T> angular_constraint_matrix(MATRIX<T,T_SPIN::dimension>::Identity_Matrix()),prismatic_constraint_matrix(MATRIX<T,TV::dimension>::Identity_Matrix());
        if(!(target_pd && joint.joint_function && joint.joint_function->active)){
            joint.Angular_Constraint_Matrix(parent.Frame(),angular_constraint_matrix);joint.Prismatic_Constraint_Matrix(parent.Frame(),prismatic_constraint_matrix);}
        RIGID_BODY<TV>::Apply_Sticking_And_Angular_Sticking_Impulse(parent,child,location,delta_relative_twist,angular_constraint_matrix,prismatic_constraint_matrix);}
    else{
        T_WORLD_SPACE_INERTIA_TENSOR I1_inv_plus_I2_inv;
        if(!parent.Has_Infinite_Inertia()) I1_inv_plus_I2_inv+=parent.World_Space_Inertia_Tensor_Inverse();
        if(!child.Has_Infinite_Inertia()) I1_inv_plus_I2_inv+=child.World_Space_Inertia_Tensor_Inverse();
        RIGID_BODY<TV>::Apply_Impulse(parent,child,TV(),TV(),I1_inv_plus_I2_inv.Inverse()*delta_relative_twist.angular);}
}
//####################################################################################
// Function Poststabilization_Projection
//####################################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Poststabilization_Projection(ARRAY_VIEW<TWIST<TV> > twist,const bool symmetric)
{
    if(!use_poststab_in_cg || !use_poststab) return;

    int iterations=poststabilization_projection_iterations;if(symmetric) iterations=(iterations+1)/2;
    for(int i=1;i<=iterations;i++) for(JOINT_ID j(1);j<=joint_mesh.Size();j++) if(joint_mesh.Is_Active(j)) Poststabilization_Projection_Joint(j,twist);
    if(symmetric) for(int i=1;i<=iterations;i++) for(JOINT_ID j=joint_mesh.Size();j>=JOINT_ID(1);j--) if(joint_mesh.Is_Active(j)) Poststabilization_Projection_Joint(j,twist);
}
//####################################################################################
// Function Poststabilization_Projection_Joint
//####################################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Poststabilization_Projection_Joint(const JOINT_ID joint_id,ARRAY_VIEW<TWIST<TV> > twist)
{
    int parent_index=Parent(joint_id)->particle_index,child_index=Child(joint_id)->particle_index;
    VECTOR<T,2*TWIST<TV>::dimension> t(twist(parent_index).Get_Vector(),twist(child_index).Get_Vector());
    t-=lambda_to_delta_v(joint_id)*(v_to_lambda(joint_id)*t);
    VECTOR<T,TWIST<TV>::dimension> tb;t.Get_Subvector(1,tb);twist(parent_index).Set_Vector(tb);t.Get_Subvector(TWIST<TV>::dimension+1,tb);twist(child_index).Set_Vector(tb);
}
//####################################################################################
// Function Initialize_Poststabilization_Projection
//####################################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Initialize_Poststabilization_Projection()
{
    if(!use_poststab_in_cg) return;
    v_to_lambda.Remove_All();lambda_to_delta_v.Remove_All();
    v_to_lambda.Resize(joint_mesh.Size());lambda_to_delta_v.Resize(joint_mesh.Size());
    for(int i=1;i<=joint_mesh.joints.m;i++){
        JOINT<TV>& joint=*joint_mesh.joints(i);JOINT_ID joint_id=joint.id_number;RIGID_BODY<TV>* rigid_bodies[]={Parent(joint_id),Child(joint_id)};
        if(rigid_bodies[0]->Has_Infinite_Inertia() && rigid_bodies[1]->Has_Infinite_Inertia()){
            lambda_to_delta_v(joint_id).Resize(2*(s+d),0);v_to_lambda(joint_id).Resize(0,2*(s+d));continue;}

        // get linear and angular constrained axes
        MATRIX_MXN<T> prismatic_constraints,angular_constraints;
        joint.Prismatic_Constraint_Matrix(rigid_bodies[0]->Frame(),prismatic_constraints);
        joint.Angular_Constraint_Matrix(rigid_bodies[0]->Frame(),angular_constraints);
        if(constrain_pd_directions && joint.joint_function && joint.joint_function->active){
            joint.joint_function->Prismatic_Constraint_Matrix(rigid_bodies[0]->Frame(),prismatic_constraints);
            joint.joint_function->Angular_Constraint_Matrix(rigid_bodies[0]->Frame(),angular_constraints);}
        int p=prismatic_constraints.Columns(),a=angular_constraints.Columns();

        MATRIX_MXN<T> R_D[2],M_inverse_R_D[2];
        TV location=joint.Location(*rigid_bodies[0],*rigid_bodies[1]);
        for(int j=0;j<=1;j++){
            MATRIX_MXN<T> r_star_p=prismatic_constraints.Cross_Product_Matrix_Times(location-rigid_bodies[j]->X());
            R_D[j].Resize(d+s,p+a);
            R_D[j].Set_Submatrix(1,1,prismatic_constraints);
            R_D[j].Set_Submatrix(d+1,1,r_star_p);
            R_D[j].Set_Submatrix(d+1,p+1,angular_constraints);
            M_inverse_R_D[j].Resize(d+s,p+a);
            if(!rigid_bodies[j]->Has_Infinite_Inertia()){
                M_inverse_R_D[j].Set_Submatrix(1,1,prismatic_constraints/rigid_bodies[j]->Mass());
                T_WORLD_SPACE_INERTIA_TENSOR inverse_inertia=rigid_bodies[j]->World_Space_Inertia_Tensor_Inverse();
                M_inverse_R_D[j].Set_Submatrix(d+1,1,inverse_inertia*r_star_p);
                M_inverse_R_D[j].Set_Submatrix(d+1,p+1,inverse_inertia*angular_constraints);}}

        lambda_to_delta_v(joint_id).Resize(2*(d+s),p+a);
        lambda_to_delta_v(joint_id).Set_Submatrix(1,1,M_inverse_R_D[0]);
        lambda_to_delta_v(joint_id).Set_Submatrix(d+s+1,1,-M_inverse_R_D[1]);

        v_to_lambda(joint_id).Resize(p+a,2*(d+s));
        MATRIX_MXN<T> DT_g_D=R_D[0].Transpose_Times(M_inverse_R_D[0])+R_D[1].Transpose_Times(M_inverse_R_D[1]);
        MATRIX_MXN<T> inverse_DT_g_D;DT_g_D.Cholesky_Inverse(inverse_DT_g_D);
        v_to_lambda(joint_id).Set_Submatrix(1,1,inverse_DT_g_D.Times_Transpose(R_D[0]));
        v_to_lambda(joint_id).Set_Submatrix(1,d+s+1,-inverse_DT_g_D.Times_Transpose(R_D[1]));}
}
//#####################################################################
// Function Generate_Process_List_Using_Contact_Graph
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Generate_Process_List_Using_Contact_Graph(const RIGID_BODY_CONTACT_GRAPH<TV>& contact_graph)
{
    process_list.Exact_Resize(contact_graph.Number_Of_Levels());
    for(int level=1;level<=process_list.m;level++) process_list(level).Remove_All();
    for(int joint=1;joint<=joint_mesh.joints.m;joint++){JOINT_ID joint_id=joint_mesh.joints(joint)->id_number;
        const int parent_id=Parent_Id(joint_id),child_id=Child_Id(joint_id);
        const int level_of_joint=max(contact_graph.directed_graph.Level_Of_Node(parent_id),contact_graph.directed_graph.Level_Of_Node(child_id));
        process_list(level_of_joint).Append(joint_id);}
}
//#####################################################################
// Function Apply_Poststabilization
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Apply_Poststabilization(bool test_system,bool print_matrix,const bool target_pd,const bool skip_global_post_stabilized_joints,const bool angular_damping_only)
{
    if(!use_poststab) return;

    if(use_krylov_poststab){
        PHYSBAM_ASSERT(!target_pd);
        Apply_Poststabilization_With_CG(0,false,test_system,print_matrix);
        return;}

    for(int k=1;k<=poststabilization_iterations;k++) for(int i=1;i<=joint_mesh.joints.m;i++){JOINT<TV>& joint=*joint_mesh.joints(i);
        if((angular_damping_only && !joint.angular_damping) || (skip_global_post_stabilized_joints && joint.global_post_stabilization)) continue;
        Apply_Poststabilization_To_Joint(joint.id_number,target_pd);}
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Read(const STREAM_TYPE stream_type,const std::string& directory,const int frame)
{
    int local_frame=frame;
    std::string arb_state_list_name=STRING_UTILITIES::string_sprintf("%s/common/arb_state_list",directory.c_str());
    if(FILE_UTILITIES::File_Exists(arb_state_list_name)){
        if(!frame_list){frame_list=new ARRAY<int>;FILE_UTILITIES::Read_From_File(stream_type,arb_state_list_name,*frame_list);}
        local_frame=(*frame_list)(frame_list->Binary_Search(frame));}
    if(last_read!=local_frame){
        std::istream* input=FILE_UTILITIES::Safe_Open_Input(STRING_UTILITIES::string_sprintf("%s/%d/arb_state",directory.c_str(),frame));
        TYPED_ISTREAM typed_input(*input,stream_type);
        joint_mesh.Read(typed_input,directory,frame);
        last_read=local_frame;delete input;}
    muscle_list->Read(stream_type,directory,frame);
    std::string muscle_activations_filename=STRING_UTILITIES::string_sprintf("%s/%d/muscle_activations",directory.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(muscle_activations_filename)) FILE_UTILITIES::Read_From_File(stream_type,muscle_activations_filename,muscle_activations);
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Write(const STREAM_TYPE stream_type,const std::string& directory,const int frame)
{
    if(joint_mesh.joints.m>0 && !(check_stale && !is_stale)){
        if(check_stale){
            if(!frame_list) frame_list=new ARRAY<int>;frame_list->Append(frame);
            FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/common/arb_state_list",directory.c_str()),*frame_list);
            is_stale=false;}
        std::ostream* output=FILE_UTILITIES::Safe_Open_Output(STRING_UTILITIES::string_sprintf("%s/%d/arb_state",directory.c_str(),frame));
        TYPED_OSTREAM typed_output(*output,stream_type);
        joint_mesh.Write(typed_output,directory,frame);
        delete output;}
    muscle_list->Write(stream_type,directory,frame);Output_Articulation_Points(stream_type,directory,frame);
    if(muscle_activations.m>0) FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/muscle_activations",directory.c_str(),frame),muscle_activations);
}
//#####################################################################
// Function Effective_Inertia_Inverse
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Effective_Inertia_Inverse(MATRIX<T,d+s>& inertia_inverse,JOINT_ID joint_id) const
{
    inertia_inverse=MATRIX<T,d+s>();
    const RIGID_BODY<TV>& parent_body=*Parent(joint_id);
    const RIGID_BODY<TV>& child_body=*Child(joint_id);
    T_WORLD_SPACE_INERTIA_TENSOR inertia_inverse_parent=parent_body.World_Space_Inertia_Tensor_Inverse();
    T_WORLD_SPACE_INERTIA_TENSOR inertia_inverse_child=child_body.World_Space_Inertia_Tensor_Inverse();
    TV location=(T).5*(parent_body.World_Space_Point(joint_mesh(joint_id)->frame_pj.t)+child_body.World_Space_Point(joint_mesh(joint_id)->frame_cj.t));
    MATRIX<T,s,d> cross_term=inertia_inverse_parent*MATRIX<T,s,d>::Cross_Product_Matrix(location-parent_body.X())+inertia_inverse_child*MATRIX<T,s,d>::Cross_Product_Matrix(location-child_body.X());
    inertia_inverse.Add_To_Submatrix(1,1,MATRIX<T,d>(parent_body.Impulse_Factor(location)+child_body.Impulse_Factor(location)));
    inertia_inverse.Add_To_Submatrix(1,d+1,cross_term.Transposed());
    inertia_inverse.Add_To_Submatrix(d+1,1,cross_term);
    inertia_inverse.Add_To_Submatrix(d+1,d+1,inertia_inverse_parent+inertia_inverse_child);
}
//#####################################################################
// Function Substitute_Joint_Parent_Body
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Substitute_Joint_Parent_Body(JOINT_ID joint_id,int new_parent,const FRAME<TV>& frame)
{
    JOINT<TV>& joint=*joint_mesh(joint_id);
    joint.Set_Joint_To_Parent_Frame(frame);
    joint_mesh.undirected_graph.Modify_Edge(new_parent,Child_Id(joint_id),joint_id);
    if(joint.joint_function) joint.joint_function->parent=Parent(joint_id);
    assert(Parent_Id(joint_id)==new_parent);
}
//#####################################################################
// Function Substitute_Joint_Parent_Body
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Substitute_Joint_Parent_Body(JOINT_ID joint_id,int new_parent)
{
    Substitute_Joint_Parent_Body(joint_id,new_parent,FRAME<TV>(rigid_body_collection.Rigid_Body(new_parent).Frame().Inverse_Times(Parent(joint_id)->Frame()*joint_mesh(joint_id)->F_pj())));
}
//#####################################################################
// Function Substitute_Joint_Child_Body
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Substitute_Joint_Child_Body(JOINT_ID joint_id,int new_child,const FRAME<TV>& frame)
{
    JOINT<TV>& joint=*joint_mesh(joint_id);
    joint.Set_Joint_To_Child_Frame(frame);
    joint_mesh.undirected_graph.Modify_Edge(Parent_Id(joint_id),new_child,joint_id);
    if(joint.joint_function) joint.joint_function->child=Child(joint_id);
    assert(Child_Id(joint_id)==new_child);
}
//#####################################################################
// Function Substitute_Joint_Child_Body
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Substitute_Joint_Child_Body(JOINT_ID joint_id,int new_child)
{
    Substitute_Joint_Child_Body(joint_id,new_child,FRAME<TV>(rigid_body_collection.Rigid_Body(new_child).Frame().Inverse_Times(Child(joint_id)->Frame()*joint_mesh(joint_id)->F_cj())));
}
//#####################################################################
// Function Resize_System_Vector
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Resize_System_Vector(ARTICULATED_VECTOR<TV>& v) const
{
    v.v.Resize(joint_mesh.Size());
}
//#####################################################################
// Function Apply_Poststabilization_With_CG
//#####################################################################
template<class TV> void ARTICULATED_RIGID_BODY_BASE<TV>::
Apply_Poststabilization_With_CG(T dt,bool correct_position,bool test_system,bool print_matrix)
{
    static int solve_id=0;
    solve_id++;
    ARTICULATED_SYSTEM<TV> system(debug_cast<ARTICULATED_RIGID_BODY<TV>&>(*this));
    system.Initialize();

    ARTICULATED_VECTOR<TV> rhs,x,q,s,r,k,z;
    Resize_System_Vector(rhs);
    Resize_System_Vector(x);
    Resize_System_Vector(q);
    Resize_System_Vector(s);
    Resize_System_Vector(r);
    Resize_System_Vector(k);
    Resize_System_Vector(z);

    if(test_system) system.Test_System(r,k,z);
    if(print_matrix){
        {std::stringstream ss;ss<<"arb solve id "<<solve_id<<std::endl;LOG::filecout(ss.str());}
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("M-%i.txt",solve_id).c_str()).Write("M",system,q,s);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("P-%i.txt",solve_id).c_str()).Write_Projection("P",system,q);}

    {ARTICULATED_SYSTEM<TV> system(debug_cast<ARTICULATED_RIGID_BODY<TV>&>(*this));
    system.break_loops=true;
    system.Initialize();
    if(print_matrix){
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("BM-%i.txt",solve_id).c_str()).Write("BM",system,q,s);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("BP-%i.txt",solve_id).c_str()).Write_Projection("BP",system,q);}}

    if(correct_position){for(JOINT_ID j(1);j<=joint_mesh.Size();j++) if(joint_mesh.Is_Active(j)) rhs.v(j)=-Joint_Error(j)/dt;}
    else{
        ARRAY<TWIST<TV> > twist;twist.Resize(rigid_body_collection.rigid_body_particle.V.Size());
        for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(rigid_body_collection.rigid_body_particle.V(i),rigid_body_collection.rigid_body_particle.angular_velocity(i));
        system.Scatter(twist,rhs.v);rhs*=-(T)1;
        for(int i=1;i<=twist.Size();i++){rigid_body_collection.rigid_body_particle.V(i)=twist(i).linear;rigid_body_collection.rigid_body_particle.angular_velocity(i)=twist(i).angular;}}
    if(print_matrix) OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("b-%i.txt",solve_id).c_str()).Write("b",rhs);

    CONJUGATE_RESIDUAL<T> cr;
    CONJUGATE_GRADIENT<T> cg;
    SYMMQMR<T> qm;
    KRYLOV_SOLVER<T>* solver=&cg;
    if(correct_position) solver=&cr;
//    if(!correct_position) system.internal_x=&x;
    solver->Solve(system,x,rhs,q,s,r,k,z,(T)1e-3,0,1000);
    if(print_matrix) OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("x-%i.txt",solve_id).c_str()).Write("x",x);

    ARRAYS_COMPUTATIONS::Fill(system.intermediate_twists,TWIST<TV>());
    system.Gather(x.v,system.intermediate_twists);
    system.Inverse_Mass(system.intermediate_twists);
    for(int i=1;i<=rigid_body_collection.rigid_body_particle.V.Size();i++){
        rigid_body_collection.rigid_body_particle.V(i)+=system.intermediate_twists(i).linear;
        rigid_body_collection.rigid_body_particle.angular_velocity(i)+=system.intermediate_twists(i).angular;}
    rigid_body_collection.Update_Angular_Momentum();
}
//#####################################################################
// Function Joint_Error
//#####################################################################
template<class TV> TWIST<TV> ARTICULATED_RIGID_BODY_BASE<TV>::
Joint_Error(JOINT_ID joint_id) const
{
    const JOINT<TV>& joint=*joint_mesh(joint_id);
    const RIGID_BODY<TV>& parent=*Parent(joint_id),&child=*Child(joint_id);
    FRAME<TV> wpj=parent.Frame()*joint.F_pj();
    FRAME<TV> wcj=child.Frame()*joint.F_cj();
    FRAME<TV> desired_frame=wpj.Inverse()*wcj;
    if(joint.Has_Prismatic_Constraint()) joint.Constrain_Prismatically(desired_frame.t);
    if(joint.Has_Angular_Constraint()){
        T_SPIN angle=desired_frame.r.Euler_Angles();
        joint.Constrain_Angles(angle);
        desired_frame.r=ROTATION<TV>::From_Euler_Angles(angle);}

    FRAME<TV> relative_world_frame=wpj*desired_frame*wcj.Inverse();
    return relative_world_frame.Delta();
}
//####################################################################################
#define INSTANTIATION_HELPER(T,d,CF) \
    template void ARTICULATED_RIGID_BODY_BASE<VECTOR<T,d> >::Compute_Constraint_Correcting_Impulse<CF<VECTOR<T,d> > >(CF<VECTOR<T,d> >&,CF<VECTOR<T,d> >::T_IMPULSE&) const;
#define INSTANTIATION_HELPER_D(T,d) \
    INSTANTIATION_HELPER(T,d,LINEAR_CONSTRAINT_FUNCTION) \
    INSTANTIATION_HELPER(T,d,ANGULAR_CONSTRAINT_FUNCTION) \
    INSTANTIATION_HELPER(T,d,LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION)
template class ARTICULATED_RIGID_BODY_BASE<VECTOR<float,2> >;
template class ARTICULATED_RIGID_BODY_BASE<VECTOR<float,3> >;
INSTANTIATION_HELPER_D(float,2);
INSTANTIATION_HELPER_D(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ARTICULATED_RIGID_BODY_BASE<VECTOR<double,2> >;
template class ARTICULATED_RIGID_BODY_BASE<VECTOR<double,3> >;
INSTANTIATION_HELPER_D(double,2);
INSTANTIATION_HELPER_D(double,3);
#endif
}
