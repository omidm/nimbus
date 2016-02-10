//#####################################################################
// Copyright 2004-2007, Kevin Der, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Optimization/LINEAR_PROGRAMMING.h>
#include <PhysBAM_Tools/Optimization/QUADRATIC_PROGRAMMING.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_FRAME.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR_ND.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/ATTACHMENT_POINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> ARTICULATED_RIGID_BODY<VECTOR<T,3> >::
ARTICULATED_RIGID_BODY(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
    :ARTICULATED_RIGID_BODY_BASE<TV>(rigid_body_collection_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> ARTICULATED_RIGID_BODY<VECTOR<T,3> >::
~ARTICULATED_RIGID_BODY()
{}
//####################################################################################
// Function Compute_Desired_PD_Velocity
//####################################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,3> >::
Compute_Desired_PD_Velocity(const T dt,const T time)
{
    for(JOINT_ID i(1);i<=joint_mesh.Size();i++)
        if(joint_mesh.Is_Active(i) && joint_mesh(i)->joint_function && joint_mesh(i)->joint_function->active && (!Parent(i)->Has_Infinite_Inertia() || !Child(i)->Has_Infinite_Inertia()))
            joint_mesh(i)->joint_function->Compute_Desired_PD_Velocity(dt,time);
}
//####################################################################################
// Function Compute_Position_Based_State
//####################################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,3> >::
Compute_Position_Based_State(const T dt,const T time)
{
    const ARRAY<int>& dynamic_rigid_body_particles=rigid_body_collection.dynamic_rigid_body_particles;
    Compute_Desired_PD_Velocity(dt,time);

    if(!global_post_stabilization) return;

    // build global poststabilization matrix
    ARRAY<ARRAY<PAIR<int,bool> >,int> joints_on_rigid_body(rigid_body_collection.rigid_body_particle.array_collection->Size()); // stores joint id and whether it's parent for joint.
    ARRAY<TV> joint_location(joint_mesh.joints.m);
    ARRAY<MATRIX_MXN<T> > joint_constraint_matrix(joint_mesh.joints.m),joint_muscle_control_matrix(joint_mesh.joints.m);
    joint_constrained_dimensions.Resize(joint_mesh.joints.m);ARRAYS_COMPUTATIONS::Fill(joint_constrained_dimensions,0);
    joint_muscle_control_dimensions.Resize(joint_mesh.joints.m);ARRAYS_COMPUTATIONS::Fill(joint_muscle_control_dimensions,0);
    joint_angular_constraint_matrix.Resize(joint_mesh.joints.m); // angular directions constrained by joint or PD
    joint_angular_muscle_control_matrix.Resize(joint_mesh.joints.m); // angular directions which are unconstrained (free for muscle to try to control)

    for(int i=1;i<=joint_mesh.joints.m;i++){
        JOINT<TV>& joint=*joint_mesh.joints(i);RIGID_BODY<TV>& parent=*Parent(joint.id_number),&child=*Child(joint.id_number);
        if(parent.Has_Infinite_Inertia() && child.Has_Infinite_Inertia()) continue;
        joints_on_rigid_body(Parent_Id(joint.id_number)).Append(PAIR<int,bool>(i,true));
        joints_on_rigid_body(Child_Id(joint.id_number)).Append(PAIR<int,bool>(i,false));

        if(!joint.global_post_stabilization){
            assert(!(joint.joint_function && joint.joint_function->active) || !joint.joint_function->muscle_control); // need global post stab for muscle control
            continue;}

        joint_location(i)=joint.Location(parent,child);

        // assumes full positional constraint and allows variable angular constraints
        // note: PD control counts as a constrained direction (it simply has a nonzero desired relative velocity for that direction)
        assert(joint.Has_Prismatic_Constraint()); // assume this for now
        joint.Angular_Constraint_Matrix(parent.Frame(),joint_angular_constraint_matrix(i),&joint_angular_muscle_control_matrix(i));
        if(joint.joint_function && joint.joint_function->active && (use_pd_actuators || (use_muscle_actuators && !joint.joint_function->muscle_control)))
            joint.joint_function->Angular_Constraint_Matrix(parent.Frame(),joint_angular_constraint_matrix(i),&joint_angular_muscle_control_matrix(i));

        // full constraint matrix (linear and angular)
        joint_constrained_dimensions(i)=3+joint_angular_constraint_matrix(i).Columns(); // includes PD dimensions
        joint_constraint_matrix(i)=MATRIX_MXN<T>(6,joint_constrained_dimensions(i));
        joint_constraint_matrix(i).Set_Submatrix(1,1,MATRIX<T,3>::Identity_Matrix());
        joint_constraint_matrix(i).Set_Submatrix(4,4,joint_angular_constraint_matrix(i));

        // full muscle control matrix
        if(use_muscle_actuators){
            joint_muscle_control_dimensions(i)=(joint.joint_function && joint.joint_function->active && joint.joint_function->muscle_control)?3-joint_angular_constraint_matrix(i).Columns():0; // muscle control dimensions
            joint_muscle_control_matrix(i)=MATRIX_MXN<T>(6,joint_muscle_control_dimensions(i));
            joint_muscle_control_matrix(i).Set_Submatrix(4,1,joint_angular_muscle_control_matrix(i));}}

    joint_offset_in_post_stabilization_matrix.Resize(joint_mesh.joints.m);
    for(int i=1;i<=joint_mesh.joints.m;i++) joint_offset_in_post_stabilization_matrix(i)=(i==1)?1:joint_offset_in_post_stabilization_matrix(i-1)+joint_constrained_dimensions(i-1);
    int total_constrained_dimensions=ARRAYS_COMPUTATIONS::Sum(joint_constrained_dimensions);
    global_post_stabilization_matrix_11=MATRIX_MXN<T>(total_constrained_dimensions,total_constrained_dimensions);
    {std::stringstream ss;ss<<"*** Total constrained dimensions "<<total_constrained_dimensions<<std::endl;LOG::filecout(ss.str());}

    if(use_muscle_actuators){
        joint_offset_in_muscle_control_matrix.Resize(joint_mesh.joints.m);
        for(int i=1;i<=joint_mesh.joints.m;i++) joint_offset_in_muscle_control_matrix(i)=(i==1)?1:joint_offset_in_muscle_control_matrix(i-1)+joint_muscle_control_dimensions(i-1);
        int total_muscle_control_dimensions=ARRAYS_COMPUTATIONS::Sum(joint_muscle_control_dimensions);assert(total_muscle_control_dimensions);
        global_post_stabilization_matrix_12=MATRIX_MXN<T>(total_constrained_dimensions,muscle_list->muscles.m);
        global_post_stabilization_matrix_21=MATRIX_MXN<T>(total_muscle_control_dimensions,total_constrained_dimensions);
        global_post_stabilization_matrix_22=MATRIX_MXN<T>(total_muscle_control_dimensions,muscle_list->muscles.m);
        {std::stringstream ss;ss<<"Muscle control: controlling "<<total_muscle_control_dimensions<<" dof with "<<muscle_list->muscles.m<<" muscles"<<std::endl;LOG::filecout(ss.str());}}

    for(int i=1;i<=dynamic_rigid_body_particles.m;i++){int p=dynamic_rigid_body_particles(i);
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(p);SYMMETRIC_MATRIX<T,3> I_inverse=rigid_body.World_Space_Inertia_Tensor_Inverse();
        int id=rigid_body.particle_index;
        for(int i=1;i<=joints_on_rigid_body(id).m;i++) for(int j=i;j<=joints_on_rigid_body(id).m;j++){
            int joint_index_1=joints_on_rigid_body(id)(i).x,joint_index_2=joints_on_rigid_body(id)(j).x;

            // joints that won't take part in global post stabilization are skipped
            // we currently assume that a connected component must either be all in or all out of global post stabilization but not mixed
            if(!joint_mesh.joints(joint_index_1)->global_post_stabilization || !joint_mesh.joints(joint_index_2)->global_post_stabilization){
                //assert(!joint_mesh.joints(joint_index_1)->global_post_stabilization && !joint_mesh.joints(joint_index_2)->global_post_stabilization);
                continue;}

            // compute 6x6 matrix assuming impulse applied at r2 and measured at r1
            TV r1=joint_location(joint_index_1)-rigid_body.X(),r2=joint_location(joint_index_2)-rigid_body.X();
            MATRIX<T,3> r_cross_1=MATRIX<T,3>::Cross_Product_Matrix(r1),r_cross_2=MATRIX<T,3>::Cross_Product_Matrix(r2);
            MATRIX<T,3> c21=I_inverse.Times_Cross_Product_Matrix(r2),c12=I_inverse.Cross_Product_Matrix_Times(-r1);
            MATRIX<T,3> c11=c21.Cross_Product_Matrix_Times(-r1)+1/rigid_body.Mass();
            SYMMETRIC_MATRIX<T,3> c22=I_inverse;
            MATRIX_MXN<T> C(6,6);C.Add_To_Submatrix(1,1,c11);C.Add_To_Submatrix(1,4,c12);C.Add_To_Submatrix(4,1,c21);C.Add_To_Submatrix(4,4,c22);
            // important: sign of block is negative if the joints at which impulse is applied and measured are opposite types (parent/child)
            if(joints_on_rigid_body(id)(i).y!=joints_on_rigid_body(id)(j).y) C*=-1;

            // compute projected matrix
            MATRIX_MXN<T> C_times_joint_constraint_matrix_2=C*joint_constraint_matrix(joint_index_2);
            MATRIX_MXN<T> projected_constrained_C=joint_constraint_matrix(joint_index_1).Transpose_Times(C_times_joint_constraint_matrix_2);

            // build global_post_stabilization_matrix_11 which relates constraint impulses (post stabilization impulses) to constrained velocity directions
            int istart=joint_offset_in_post_stabilization_matrix(joint_index_1),jstart=joint_offset_in_post_stabilization_matrix(joint_index_2);
            global_post_stabilization_matrix_11.Add_To_Submatrix(istart,jstart,projected_constrained_C);
            if(joint_index_1!=joint_index_2) global_post_stabilization_matrix_11.Add_To_Submatrix(jstart,istart,projected_constrained_C.Transposed());

            if(use_muscle_actuators){ // build global_post_stabilization_matrix_21 which relates constraint impulses to muscle controlled velocity directions
                global_post_stabilization_matrix_21.Add_To_Submatrix(joint_offset_in_muscle_control_matrix(joint_index_1),jstart,
                    joint_muscle_control_matrix(joint_index_1).Transpose_Times(C_times_joint_constraint_matrix_2));
                if(joint_index_1!=joint_index_2)
                    global_post_stabilization_matrix_21.Add_To_Submatrix(joint_offset_in_muscle_control_matrix(joint_index_2),istart,
                        joint_muscle_control_matrix(joint_index_2).Transpose_Times(C.Transpose_Times(joint_constraint_matrix(joint_index_1))));}}}

    if(use_muscle_actuators){
        muscle_list->Initialize_Muscle_Attachments_On_Rigid_Body();
        for(int i=1;i<=joint_mesh.joints.m;i++) if(joint_mesh.joints(i)->global_post_stabilization){
            int parent_id=Parent_Id(joint_mesh.joints(i)->id_number),child_id=Child_Id(joint_mesh.joints(i)->id_number);
            // joint should not have any muscle control dimensions if both bodies have infinite inertia
            assert(!rigid_body_collection.Rigid_Body(parent_id).Has_Infinite_Inertia() || !rigid_body_collection.Rigid_Body(child_id).Has_Infinite_Inertia());
            int body_indices[]={parent_id,child_id};
            for(int t=0;t<2;t++){ // loop over parent and child
                int body_id=body_indices[t];RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(body_id);
                if(!body.Has_Infinite_Inertia() && muscle_list->muscle_attachments_on_rigid_body(body_id).m){
                    TV r_joint=joint_location(i)-body.X();
                    SYMMETRIC_MATRIX<T,3> I_inverse=body.World_Space_Inertia_Tensor_Inverse();
                    for(int j=1;j<=muscle_list->muscle_attachments_on_rigid_body(body_id).m;j++){
                        TRIPLE<int,ATTACHMENT_POINT<TV>*,ATTACHMENT_POINT<TV>*>& muscle_attachments=muscle_list->muscle_attachments_on_rigid_body(body_id)(j);
                        int muscle_index=muscle_attachments.x;
                        // compute 6x1 vector which is the change in linear and angular velocity in response to unit muscle impulse along its line of action
                        TV r_attach=muscle_attachments.y->Embedded_Position()-body.X();
                        TV direction=(muscle_attachments.z->Embedded_Position()-muscle_attachments.y->Embedded_Position()).Normalized();
                        TV c21_along_direction=I_inverse*TV::Cross_Product(r_attach,direction);
                        TV c11_along_direction=TV::Cross_Product(c21_along_direction,r_joint)+direction/body.Mass();
                        VECTOR_ND<T> C_along_direction(6);C_along_direction.Set_Subvector(1,c11_along_direction);C_along_direction.Set_Subvector(4,c21_along_direction);
                        if(t==1) C_along_direction=-C_along_direction; // child gets negative sign because the relative velocity measures parent vel - child vel
                        global_post_stabilization_matrix_12.Add_To_Submatrix(joint_offset_in_post_stabilization_matrix(i),muscle_index,
                            joint_constraint_matrix(i).Transpose_Times(C_along_direction));
                        global_post_stabilization_matrix_22.Add_To_Submatrix(joint_offset_in_muscle_control_matrix(i),muscle_index,
                            joint_muscle_control_matrix(i).Transpose_Times(C_along_direction));}}}}}

    for(int i=1;i<=muscle_list->muscles.m;i++) muscle_list->muscles(i)->Update_Segments();
}
//####################################################################################
// Function Solve_For_Muscle_Control
//####################################################################################
// Uses least squares
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,3> >::
Solve_For_Muscle_Control(MATRIX_MXN<T>& A,const VECTOR_ND<T>& b,VECTOR_ND<T>& x,const T dt)
{
    T min_activation_penalty=1; // weight used in minimization of sum of activations squared
    int number_of_muscles=A.n;

    if(enforce_nonnegative_activations){
        T tolerance=(T)1e-8,step_tolerance=(T)1e-6;

        if(verbose){
            {std::stringstream ss;ss<<"######################### NEW SOLVE ####################################"<<std::endl;LOG::filecout(ss.str());}
            {std::stringstream ss;ss<<"*** MUSCLE SYSTEM BEFORE TRANSFORM: "<<std::endl;LOG::filecout(ss.str());}
            {std::stringstream ss;ss<<A<<std::endl<<b<<std::endl;LOG::filecout(ss.str());}}

        // get initial factorization
        MATRIX_MXN<T>& A_after_QR(A);VECTOR_ND<T> b_after_QR(b);VECTOR_ND<int> permute(A.n);

        // set up muscle bounds
        // NOTE: we normalize x to be between 0 to 1 so our optimization is hopefully better conditioned
        // (this assumes x_min=0).  We multiply x back with the peak force before returning below.
        ARRAY<PAIR<bool,T> > x_min(number_of_muscles),x_max(number_of_muscles);
        for(int j=1;j<=number_of_muscles;j++){
            MUSCLE<TV>* muscle=muscle_list->muscles(j);
            // TODO: include tendon length 
            // TODO: is this min necessary?
            T max_force=muscle->Force(1);T actuation_min=min(muscle->Passive_Force(muscle->Total_Length())/max_force,(T)1);
            x_min(j)=PAIR<bool,T>(true,actuation_min);x_max(j)=PAIR<bool,T>(true,1);
            for(int i=1;i<=A.m;i++) A_after_QR(i,j)*=max_force*dt;} // multiply by the maximum value x could have had with no normalization

        A_after_QR.In_Place_Robust_Householder_QR_Solve(b_after_QR,permute);
        int equations_to_keep=A_after_QR.Number_Of_Nonzero_Rows(zero_row_tolerance);

        if(verbose){
            {std::stringstream ss;ss<<"*** MUSCLE SYSTEM AFTER TRANSFORM: "<<std::endl;LOG::filecout(ss.str());}
            {std::stringstream ss;ss<<A_after_QR<<std::endl<<b_after_QR<<std::endl<<permute<<std::endl;LOG::filecout(ss.str());}
            {std::stringstream ss;ss<<"equations_to_keep = "<<equations_to_keep<<std::endl;LOG::filecout(ss.str());}}

        MATRIX_MXN<T> B(equations_to_keep);MATRIX_MXN<T> S(B.n,0),N(equations_to_keep,A.n-equations_to_keep),epsilon_hat(A.m-equations_to_keep,A.n);
        A_after_QR.Get_Submatrix(1,1,B);A_after_QR.Get_Submatrix(1,equations_to_keep+1,N);A_after_QR.Get_Submatrix(equations_to_keep+1,1,epsilon_hat);
        MATRIX_MXN<T> epsilon_hat_unpermuted=epsilon_hat.Unpermute_Columns(permute); // store as unpermuted quantities

        MATRIX_MXN<T> D(number_of_muscles); // this is actually the square root of the weight
        T sqrt_min_activation_penalty=sqrt(min_activation_penalty);for(int i=1;i<=number_of_muscles;i++) D(i,i)=sqrt_min_activation_penalty;

        VECTOR_ND<int> permute_B(B.n),permute_S,permute_N(N.n);permute.Get_Subvector(1,permute_B);permute.Get_Subvector(B.n+1,permute_N);
        VECTOR_ND<T> b(B.m),f_hat(epsilon_hat.m),x_B(B.n),x_S,x_N(N.n);
        b_after_QR.Get_Subvector(1,b);b_after_QR.Get_Subvector(equations_to_keep+1,f_hat);

        // Solve for x_B given previous guess for remainder of vector.  If x_B is within bounds then we have a good initial guess which we can pass to QP.
        bool last_values_are_feasible=false;
        if(last_muscle_actuations.n && last_muscle_actuations.n==number_of_muscles){ // TODO: this assumes muscle list doesn't change dynamically!
            for(int i=1;i<=x_N.n;i++){int p=permute_N(i);x_N(i)=last_muscle_actuations(p);if(x_min(p).x) x_N(i)=max(x_N(i),x_min(p).y);if(x_max(p).x) x_N(i)=min(x_N(i),x_max(p).y);}
            x_B=B.Upper_Triangular_Solve(b-N*x_N);
            last_values_are_feasible=true;
            for(int i=1;i<=x_B.n;i++){
                if((x_min(permute_B(i)).x && x_B(i)<x_min(permute_B(i)).y) || (x_max(permute_B(i)).x && x_B(i)>x_max(permute_B(i)).y)){last_values_are_feasible=false;break;}}
            if(last_values_are_feasible){
                {std::stringstream ss;ss<<"Using last activations as feasible initial guess"<<std::endl;LOG::filecout(ss.str());}
                if(verbose) {std::stringstream ss;ss<<"Before:\nN:\n"<<N<<"\npermute_N:\n"<<permute_N<<"\nx_N:\n"<<x_N<<std::endl;LOG::filecout(ss.str());}
                ARRAY<bool> move_to_S(N.n);
                for(int i=1;i<=x_N.n;i++) if((!x_min(permute_N(i)).x || x_N(i)>x_min(permute_N(i)).y) && (!x_max(permute_N(i)).x || x_N(i)<x_max(permute_N(i)).y)) move_to_S(i)=true;
                S.Resize(B.n,move_to_S.Number_True());permute_S.Resize(S.n);x_S.Resize(S.n);
                MATRIX_MXN<T> new_N(B.n,move_to_S.Number_False());VECTOR_ND<int> new_permute_N(new_N.n);VECTOR_ND<T> new_x_N(new_N.n);
                int Sj=0,Nj=0;
                for(int j=1;j<=move_to_S.m;j++){
                    if(move_to_S(j)){Sj++;permute_S(Sj)=permute_N(j);x_S(Sj)=x_N(j);for(int i=1;i<=S.m;i++) S(i,Sj)=N(i,j);}
                    else{Nj++;new_permute_N(Nj)=permute_N(j);new_x_N(Nj)=x_N(j);for(int i=1;i<=N.m;i++) new_N(i,Nj)=N(i,j);}}
                N=new_N;permute_N=new_permute_N;x_N=new_x_N;
                if(verbose) {std::stringstream ss;ss<<"After:\nN:\n"<<N<<"\nS:\n"<<S<<"\npermute_N:\n"<<permute_N<<"\npermute_S:\n"<<permute_S<<"\nx_N:\n"<<x_N<<"\nx_S:\n"<<x_S<<std::endl;LOG::filecout(ss.str());}}}

        if(!last_values_are_feasible) LINEAR_PROGRAMMING<T>::Find_Feasible_Solution(B,N,x_B,b,x_N,permute_B,permute_N,x,x_min,x_max,tolerance,step_tolerance,verbose);
        else{
            x.Set_Subvector(1,x_B);x.Set_Subvector(B.n+1,x_S);x.Set_Subvector(B.n+S.n+1,x_N);
            permute.Set_Subvector(1,permute_B);permute.Set_Subvector(B.n+1,permute_S);permute.Set_Subvector(B.n+S.n+1,permute_N);
            x=x.Unpermute(permute);}

        VECTOR_ND<T> x_lp=x;

        QUADRATIC_PROGRAMMING<T>::Find_Optimal_Solution(B,S,N,x_B,x_S,b,x_N,permute_B,permute_S,permute_N,D,epsilon_hat_unpermuted,f_hat,x,x_min,x_max,tolerance,step_tolerance,verbose);

        last_muscle_actuations=x;
        for(int i=1;i<=x.n;i++){T max_force=muscle_list->muscles(i)->Force(1);x_lp(i)*=max_force*dt;x(i)*=max_force*dt;}

        if(verbose){
            if(last_values_are_feasible) {std::stringstream ss;ss<<"===> Initial feasible result x (using last activations):"<<std::endl<<x_lp<<std::endl;LOG::filecout(ss.str());}
            else {std::stringstream ss;ss<<"===> LP result x:"<<std::endl<<x_lp<<std::endl;LOG::filecout(ss.str());}
            {std::stringstream ss;ss<<"===> QP result x:"<<std::endl<<x<<std::endl;LOG::filecout(ss.str());}}
        {std::stringstream ss;ss<<"===> LP: muscle activation sum is "<<(D*x_lp).Magnitude()<<std::endl;LOG::filecout(ss.str());}
        {std::stringstream ss;ss<<"===> QP: muscle activation sum is "<<(D*x).Magnitude()<<std::endl;LOG::filecout(ss.str());}}
    else{
        bool debug=(getenv("DEBUG_QR")!=0);
        T negative_penalty_factor=(T)10,min_negative_penalty=(T).1,max_negative_penalty=(T)1e6;
        int iterations=enforce_nonnegative_activations?activation_optimization_iterations:1;

        // OLD METHOD
        if(debug || getenv("OLD_LS")){
            VECTOR_ND<T> negative_penalty(number_of_muscles);
            MATRIX_MXN<T> A_transpose_A=A.Normal_Equations_Matrix();
            int iterations=enforce_nonnegative_activations?activation_optimization_iterations:1;
            for(int optimization_iteration=1;optimization_iteration<=iterations;optimization_iteration++){

                //if(optimization_iteration==2) for(int i=1;i<=number_of_muscles;i++){
                //    if(x(i)<=0) negative_penalty(i)=max_negative_penalty;}
                if(optimization_iteration>1) for(int i=1;i<=number_of_muscles;i++){
                    if(x(i)>0) negative_penalty(i)=0;
                    else negative_penalty(i)=clamp(negative_penalty_factor*negative_penalty(i),min_negative_penalty,max_negative_penalty);}

                // assume all relative weights are 1
                for(int i=1;i<=number_of_muscles;i++) A_transpose_A(i,i)+=min_activation_penalty+negative_penalty(i);
                VECTOR_ND<T> A_transpose_b=A.Transpose_Times(b);
                x=A_transpose_A.Cholesky_Solve(A_transpose_b);}}
        if(!getenv("OLD_LS")){
            VECTOR_ND<T> negative_penalty(number_of_muscles);
            MATRIX_MXN<T> transformed_A(A);VECTOR_ND<T> transformed_b(b);VECTOR_ND<int> permute(A.n);
            transformed_A.In_Place_Robust_Householder_QR_Solve(transformed_b,permute);

            VECTOR_ND<T> old_method_x;
            if(debug){
                old_method_x=x;
                {std::stringstream ss;ss<<"*** MUSCLE SYSTEM BEFORE TRANSFORM: "<<std::endl;LOG::filecout(ss.str());}
                {std::stringstream ss;ss<<A<<std::endl;LOG::filecout(ss.str());}
                {std::stringstream ss;ss<<b<<std::endl;LOG::filecout(ss.str());}
                {std::stringstream ss;ss<<"*** MUSCLE SYSTEM AFTER TRANSFORM: "<<std::endl;LOG::filecout(ss.str());}
                {std::stringstream ss;ss<<transformed_A<<std::endl;LOG::filecout(ss.str());}
                {std::stringstream ss;ss<<transformed_b<<std::endl;LOG::filecout(ss.str());}
                {std::stringstream ss;ss<<permute<<std::endl;LOG::filecout(ss.str());}}

            int equations_to_keep=transformed_A.Number_Of_Nonzero_Rows((T)1e-6);
            MATRIX_MXN<T> R(equations_to_keep,equations_to_keep);MATRIX_MXN<T> U(equations_to_keep,A.n-equations_to_keep),G(A.m-equations_to_keep,A.n-equations_to_keep);
            transformed_A.Get_Submatrix(1,1,R);transformed_A.Get_Submatrix(1,equations_to_keep+1,U);transformed_A.Get_Submatrix(equations_to_keep+1,equations_to_keep+1,G);
            VECTOR_ND<T> b1(R.n),b2(G.m),x1(R.n),x2(U.n);
            transformed_b.Get_Subvector(1,b1);transformed_b.Get_Subvector(b1.n+1,b2);

            // Z is
            // ( -R^-1 * U )
            // (     I     )
            //
            // Z_rhs is
            // ( -R^-1 * b1 )
            // (     0      )
            MATRIX_MXN<T> negative_R_inverse_U(U.m,U.n);VECTOR_ND<T> column(U.m),u(U.m);
            for(int j=1;j<=U.n;j++){
                for(int i=1;i<=U.m;i++) u(i)=U(i,j);
                column=R.Upper_Triangular_Solve(u);
                for(int i=1;i<=U.m;i++) negative_R_inverse_U(i,j)=-column(i);}
            VECTOR_ND<T> negative_R_inverse_b1(-R.Upper_Triangular_Solve(b1));
            MATRIX_MXN<T> Z(A.n,U.n);Z.Add_To_Submatrix(1,1,negative_R_inverse_U);Z.Add_To_Submatrix(equations_to_keep+1,1,MATRIX_MXN<T>::Identity_Matrix(U.n));
            VECTOR_ND<T> Z_rhs(A.n);Z_rhs.Set_Subvector(1,negative_R_inverse_b1);

            MATRIX_MXN<T> sqrt_D(number_of_muscles,number_of_muscles);
            MATRIX_MXN<T> ls_A(G.m+A.n,G.n);VECTOR_ND<T> ls_b(G.m+A.n);
            MATRIX_MXN<T> A_transpose_A=A.Normal_Equations_Matrix();
            for(int optimization_iteration=1;optimization_iteration<=iterations;optimization_iteration++){

                //if(optimization_iteration==2) for(int i=1;i<=number_of_muscles;i++){
                //    if(x(i)<=0) negative_penalty(i)=max_negative_penalty;}
                if(optimization_iteration>1) for(int i=1;i<=number_of_muscles;i++){
                    if(x(i)>0) negative_penalty(i)=0;
                    else negative_penalty(i)=clamp(negative_penalty_factor*negative_penalty(i),min_negative_penalty,max_negative_penalty);}

                for(int i=1;i<=number_of_muscles;i++) sqrt_D(i,i)=sqrt(min_activation_penalty+negative_penalty(i));

                // least squares matrix is
                // ( G         )
                // ( D^1/2 * Z )
                ls_A.Set_Zero_Matrix();ls_A.Add_To_Submatrix(1,1,G);ls_A.Add_To_Submatrix(G.m+1,1,sqrt_D*Z);
                ls_b.Set_Subvector(1,b2);ls_b.Set_Subvector(G.m+1,sqrt_D*Z_rhs);

                x2=ls_A.Normal_Equations_Solve(ls_b);
                x1=negative_R_inverse_U*x2-negative_R_inverse_b1;
                x.Set_Subvector(1,x1);x.Set_Subvector(x1.n+1,x2);}

            x=x.Unpermute(permute);

            if(debug){
                {std::stringstream ss;ss<<"==> R "<<std::endl<<R<<std::endl;LOG::filecout(ss.str());}
                {std::stringstream ss;ss<<"==> U "<<std::endl<<U<<std::endl;LOG::filecout(ss.str());}
                {std::stringstream ss;ss<<"==> G "<<std::endl<<G<<std::endl;LOG::filecout(ss.str());}
                {std::stringstream ss;ss<<"==> b1 "<<std::endl<<b1<<std::endl;LOG::filecout(ss.str());}
                {std::stringstream ss;ss<<"==> b2 "<<std::endl<<b2<<std::endl;LOG::filecout(ss.str());}
                {std::stringstream ss;ss<<"==> Z "<<std::endl<<Z<<std::endl;LOG::filecout(ss.str());}
                {std::stringstream ss;ss<<"==> Z_rhs "<<std::endl<<Z_rhs<<std::endl;LOG::filecout(ss.str());}
                {std::stringstream ss;ss<<"OLD METHOD X:"<<std::endl<<old_method_x<<std::endl;LOG::filecout(ss.str());}
                {std::stringstream ss;ss<<"RESIDUAL: "<<(A*old_method_x-b).Magnitude()<<std::endl;LOG::filecout(ss.str());}
                {std::stringstream ss;ss<<"NEW METHOD X:"<<std::endl<<x<<std::endl;LOG::filecout(ss.str());}
                {std::stringstream ss;ss<<"RESIDUAL: "<<(A*x-b).Magnitude()<<std::endl;LOG::filecout(ss.str());}}}}

    if(clamp_negative_activations) for(int i=1;i<=x.n;i++) if(x(i)<0) x(i)=0;
}
//####################################################################################
// Function Solve_Minimum_Norm_Solution_For_Linear_Constraints
//####################################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,3> >::
Solve_Minimum_Norm_Solution_For_Linear_Constraints(MATRIX_MXN<T>& A,const VECTOR_ND<T>& b,VECTOR_ND<T>& x,const T zero_row_tolerance,const bool verbose)
{
    assert(A.m==A.n);
    if(verbose){
        {std::stringstream ss;ss<<"*** SYSTEM BEFORE TRANSFORM: "<<std::endl;LOG::filecout(ss.str());}
        {std::stringstream ss;ss<<A<<std::endl<<b<<std::endl;LOG::filecout(ss.str());}}

    // get initial factorization
    MATRIX_MXN<T> A_after_QR(A);VECTOR_ND<T> b_after_QR(b);VECTOR_ND<int> permute(A.n);
    A_after_QR.In_Place_Robust_Householder_QR_Solve(b_after_QR,permute);
    int equations_to_keep=A_after_QR.Number_Of_Nonzero_Rows(zero_row_tolerance);

    if(verbose){
        {std::stringstream ss;ss<<"*** SYSTEM AFTER TRANSFORM: "<<std::endl;LOG::filecout(ss.str());}
        {std::stringstream ss;ss<<A_after_QR<<std::endl<<b_after_QR<<std::endl<<permute<<std::endl;LOG::filecout(ss.str());}}
    {std::stringstream ss;ss<<"Global poststabilization: got "<<equations_to_keep<<" dof"<<std::endl;LOG::filecout(ss.str());}

    MATRIX_MXN<T> B(equations_to_keep);A_after_QR.Get_Submatrix(1,1,B);

    if(equations_to_keep==A.n) x=B.Upper_Triangular_Solve(b_after_QR).Unpermute(permute); // unique solution
    else{
        // extract submatrices and subvectors
        int epsilon_size=A.n-equations_to_keep;
        MATRIX_MXN<T> epsilon(epsilon_size);A_after_QR.Get_Submatrix(equations_to_keep+1,equations_to_keep+1,epsilon);
        MATRIX_MXN<T> S(equations_to_keep,epsilon_size);A_after_QR.Get_Submatrix(1,equations_to_keep+1,S);
        VECTOR_ND<T> b_B(equations_to_keep),b_epsilon(epsilon_size);b_after_QR.Get_Subvector(1,b_B);b_after_QR.Get_Subvector(equations_to_keep+1,b_epsilon);
        // construct the matrix
        MATRIX_MXN<T> B_inverse_S=B.Upper_Triangular_Solve(S);
        MATRIX_MXN<T> C(epsilon_size);MATRIX_MXN<T>::Add_Transpose_Times(B_inverse_S,B_inverse_S,C);MATRIX_MXN<T>::Add_Transpose_Times(epsilon,epsilon,C);C.Add_Identity_Matrix();
        // construct the rhs
        VECTOR_ND<T> x_B(B.Upper_Triangular_Solve(b_B));
        VECTOR_ND<T> rhs(epsilon_size);MATRIX_MXN<T>::Add_Transpose_Times(epsilon,b_epsilon,rhs);MATRIX_MXN<T>::Add_Transpose_Times(B_inverse_S,x_B,rhs);
        // compute solution vector 
        VECTOR_ND<T> x_S=C.Cholesky_Solve(rhs);
        MATRIX_MXN<T>::Subtract_Times(S,x_S,b_B);x_B=B.Upper_Triangular_Solve(b_B);
        x.Resize(A.n);x.Set_Subvector(1,x_B);x.Set_Subvector(B.n+1,x_S);x=x.Unpermute(permute);}

    if(verbose) {std::stringstream ss;ss<<"Result of linearly constrained optimization:\n"<<x<<std::endl;LOG::filecout(ss.str());}
}
//####################################################################################
// Function Solve_Velocities_for_PD
//####################################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,3> >::
Solve_Velocities_for_PD(const T time,const T dt,bool test_system,bool print_matrix)
{
    assert(global_post_stabilization || (!use_pd_actuators && !use_muscle_actuators)); // must use global post stabilization for pd/muscle actuation
    

    if(global_post_stabilization && global_post_stabilization_matrix_11.n+global_post_stabilization_matrix_21.m){
        // build constrained and muscle controlled delta relative joint velocities (also get joint locations)
        VECTOR_ND<T> constrained_delta_relative_joint_velocities(global_post_stabilization_matrix_11.n),joint_constraint_impulses,muscle_controlled_delta_relative_joint_velocities;
        if(use_muscle_actuators){muscle_activations.Resize(muscle_list->muscles.m);muscle_controlled_delta_relative_joint_velocities.Resize(global_post_stabilization_matrix_21.m);}
        ARRAY<TV> joint_locations(joint_mesh.joints.m);
        for(int i=1;i<=joint_mesh.joints.m;i++) if(joint_constrained_dimensions(i) || joint_muscle_control_dimensions(i)){
            TWIST<TV> delta_relative_twist;Delta_Relative_Twist(joint_mesh.joints(i)->id_number,true,joint_locations(i),delta_relative_twist);
            if(joint_constrained_dimensions(i)){
                assert(joint_constrained_dimensions(i)>=3); // assume unconstrained linear directions, constrained angular directions
                int istart=joint_offset_in_post_stabilization_matrix(i);
                constrained_delta_relative_joint_velocities.Set_Subvector(istart,delta_relative_twist.linear);
                for(int ii=1;ii<=joint_constrained_dimensions(i)-3;ii++){
                    TV u;joint_angular_constraint_matrix(i).Get_Column(ii,u);constrained_delta_relative_joint_velocities(istart+ii+2)=TV::Dot_Product(u,delta_relative_twist.angular);}}
            if(joint_muscle_control_dimensions(i)){ // assume muscle control only affects angular dof's
                int istart=joint_offset_in_muscle_control_matrix(i);
                for(int ii=1;ii<=joint_muscle_control_dimensions(i);ii++){
                    TV u;joint_angular_muscle_control_matrix(i).Get_Column(ii,u);muscle_controlled_delta_relative_joint_velocities(istart+ii-1)=TV::Dot_Product(u,delta_relative_twist.angular);}}}

        if(!use_muscle_actuators || !muscle_list->muscles.m){ // With no muscles we can solve the system directly
            Solve_Minimum_Norm_Solution_For_Linear_Constraints(global_post_stabilization_matrix_11,constrained_delta_relative_joint_velocities,
                joint_constraint_impulses,zero_row_tolerance,verbose);}
        else{
            // With muscles we need to do a least squares solve
            MATRIX_MXN<T> global_post_stabilization_matrix_11_inverse;global_post_stabilization_matrix_11.Cholesky_Inverse(global_post_stabilization_matrix_11_inverse);

            MATRIX_MXN<T> matrix_21_times_11_inverse=global_post_stabilization_matrix_21*global_post_stabilization_matrix_11_inverse;
            MATRIX_MXN<T> A=matrix_21_times_11_inverse*global_post_stabilization_matrix_12-global_post_stabilization_matrix_22;
            VECTOR_ND<T> b(matrix_21_times_11_inverse*constrained_delta_relative_joint_velocities-muscle_controlled_delta_relative_joint_velocities);
            VECTOR_ND<T> muscle_impulse_magnitudes(muscle_list->muscles.m);

            Solve_For_Muscle_Control(A,b,muscle_impulse_magnitudes,dt);

            // apply and increment muscle impulses
            for(int i=1;i<=muscle_list->muscles.m;i++){
                muscle_activations(i)=(muscle_impulse_magnitudes(i)/dt-muscle_list->muscles(i)->Passive_Force(muscle_list->muscles(i)->Total_Length()))/muscle_list->muscles(i)->Force(1);
                muscle_list->muscles(i)->Apply_Fixed_Impulse_At_All_Points(muscle_impulse_magnitudes(i));}

            joint_constraint_impulses=global_post_stabilization_matrix_11_inverse*(constrained_delta_relative_joint_velocities-global_post_stabilization_matrix_12*muscle_impulse_magnitudes);}

        // apply poststabilization impulses
        for(int i=1;i<=joint_mesh.joints.m;i++) if(joint_constrained_dimensions(i)){
            assert(joint_constrained_dimensions(i)>=3); // unconstrained linear directions, constrained angular directions
            int istart=joint_offset_in_post_stabilization_matrix(i);
            TV impulse;joint_constraint_impulses.Get_Subvector(istart,impulse);
            TV angular_impulse;for(int ii=1;ii<=joint_constrained_dimensions(i)-3;ii++){
                TV u;joint_angular_constraint_matrix(i).Get_Column(ii,u);angular_impulse+=joint_constraint_impulses(istart+ii+2)*u;}
            if(joint_mesh.joints(i)->impulse_accumulator) joint_mesh.joints(i)->impulse_accumulator->Add_Impulse(joint_locations(i),TWIST<TV>(impulse,angular_impulse));
            RIGID_BODY<TV>::Apply_Impulse(*Parent(joint_mesh.joints(i)->id_number),*Child(joint_mesh.joints(i)->id_number),joint_locations(i),impulse,angular_impulse);}}

    // the second flag will cause poststabilization to be applied to joints not processed by global post stabilization
    Apply_Poststabilization(use_pd_actuators,global_post_stabilization);

    if(use_angular_damping){
        // precompute damped stuff
        precomputed_desired_damped_angular_velocities.Resize(joint_mesh.joints.m);
        int njoints=0;
        for(int i=1;i<=joint_mesh.joints.m;i++) if(joint_mesh.joints(i)->angular_damping){
            njoints++;JOINT<TV>& joint=*joint_mesh.joints(i);RIGID_BODY<TV> &parent=*Parent(joint.id_number),&child=*Child(joint.id_number);
            parent.Update_Angular_Velocity();child.Update_Angular_Velocity();
            precomputed_desired_damped_angular_velocities(i)=(std::exp(-joint.angular_damping*dt)-1)*RIGID_BODY<TV>::Relative_Angular_Velocity(parent,child);}
            {std::stringstream ss;ss<<"Angular damping of "<<njoints<<" joints"<<std::endl;LOG::filecout(ss.str());}
        Apply_Poststabilization(use_pd_actuators,false,true);}
}
//####################################################################################
// Function Create_Joint_Function
//####################################################################################
template<class T> JOINT_FUNCTION<VECTOR<T,3> >* ARTICULATED_RIGID_BODY<VECTOR<T,3> >::
Create_Joint_Function(const JOINT_ID joint_id)
{
    joint_mesh(joint_id)->Set_Joint_Function(new JOINT_FUNCTION<TV>(*joint_mesh(joint_id),*Parent(joint_id),*Child(joint_id)));
    return joint_mesh(joint_id)->joint_function;
}
//####################################################################################
// Function Output_Articulation_Points
//####################################################################################
template<class T> void ARTICULATED_RIGID_BODY<VECTOR<T,3> >::
Output_Articulation_Points(const STREAM_TYPE stream_type,const std::string& output_directory,const int frame) const
{
    if(joint_mesh.joints.m==0) return;
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(STRING_UTILITIES::string_sprintf("%s/%d/arb_info",output_directory.c_str(),frame));
    TYPED_OSTREAM typed_output(*output,stream_type);
    Write_Binary(typed_output,joint_mesh.joints.m*2);
    for(int i=1;i<=joint_mesh.joints.m;i++){
        JOINT<TV>& joint=*joint_mesh.joints(i);
        const RIGID_BODY<TV> &parent=*Parent(joint.id_number),&child=*Child(joint.id_number);
        TV ap1=parent.World_Space_Point(joint.F_pj().t),ap2=child.World_Space_Point(joint.F_cj().t);
        Write_Binary(typed_output,ap1,parent.Frame()*joint.F_pj(),ap2,child.Frame()*joint.F_pj());}
    delete output;
}
//####################################################################################
template class ARTICULATED_RIGID_BODY<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ARTICULATED_RIGID_BODY<VECTOR<double,3> >;
#endif
