//#####################################################################
// Copyright 2007-2008, Nipun Kwatra, Craig Schroeder, Tamar Shinar, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARTICULATED_RIGID_BODY_1D
//#####################################################################
#ifndef __ARTICULATED_RIGID_BODY_1D__
#define __ARTICULATED_RIGID_BODY_1D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_ID.h>
namespace PhysBAM{
template<class TV> class RIGID_BODY_CONTACT_GRAPH;
template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> class ARTICULATED_RIGID_BODY;
template<class TV> class MUSCLE_LIST;
template<class TV> class RIGID_BODY;
template<class TV> class JOINT_MESH;

template<class T>
class ARTICULATED_RIGID_BODY<VECTOR<T,1> >:public NONCOPYABLE
{
    typedef VECTOR<T,1> TV;

public:
    ARRAY<ARRAY<JOINT_ID> > process_list;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    JOINT_MESH<TV> &joint_mesh;
    MUSCLE_LIST<TV>* muscle_list;
    ARRAY<T> muscle_activations;
    bool use_epsilon_scale;
    int contact_level_iterations,shock_propagation_level_iterations,actuation_iterations;
    bool use_shock_propagation,do_final_pass;
    bool use_pd_actuators,use_muscle_actuators;
    bool constrain_pd_directions;
    bool use_krylov_prestab;
    int max_iterations;

    ARTICULATED_RIGID_BODY(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input);

    ~ARTICULATED_RIGID_BODY();

    bool Has_Actuators() const
    {return use_muscle_actuators || use_pd_actuators;}

//#####################################################################
    void Read(const STREAM_TYPE stream_type,const std::string& directory,const int frame);
    void Write(const STREAM_TYPE stream_type,const std::string& directory,const int frame) const;
    int Parent_Id(JOINT_ID joint_id) const;
    int Child_Id(JOINT_ID joint_id) const;
    RIGID_BODY<TV>* Parent(JOINT_ID joint_id);
    const RIGID_BODY<TV>* Parent(JOINT_ID joint_id) const;
    RIGID_BODY<TV>* Child(JOINT_ID joint_id);
    const RIGID_BODY<TV>* Child(JOINT_ID joint_id) const;
    void Apply_Poststabilization(bool test_system,bool print_matrix,const bool target_pd=false,const bool no_global_post_stabilization_only=false,const bool angular_damping_only=false);
    void Store_Velocities_And_Momenta();
    void Restore_Velocities_And_Momenta();
    void Compute_Position_Based_State(const T dt,const T time);
    void Apply_Prestabilization_To_Joint(const JOINT_ID joint_id,const T dt,const T epsilon_scale=1);
    void Compute_Desired_PD_Velocity(const T dt,const T time);
    void Solve_Velocities_for_PD(const T time,const T dt,bool test_system,bool print_matrix);
    void Initialize_Poststabilization_Projection();
    void Poststabilization_Projection(ARRAY_VIEW<TWIST<TV> > twist,const bool symmetric=false);
    void Generate_Process_List_Using_Contact_Graph(const RIGID_BODY_CONTACT_GRAPH<TV>& contact_graph);
    void Substitute_Joint_Parent_Body(JOINT_ID joint_id,int new_parent,const FRAME<TV>& frame);
    void Substitute_Joint_Child_Body(JOINT_ID joint_id,int new_child,const FRAME<TV>& frame);
    void Substitute_Joint_Parent_Body(JOINT_ID joint_id,int new_parent);
    void Substitute_Joint_Child_Body(JOINT_ID joint_id,int new_child);
    void Apply_Poststabilization_With_CG(T dt,bool correct_position,bool test_system,bool print_matrix){}
//#####################################################################
};
}
#endif
