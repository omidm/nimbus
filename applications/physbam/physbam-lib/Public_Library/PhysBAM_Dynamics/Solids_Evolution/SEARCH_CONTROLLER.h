//#####################################################################
// Copyright 2008-2009, Jon Gretarsson, Michael Lentine, Craig Schroeder, Bridget Vuong.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEARCH_CONTROLLER
//#####################################################################
#ifndef __SEARCH_CONTROLLER__
#define __SEARCH_CONTROLLER__

#define NOFUNCTION 0
#define FORCE 1
#define DRAG 2

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/FLUIDS_PARAMETERS_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER.h>
namespace PhysBAM{

template<class TV> class JOINT_FUNCTION;
template<class T_GRID>
class ENVIRONMENTAL_STATE
{
    typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
  public:
    typedef int HAS_TYPED_READ_WRITE;
    
    bool incorporate_fluids;
    T time, force_to_stay;
    ARRAY<ROTATION<TV>,JOINT_ID> angles;
    
    T_GRID grid;
    T_ARRAYS_SCALAR pressures;

    T external_force_mag;
    TV external_force_dir;

    ENVIRONMENTAL_STATE()
        :incorporate_fluids(false){}

    ENVIRONMENTAL_STATE(bool incorporate_fluids,T_GRID* grid)
    {
        this->incorporate_fluids=incorporate_fluids;
        if(incorporate_fluids) this->grid=*grid;
    }

    bool operator==(const ENVIRONMENTAL_STATE<T_GRID>& state) const
    {
        T threshold=(T)1e-5;
        if(incorporate_fluids) for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()) if(abs(this->pressures(iterator.Cell_Index())-state.pressures(iterator.Cell_Index()))>threshold) return false;
        else if(abs(this->external_force_mag-state.external_force_mag)>1e-5||(this->external_force_dir-state.external_force_dir).Magnitude()>threshold) return false;
        TV vector;vector.Fill(0);vector(1)=1;for(JOINT_ID i(1);i<=this->angles.Size();i++){if(((this->angles(i).Inverse_Rotate(state.angles(i).Rotate(vector)))-vector).Magnitude()>threshold) return false;}
        return true;
    }
 
    void Read(TYPED_ISTREAM& typed_input)
    {
        Read_Binary(typed_input,time);
        Read_Binary(typed_input,incorporate_fluids);
        Read_Binary(typed_input,force_to_stay);
        Read_Binary(typed_input,angles);
        if(incorporate_fluids){
            Read_Binary(typed_input,pressures);
            Read_Binary(typed_input,grid);}
        else{
            Read_Binary(typed_input,external_force_mag);
            Read_Binary(typed_input,external_force_dir);}
    }

    void Write(TYPED_OSTREAM& typed_output) const
    {
        Write_Binary(typed_output,time);
        Write_Binary(typed_output,incorporate_fluids);
        Write_Binary(typed_output,force_to_stay);
        Write_Binary(typed_output,angles);
        if(incorporate_fluids){
            Write_Binary(typed_output,pressures);
            Write_Binary(typed_output,grid);}
        else{
            Write_Binary(typed_output,external_force_mag);
            Write_Binary(typed_output,external_force_dir);}
    }
};

class RIGID_BODY_COLLISION_MANAGER_HASH;

template<class T_GRID>
class SEARCH_CONTROLLER:public NONCOPYABLE
{
    typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::VECTOR_T TV;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename RIGID_BODY_POLICY<TV>::WORLD_SPACE_INERTIA_TENSOR T_WORLD_SPACE_INERTIA_TENSOR;
    typedef typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR T_INERTIA_TENSOR;
    typedef typename TV::SPIN T_SPIN;
    typedef typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER T_CLUSTER;
public:
    typedef int HAS_TYPED_READ_WRITE;

    struct JOINT_STATE
    {
        int parent,child;
        JOINT<TV>* joint;
    };

    SOLIDS_FLUIDS_DRIVER<TV>* driver;
    JOINT_MESH<TV>& joint_mesh;
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    HASHTABLE<int> not_affected_by_fluid; // rigid body index, not present = affected by fluid

    ARRAY<ENVIRONMENTAL_STATE<T_GRID>*> states;
    bool use_random,collecting_data;
    ARRAY<ROTATION<TV>,JOINT_ID>* target_angles;
    DIRECTED_GRAPH<int>* states_graph;
    ARRAY<ENVIRONMENTAL_STATE<T_GRID>*> graph_index_to_state_map;

private: 
    ARRAY<bool,JOINT_ID> pd_state_save,joint_function_active_save;
    ARRAY<T,JOINT_ID> k_p_save;
    ARRAY<TV> X_save,V_save,rigid_velocity_save,rigid_X_save;
    ARRAY<T_SPIN> rigid_angular_momentum_save;
    ARRAY<ROTATION<TV> > rigid_rotation_save;
    ARRAY<TV> V_pd,rigid_velocity_pd;
    ARRAY<T_SPIN> rigid_angular_momentum_pd;
    bool project_at_frame_boundaries_save;
    bool fluid_affects_solid_save;
    bool constrain_pd_save;
    T time_save;

    T_FACE_ARRAYS_SCALAR fluid_velocity_save,fluid_velocity_save_for_projection_save,beta_face_save;
    T_FACE_ARRAYS_SCALAR fluid_velocity_pd_save;
    T_FACE_ARRAYS_BOOL psi_N_save;
    T_ARRAYS_SCALAR fluid_pressure_save,fluid_pressure_save_for_projection_save;
    T_ARRAYS_SCALAR fluid_pressure_pd_save;
    T_ARRAYS_BOOL psi_D_save;
    T_ARRAYS_INT filled_region_colors_save;

public:
    FLUIDS_PARAMETERS_UNIFORM<T_GRID>* fluids_parameters;
    VECTOR_ND<T> negative_gradient;

    ARRAY<int> global_clusters;
    ARRAY<PAIR<int,int>,JOINT_ID> joint_clusters;
    ARRAY<ARRAY<JOINT_STATE>,JOINT_ID> removed_joints;
    ARRAY<ARRAY<PAIR<JOINT_STATE,JOINT_STATE> >,JOINT_ID> swap_joints;
    ARRAY<JOINT_STATE> removed_joints_global;
    ARRAY<PAIR<JOINT_STATE,JOINT_STATE> > swap_joints_global;
    ARRAY<RIGID_BODY<TV>*> stored_bodies;

    ARRAY<int,JOINT_ID> objective;
    bool minimize,solve_minimization,line_search,use_drag_direction,steady_state_drag;
    T steady_state_drag_time;

    TV drag_direction;

    ARRAY<ARRAY<int> > bone_hierarchy;
    int root_particle_index;

    ARRAY<ROTATION<TV>,JOINT_ID> current_angle;
    ARRAY<PAIR<T,T> > dF_array_multipliers;
    T real_dx,dx,dt_hyp,dt_per_search_step,last_search_time,min_multiplier;
    bool hypothetical_step,pd_step,drag_step;
    mutable bool initialized;
    int max_iterations;
    T threshold;
    int rigid_body_particles_number;
    bool incorporate_fluids,use_projection,use_clusters;
    ARRAY<JOINT_ID> current_joints,strain_joints;
    T old_time;

    SEARCH_CONTROLLER(SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,SOLIDS_FLUIDS_DRIVER<TV>* driver_input)
        :driver(driver_input),joint_mesh(solid_body_collection_input.rigid_body_collection.articulated_rigid_body.joint_mesh),solid_body_collection(solid_body_collection_input),
        use_random(false),collecting_data(false),states_graph(0),fluids_parameters(0),minimize(true),solve_minimization(false),
        line_search(false),use_drag_direction(true),steady_state_drag(false),steady_state_drag_time((T)4),root_particle_index(-1),real_dx(15*(T)pi/180),dx(15*(T)pi/180),dt_hyp(0),
        dt_per_search_step(0),last_search_time(-1),min_multiplier((T)1e-4),hypothetical_step(false),pd_step(false),drag_step(false),initialized(false),max_iterations(100),threshold((T)1e-1),
        rigid_body_particles_number(0),incorporate_fluids(false),use_projection(false),use_clusters(false),old_time((T)-1)
    {}

    ~SEARCH_CONTROLLER()
    {Remove_Clusters();}

    JOINT<TV>* Joint(int parent_id,int child_id) const
    {
        for(int i=1;i<=joint_mesh.joints.m;i++){
            JOINT_ID joint_id=joint_mesh.joints(i)->id_number;
            if(Parent_Id(joint_id)==parent_id && Child_Id(joint_id)==child_id) return joint_mesh.joints(i);}
        return 0;
    }
    
    int Parent_Id(JOINT_ID joint_id) const
    {return joint_mesh.Parent_Id(joint_id);}

    int Child_Id(JOINT_ID joint_id) const
    {return joint_mesh.Child_Id(joint_id);}

    RIGID_BODY<TV>* Parent(JOINT_ID joint_id)
    {return &solid_body_collection.rigid_body_collection.Rigid_Body(Parent_Id(joint_id));}

    const RIGID_BODY<TV>* Parent(JOINT_ID joint_id) const
    {return &solid_body_collection.rigid_body_collection.Rigid_Body(Parent_Id(joint_id));}

    RIGID_BODY<TV>* Child(JOINT_ID joint_id)
    {return &solid_body_collection.rigid_body_collection.Rigid_Body(Child_Id(joint_id));}

    const RIGID_BODY<TV>* Child(JOINT_ID joint_id) const
    {return &solid_body_collection.rigid_body_collection.Rigid_Body(Child_Id(joint_id));}

    bool Fluid_Node() const
    {SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters=driver->example.solids_fluids_parameters;
    return (!solids_fluids_parameters.mpi_solid_fluid || solids_fluids_parameters.mpi_solid_fluid->Fluid_Node());}

    bool Simulate_Solids() const
    {SOLIDS_PARAMETERS<TV>& solids_parameters=driver->example.solids_parameters;SOLID_BODY_COLLECTION<TV>& solid_body_collection=driver->example.solid_body_collection;
    SOLIDS_FLUIDS_PARAMETERS<TV>& solids_fluids_parameters=driver->example.solids_fluids_parameters;
    return ((solids_fluids_parameters.mpi_solid_fluid || solid_body_collection.deformable_body_collection.simulate) && solid_body_collection.deformable_body_collection.particles.array_collection->Size()) || ((solids_fluids_parameters.mpi_solid_fluid || solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies) && solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size());}

//#####################################################################
    void Update_Position_Based_State(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    JOINT_FUNCTION<TV>* Create_Joint_Function(const JOINT_ID joint_id);
    void Save_Position(ARRAY<TV>& X,ARRAY<TV>& rigid_X,ARRAY<ROTATION<TV> >& rigid_rotation);
    void Restore_Position(ARRAY_VIEW<const TV> X,ARRAY_VIEW<const TV> rigid_X,ARRAY_VIEW<const ROTATION<TV> > rigid_rotation);
    void Save_Velocity(ARRAY<TV>& V,ARRAY<TV>& rigid_velocity,ARRAY<T_SPIN>& rigid_angular_momentum);
    void Restore_Velocity(ARRAY<TV>& V,ARRAY<TV>& rigid_velocity,ARRAY<T_SPIN>& rigid_angular_momentum);
    void Save_Nested_State(T_FACE_ARRAYS_SCALAR& face_velocities);
    void Restore_Nested_State(T_FACE_ARRAYS_SCALAR& face_velocities);
    void Save_State(T_FACE_ARRAYS_SCALAR& face_velocities);
    void Restore_State(T_FACE_ARRAYS_SCALAR& face_velocities);
    void Save_PD_State();
    void Restore_PD_State();
    void Propogate_Solid_Helper(ARRAY<int>& cluster_bodies,TV& cluster_translation,TWIST<TV>& cluster_twist);
    void Make_Cluster_List(int parent,ARRAY<int>& child_bodies);
    void Project_Solid_Velocities(const T time);
    void Project_Velocities(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Set_PD_Targets();
    void Zero_PD_Targets();
    void Create_Clusters_From_Joint_List(const ARRAY<bool,JOINT_ID>& blocking_joint,ARRAY<ARRAY<int> >& body_lists,ARRAY<VECTOR<int,2>,JOINT_ID>& adjacent_lists);
    void Create_All_Clusters(RIGID_BODY_COLLISION_MANAGER_HASH* collision_manager);
    void Remove_Clusters();
    T Evaluate_Force(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const T dx,const T magnitude,const VECTOR_ND<T>* grad);
    PAIR<T,T_SPIN> Evaluate_Force_For_Joint(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const JOINT_ID joint,const int dimension,const T dx_multiplier);
    T Evaluate_Force_To_Stay(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,int hack=false);
    T Evaluate_Force_To_Stay_For_Joint(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,ARRAY<JOINT_ID>& joints);
    void Add_Solid_Drag(const T time,ARRAY<TV> &F,ARRAY<TWIST<TV> > &rigid_F);
    void Add_Fluid_Drag(const T dt,const T time,ARRAY<TV> &F,ARRAY<TWIST<TV> > &rigid_F);
    T Evaluate_Drag(const T dt,const T time);
    T Evaluate_Drag_For_Joint(const T dt,const T time,ARRAY<JOINT_ID>& joints);
    void Take_Hypothetical_Step(T_FACE_ARRAYS_SCALAR& face_velocities,const T time,bool simulate_fluids);
    void Take_Hypothetical_Fluids_Step(const T time);
    bool Compute_Gradient_From_Graph(ENVIRONMENTAL_STATE<T_GRID>* current_state);
    void Compute_Gradient(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Normalize_Gradient(VECTOR_ND<T> &grad);
    VECTOR_ND<T> Normalized_Gradient(const VECTOR_ND<T> &grad);
    T Golden_Section_Search(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,T a,T b,const T dx,T* F=NULL);
    T Steepest_Descent(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Advance_One_Time_Step_Position(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Preprocess_Drag_Substep(const T time);
    void Read(const STREAM_TYPE stream_type,const std::string& directory,const int frame);
    void Read(TYPED_ISTREAM& typed_input);
    void Write(const STREAM_TYPE stream_type,const std::string& directory,const int frame) const;
    void Write(TYPED_OSTREAM& typed_output) const;
//#####################################################################
};

}
#endif
