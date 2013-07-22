//#####################################################################
// Copyright 2008, Nipun Kwatra, Craig Schroeder, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ASYNCHRONOUS_EVOLUTION
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Topology/TETRAHEDRON_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SCALED_DEFORMABLES_FORCES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/ASYNCHRONOUS_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> ASYNCHRONOUS_EVOLUTION<TV>::
ASYNCHRONOUS_EVOLUTION(SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,SOLIDS_EVOLUTION<TV>* solids_evolution_input,const T cfl,const bool use_projection_input,
    const T projection_rigidity_input)
    :solid_body_collection(solid_body_collection_input),solids_evolution(solids_evolution_input),
    asynchronous_mode(COARSE_SCALE),time_n(0),time_np1(0),dt_epsilon(0),initialized(false),
    use_velocity_averaging(true),use_adaptive(false),treat_mixed_particles_kinematic_in_finescale(false),
    use_squared_ratio_for_implicit_velocity_independent_scale(true),
    use_projection(use_projection_input),projection_rigidity(projection_rigidity_input)
{
    solid_body_collection.Set_CFL_Number(cfl);
    dynamic_cast<NEWMARK_EVOLUTION<TV>*>(solids_evolution)->asynchronous_evolution=this;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> ASYNCHRONOUS_EVOLUTION<TV>::
~ASYNCHRONOUS_EVOLUTION()
{
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Initialize()
{
    if(!use_projection && (finescale_forces_indices.m+scaled_coarsescale_forces_indices.m)!=solid_body_collection.deformable_body_collection.deformables_forces.m){
        // TODO: remove allowance for the case when use_projection=true 
        PHYSBAM_FATAL_ERROR("Forces should only be added through ASYNCHRONOUS_EVOLUTION");}
    Store_State((T)0);
    Setup_Particle_Indices_From_Maps();
    Setup_Scaled_Forces_Rewind_Particles();
    initialized=true;
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Preprocess_Frame(const int frame,const T time_at_frame,const T time)
{
    if(use_adaptive){
        Set_Asynchronous_Within_Distance((T)1);
        Setup_Scaled_Forces_Rewind_Particles();}
    time_np1=time_at_frame;
    Store_State(time);
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Preprocess_Substep(const T dt,const T time)
{
    if(time+dt>=time_np1) Setup_Coarsescale_Simulation(dt,time);
    else Setup_Finescale_Simulation();
}
//#####################################################################
// Function Postprocess_Substep
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Postprocess_Substep(const T dt,const T time)
{
    if(treat_mixed_particles_kinematic_in_finescale && asynchronous_mode==FINE_SCALE) Accumulate_Finescale_Impulse_On_Mixed_Particles(dt,time);
}
//#####################################################################
// Function Setup_Finescale_Simulation
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Setup_Finescale_Simulation()
{
    if(asynchronous_mode==FINE_SCALE) return;
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Setup_Finescale_Simulation (ASYNCHRONOUS_EVOLUTION)",0,0);
    Disable_Forces(scaled_coarsescale_forces_indices,saved_force_settings);
    asynchronous_mode=FINE_SCALE;
}
//#####################################################################
// Function Setup_Coarsescale_Simulation
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Setup_Coarsescale_Simulation(const T dt,const T time)
{
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Setup_Coarsescale_Simulation (ASYNCHRONOUS_EVOLUTION)",0,0);
    dt_epsilon=dt;
    Store_Time_Last_Finescale_State();
    if(treat_mixed_particles_kinematic_in_finescale) Apply_Accumulated_Finescale_Impulse_To_Mixed_Particles();
    Enable_Forces(scaled_coarsescale_forces_indices,saved_force_settings);
    Set_Scaled_Coarsescale_Forces_Scale(dt_epsilon,time_np1-time_n,true);
    Set_Scaled_Coarsescale_Forces_Rewind_State(true,true);
    asynchronous_mode=COARSE_SCALE;
}
//#####################################################################
// Function Disable_Forces
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Disable_Forces(const ARRAY<int>& disable_forces_indices,ARRAY<VECTOR<bool,3> >& saved_force_settings)
{
    for(int i=1;i<=disable_forces_indices.m;i++){int index=disable_forces_indices(i);
        DEFORMABLES_FORCES<TV>& force=*solid_body_collection.deformable_body_collection.deformables_forces(index);
        // Don't disable if already disabled. Prevents saved_force_settings to be overwritten for consecutive calls.
        if(!(force.use_velocity_independent_forces || force.use_velocity_dependent_forces || force.use_implicit_velocity_independent_forces)) continue;

        saved_force_settings(index)[1]=force.use_velocity_independent_forces;
        saved_force_settings(index)[2]=force.use_velocity_dependent_forces;
        saved_force_settings(index)[3]=force.use_implicit_velocity_independent_forces;

        force.use_velocity_independent_forces=false;
        force.use_velocity_dependent_forces=false;
        force.use_implicit_velocity_independent_forces=false;}
}
//#####################################################################
// Function Enable_Forces
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Enable_Forces(const ARRAY<int>& enable_forces_indices,const ARRAY<VECTOR<bool,3> >& saved_force_settings)
{
    for(int i=1;i<=enable_forces_indices.m;i++){int index=enable_forces_indices(i);
        solid_body_collection.deformable_body_collection.deformables_forces(index)->use_velocity_independent_forces=saved_force_settings(index)[1];
        solid_body_collection.deformable_body_collection.deformables_forces(index)->use_velocity_dependent_forces=saved_force_settings(index)[2];
        solid_body_collection.deformable_body_collection.deformables_forces(index)->use_implicit_velocity_independent_forces=saved_force_settings(index)[3];}
}
//#####################################################################
// Function Rewind_Coarsescale_Particle_Positions
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Rewind_Coarsescale_Particle_Positions()
{
    for(int i=1;i<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();i++) if (!finescale_forces_particles_map.Contains(i)){
        solid_body_collection.deformable_body_collection.particles.X(i)=X_n(i);}
    for(int i=1;i<=solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();i++) if (!finescale_forces_rigid_body_particles_map.Contains(i)){
        solid_body_collection.rigid_body_collection.rigid_body_particle.X(i)=rigid_X_n(i);
        solid_body_collection.rigid_body_collection.rigid_body_particle.rotation(i)=rigid_rotation_n(i);}
}
//#####################################################################
// Function Store_State
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Store_State(const T time)
{
    X_n=solid_body_collection.deformable_body_collection.particles.X;
    V_n=solid_body_collection.deformable_body_collection.particles.V;
    rigid_V_n=solid_body_collection.rigid_body_collection.rigid_body_particle.V;
    rigid_angular_velocity_n=solid_body_collection.rigid_body_collection.rigid_body_particle.angular_velocity;
    rigid_X_n=solid_body_collection.rigid_body_collection.rigid_body_particle.X;
    rigid_rotation_n=solid_body_collection.rigid_body_collection.rigid_body_particle.rotation;
    angular_momentum_n=solid_body_collection.rigid_body_collection.rigid_body_particle.angular_momentum;
    time_n=time;
}
//#####################################################################
// Function Store_Time_np1_State
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Store_Time_np1_State()
{
    X_np1=solid_body_collection.deformable_body_collection.particles.X;
    V_np1=solid_body_collection.deformable_body_collection.particles.V;
    rigid_V_np1=solid_body_collection.rigid_body_collection.rigid_body_particle.V;
    rigid_angular_velocity_np1=solid_body_collection.rigid_body_collection.rigid_body_particle.angular_velocity;
    rigid_X_np1=solid_body_collection.rigid_body_collection.rigid_body_particle.X;
    rigid_rotation_np1=solid_body_collection.rigid_body_collection.rigid_body_particle.rotation;
    angular_momentum_np1=solid_body_collection.rigid_body_collection.rigid_body_particle.angular_momentum;
}
//#####################################################################
// Function Store_Time_Last_Finescale_State
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Store_Time_Last_Finescale_State()
{
    X_last_finescale=solid_body_collection.deformable_body_collection.particles.X;
    V_last_finescale=solid_body_collection.deformable_body_collection.particles.V;
    rigid_V_last_finescale=solid_body_collection.rigid_body_collection.rigid_body_particle.V;
    rigid_angular_velocity_last_finescale=solid_body_collection.rigid_body_collection.rigid_body_particle.angular_velocity;
    rigid_X_last_finescale=solid_body_collection.rigid_body_collection.rigid_body_particle.X;
    rigid_rotation_last_finescale=solid_body_collection.rigid_body_collection.rigid_body_particle.rotation;
    angular_momentum_last_finescale=solid_body_collection.rigid_body_collection.rigid_body_particle.angular_momentum;
}
//#####################################################################
// Function Add_Finescale_Force
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Add_Finescale_Force(DEFORMABLES_FORCES<TV>* force,const ARRAY<int>& affected_particle_indices_list,const ARRAY<int>& affected_rigid_body_particle_indices_list,const bool implicit_velocity_independent)
{
    force->use_implicit_velocity_independent_forces=implicit_velocity_independent;
    finescale_forces_indices.Append(solid_body_collection.Add_Force(force));
    Add_Force_Settings(force);

    finescale_forces_particles_map.Set_All(affected_particle_indices_list);
    finescale_forces_rigid_body_particles_map.Set_All(affected_rigid_body_particle_indices_list);

    if(initialized){Setup_Particle_Indices_From_Maps();Setup_Scaled_Forces_Rewind_Particles();}
}
//#####################################################################
// Function Add_Coarsescale_Force
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Add_Coarsescale_Force(DEFORMABLES_FORCES<TV>* force,const ARRAY<int>& affected_particle_indices_list,const ARRAY<int>& affected_rigid_body_particle_indices_list,const bool implicit_velocity_independent,const bool add_to_blobs)
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;

    force->use_implicit_velocity_independent_forces=implicit_velocity_independent;
    SCALED_DEFORMABLES_FORCES<TV>* scaled_coarsescale_force=new SCALED_DEFORMABLES_FORCES<TV>(force,particles,X_n,V_n);
    scaled_coarsescale_forces_indices.Append(solid_body_collection.Add_Force(scaled_coarsescale_force));
    Add_Force_Settings(scaled_coarsescale_force);
    ARRAY<int> current_force_index;current_force_index.Append(scaled_coarsescale_forces_indices.Last());    
    Disable_Forces(current_force_index,saved_force_settings);
    
    coarsescale_forces_particles_map.Set_All(affected_particle_indices_list);
    coarsescale_forces_rigid_body_particles_map.Set_All(affected_rigid_body_particle_indices_list);

    if(use_projection && add_to_blobs) Add_Blob_From_Particles(affected_particle_indices_list);
    
    if(initialized){Setup_Particle_Indices_From_Maps();Setup_Scaled_Forces_Rewind_Particles();}
}
//#####################################################################
// Function Add_Force_Settings
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Add_Force_Settings(DEFORMABLES_FORCES<TV>* force)
{
    VECTOR<bool,3> force_setting(force->use_velocity_independent_forces,force->use_velocity_dependent_forces,force->use_implicit_velocity_independent_forces);
    saved_force_settings.Append(force_setting);
}
//#####################################################################
// Function Setup_Particle_Indices_From_Maps
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Setup_Particle_Indices_From_Maps()
{
    // coarsescale_particles are particles not affected by the finescale forces
    finescale_forces_particles_map.Get_Complementary_Keys(IDENTITY_ARRAY<int>(solid_body_collection.deformable_body_collection.particles.array_collection->Size()),coarsescale_particle_indices);
    finescale_forces_rigid_body_particles_map.Get_Complementary_Keys(IDENTITY_ARRAY<int>(solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size()),coarsescale_rigid_body_particle_indices);

    // both_forces_particles are particles affected by both finescale and coarsescale forces
    both_forces_particles_indices.Remove_All();both_forces_rigid_body_particles_indices.Remove_All();
    for(HASHTABLE_ITERATOR<int> it(finescale_forces_particles_map);it.Valid();it.Next()){
        if(coarsescale_forces_particles_map.Contains(it.Key())) both_forces_particles_indices.Append(it.Key());}
    for(HASHTABLE_ITERATOR<int> it(finescale_forces_rigid_body_particles_map);it.Valid();it.Next()){
        if(coarsescale_forces_rigid_body_particles_map.Contains(it.Key())) both_forces_rigid_body_particles_indices.Append(it.Key());}

    // initialize auxiliary structures
    if(treat_mixed_particles_kinematic_in_finescale){
        accumulated_finescale_impulse_on_mixed_particles.Resize(both_forces_particles_indices.m,true,false);}
}
//#####################################################################
// Function Setup_Scaled_Forces_Rewind_Particles
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Setup_Scaled_Forces_Rewind_Particles()
{
    // rewind coarsescale and mixed particles
    ARRAY<int> rewind_to_time_n_particle_indices=both_forces_particles_indices;
    rewind_to_time_n_particle_indices.Append_Elements(coarsescale_particle_indices);

    for(int i=1;i<=scaled_coarsescale_forces_indices.m;i++){int index=scaled_coarsescale_forces_indices(i);
        SCALED_DEFORMABLES_FORCES<TV>* scaled_force=dynamic_cast<SCALED_DEFORMABLES_FORCES<TV>* >(solid_body_collection.deformable_body_collection.deformables_forces(index));
        scaled_force->Set_Rewind_To_Time_n_Particle_Indices(rewind_to_time_n_particle_indices);}
}
//#####################################################################
// Function Set_Scaled_Coarsescale_Forces_Scale
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Set_Scaled_Coarsescale_Forces_Scale(const T dt_partial,const T dt_full,const bool position_update)
{
    T ratio=dt_full/dt_partial;
    T scale=ratio;
    if(position_update) scale=(T)2.*ratio-(T)1.;
    T scale_implicit=use_squared_ratio_for_implicit_velocity_independent_scale?scale*scale:scale;
    for(int i=1;i<=scaled_coarsescale_forces_indices.m;i++){int index=scaled_coarsescale_forces_indices(i);
        dynamic_cast<SCALED_DEFORMABLES_FORCES<TV>*>(solid_body_collection.deformable_body_collection.deformables_forces(index))->Set_Scale(scale,scale,scale_implicit);}
}
//#####################################################################
// Function Set_Scaled_Coarsescale_Forces_Rewind_State
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Set_Scaled_Coarsescale_Forces_Rewind_State(const bool rewind_positions,const bool rewind_velocities)
{
    for(int i=1;i<=scaled_coarsescale_forces_indices.m;i++){int index=scaled_coarsescale_forces_indices(i);
        dynamic_cast<SCALED_DEFORMABLES_FORCES<TV>* >(solid_body_collection.deformable_body_collection.deformables_forces(index))->Set_Rewind_State(rewind_positions,rewind_velocities);}
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Update_Time_Varying_Material_Properties(const T time)
{
    if(use_projection) Compute_Blob_Properties();

    if((asynchronous_mode==FINE_SCALE) || (dt_epsilon==0)) return;
    if(time>=time_np1-dt_epsilon/2){
        Set_Scaled_Coarsescale_Forces_Scale(dt_epsilon,time_np1-time_n,false);
        Set_Scaled_Coarsescale_Forces_Rewind_State(false,true);}
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Set_External_Positions(ARRAY_VIEW<FRAME<TV> > X,const T time)
{
    if(asynchronous_mode==FINE_SCALE || time!=time_np1) return;
    PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Set_External_Positions(ARRAY_VIEW<TV> X,const T time)
{
    if(use_projection || !use_velocity_averaging || asynchronous_mode==FINE_SCALE || time!=time_np1) return;
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Set_External_Positions (before). time=%f, time_np1=%f, time_n=%f (ASYNCHRONOUS_EVOLUTION)",time,time_np1,time_n),0,1);

    for(int i=1;i<=coarsescale_particle_indices.m;i++){int p=coarsescale_particle_indices(i);
        solid_body_collection.deformable_body_collection.particles.X(p)=X_n(p);}
    solid_body_collection.deformable_body_collection.particles.Euler_Step_Position(coarsescale_particle_indices,time_np1-time_n);

    PHYSBAM_DEBUG_WRITE_SUBSTEP("Set_External_Positions (after) (ASYNCHRONOUS_EVOLUTION)",0,1);
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time)
{
    if(!treat_mixed_particles_kinematic_in_finescale || asynchronous_mode==COARSE_SCALE) return;
    for(int i=1;i<=both_forces_particles_indices.m;i++){int p=both_forces_particles_indices(i);
        V(p)=V_n(p);}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time)
{
    if(!treat_mixed_particles_kinematic_in_finescale || asynchronous_mode==COARSE_SCALE) return;
    for(int i=1;i<=both_forces_particles_indices.m;i++){int p=both_forces_particles_indices(i);
        V(p)=TV();}
}
//#####################################################################
// Function Limit_Solids_Dt
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Limit_Solids_Dt(T& dt,const T time)
{
    // if(asynchronous_mode==COARSE_SCALE && fixed_dt) dt=fixed_dt;
}
//#####################################################################
// Function Get_Active_Forces_Contributions
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Get_Active_Forces_Contributions(const T dt,const T time,ARRAY<TV>& force_on_particles,ARRAY<TWIST<TV> >& force_on_rigid_body_particles)
{
    rigid_B_full.Resize(solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size(),true,false);
    rigid_F_full.Resize(solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size(),true,false);
    B_full.Resize(solid_body_collection.deformable_body_collection.particles.array_collection->Size(),true,false);
    F_full.Resize(solid_body_collection.deformable_body_collection.particles.array_collection->Size(),true,false);
    ARRAYS_COMPUTATIONS::Fill(rigid_B_full,TWIST<TV>());
    ARRAYS_COMPUTATIONS::Fill(rigid_F_full,TWIST<TV>());
    ARRAYS_COMPUTATIONS::Fill(B_full,TV());
    ARRAYS_COMPUTATIONS::Fill(F_full,TV());

    solid_body_collection.Add_Velocity_Independent_Forces(B_full,rigid_B_full,time+dt);

    ARRAY<TWIST<TV> > twist;twist.Resize(solid_body_collection.rigid_body_collection.rigid_body_particle.V.Size());
    for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(solid_body_collection.rigid_body_collection.rigid_body_particle.V(i),solid_body_collection.rigid_body_collection.rigid_body_particle.angular_velocity(i));
    GENERALIZED_VELOCITY<TV> F(F_full,rigid_F_full,solid_body_collection),B(B_full,rigid_B_full,solid_body_collection),
        V(solid_body_collection.deformable_body_collection.particles.V,twist,solid_body_collection);

    solid_body_collection.Implicit_Velocity_Independent_Forces(V.V.array,V.rigid_V.array,F.V.array,F.rigid_V.array,dt,time+dt);
    solid_body_collection.Add_Velocity_Dependent_Forces(V.V.array,V.rigid_V.array,F.V.array,F.rigid_V.array,time+dt);
        
    for(int i=1;i<=F.V.Size();i++){int p=F.V.indices(i);
        force_on_particles(p)=F.V(i)+B.V(i);}
    for(int i=1;i<=F.rigid_V.Size();i++){int p=F.rigid_V.indices(i);
        force_on_rigid_body_particles(p)=F.rigid_V(i)+B.rigid_V(i);}
    for(int i=1;i<twist.Size();i++){
        solid_body_collection.rigid_body_collection.rigid_body_particle.V(i)=twist(i).linear;solid_body_collection.rigid_body_collection.rigid_body_particle.angular_velocity(i)=twist(i).angular;}
}
//#####################################################################
// Function Average_Velocities
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Average_Velocities()
{
    ARRAY<TWIST<TV> > twist;twist.Resize(solid_body_collection.rigid_body_collection.rigid_body_particle.V.Size());
    for(int i=1;i<=twist.Size();i++) twist(i)=TWIST<TV>(solid_body_collection.rigid_body_collection.rigid_body_particle.V(i),solid_body_collection.rigid_body_collection.rigid_body_particle.angular_velocity(i));
    ARRAY<TWIST<TV> > twist_last_finescale;twist_last_finescale.Resize(rigid_V_last_finescale.Size());
    for(int i=1;i<=twist.Size();i++) twist_last_finescale(i)=TWIST<TV>(rigid_V_last_finescale(i),rigid_angular_velocity_last_finescale(i));
    GENERALIZED_VELOCITY<TV> GV_last_finescale(V_last_finescale,twist_last_finescale,solid_body_collection),
        V(solid_body_collection.deformable_body_collection.particles.V,twist,solid_body_collection);
    V+=GV_last_finescale;
    V*=(T).5;
    for(int i=1;i<twist.Size();i++){
        solid_body_collection.rigid_body_collection.rigid_body_particle.V(i)=twist(i).linear;solid_body_collection.rigid_body_collection.rigid_body_particle.angular_velocity(i)=twist(i).angular;
        rigid_V_last_finescale(i)=twist_last_finescale(i).linear;rigid_angular_velocity_last_finescale(i)=twist_last_finescale(i).angular;}
}
//#####################################################################
// Function Position_Velocity_Update
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Position_Velocity_Update(const T dt,const T time)
{
    assert(dt==dt_epsilon);
    if(use_projection) PHYSBAM_FATAL_ERROR();
    
    if(use_velocity_averaging){
        Average_Velocities();
        return;}

    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particle;

    finescale_force_on_particles.Resize(particles.array_collection->Size(),false,false);
    coarsescale_force_on_particles.Resize(particles.array_collection->Size(),false,false);
    finescale_force_on_rigid_body_particles.Resize(rigid_body_particles.array_collection->Size(),false,false);
    coarsescale_force_on_rigid_body_particles.Resize(rigid_body_particles.array_collection->Size(),false,false);

    // Get finescale forces
    Disable_Forces(scaled_coarsescale_forces_indices,saved_force_settings);
    Get_Active_Forces_Contributions(dt,time,finescale_force_on_particles,finescale_force_on_rigid_body_particles);
    Enable_Forces(scaled_coarsescale_forces_indices,saved_force_settings);
    // Get coarsescale forces
    Disable_Forces(finescale_forces_indices,saved_force_settings);
    Get_Active_Forces_Contributions(dt,time,coarsescale_force_on_particles,coarsescale_force_on_rigid_body_particles);
    Enable_Forces(finescale_forces_indices,saved_force_settings);

    T ratio=(time_np1-time_n)/dt;
    for(int p=1;p<=particles.array_collection->Size();p++){T one_over_mass=particles.one_over_mass(p);
        particles.V(p)=V_last_finescale(p)+(T).5*dt*one_over_mass*(finescale_force_on_particles(p)+ratio*coarsescale_force_on_particles(p));}
    for(int p=1;p<=rigid_body_particles.array_collection->Size();p++){
        RIGID_BODY_MASS<TV,true>world_space_rigid_mass_inverse=solids_evolution->world_space_rigid_mass_inverse(p);
        TWIST<TV> twist_last_finescale;twist_last_finescale.linear=rigid_V_last_finescale(p);twist_last_finescale.angular=rigid_angular_velocity_last_finescale(p);
        TWIST<TV> twist=twist_last_finescale+
            (T).5*dt*(world_space_rigid_mass_inverse*(finescale_force_on_rigid_body_particles(p)+ratio*coarsescale_force_on_rigid_body_particles(p)));
        rigid_body_particles.V(p)=twist.linear;rigid_body_particles.angular_velocity(p)=twist.angular;}
}
//#####################################################################
// Function Accumulate_Finescale_Impulse_On_Mixed_Particles
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Accumulate_Finescale_Impulse_On_Mixed_Particles(const T dt, const T time)
{
    finescale_force_on_particles.Resize(solid_body_collection.deformable_body_collection.particles.array_collection->Size(),false,false);
    finescale_force_on_rigid_body_particles.Resize(solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size(),false,false);
    Get_Active_Forces_Contributions(dt,time,finescale_force_on_particles,finescale_force_on_rigid_body_particles);
    accumulated_finescale_impulse_on_mixed_particles+=finescale_force_on_particles.Subset(both_forces_particles_indices)*dt;
}
//#####################################################################
// Function Apply_Accumulated_Finescale_Impulse_To_Mixed_Particles
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Apply_Accumulated_Finescale_Impulse_To_Mixed_Particles()
{
    for(int i=1;i<=both_forces_particles_indices.m;i++){int p=both_forces_particles_indices(i);
        solid_body_collection.deformable_body_collection.particles.V(p)+=accumulated_finescale_impulse_on_mixed_particles(i);}
    ARRAYS_COMPUTATIONS::Fill(accumulated_finescale_impulse_on_mixed_particles,TV());
}
//#####################################################################
// Function Set_Force_Active_Particles
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Set_Force_Active_Particles(DEFORMABLES_FORCES<TV>* force,const ARRAY<bool>& fine_list,bool is_fine)
{
    ARRAY<int> list;
    FORCE_ELEMENTS* fe=0;
    if(LINEAR_SPRINGS<TV>* spring=dynamic_cast<LINEAR_SPRINGS<TV>*>(force)){
        for(int i=1;i<=spring->segment_mesh.elements.m;i++) if(fine_list.Subset(spring->segment_mesh.elements(i)).Contains(false)!=is_fine) list.Append(i);
        fe=&spring->force_segments;
        if(!is_fine) for(int i=1;i<=list.m;i++) for(int j=1;j<=2;j++) coarsescale_forces_particles_map.Set(spring->segment_mesh.elements(list(i))(j));}
    else if(GRAVITY<TV>* gravity=dynamic_cast<GRAVITY<TV>*>(force)){
        for(int i=1;i<=gravity->influenced_particles->m;i++) if(fine_list((*gravity->influenced_particles)(i))==is_fine) list.Append((*gravity->influenced_particles)(i));
        fe=&gravity->force_particles;
        if(!is_fine) coarsescale_forces_particles_map.Set_All(list);}
    else if(LINEAR_ALTITUDE_SPRINGS<TV,3>* spring=dynamic_cast<LINEAR_ALTITUDE_SPRINGS<TV,3>*>(force)){
        for(int i=1;i<=spring->mesh.elements.m;i++) if(fine_list.Subset(spring->mesh.elements(i)).Contains(false)!=is_fine) list.Append(i);
        fe=&spring->force_elements;
        if(!is_fine) for(int i=1;i<=list.m;i++) for(int j=1;j<=4;j++) coarsescale_forces_particles_map.Set(spring->mesh.elements(list(i))(j));}
    else if(SCALED_DEFORMABLES_FORCES<TV>* scaled_force=dynamic_cast<SCALED_DEFORMABLES_FORCES<TV>*>(force)){
        Set_Force_Active_Particles(scaled_force->base_force,fine_list,is_fine);}
    else PHYSBAM_FATAL_ERROR("Don't know how to handle this force");
    if(fe) PHYSBAM_FATAL_ERROR("This code got removed.");
}
//#####################################################################
// Function Set_Asynchronous_Particles
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Set_Asynchronous_Particles(const ARRAY<int>& async_list)
{
    ARRAY<bool> fine_list(solid_body_collection.deformable_body_collection.particles.array_collection->Size());
    INDIRECT_ARRAY<ARRAY<bool>,ARRAY<int>&> fine_subset=fine_list.Subset(async_list);
    ARRAYS_COMPUTATIONS::Fill(fine_subset,true);
    coarsescale_forces_particles_map.Remove_All();
    for(int i=1;i<=scaled_coarsescale_forces_indices.m;i++) Set_Force_Active_Particles(solid_body_collection.deformable_body_collection.deformables_forces(scaled_coarsescale_forces_indices(i)),fine_list,false);
    for(int i=1;i<=finescale_forces_indices.m;i++) Set_Force_Active_Particles(solid_body_collection.deformable_body_collection.deformables_forces(finescale_forces_indices(i)),fine_list,true);

    // update particle maps and lists
    finescale_forces_particles_map.Remove_All();
    finescale_forces_particles_map.Set_All(async_list);
    finescale_forces_particles_map.Get_Complementary_Keys(IDENTITY_ARRAY<int>(solid_body_collection.deformable_body_collection.particles.array_collection->Size()),coarsescale_particle_indices);
    both_forces_particles_indices.Remove_All();
    for(HASHTABLE_ITERATOR<int> it(finescale_forces_particles_map);it.Valid();it.Next()){
        if(coarsescale_forces_particles_map.Contains(it.Key())) both_forces_particles_indices.Append(it.Key());}
}
//#####################################################################
// Function Set_Asynchronous_Within_Distance
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Set_Asynchronous_Within_Distance(const T distance)
{
    HASHTABLE<int> async_map;
    ARRAY<int> async_list;
    for(int j=1;j<=solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();j++){
        RIGID_BODY<TV>& body=solid_body_collection.rigid_body_collection.Rigid_Body(j);
        // choose from all particles
        for(int i=1;i<=solid_body_collection.deformable_body_collection.particles.array_collection->Size();i++)
            if(body.Implicit_Geometry_Lazy_Inside(solid_body_collection.deformable_body_collection.particles.X(i),distance))
            async_map.Set(i);}
        /* // choose from boundary-one-ring particles
        for(int i=1;i<=coarsescale_particles_pool.m;i++){int p=coarsescale_particles_pool(i);
            if(body.Implicit_Geometry_Lazy_Inside(solid_body_collection.deformable_body_collection.particles.X(p),distance))
            async_map.Set(p);}}*/
    // add boundary-one-ring particles if adaptive
    if(async_map.Size()) async_map.Set_All(coarsescale_particles_pool);
    async_map.Get_Keys(async_list);
    Set_Asynchronous_Particles(async_list);
    {std::stringstream ss;ss<<"Updated adaptive asynchronous particles number: "<<async_list.m<<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
// Function Add_Blob_From_Particles
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Add_Blob_From_Particles(const ARRAY<int>& blob_particles)
{
    {std::stringstream ss;ss<<"Adding blob of "<<blob_particles.m<<" particles"<<std::endl;LOG::filecout(ss.str());}
    blobs.Append(BLOB<TV>());
    blobs.Last().blob_particles=blob_particles;
}
//#####################################################################
// Function Compute_Blob_Properties
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Compute_Blob_Properties()
{
    for(int i=1;i<=blobs.m;i++){
        blobs(i).mass=0;
        blobs(i).inertia=T_SYMMETRIC_MATRIX();
        blobs(i).center=TV();
        for(int j=1;j<=blobs(i).blob_particles.m;j++){int p=blobs(i).blob_particles(j);
            blobs(i).mass+=solid_body_collection.deformable_body_collection.particles.mass(p);
            blobs(i).center+=solid_body_collection.deformable_body_collection.particles.mass(p)*solid_body_collection.deformable_body_collection.particles.X(p);}

        blobs(i).center/=blobs(i).mass;
        for(int j=1;j<=blobs(i).blob_particles.m;j++){int p=blobs(i).blob_particles(j);
            TV r=solid_body_collection.deformable_body_collection.particles.X(p)-blobs(i).center;
            blobs(i).inertia+=solid_body_collection.deformable_body_collection.particles.mass(p)*
                (MATRIX_POLICY<TV>::CROSS_PRODUCT_MATRIX::Cross_Product_Matrix(r).Times_Transpose(MATRIX_POLICY<TV>::CROSS_PRODUCT_MATRIX::Cross_Product_Matrix(r))).Symmetric_Part();}}
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Project(GENERALIZED_VELOCITY<TV>& V) const
{
    if(!use_projection || !asynchronous_mode==FINE_SCALE) return;
    for(int i=1;i<=blobs.m;i++){
        TV impulse;T_SPIN torque;
        for(int j=1;j<=blobs(i).blob_particles.m;j++){int p=blobs(i).blob_particles(j);
            impulse+=V.V.array(p);
            torque+=TV::Cross_Product(solid_body_collection.deformable_body_collection.particles.X(p)-blobs(i).center,V.V.array(p));}

        impulse/=blobs(i).mass;
        torque=blobs(i).inertia.Solve_Linear_System(torque);

        for(int j=1;j<=blobs(i).blob_particles.m;j++){int p=blobs(i).blob_particles(j);
            TV v_rigid_projected=solid_body_collection.deformable_body_collection.particles.mass(p)*(impulse+TV::Cross_Product(torque,solid_body_collection.deformable_body_collection.particles.X(p)-blobs(i).center));
            V.V.array(p)=((T)1-projection_rigidity)*V.V.array(p)+projection_rigidity*v_rigid_projected;}}
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void ASYNCHRONOUS_EVOLUTION<TV>::
Boundary_Conditions(GENERALIZED_VELOCITY<TV>& V) const
{
    if(!use_projection || !asynchronous_mode==FINE_SCALE) return;
    for(int i=1;i<=V.V.Size();i++) V.V(i)*=solid_body_collection.deformable_body_collection.particles.mass(i);
    Project(V);
    for(int i=1;i<=V.V.Size();i++) V.V(i)*=solid_body_collection.deformable_body_collection.particles.one_over_mass(i);
}
//#####################################################################
template class ASYNCHRONOUS_EVOLUTION<VECTOR<float,1> >;
template class ASYNCHRONOUS_EVOLUTION<VECTOR<float,2> >;
template class ASYNCHRONOUS_EVOLUTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ASYNCHRONOUS_EVOLUTION<VECTOR<double,1> >;
template class ASYNCHRONOUS_EVOLUTION<VECTOR<double,2> >;
template class ASYNCHRONOUS_EVOLUTION<VECTOR<double,3> >;
#endif
