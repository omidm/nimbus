//#####################################################################
// Copyright 2006-2009, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_MATRIX_0X0.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_MATRIX_1X1.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_1D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGIDS_FORCES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <climits>
namespace PhysBAM{
template<class TV>
struct ALLOCATE_BODY_HELPER:public ALLOCATE_HELPER<TV>
{
    RIGID_BODY_COLLECTION<TV>& collection;
    ALLOCATE_BODY_HELPER(RIGID_BODY_COLLECTION<TV>& collection_input):collection(collection_input) {}
    RIGID_BODY<TV>* Create(int index=0) PHYSBAM_OVERRIDE {return new RIGID_BODY<TV>(collection,true,index);}
    virtual ~ALLOCATE_BODY_HELPER(){}
};
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_BODY_COLLECTION<TV>::
RIGID_BODY_COLLECTION(RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>* rigids_example_forces_and_velocities_input,COLLISION_GEOMETRY_COLLECTION<TV>* collision_body_list_input,
    ARRAY_COLLECTION* array_collection)
    :rigid_body_particle(array_collection?array_collection:new ARRAY_COLLECTION),
    rigid_geometry_collection(rigid_body_particle,rigids_example_forces_and_velocities_input,collision_body_list_input,new ALLOCATE_BODY_HELPER<TV>(*this)),
    articulated_rigid_body(*new ARTICULATED_RIGID_BODY<TV>(*this)),rigid_body_cluster_bindings(*new RIGID_BODY_CLUSTER_BINDINGS<TV>(*this,articulated_rigid_body)),
    rigids_example_forces_and_velocities(rigids_example_forces_and_velocities_input),dynamic_rigid_body_particles(0),static_rigid_bodies(rigid_geometry_collection.static_rigid_geometry),
    kinematic_rigid_bodies(rigid_geometry_collection.kinematic_rigid_geometry),static_and_kinematic_rigid_bodies(0),print_diagnostics(false),print_residuals(false),print_energy(false),iterations_used_diagnostic(0)
{rigid_body_particle.Store_Velocity();}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_BODY_COLLECTION<TV>::
~RIGID_BODY_COLLECTION()
{
    rigids_forces.Delete_Pointers_And_Clean_Memory();
    delete &articulated_rigid_body;
    delete &rigid_body_cluster_bindings;
}
//#####################################################################
// Function Rigid_Body
//#####################################################################
template<class TV> RIGID_BODY<TV>& RIGID_BODY_COLLECTION<TV>::
Rigid_Body(const int particle_index)
{
    return *debug_cast<RIGID_BODY<TV>*>(rigid_body_particle.rigid_geometry(particle_index));
}
//#####################################################################
// Function Rigid_Body
//#####################################################################
template<class TV> const RIGID_BODY<TV>& RIGID_BODY_COLLECTION<TV>::
Rigid_Body(const int particle_index) const
{
    return *debug_cast<RIGID_BODY<TV>*>(rigid_body_particle.rigid_geometry(particle_index));
}
//#####################################################################
// Function Add_Rigid_Body
//#####################################################################
// structures already added to their respective lists
template<class TV> int RIGID_BODY_COLLECTION<TV>::
Add_Rigid_Body(RIGID_BODY<TV>* rigid_body,const int simplicial_boundary_id,const int implicit_object_id,const int simplicial_interior_id)
{
    int id=rigid_body->particle_index;
    if(simplicial_boundary_id) rigid_body_particle.structure_ids(id)(1)=simplicial_boundary_id;
    if(implicit_object_id) rigid_body_particle.structure_ids(id)(2)=implicit_object_id;
    if(simplicial_interior_id) rigid_body_particle.structure_ids(id)(3)=simplicial_interior_id;
    for(int i=1;i<=rigid_body_particle.structure_ids(id).m;i++) if(rigid_body_particle.structure_ids(id)(i) && !rigid_geometry_collection.structure_list.Element(rigid_body_particle.structure_ids(id)(i))) PHYSBAM_FATAL_ERROR();
    return id;
}
//#####################################################################
// Function Add_Rigid_Body
//#####################################################################
// adds body's segmented curve and implicit curve to their respective lists
template<class TV> int RIGID_BODY_COLLECTION<TV>::
Add_Rigid_Body_And_Geometry(RIGID_BODY<TV>* rigid_body)
{
    return rigid_geometry_collection.Add_Rigid_Geometry(rigid_body);
}
//#####################################################################
// Function Add_Rigid_Body
//#####################################################################
template<class TV> int RIGID_BODY_COLLECTION<TV>::
Add_Rigid_Body(const STREAM_TYPE stream_type,const std::string& basename,const T scaling_factor,const bool read_simplicial_boundary,const bool read_implicit_object,
    const bool read_simplicial_interior,const bool read_rgd_file)
{
    return Add_Rigid_Body(stream_type,false,basename,scaling_factor,read_simplicial_boundary,read_implicit_object,read_simplicial_interior,read_rgd_file);
}
//#####################################################################
// Function Add_Rigid_Body
//#####################################################################
template<class TV> int RIGID_BODY_COLLECTION<TV>::
Add_Rigid_Body(const STREAM_TYPE stream_type,const bool thin_shell,const std::string& basename,const T scaling_factor,const bool read_simplicial_boundary,const bool read_implicit_object,
    const bool read_simplicial_interior,const bool read_rgd_file)
{
    RIGID_BODY<TV>* rigid_body=new RIGID_BODY<TV>(*this,true);
    rigid_body->thin_shell=thin_shell;

    // rigid body
    std::string rgd=TV::dimension==2?"rgd2d":"rgd";
    if(read_rgd_file){
        try{FILE_UTILITIES::Read_From_File(stream_type,basename+"."+rgd,rigid_body->Mass(),rigid_body->Inertia_Tensor(),rigid_body->X(),rigid_body->Rotation());}
        catch(FILESYSTEM_ERROR&){{std::stringstream ss;ss<<"Note: No "<<rgd<<" file for "<<basename<<" (using default values)"<<std::endl;LOG::filecout(ss.str());}}}
    if(scaling_factor!=1) rigid_body->Rescale(scaling_factor);
    rigid_body->Update_Angular_Velocity();

    int id=rigid_geometry_collection.Add_Rigid_Geometry(rigid_body,stream_type,basename,scaling_factor,read_simplicial_boundary,read_implicit_object,read_simplicial_interior,read_rgd_file);
    return id;
}
//#####################################################################
// Function Reset_Impulse_Accumulators
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Reset_Impulse_Accumulators()
{
    for(int i=1;i<=rigid_body_particle.array_collection->Size();i++)
        if(Is_Active(i) && rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Contains(i) && rigid_geometry_collection.collision_body_list->bodies(rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(i))->impulse_accumulator)
            rigid_geometry_collection.collision_body_list->bodies(rigid_geometry_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(i))->impulse_accumulator->Reset();
}
//#####################################################################
// Function Update_Angular_Velocity
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Update_Angular_Velocity()
{
    for(int p=1;p<=rigid_body_particle.array_collection->Size();p++) if(Is_Active(p)) Rigid_Body(p).Update_Angular_Velocity();
}
//#####################################################################
// Function Update_Angular_Momentum
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Update_Angular_Momentum()
{
    for(int p=1;p<=rigid_body_particle.array_collection->Size();p++) if(Is_Active(p)) Rigid_Body(p).Update_Angular_Momentum();
}
//#####################################################################
// Function Update_Angular_Velocity
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Update_Angular_Velocity(const ARRAY<int>& rigid_body_particles)
{
    for(int i=1;i<=rigid_body_particles.m;i++) Rigid_Body(rigid_body_particles(i)).Update_Angular_Velocity();
}
//#####################################################################
// Function Update_Angular_Momentum
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Update_Angular_Momentum(const ARRAY<int>& rigid_body_particles)
{
    for(int i=1;i<=rigid_body_particles.m;i++) Rigid_Body(rigid_body_particles(i)).Update_Angular_Momentum();
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Read(const STREAM_TYPE stream_type,const std::string& directory,const int frame,ARRAY<int>* needs_init,ARRAY<int>* needs_destroy)
{
    if(!stream_type.use_doubles)
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,float>::Read(stream_type,directory,frame,rigid_geometry_collection,needs_init,needs_destroy);
    else
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,double>::Read(stream_type,directory,frame,rigid_geometry_collection,needs_init,needs_destroy);
#else
        PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif

}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Write(const STREAM_TYPE stream_type,const std::string& directory,const int frame) const
{
    if(!stream_type.use_doubles)
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,float>::Write(stream_type,directory,frame,rigid_geometry_collection);
    else
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,double>::Write(stream_type,directory,frame,rigid_geometry_collection);
#else
        PHYSBAM_FATAL_ERROR("Cannot read doubles");
#endif
    articulated_rigid_body.Write(stream_type,directory,frame);
}
//#####################################################################
// Function Update_Simulated_Particles
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Update_Simulated_Particles()
{
    int rigid_particles_number=rigid_body_particle.array_collection->Size();

    ARRAY<bool> particle_is_simulated(rigid_particles_number);
    INDIRECT_ARRAY<ARRAY<bool>,ARRAY<int>&> simulated_subset=particle_is_simulated.Subset(rigid_body_particle.array_collection->deletion_list);
    ARRAYS_COMPUTATIONS::Fill(simulated_subset,false);
    for(int i=1;i<=rigid_particles_number;i++)
        if(Is_Active(i) && Rigid_Body(i).Is_Simulated()) // TODO: Can't everything be defaulted to true?
            particle_is_simulated(i)=true;
    rigids_example_forces_and_velocities->Set_Rigid_Particle_Is_Simulated(particle_is_simulated);
    for(int i=1;i<=rigid_particles_number;i++)
        if(!Is_Active(i) || !Rigid_Body(i).Is_Simulated())
            particle_is_simulated(i)=false;

    simulated_rigid_body_particles.Remove_All();
    dynamic_rigid_body_particles.Remove_All();

    for(int p=1;p<=rigid_particles_number;p++) if(particle_is_simulated(p)) simulated_rigid_body_particles.Append(p);

    rigid_body_cluster_bindings.Clear_Hard_Bound_Particles(particle_is_simulated);

    for(int p=1;p<=rigid_particles_number;p++) if(particle_is_simulated(p)) dynamic_rigid_body_particles.Append(p);

    static_rigid_bodies.Remove_All();kinematic_rigid_bodies.Remove_All();static_and_kinematic_rigid_bodies.Remove_All();
    for(int p=1;p<=rigid_particles_number;p++) if(Is_Active(p)){RIGID_BODY<TV>& rigid_body=Rigid_Body(p);
        if(rigid_body.is_static){static_rigid_bodies.Append(p);static_and_kinematic_rigid_bodies.Append(p);}
        if(rigid_body_particle.kinematic(p)){kinematic_rigid_bodies.Append(p);static_and_kinematic_rigid_bodies.Append(p);}}

    ARRAY<bool> rigid_particle_is_simulated(rigid_particles_number);
    INDIRECT_ARRAY<ARRAY<bool>,ARRAY<int>&> simulated_particle_subset=rigid_particle_is_simulated.Subset(simulated_rigid_body_particles);
    ARRAYS_COMPUTATIONS::Fill(simulated_particle_subset,true);
    for(int i=1;i<=rigids_forces.m;i++) rigids_forces(i)->Update_Mpi(rigid_particle_is_simulated);
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const
{
    for(int k=1;k<=rigids_forces.m;k++)
        if(rigids_forces(k)->use_velocity_independent_forces) rigids_forces(k)->Add_Velocity_Independent_Forces(rigid_F_full,time);
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
// can depend on position too
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T time) const
{
    for(int k=1;k<=rigids_forces.m;k++)
        if(rigids_forces(k)->use_velocity_dependent_forces) rigids_forces(k)->Add_Velocity_Dependent_Forces(rigid_V_full,rigid_F_full,time);
}
//#####################################################################
// Function Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V_full,ARRAY_VIEW<TWIST<TV> > rigid_F_full,const T scale,const T time) const
{
    assert(rigid_F_full.Size()==rigid_body_particle.array_collection->Size());
    INDIRECT_ARRAY<ARRAY_VIEW<TWIST<TV> > > rigid_F_subset(rigid_F_full.Subset(dynamic_rigid_body_particles));
    ARRAYS_COMPUTATIONS::Fill(rigid_F_subset,TWIST<TV>()); // note we zero here because we will scale the forces below
    bool added=false;
    for(int k=1;k<=rigids_forces.m;k++) if(rigids_forces(k)->use_implicit_velocity_independent_forces){
        rigids_forces(k)->Add_Implicit_Velocity_Independent_Forces(rigid_V_full,rigid_F_full,time);added=true;}
    if(added) rigid_F_full.Subset(simulated_rigid_body_particles)*=scale;
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Update_Position_Based_State(const T time)
{
    for(int k=1;k<=rigids_forces.m;k++)
        if(rigids_forces(k)->use_position_based_state) rigids_forces(k)->Update_Position_Based_State(time);
}
//#####################################################################
// Function Compute_Energy
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Compute_Energy(const T time,T& kinetic_energy,T& potential_energy) const
{
    potential_energy=0;
    kinetic_energy=0;
    for(int i=1;i<=rigids_forces.m;i++) potential_energy+=rigids_forces(i)->Potential_Energy(time);
    for(int i=1;i<=dynamic_rigid_body_particles.m;i++){int p=dynamic_rigid_body_particles(i);
        kinetic_energy+=Rigid_Body(p).Kinetic_Energy();}
}
//#####################################################################
// Function Print_Energy
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Print_Energy(const T time,const int step) const
{
    if(print_energy){
        T potential_energy=0,kinetic_energy=0;
        Compute_Energy(time,kinetic_energy,potential_energy);
        {std::stringstream ss;ss<<"total energy = "<<(potential_energy+kinetic_energy)<<"    (KE = "<<kinetic_energy<<"   PE = "<<potential_energy<<")  Step "<<step<<std::endl;LOG::filecout(ss.str());}}
}
//#####################################################################
// Function CFL_Rigid
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_BODY_COLLECTION<TV>::
CFL_Rigid(const RIGID_BODY_EVOLUTION_PARAMETERS<TV>& rigid_body_evolution_parameters,const bool verbose_dt)
{
    static T static_min_bounding_box_width=FLT_MAX;
    T min_bounding_box_width=FLT_MAX;
    for(int i(1);i<=rigid_body_particle.array_collection->Size();i++) if(Is_Active(i)){
            const RANGE<TV>& box=Rigid_Body(i).Object_Space_Bounding_Box();
            TV edge_lengths=box.Edge_Lengths();min_bounding_box_width=min(min_bounding_box_width,edge_lengths.Min());}
    if(min_bounding_box_width!=static_min_bounding_box_width){
        static_min_bounding_box_width=min_bounding_box_width;
        LOG::Stat("minimum rigid body bounding box width",min_bounding_box_width);}

    T max_distance_per_time_step=rigid_body_evolution_parameters.max_rigid_body_linear_movement_fraction_per_time_step*min_bounding_box_width;
    T dt=FLT_MAX;
    bool no_active_bodies=true;
    for(int p=1;p<=rigid_body_particle.array_collection->Size();p++) if(Is_Active(p)){
        dt=min(dt,Rigid_Body(p).CFL(max_distance_per_time_step,rigid_body_evolution_parameters.max_rigid_body_rotation_per_time_step,verbose_dt));
        no_active_bodies=false;}
    if(no_active_bodies) return FLT_MAX; // don't apply rigid dt bounds if there aren't any active rigid bodies
    dt=Robust_Multiply(rigid_body_evolution_parameters.rigid_cfl,dt);
    T dt_clamped=clamp(dt,rigid_body_evolution_parameters.rigid_minimum_dt,rigid_body_evolution_parameters.rigid_maximum_dt);
    if(dt_clamped>dt && verbose_dt) {std::stringstream ss;ss<<"Warning: taking larger time step ("<<dt_clamped<<") than CFL dt ("<<dt<<")"<<std::endl;LOG::filecout(ss.str());}
    return dt_clamped;
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class TV> int RIGID_BODY_COLLECTION<TV>::
Add_Force(RIGIDS_FORCES<TV>* force)
{
    rigids_forces.Append(force);
    force->Set_CFL_Number((T).5);
    return rigids_forces.m;
}
//#####################################################################
// Function Save_Potential_Energy
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Save_Potential_Energy(const T time)
{
    for(int i=1;i<=rigids_forces.m;i++) rigids_forces(i)->Save_Potential_Energy(time);
}
//#####################################################################
// Function Compute_Energy_Error
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Compute_Energy_Error(ARRAY_VIEW<const TWIST<TV> > rigid_velocity_save,const T time,const T dt)
{
    for(int i=1;i<=rigids_forces.m;i++) rigids_forces(i)->Compute_Energy_Error(rigid_velocity_save,time,dt);
}
//#####################################################################
// Function Add_Energy_Correction_Force
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Add_Energy_Correction_Force(ARRAY_VIEW<const TWIST<TV> > rigid_velocity_save,const T time,const T dt)
{
    for(int i=1;i<=rigids_forces.m;i++) rigids_forces(i)->Add_Energy_Correction_Force(rigid_velocity_save,time,dt);
}
//#####################################################################
// Function Compute_Previously_Applied_Forces
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Compute_Previously_Applied_Forces()
{
    for(int i=1;i<=rigids_forces.m;i++) rigids_forces(i)->Compute_Previously_Applied_Forces();
}
//#####################################################################
// Function Store_Velocities
//#####################################################################
template<class TV> void RIGID_BODY_COLLECTION<TV>::
Store_Velocities()
{
    for(int i=1;i<=rigids_forces.m;i++) rigids_forces(i)->Store_Velocities();
}
//#####################################################################
template class RIGID_BODY_COLLECTION<VECTOR<float,1> >;
template class RIGID_BODY_COLLECTION<VECTOR<float,2> >;
template class RIGID_BODY_COLLECTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_BODY_COLLECTION<VECTOR<double,1> >;
template class RIGID_BODY_COLLECTION<VECTOR<double,2> >;
template class RIGID_BODY_COLLECTION<VECTOR<double,3> >;
#endif
}
