//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BINDING_LIST
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/PARTICLE_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BINDING_LIST<TV>::
BINDING_LIST(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection)
    :particles(deformable_body_collection.particles),deformable_body_collection(&deformable_body_collection),last_read(-1),is_stale(true),frame_list(0)
{}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BINDING_LIST<TV>::
BINDING_LIST(PARTICLES<TV>& particles)
    :particles(particles),deformable_body_collection(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BINDING_LIST<TV>::
~BINDING_LIST()
{
    bindings.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Clean_Memory()
{
    bindings.Delete_Pointers_And_Clean_Memory();binding_index_from_particle_index.Clean_Memory();
}
//#####################################################################
// Function Add_Binding
//#####################################################################
template<class TV> int BINDING_LIST<TV>::
Add_Binding(BINDING<TV>* binding)
{
    bindings.Append(binding);
    if(binding_index_from_particle_index.m<binding->particle_index) binding_index_from_particle_index.Resize(binding->particle_index);
    binding_index_from_particle_index(binding->particle_index)=bindings.m;
    return bindings.m;
}
//#####################################################################
// Function Update_Binding_Index_From_Particle_Index
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Update_Binding_Index_From_Particle_Index()
{
    binding_index_from_particle_index.Clean_Memory();
    int max_particle_index=0;for(int b=1;b<=bindings.m;b++) max_particle_index=max(max_particle_index,bindings(b)->particle_index);
    binding_index_from_particle_index.Resize(max_particle_index);
    for(int b=1;b<=bindings.m;b++){assert(!binding_index_from_particle_index(bindings(b)->particle_index));binding_index_from_particle_index(bindings(b)->particle_index)=b;}
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    for(int b=1;b<=bindings.m;b++) bindings(b)->Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Compute_Dependency_Closure_Based_On_Embedding
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Compute_Dependency_Closure_Based_On_Embedding(SEGMENT_MESH& dependency_mesh) const
{
    ARRAY<VECTOR<int,2> > dependencies;dependency_mesh.elements.Exchange(dependencies);dependency_mesh.Refresh_Auxiliary_Structures();
    for(int i=1;i<=dependencies.m;i++){VECTOR<int,2>& dependency=dependencies(i);
        ARRAY<int> parents1(Dynamic_Parents(dependency.x)),parents2(Dynamic_Parents(dependency.y));
        for(int i=1;i<=parents1.m;i++) for(int j=1;j<=parents2.m;j++) dependency_mesh.Add_Element_If_Not_Already_There(VECTOR<int,2>(parents1(i),parents2(j)));}
}
//#####################################################################
// Function Compute_Particle_Closure_Based_On_Embedding
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Compute_Particle_Closure_Based_On_Embedding(ARRAY<int>& particle_set) const
{
    ARRAY<bool> particle_is_present(particles.array_collection->Size());
    INDIRECT_ARRAY<ARRAY<bool>,ARRAY<int>&> subset=particle_is_present.Subset(particle_set);
    ARRAYS_COMPUTATIONS::Fill(subset,true);
    for(int b=1;b<=bindings.m;b++){BINDING<TV>& binding=*bindings(b);
        if(!particle_is_present.Subset(binding.Parents()).Contains(false) && !particle_is_present(binding.particle_index)) particle_set.Append(binding.particle_index);}
}
//#####################################################################
// Function Clamp_Particles_To_Embedded_Positions
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Clamp_Particles_To_Embedded_Positions() const
{
    for(int i=1;i<=bindings.m;i++) bindings(i)->Clamp_To_Embedded_Position();
}
//#####################################################################
// Function Clamp_Particles_To_Embedded_Velocities
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Clamp_Particles_To_Embedded_Velocities() const
{
    for(int i=1;i<=bindings.m;i++) bindings(i)->Clamp_To_Embedded_Velocity();
}
//#####################################################################
// Function Clamp_Particles_To_Embedded_Positions
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Clamp_Particles_To_Embedded_Positions(ARRAY_VIEW<TV> X) const
{
    for(int i=1;i<=bindings.m;i++) X(bindings(i)->particle_index)=bindings(i)->Embedded_Position(X);
}
//#####################################################################
// Function Clamp_Particles_To_Embedded_Velocities
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Clamp_Particles_To_Embedded_Velocities(ARRAY_VIEW<TV> V) const
{
    for(int i=1;i<=bindings.m;i++) V(bindings(i)->particle_index)=bindings(i)->Embedded_Velocity(V);
}
//#####################################################################
// Function Clamp_Particles_To_Embedded_Velocities
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Clamp_Particles_To_Embedded_Velocities(ARRAY_VIEW<TV> V,ARRAY_VIEW<const TWIST<TV> > twist) const
{
    for(int i=1;i<=bindings.m;i++) V(bindings(i)->particle_index)=bindings(i)->Embedded_Velocity(V,twist);
}
//#####################################################################
// Function Distribute_Force_To_Parents
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full) const
{
    for(int i=1;i<=bindings.m;i++) bindings(i)->Distribute_Force_To_Parents(F_full,F_full(bindings(i)->particle_index));
}
//#####################################################################
// Function Distribute_Force_To_Parents
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,ARRAY_VIEW<TWIST<TV> > wrench_full) const
{
    for(int i=1;i<=bindings.m;i++) bindings(i)->Distribute_Force_To_Parents(F_full,wrench_full,F_full(bindings(i)->particle_index));
}
//#####################################################################
// Function Distribute_Mass_To_Parents
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Distribute_Mass_To_Parents() const
{
    for(int b=1;b<=bindings.m;b++) bindings(b)->Distribute_Mass_To_Parents(particles.mass);
}
//#####################################################################
// Function Adjust_Parents_For_Changes_In_Surface_Children
//#####################################################################
// TODO: get rid of this function by completely disallowing embedded particle modification
template<class TV> int BINDING_LIST<TV>::
Adjust_Parents_For_Changes_In_Surface_Children(const ARRAY<bool>& particle_on_surface)
{
    int interactions=0;
    HASHTABLE<int,TV> total_delta_X,total_delta_V;
    HASHTABLE<int> self_collision;
    HASHTABLE<int,int> number_of_children;

    PHYSBAM_ASSERT(deformable_body_collection);

    // figure out which embedded points are not where they should be or have differing velocities
    for(int i=1;i<=bindings.m;i++){
        BINDING<TV>& binding=*bindings(i);
        int p=binding.particle_index;
        const ARRAY<int>& parents=binding.Parents();
        if(particle_on_surface(p)){
            TV Xp=binding.Embedded_Position(),Vp=binding.Embedded_Velocity();
            if(particles.X(p)!=Xp || particles.V(p)!=Vp){interactions++;
                TV delta_X=particles.X(p)-Xp,delta_V=particles.V(p)-Vp;
                for(int j=1;j<=parents.m;j++) if(!particle_on_surface(parents(j))){
                        self_collision.Set(parents(j));total_delta_X.Get_Or_Insert(parents(j))+=delta_X;total_delta_V.Get_Or_Insert(parents(j))+=delta_V;}}}
        for(int j=1;j<=parents.m;j++) number_of_children.Get_Or_Insert(parents(j))++;} // TODO: consider only counting children that are on the surface

    // adjust parents to match changes in children, skipping parents on the surface since they may also be self-colliding
    for(typename HASHTABLE<int>::ITERATOR it(self_collision);it.Valid();it.Next()){int p=it.Key();
        T one_over_number_of_children=(T)1/number_of_children.Get(p);
        particles.X(p)+=total_delta_X.Get(p)*one_over_number_of_children;
        particles.V(p)+=total_delta_V.Get(p)*one_over_number_of_children;}

    return interactions; // self-collision/etc. may require bindings to drift, so do not resynchronize
}
//#####################################################################
// Function Adjust_Parents_For_Changes_In_Surface_Children_Velocities
//#####################################################################
// TODO: get rid of this function by completely disallowing embedded particle modification
template<class TV> int BINDING_LIST<TV>::
Adjust_Parents_For_Changes_In_Surface_Children_Velocities(const ARRAY<bool>& particle_on_surface)
{
    int interactions=0;
    HASHTABLE<int,TV> total_delta_V;
    HASHTABLE<int> self_collision;
    HASHTABLE<int,int> number_of_children;

    PHYSBAM_ASSERT(deformable_body_collection);

    // figure out which embedded points are not where they should be or have differing velocities
    for(int i=1;i<=bindings.m;i++){
        BINDING<TV>& binding=*bindings(i);
        int p=binding.particle_index;
        const ARRAY<int>& parents=binding.Parents();
        if(particle_on_surface(p)){
            TV Vp=binding.Embedded_Velocity();
            if(particles.V(p)!=Vp){interactions++;
                TV delta_V=particles.V(p)-Vp;
                for(int j=1;j<=parents.m;j++) if(!particle_on_surface(parents(j))){
                        self_collision.Set(parents(j));total_delta_V.Get_Or_Insert(parents(j))+=delta_V;}}}
        for(int j=1;j<=parents.m;j++) number_of_children.Get_Or_Insert(parents(j))++;} // TODO: consider only counting children that are on the surface

    // adjust parents to match changes in children, skipping parents on the surface since they may also be self-colliding
    for(typename HASHTABLE<int>::ITERATOR it(self_collision);it.Valid();it.Next()){int p=it.Key();
        T one_over_number_of_children=(T)1/number_of_children.Get(p);
        particles.V(p)+=total_delta_V.Get(p)*one_over_number_of_children;}

    return interactions; // self-collision/etc. may require bindings to drift, so do not resynchronize
}
//#####################################################################
// Function Read
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Read(TYPED_ISTREAM& input)
{
    Clean_Memory();
    int m;Read_Binary(input,m);bindings.Preallocate(m);
    for(int k=1;k<=m;k++){BINDING<TV>* binding=BINDING<TV>::Create(input,particles);
        Add_Binding(binding);}
}
//#####################################################################
// Function Write
//#####################################################################
template<class TV> void BINDING_LIST<TV>::
Write(TYPED_OSTREAM& output) const
{
    Write_Binary(output,bindings.m);
    for(int k=1;k<=bindings.m;k++) bindings(k)->Write(output);
}
//#####################################################################
template class BINDING_LIST<VECTOR<float,1> >;
template class BINDING_LIST<VECTOR<float,2> >;
template class BINDING_LIST<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BINDING_LIST<VECTOR<double,1> >;
template class BINDING_LIST<VECTOR<double,2> >;
template class BINDING_LIST<VECTOR<double,3> >;
#endif
