//#####################################################################
// Copyright 2006-2008, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOFT_BINDINGS
//#####################################################################
#ifndef __SOFT_BINDINGS__
#define __SOFT_BINDINGS__

#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
namespace PhysBAM{

template<class TV> class BINDING;
template<class TV> class PARTICLES;

template<class TV>
class SOFT_BINDINGS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    typedef int HAS_TYPED_READ_WRITE;

    BINDING_LIST<TV>& binding_list;
    PARTICLES<TV>& particles;
    ARRAY<VECTOR<int,2> > bindings; // (particle_index,parent) pairs
    ARRAY<bool> use_impulses_for_collisions;
    ARRAY<int> bindings_using_impulses_for_collisions;
    ARRAY<int> bindings_using_forces_for_collisions;
    ARRAY<int> binding_index_from_particle_index;
    SEGMENT_MESH* binding_mesh;
    bool use_gauss_seidel_for_impulse_based_collisions; // TODO: set per binding and add to I/O
    int last_read;
    mutable bool is_stale;
    mutable ARRAY<int>* frame_list;
public:
    SOFT_BINDINGS(BINDING_LIST<TV>& binding_list_input);

    virtual ~SOFT_BINDINGS();

    void Clean_Memory()
    {bindings.Clean_Memory();use_impulses_for_collisions.Clean_Memory();
    bindings_using_forces_for_collisions.Clean_Memory();bindings_using_impulses_for_collisions.Clean_Memory();
    binding_index_from_particle_index.Clean_Memory();delete binding_mesh;binding_mesh=0;}

    int Add_Binding(const VECTOR<int,2>& binding,const bool use_impulses_for_collisions_input)
    {bindings.Append(binding);use_impulses_for_collisions.Append(use_impulses_for_collisions_input);
    if(binding_index_from_particle_index.m<binding.x) binding_index_from_particle_index.Resize(binding.x);
    binding_index_from_particle_index(binding.x)=bindings.m;
    return bindings.m;}

    bool Particle_Is_Bound(const int particle_index) const
    {return particle_index<=binding_index_from_particle_index.m && binding_index_from_particle_index(particle_index);}

    T One_Over_Effective_Mass(const int particle_index) const
    {if(!Particle_Is_Bound(particle_index)) return binding_list.One_Over_Effective_Mass(particle_index);
    return binding_list.One_Over_Effective_Mass(bindings(binding_index_from_particle_index(particle_index)).y);}

    int Parent(const int particle_index) const
    {return Particle_Is_Bound(particle_index)?bindings(binding_index_from_particle_index(particle_index)).y:0;}

    int Driftless_Particle(const int particle_index) const
    {return Particle_Is_Bound(particle_index)?bindings(binding_index_from_particle_index(particle_index)).y:particle_index;}

    int Soft_Binding(const int particle_index) const
    {return particle_index<=binding_index_from_particle_index.m?binding_index_from_particle_index(particle_index):0;}

    BINDING<TV>* Hard_Binding(const int particle_index) const
    {return Particle_Is_Bound(particle_index)?binding_list.Binding(Parent(particle_index)):0;}

    ARRAY<int> Parents(const int particle_index) const
    {if(!Particle_Is_Bound(particle_index)) return ARRAY<int>();
    int parent=bindings(binding_index_from_particle_index(particle_index)).y;
    ARRAY<int> parents=binding_list.Parents(parent);parents.Append(parent);
    return parents;}

//#####################################################################
    void Initialize_Binding_Mesh(const bool exclude_particles_using_impulses_for_collisions=false);
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const;
    void Set_Mass_From_Effective_Mass();
    void Remove_Soft_Bound_Particles(ARRAY<int>& particles) const;
    int Adjust_Parents_For_Changes_In_Surface_Children(const ARRAY<bool>& particle_on_surface); // TODO: names are needlessly collision specific
    int Adjust_Parents_For_Changes_In_Surface_Children_Velocities(const ARRAY<bool>& particle_on_surface); // same as above, but touches only velocities
    bool Need_Bindings_Mapped() const;
    void Map_Forces_From_Parents(ARRAY_VIEW<TV> F_full,ARRAY_VIEW<const TWIST<TV> > wrench_full) const;
    void Clamp_Particles_To_Embedded_Positions(const bool bindings_using_impulses_for_collisions_only=false) const;
    void Clamp_Particles_To_Embedded_Velocities(const bool bindings_using_impulses_for_collisions_only=false) const;
    void Update_Binding_Index_From_Particle_Index();
    void Read(TYPED_ISTREAM& input);
    void Write(TYPED_OSTREAM& output) const;
//#####################################################################
};
}
#endif
