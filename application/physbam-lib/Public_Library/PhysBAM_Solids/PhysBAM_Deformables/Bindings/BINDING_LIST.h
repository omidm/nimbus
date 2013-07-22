//#####################################################################
// Copyright 2006-2008, Ron Fedkiw, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BINDING_LIST
//#####################################################################
#ifndef __BINDING_LIST__
#define __BINDING_LIST__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLES_PARTICLES_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
namespace PhysBAM{

template<class TV> class DEFORMABLE_BODY_COLLECTION;

template<class TV>
class BINDING_LIST:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    typedef int HAS_TYPED_READ_WRITE;

    PARTICLES<TV>& particles;
    ARRAY<BINDING<TV>*> bindings;
    ARRAY<int> binding_index_from_particle_index;
    DEFORMABLE_BODY_COLLECTION<TV>* deformable_body_collection;
    int last_read;
    mutable bool is_stale;
    mutable ARRAY<int>* frame_list;

    BINDING_LIST(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection);
    BINDING_LIST(PARTICLES<TV>& particles); // use of this constructor disables some features
    virtual ~BINDING_LIST();

    int Binding_Index_From_Particle_Index(const int particle) const
    {assert(particle<=particles.array_collection->Size());return binding_index_from_particle_index.m<particle?0:binding_index_from_particle_index(particle);}

    TV Embedded_Position(const int particle_index) const
    {return bindings(binding_index_from_particle_index(particle_index))->Embedded_Position();}

    TV Embedded_Position(const int particle_index,ARRAY_VIEW<const TV> X_input) const
    {return bindings(binding_index_from_particle_index(particle_index))->Embedded_Position(X_input);}

    TV Embedded_Velocity(const int particle_index) const
    {return bindings(binding_index_from_particle_index(particle_index))->Embedded_Velocity();}

    BINDING<TV>* Binding(const int particle_index) const
    {if(binding_index_from_particle_index.m<particle_index || !binding_index_from_particle_index(particle_index)) return 0;
    return bindings(binding_index_from_particle_index(particle_index));}

    TV V(const int particle_index) const
    {if(BINDING<TV>* binding=Binding(particle_index)) return binding->Embedded_Velocity();return particles.V(particle_index);}

    void Apply_Impulse(const int particle_index,const TV& impulse)
    {if(BINDING<TV>* binding=Binding(particle_index)) binding->Apply_Impulse(impulse);else particles.V(particle_index)+=particles.one_over_mass(particle_index)*impulse;}

    T One_Over_Effective_Mass(const int particle_index) const
    {if(BINDING<TV>* binding=Binding(particle_index)) return binding->One_Over_Effective_Mass();return particles.one_over_mass(particle_index);}

    T One_Over_Effective_Mass(const int particle_index,const TV& direction) const
    {if(BINDING<TV>* binding=Binding(particle_index)) return binding->One_Over_Effective_Mass(direction);return particles.one_over_mass(particle_index);}

    ARRAY<int> Parents(const int particle_index) const
    {if(BINDING<TV>* binding=Binding(particle_index)) return binding->Parents();return ARRAY<int>();}

    ARRAY<int> Dynamic_Parents(const int particle_index) const
    {ARRAY<int> parents(Parents(particle_index));if(!parents.m) parents.Append(particle_index);return parents;}

    template<class T_ARRAY>
    void Clear_Hard_Bound_Particles(T_ARRAY& array) const
    {for(int b=1;b<=bindings.m;b++) array(bindings(b)->particle_index)=typename T_ARRAY::ELEMENT();}

//#####################################################################
    void Clean_Memory();
    int Add_Binding(BINDING<TV>* binding);
    void Update_Binding_Index_From_Particle_Index();
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const;
    void Compute_Dependency_Closure_Based_On_Embedding(SEGMENT_MESH& dependency_mesh) const;
    void Compute_Particle_Closure_Based_On_Embedding(ARRAY<int>& particle_set) const;
    void Clamp_Particles_To_Embedded_Positions() const;
    void Clamp_Particles_To_Embedded_Velocities() const;
    void Clamp_Particles_To_Embedded_Positions(ARRAY_VIEW<TV> X) const;
    void Clamp_Particles_To_Embedded_Velocities(ARRAY_VIEW<TV> V) const;
    void Clamp_Particles_To_Embedded_Velocities(ARRAY_VIEW<TV> V,ARRAY_VIEW<const TWIST<TV> > twist) const;
    void Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full) const;
    void Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,ARRAY_VIEW<TWIST<TV> > wrench_full) const;
    void Distribute_Mass_To_Parents() const;
    int Adjust_Parents_For_Changes_In_Surface_Children(const ARRAY<bool>& particle_on_surface); // TODO: names are needlessly collision specific
    int Adjust_Parents_For_Changes_In_Surface_Children_Velocities(const ARRAY<bool>& particle_on_surface); // same as above, but touches only velocities
    void Read(TYPED_ISTREAM& input);
    void Write(TYPED_OSTREAM& output) const;
//#####################################################################
};
}
#endif
