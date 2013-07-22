//#####################################################################
// Copyright 2006, Ronald Fedkiw, Geoffrey Irving, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BINDING
//##################################################################### 
#ifndef __BINDING__
#define __BINDING__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Geometry/Registry/REGISTRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
namespace PhysBAM{

template<class TV> class PARTICLES;
class SEGMENT_MESH;
template<class ID> class SPARSE_UNION_FIND;

template<class TV>
class BINDING:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
public:
    typedef TV VECTOR_T;

    PARTICLES<TV>& particles;
    int particle_index;

    BINDING(PARTICLES<TV>& particles_input)
        :particles(particles_input),particle_index(0)
    {}

    BINDING(PARTICLES<TV>& particles_input,const int particle_index_input)
        :particles(particles_input),particle_index(particle_index_input)
    {}

    virtual ~BINDING()
    {}

    virtual T One_Over_Effective_Mass(const TV& direction) const // by default return the direction-independent effective mass
    {return One_Over_Effective_Mass();}

    void Clamp_To_Embedded_Position(ARRAY_VIEW<TV> X)
    {X(particle_index)=Embedded_Position(X);}

    void Clamp_To_Embedded_Position()
    {particles.X(particle_index)=Embedded_Position();}

    void Clamp_To_Embedded_Velocity()
    {particles.V(particle_index)=Embedded_Velocity();}

protected:    
    virtual void Read_Helper(TYPED_ISTREAM& input)
    {Read_Binary(input,particle_index);}

    virtual void Write_Helper(TYPED_OSTREAM& output) const
    {Write_Binary(output,particle_index);}
public:

//#####################################################################
    virtual TV Embedded_Position() const=0;
    virtual TV Embedded_Position(ARRAY_VIEW<const TV> X) const=0;
    virtual TV Embedded_Velocity() const=0;
    virtual TV Embedded_Velocity(ARRAY_VIEW<const TV> V) const=0;
    virtual TV Embedded_Velocity(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > twist) const=0;
    virtual TV Embedded_Acceleration(ARRAY_VIEW<const TV> F,ARRAY_VIEW<const TWIST<TV> > wrench) const=0;
    virtual T One_Over_Effective_Mass() const=0; // returns a lower bound when effective mass is direction dependent
    virtual void Apply_Impulse(const TV& impulse)=0;
    virtual void Apply_Displacement_To_Parents_Based_On_Embedding(const TV& dX,const ARRAY<bool>* skip_particle)=0; // skip_particle.m == particles.array_collection->Size()
    virtual void Apply_Velocity_Change_To_Parents_Based_On_Embedding(const TV& dV,const ARRAY<bool>* skip_particle)=0; // skip_particle.m == particles.array_collection->Size()
    virtual void Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,const TV& force) const=0;
    virtual void Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,ARRAY_VIEW<TWIST<TV> > wrench_full,const TV& force) const=0;
    virtual void Distribute_Mass_To_Parents(ARRAY_VIEW<T> mass_full) const=0;
    virtual ARRAY<int> Parents() const=0;
    virtual ARRAY<T> Weights() const=0;
    virtual void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const;
    // read/write
    virtual int Name() const {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return Static_Name();}
    static int Static_Name() {return -1;}
    virtual std::string Extension() const {PHYSBAM_WARN_IF_NOT_OVERRIDDEN();return Static_Extension();}
    static std::string Static_Extension() {return "";}
    static BINDING* Create(TYPED_ISTREAM& input,PARTICLES<TV>& particles);
    void Write(TYPED_OSTREAM& output) const;
private:
    static BINDING* Create_From_Name(const int name,PARTICLES<TV>& particles);
//#####################################################################
};

//#####################################################################
// Class BINDING_REGISTRY
//#####################################################################
template<class TV> class BINDING_REGISTRY;

template<class TV>
class BINDING_REGISTRY:public REGISTRY<BINDING<TV>,int,BINDING_REGISTRY<TV> >
{
public:
    template<class T_OBJECT> static T_OBJECT* Create_Representative()
    {static PARTICLES<TV> particles;return T_OBJECT::Create(particles);}
};
}
#endif
