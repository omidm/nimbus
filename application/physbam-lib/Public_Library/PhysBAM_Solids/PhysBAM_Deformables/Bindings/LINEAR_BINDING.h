//#####################################################################
// Copyright 2006-2008, Ronald Fedkiw, Geoffrey Irving, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_BINDING
//##################################################################### 
#ifndef __LINEAR_BINDING__
#define __LINEAR_BINDING__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING.h>
namespace PhysBAM{

template<class TV,int d>
class LINEAR_BINDING:public BINDING<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef BINDING<TV> BASE;
    using BASE::particles;using BASE::particle_index;
    using BASE::One_Over_Effective_Mass; // silence -Woverloaded-virtual

    VECTOR<int,d> parents;
    VECTOR<T,d> weights; // weights should sum to 1

    LINEAR_BINDING(PARTICLES<TV>& particles_input)
        :BINDING<TV>(particles_input)
    {}
    
    LINEAR_BINDING(PARTICLES<TV>& particles_input,const int particle_index_input,const VECTOR<int,d>& parents_input,const VECTOR<T,d>& weights_input)
        :BINDING<TV>(particles_input,particle_index_input),parents(parents_input),weights(weights_input)
    {}

    LINEAR_BINDING(PARTICLES<TV>& particles_input,const int particle_index_input,const VECTOR<int,d>& parents_input,const VECTOR<T,d-1>& weights_input)
        :BINDING<TV>(particles_input,particle_index_input),parents(parents_input),weights(weights_input)
    {
        STATIC_ASSERT(d>2); // this would be confusing in the segment case because interpolation fractions refer to the other vertex
        weights[d]=1-weights_input.Sum(); 
    }

    static LINEAR_BINDING* Create(GEOMETRY_PARTICLES<TV>& particles)
    {return new LINEAR_BINDING(dynamic_cast<PARTICLES<TV>&>(particles));}

    virtual int Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    static int Static_Name()
    {return 2+(d<<16);}

    TV Embedded_Position() const PHYSBAM_OVERRIDE
    {TV result;for(int i=1;i<=d;i++) result+=weights[i]*particles.X(parents[i]);
    return result;}

    TV Embedded_Position(ARRAY_VIEW<const TV> X) const PHYSBAM_OVERRIDE
    {TV result;for(int i=1;i<=d;i++) result+=weights[i]*X(parents[i]);
    return result;}

    TV Embedded_Velocity() const  PHYSBAM_OVERRIDE
    {TV result;for(int i=1;i<=d;i++) result+=weights[i]*particles.V(parents[i]);
    return result;}

    TV Embedded_Velocity(ARRAY_VIEW<const TV> V) const  PHYSBAM_OVERRIDE
    {TV result;for(int i=1;i<=d;i++) result+=weights[i]*V(parents[i]);
    return result;}

    TV Embedded_Velocity(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > twist) const  PHYSBAM_OVERRIDE
    {TV result;for(int i=1;i<=d;i++) result+=weights[i]*V(parents[i]);
    return result;}

    TV Embedded_Acceleration(ARRAY_VIEW<const TV> F,ARRAY_VIEW<const TWIST<TV> > wrench) const PHYSBAM_OVERRIDE
    {TV a;for(int i=1;i<=d;i++) a+=weights[i]*particles.one_over_mass(parents[i])*F(parents[i]);
    return a;}

    template<class T_ARRAY_T2>
    typename T_ARRAY_T2::ELEMENT Embedded_Value(const T_ARRAY_T2& u) const
    {assert(u.Size()==particles.array_collection->Size());typedef typename T_ARRAY_T2::ELEMENT T2;
    T2 result=T2();for(int i=1;i<=d;i++) result+=weights[i]*u(parents[i]);
    return result;}

    T One_Over_Effective_Mass() const PHYSBAM_OVERRIDE
    {T result=0;for(int i=1;i<=d;i++) result+=sqr(weights[i])*particles.one_over_mass(parents[i]);
    return result;}

    void Apply_Impulse(const TV& impulse) PHYSBAM_OVERRIDE
    {for(int i=1;i<=d;i++) particles.V(parents[i])+=particles.one_over_mass(parents[i])*weights[i]*impulse;}

    void Apply_Push(const TV& impulse)
    {for(int i=1;i<=d;i++) particles.X(parents[i])+=particles.one_over_mass(parents[i])*weights[i]*impulse;}

    void Apply_Displacement_To_Parents_Based_On_Embedding(const TV& dX,const ARRAY<bool>* skip_particle) PHYSBAM_OVERRIDE
    {VECTOR<T,d> c=weights/weights.Magnitude_Squared();
    for(int i=1;i<=d;i++) if(!skip_particle || !(*skip_particle)(parents[i])) particles.X(parents[i])+=c[i]*dX;}

    void Apply_Velocity_Change_To_Parents_Based_On_Embedding(const TV& dV,const ARRAY<bool>* skip_particle) PHYSBAM_OVERRIDE
    {VECTOR<T,d> c=weights/weights.Magnitude_Squared();
    for(int i=1;i<=d;i++) if(!skip_particle || !(*skip_particle)(parents[i])) particles.V(parents[i])+=c[i]*dV;}

    void Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,const TV& force) const PHYSBAM_OVERRIDE
    {for(int i=1;i<=d;i++) F_full(parents[i])+=weights[i]*force;}

    void Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,ARRAY_VIEW<TWIST<TV> > wrench_full,const TV& force) const PHYSBAM_OVERRIDE
    {for(int i=1;i<=d;i++) F_full(parents[i])+=weights[i]*force;}

    void Distribute_Mass_To_Parents(ARRAY_VIEW<T> mass_full) const PHYSBAM_OVERRIDE
    {for(int i=1;i<=d;i++) mass_full(parents[i])+=weights[i]*particles.mass(particle_index);}

    ARRAY<int> Parents() const PHYSBAM_OVERRIDE
    {ARRAY<int> parents_array(d);for(int i=1;i<=d;i++) parents_array(i)=parents[i];return parents_array;}

    ARRAY<T> Weights() const PHYSBAM_OVERRIDE
    {ARRAY<T> weights_array(d);for(int i=1;i<=d;i++) weights_array(i)=weights[i];return weights_array;}

private:
    void Read_Helper(TYPED_ISTREAM& input) PHYSBAM_OVERRIDE
    {BINDING<TV>::Read_Helper(input);Read_Binary(input,parents,weights);}

    void Write_Helper(TYPED_OSTREAM& output) const PHYSBAM_OVERRIDE
    {BINDING<TV>::Write_Helper(output);Write_Binary(output,parents,weights);}

//#####################################################################
}; 
}
#endif
