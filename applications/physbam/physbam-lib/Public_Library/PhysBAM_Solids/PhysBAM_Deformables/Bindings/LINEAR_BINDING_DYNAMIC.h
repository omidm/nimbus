//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_BINDING_DYNAMIC
//##################################################################### 
#ifndef __LINEAR_BINDING_DYNAMIC__
#define __LINEAR_BINDING_DYNAMIC__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING.h>
namespace PhysBAM{

template<class TV>
class LINEAR_BINDING_DYNAMIC:public BINDING<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef BINDING<TV> BASE;
    using BASE::particles;using BASE::particle_index;
    using BASE::One_Over_Effective_Mass; // silence -Woverloaded-virtual

    ARRAY<int> parents;
    ARRAY<T> weights; // weights should sum to 1

    LINEAR_BINDING_DYNAMIC(PARTICLES<TV>& particles_input)
        :BINDING<TV>(particles_input)
    {}

    LINEAR_BINDING_DYNAMIC(PARTICLES<TV>& particles_input,const int particle_index_input,const int number_of_parents)
        :BINDING<TV>(particles_input,particle_index_input),parents(number_of_parents),weights(number_of_parents)
    {}

    static LINEAR_BINDING_DYNAMIC* Create(GEOMETRY_PARTICLES<TV>& particles)
    {return new LINEAR_BINDING_DYNAMIC(dynamic_cast<PARTICLES<TV>&>(particles));}

    virtual int Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    static int Static_Name()
    {return 3;}

    TV Embedded_Position() const PHYSBAM_OVERRIDE
    {TV result;for(int i=1;i<=parents.m;i++) result+=weights(i)*particles.X(parents(i));
    return result;}

    TV Embedded_Position(ARRAY_VIEW<const TV> X) const PHYSBAM_OVERRIDE
    {TV result;for(int i=1;i<=parents.m;i++) result+=weights(i)*X(parents(i));
    return result;}

    TV Embedded_Velocity() const  PHYSBAM_OVERRIDE
    {TV result;for(int i=1;i<=parents.m;i++) result+=weights(i)*particles.V(parents(i));
    return result;}

    TV Embedded_Velocity(ARRAY_VIEW<const TV> V) const  PHYSBAM_OVERRIDE
    {TV result;for(int i=1;i<=parents.m;i++) result+=weights(i)*V(parents(i));
    return result;}

    TV Embedded_Velocity(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > twist) const  PHYSBAM_OVERRIDE
    {TV result;for(int i=1;i<=parents.m;i++) result+=weights(i)*V(parents(i));
    return result;}

    TV Embedded_Acceleration(ARRAY_VIEW<const TV> F,ARRAY_VIEW<const TWIST<TV> > wrench) const PHYSBAM_OVERRIDE
    {TV a;for(int i=1;i<=parents.m;i++) a+=weights(i)*particles.one_over_mass(parents(i))*F(parents(i));
    return a;}

    T One_Over_Effective_Mass() const PHYSBAM_OVERRIDE
    {T result=0;for(int i=1;i<=parents.m;i++) result+=sqr(weights(i))*particles.one_over_mass(parents(i));
    return result;}

    void Apply_Impulse(const TV& impulse) PHYSBAM_OVERRIDE
    {for(int i=1;i<=parents.m;i++) particles.V(parents(i))+=particles.one_over_mass(parents(i))*weights(i)*impulse;}

    void Apply_Displacement_To_Parents_Based_On_Embedding(const TV& dX,const ARRAY<bool>* skip_particle) PHYSBAM_OVERRIDE
    {T one_over_weights_squared=1/ARRAYS_COMPUTATIONS::Magnitude_Squared(weights);
    for(int i=1;i<=parents.m;i++) if(!skip_particle || !(*skip_particle)(parents(i))) particles.X(parents(i))+=one_over_weights_squared*weights(i)*dX;}

    void Apply_Velocity_Change_To_Parents_Based_On_Embedding(const TV& dV,const ARRAY<bool>* skip_particle) PHYSBAM_OVERRIDE
    {T one_over_weights_squared=1/ARRAYS_COMPUTATIONS::Magnitude_Squared(weights);
    for(int i=1;i<=parents.m;i++) if(!skip_particle || !(*skip_particle)(parents(i))) particles.V(parents(i))+=one_over_weights_squared*weights(i)*dV;}

    void Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,const TV& force) const PHYSBAM_OVERRIDE
    {for(int i=1;i<=parents.m;i++) F_full(parents(i))+=weights(i)*force;}

    void Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,ARRAY_VIEW<TWIST<TV> > wrench_full,const TV& force) const PHYSBAM_OVERRIDE
    {for(int i=1;i<=parents.m;i++) F_full(parents(i))+=weights(i)*force;}

    void Distribute_Mass_To_Parents(ARRAY_VIEW<T> mass_full) const PHYSBAM_OVERRIDE
    {for(int i=1;i<=parents.m;i++) mass_full(parents(i))+=weights(i)*particles.mass(particle_index);}

    ARRAY<int> Parents() const PHYSBAM_OVERRIDE
    {return parents;}

    ARRAY<T> Weights() const PHYSBAM_OVERRIDE
    {return weights;}

private:
    void Read_Helper(TYPED_ISTREAM& input) PHYSBAM_OVERRIDE
    {BINDING<TV>::Read_Helper(input);Read_Binary(input,parents,weights);}

    void Write_Helper(TYPED_OSTREAM& output) const PHYSBAM_OVERRIDE
    {BINDING<TV>::Write_Helper(output);Write_Binary(output,parents,weights);}

//#####################################################################
}; 
}
#endif
