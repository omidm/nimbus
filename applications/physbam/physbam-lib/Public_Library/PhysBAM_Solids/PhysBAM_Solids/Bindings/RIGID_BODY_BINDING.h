//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Sergey Levine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_BINDING
//##################################################################### 
#ifndef __RIGID_BODY_BINDING__
#define __RIGID_BODY_BINDING__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
namespace PhysBAM{

template<class TV>
class RIGID_BODY_BINDING:public BINDING<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
public:
    typedef BINDING<TV> BASE;
    using BASE::particles;using BASE::particle_index;

    RIGID_BODY_COLLECTION<TV>* rigid_body_collection;
    int rigid_body_particle_index;
    TV object_space_position;

    RIGID_BODY_BINDING(PARTICLES<TV>& particles_input)
        :BINDING<TV>(particles_input),rigid_body_collection(0),rigid_body_particle_index(0)
    {}

    RIGID_BODY_BINDING(PARTICLES<TV>& particles_input,const int particle_index_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,
        const int rigid_body_particle_index_input,const TV& object_space_position_input)
        :BINDING<TV>(particles_input,particle_index_input),rigid_body_collection(&rigid_body_collection_input),rigid_body_particle_index(rigid_body_particle_index_input),
        object_space_position(object_space_position_input)
    {}

    static RIGID_BODY_BINDING* Create(GEOMETRY_PARTICLES<TV>& particles)
    {return new RIGID_BODY_BINDING(dynamic_cast<PARTICLES<TV>&>(particles));}

    virtual int Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    static int Static_Name()
    {return 4;}

    // TODO: This may not handle kinematic/static bodies
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE
    {dependency_mesh.Add_Element_If_Not_Already_There(VECTOR<int,2>(particles.array_collection->Size()+rigid_body_particle_index,particle_index));}

    RIGID_BODY<TV>& Rigid_Body() const
    {return rigid_body_collection->Rigid_Body(rigid_body_particle_index);}

    int Id_Number() const
    {return rigid_body_particle_index;}

    TV Embedded_Position() const PHYSBAM_OVERRIDE
    {return Rigid_Body().World_Space_Point(object_space_position);}

    TV Embedded_Position(ARRAY_VIEW<const TV> X) const PHYSBAM_OVERRIDE
    {PHYSBAM_NOT_IMPLEMENTED();}

    TV Embedded_Velocity() const PHYSBAM_OVERRIDE
    {return Rigid_Body().Pointwise_Object_Velocity(Embedded_Position());}

    TV Embedded_Velocity(ARRAY_VIEW<const TV> V) const PHYSBAM_OVERRIDE
    {PHYSBAM_FATAL_ERROR("Rigid body velocities required");}

    TV Embedded_Velocity(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > twist) const PHYSBAM_OVERRIDE
    {return RIGID_BODY_STATE<TV>(rigid_body_collection->Rigid_Body(rigid_body_particle_index).Frame(),twist(rigid_body_particle_index)).Pointwise_Object_Velocity(Embedded_Position());}

    TV Embedded_Acceleration(ARRAY_VIEW<const TV> F,ARRAY_VIEW<const TWIST<TV> > wrench_input) const PHYSBAM_OVERRIDE// TODO: consider including centrifugal acceleration
    {RIGID_BODY<TV>& rigid_body=Rigid_Body();
    if(rigid_body.Has_Infinite_Inertia()) return TV();
    const TWIST<TV>& wrench=wrench_input(rigid_body_particle_index);
    return wrench.linear/rigid_body.Mass()+TV::Cross_Product(rigid_body.World_Space_Inertia_Tensor_Inverse()*wrench.angular,rigid_body.World_Space_Vector(object_space_position));}

    T One_Over_Effective_Mass(const TV& direction) const  PHYSBAM_OVERRIDE// assumes direction is normalized
    {RIGID_BODY<TV>& rigid_body=Rigid_Body();
    if(rigid_body.Has_Infinite_Inertia()) return 0;
    TV object_space_direction=rigid_body.Object_Space_Vector(direction);
    return TV::Dot_Product(object_space_direction,rigid_body.Object_Space_Impulse_Factor(object_space_position)*object_space_direction);}

    T One_Over_Effective_Mass() const  PHYSBAM_OVERRIDE// return a lower bound for effective mass over all directions
    {RIGID_BODY<TV>& rigid_body=Rigid_Body();if(rigid_body.Has_Infinite_Inertia()) return 0;return rigid_body.Object_Space_Impulse_Factor(object_space_position).Fast_Eigenvalues().Max();}

    void Apply_Impulse(const TV& impulse) PHYSBAM_OVERRIDE
    {RIGID_BODY<TV>& rigid_body=Rigid_Body();
    if(rigid_body.Has_Infinite_Inertia()) return;
    T& mass=rigid_body_collection->rigid_body_particle.mass(rigid_body_particle_index);
    rigid_body_collection->rigid_body_particle.V(rigid_body_particle_index)+=impulse/mass;
    rigid_body.Update_Angular_Momentum();
    rigid_body.Angular_Momentum()+=TV::Cross_Product(rigid_body.World_Space_Vector(object_space_position),impulse);
    rigid_body.Update_Angular_Velocity();}

    void Apply_Displacement_To_Parents_Based_On_Embedding(const TV& dX,const ARRAY<bool>* skip_particle) PHYSBAM_OVERRIDE
    {PHYSBAM_NOT_IMPLEMENTED();}

    void Apply_Velocity_Change_To_Parents_Based_On_Embedding(const TV& dV,const ARRAY<bool>* skip_particle) PHYSBAM_OVERRIDE
    {PHYSBAM_NOT_IMPLEMENTED();}

    void Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,const TV& force) const PHYSBAM_OVERRIDE
    {PHYSBAM_FATAL_ERROR("Wrench required");}

    void Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,ARRAY_VIEW<TWIST<TV> > wrench_full,const TV& force) const PHYSBAM_OVERRIDE
    {if(Rigid_Body().Has_Infinite_Inertia()) return;
    wrench_full(rigid_body_particle_index).linear+=force;
    wrench_full(rigid_body_particle_index).angular+=TV::Cross_Product(Rigid_Body().World_Space_Vector(object_space_position),force);}

    void Distribute_Mass_To_Parents(ARRAY_VIEW<T> mass_full) const PHYSBAM_OVERRIDE
    {/*if(particles.mass(particle_index)) PHYSBAM_NOT_IMPLEMENTED();*/}

    ARRAY<int> Parents() const PHYSBAM_OVERRIDE
    {PHYSBAM_NOT_IMPLEMENTED();return ARRAY<int>();} // TODO: consider failing instead

    ARRAY<T> Weights() const PHYSBAM_OVERRIDE
    {PHYSBAM_NOT_IMPLEMENTED();return ARRAY<T>();} // TODO: consider failing instead

private:
    void Read_Helper(TYPED_ISTREAM& input) PHYSBAM_OVERRIDE
    {BINDING<TV>::Read_Helper(input);Read_Binary(input,rigid_body_particle_index);Read_Binary(input,object_space_position);}

    void Write_Helper(TYPED_OSTREAM& output) const PHYSBAM_OVERRIDE
    {BINDING<TV>::Write_Helper(output);Write_Binary(output,rigid_body_particle_index);Write_Binary(output,object_space_position);}

//#####################################################################
}; 
}
#endif
