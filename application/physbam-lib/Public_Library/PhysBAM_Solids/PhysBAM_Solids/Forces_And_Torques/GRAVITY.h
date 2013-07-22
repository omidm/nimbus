//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRAVITY
//#####################################################################
#ifndef __GRAVITY__
#define __GRAVITY__

#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/POINTWISE_FORCE.h>
namespace PhysBAM{

template<class TV>
class GRAVITY:public POINTWISE_FORCE<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef POINTWISE_FORCE<TV> BASE;
    using BASE::particles;using BASE::rigid_body_collection;using BASE::influenced_particles;
    using BASE::mpi_solids;using BASE::force_particles;using BASE::force_rigid_body_particles;
    typedef typename FORCE_ELEMENTS::ITERATOR ELEMENT_ITERATOR;

    T gravity;
    TV downward_direction;
public:

    GRAVITY(PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,ARRAY<int>* influenced_particles_input,
        ARRAY<int>* influenced_rigid_body_particles_input,const T gravity_input=9.8,const TV& downward_direction_input=-TV::Axis_Vector(2-(TV::dimension==1)))
        :POINTWISE_FORCE<TV>(particles_input,rigid_body_collection_input,influenced_particles_input,influenced_rigid_body_particles_input),gravity(gravity_input),
        downward_direction(downward_direction_input)
    {}

    GRAVITY(PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const bool influence_all_particles_input,
        const bool influence_all_rigid_body_particles_input,const T gravity_input=9.8,const TV& downward_direction_input=-TV::Axis_Vector(2-(TV::dimension==1)))
        :POINTWISE_FORCE<TV>(particles_input,rigid_body_collection_input,influence_all_particles_input,influence_all_rigid_body_particles_input),
        gravity(gravity_input),downward_direction(downward_direction_input)
    {}

    template<class T_MESH>
    GRAVITY(PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const T_MESH& mesh,ARRAY<int>* influenced_rigid_body_particles_input,
        const T gravity_input=9.8,const TV& downward_direction_input=((TV::dimension==1)?((T)-1*TV::All_Ones_Vector()):TV(VECTOR<T,2>(0,-1))))
        :POINTWISE_FORCE<TV>(particles_input,rigid_body_collection_input,mesh,influenced_rigid_body_particles_input),gravity(gravity_input),downward_direction(downward_direction_input)
    {
        mesh.elements.Flattened().Get_Unique(*influenced_particles);
    }

    virtual ~GRAVITY()
    {}

    void Set_Gravity(const TV& gravity_vector)
    {downward_direction=gravity_vector;gravity=downward_direction.Normalize();}

    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE
    {}

    void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const PHYSBAM_OVERRIDE
    {}

    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE
    {}

    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE
    {}

    int Velocity_Dependent_Forces_Size() const PHYSBAM_OVERRIDE
    {return 0;}

    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<T> aggregate,const T time) const PHYSBAM_OVERRIDE
    {}

    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE
    {}

//#####################################################################
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE;
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
