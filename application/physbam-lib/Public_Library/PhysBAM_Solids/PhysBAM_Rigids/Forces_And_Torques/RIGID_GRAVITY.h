//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_GRAVITY
//#####################################################################
#ifndef __RIGID_GRAVITY__
#define __RIGID_GRAVITY__

#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGID_POINTWISE_FORCE.h>
namespace PhysBAM{

template<class TV>
class RIGID_GRAVITY:public RIGID_POINTWISE_FORCE<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef RIGID_POINTWISE_FORCE<TV> BASE;
    typedef typename FORCE_ELEMENTS::ITERATOR ELEMENT_ITERATOR;
    T gravity;
    TV downward_direction;
    
    using BASE::force_rigid_body_particles;using BASE::rigid_body_collection;
public:

    RIGID_GRAVITY(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,ARRAY<int>* influenced_rigid_body_particles_input,
        const T gravity_input=9.8,const TV& downward_direction_input=-TV::Axis_Vector(2-(TV::dimension==1)))
        :BASE(rigid_body_collection_input,influenced_rigid_body_particles_input),gravity(gravity_input),downward_direction(downward_direction_input)
    {}

    RIGID_GRAVITY(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const bool influence_all_rigid_body_particles_input,
        const T gravity_input=9.8,const TV& downward_direction_input=-TV::Axis_Vector(2-(TV::dimension==1)))
        :BASE(rigid_body_collection_input,influence_all_rigid_body_particles_input),gravity(gravity_input),downward_direction(downward_direction_input)
    {}

    virtual ~RIGID_GRAVITY()
    {}

    void Set_Gravity(const TV& gravity_vector)
    {downward_direction=gravity_vector;gravity=downward_direction.Normalize();}

    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE
    {}

    int Velocity_Dependent_Forces_Size() const PHYSBAM_OVERRIDE
    {return 0;}

    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<T> aggregate,const T time) const PHYSBAM_OVERRIDE
    {}

    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE
    {}

//#####################################################################
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE;
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
