//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONTACT
//##################################################################### 
#ifndef __CONTACT__
#define __CONTACT__

#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>

namespace PhysBAM{

//#####################################################################
// Class CONTACT
//#####################################################################
template<class TV> class CONTACT
{
public:
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    enum WORKAROUND {d=TV::dimension,s=T_SPIN::dimension};

    VECTOR<int,2> id;

    TV location;
    TV normal;
    VECTOR<TV,d-1> tangent;
    T distance;
    T coefficient_of_friction;

    VECTOR<TWIST<TV>,2> normal_constraint;
    VECTOR<TWIST<TV>,2> inverse_mass_times_normal_constraint;
    T normal_relative_velocity;

    VECTOR<VECTOR<TWIST<TV>,d-1>,2> tangent_constraint;
    VECTOR<VECTOR<TWIST<TV>,d-1>,2> inverse_mass_times_tangent_constraint;
    
    T normal_diagonal;
    VECTOR<T,d-1> tangent_diagonal;
    
    CONTACT():
        id(0,0)
    {}

    CONTACT(RIGID_BODY<TV>& body_1,RIGID_BODY<TV>& body_2,TV& _location,TV& _normal,T _distance,T dt):
        location(_location),normal(_normal),distance(_distance),coefficient_of_friction(RIGID_BODY<TV>::Coefficient_Of_Friction(body_1,body_2))
    {
        id(1)=body_1.particle_index;
        id(2)=body_2.particle_index;

        TV r_1=location-body_1.X();
        TV r_2=location-body_2.X();

        normal_constraint(1).linear=-normal;
        normal_constraint(1).angular=TV::Cross_Product(normal,r_1);

        normal_constraint(2).linear=normal;
        normal_constraint(2).angular=-TV::Cross_Product(normal,r_2);

        inverse_mass_times_normal_constraint(1)=body_1.Inertia_Inverse_Times(normal_constraint(1));
        inverse_mass_times_normal_constraint(2)=body_2.Inertia_Inverse_Times(normal_constraint(2));

        normal_diagonal=0;
        if(!body_1.Has_Infinite_Inertia())
            normal_diagonal+=TWIST<TV>::Dot_Product(normal_constraint(1),inverse_mass_times_normal_constraint(1));
        if(!body_2.Has_Infinite_Inertia())
            normal_diagonal+=TWIST<TV>::Dot_Product(normal_constraint(2),inverse_mass_times_normal_constraint(2));
        normal_relative_velocity=-distance/dt;

        Compute_Tangent_Helper(*this);
        
        for(int i=1;i<d;i++)
        {
            tangent_constraint(1)(i).linear=-tangent(i);
            tangent_constraint(1)(i).angular=TV::Cross_Product(tangent(i),r_1);

            tangent_constraint(2)(i).linear=tangent(i);
            tangent_constraint(2)(i).angular=-TV::Cross_Product(tangent(i),r_2);

            inverse_mass_times_tangent_constraint(1)(i)=body_1.Inertia_Inverse_Times(tangent_constraint(1)(i));
            inverse_mass_times_tangent_constraint(2)(i)=body_2.Inertia_Inverse_Times(tangent_constraint(2)(i));

            tangent_diagonal(i)=0;
            if(!body_1.Has_Infinite_Inertia())
                tangent_diagonal(i)+=TWIST<TV>::Dot_Product(tangent_constraint(1)(i),inverse_mass_times_tangent_constraint(1)(i));
            if(!body_2.Has_Infinite_Inertia())
                tangent_diagonal(i)+=TWIST<TV>::Dot_Product(tangent_constraint(2)(i),inverse_mass_times_tangent_constraint(2)(i));
        }
    }

    template<class T> void Compute_Tangent_Helper(CONTACT<VECTOR<T,1> >& contact) {}
    template<class T> void Compute_Tangent_Helper(CONTACT<VECTOR<T,2> >& contact)
    {
        contact.tangent(1)=contact.normal.Unit_Orthogonal_Vector();
    }
    template<class T> void Compute_Tangent_Helper(CONTACT<VECTOR<T,3> >& contact)
    {
        contact.tangent(1)=contact.normal.Unit_Orthogonal_Vector();
        contact.tangent(2)=VECTOR<T,3>::Cross_Product(contact.normal,contact.tangent(1)).Normalized();
    }
};

}
#endif
