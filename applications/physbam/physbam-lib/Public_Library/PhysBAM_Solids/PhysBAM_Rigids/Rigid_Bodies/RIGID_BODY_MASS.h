//#####################################################################
// Copyright 2006, Geoffrey Irving, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_MASS
//#####################################################################
#ifndef __RIGID_BODY_MASS__
#define __RIGID_BODY_MASS__

#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_0X0.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_POLICY.h>
namespace PhysBAM{

template<class TV,bool world_space> // world_space=false
class RIGID_BODY_MASS
{
    typedef typename TV::SCALAR T;
    typedef typename IF<world_space,typename RIGID_BODY_POLICY<TV>::WORLD_SPACE_INERTIA_TENSOR,typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR>::TYPE T_INERTIA_TENSOR;
public:
    T mass;
    T_INERTIA_TENSOR inertia_tensor;

    RIGID_BODY_MASS()
        :mass(),inertia_tensor()
    {}

    RIGID_BODY_MASS(const T& mass_input,const T_INERTIA_TENSOR& inertia_tensor_input)
        :mass(mass_input),inertia_tensor(inertia_tensor_input)
    {
        assert(Valid());
    }

    bool Valid() const
    {return mass>=0 && inertia_tensor.Positive_Semidefinite();}

    bool operator==(const RIGID_BODY_MASS& other) const
    {return mass==other.mass && inertia_tensor==other.inertia_tensor;}

    bool operator!=(const RIGID_BODY_MASS& other) const
    {return !(*this==other);}

    static RIGID_BODY_MASS Identity_Mass()
    {RIGID_BODY_MASS mass;mass.mass=1;Identity_Inertia_Tensor(mass.inertia_tensor);return mass;}

private:
    template<class T_INERTIA_TENSOR_2> static void Identity_Inertia_Tensor(T_INERTIA_TENSOR_2& inertia_tensor_2)
    {inertia_tensor_2+=1;}

    static void Identity_Inertia_Tensor(MATRIX<T,0>& zero)
    {}
public:

    RIGID_BODY_MASS& operator*=(const T a)
    {assert(a>=0);mass*=a;inertia_tensor*=a;return *this;}

    RIGID_BODY_MASS operator*(const T a) const
    {assert(a>=0);return RIGID_BODY_MASS<TV>(mass*a,inertia_tensor*a);}
    
    TWIST<TV> operator*(const TWIST<TV>& v) const
    {STATIC_ASSERT(world_space);return TWIST<TV>(mass*v.linear,inertia_tensor*v.angular);}

    T Inner_Product(const TWIST<TV>& v1,const TWIST<TV>& v2) const
    {STATIC_ASSERT(world_space);return mass*Dot_Product(v1.linear,v2.linear)+Dot_Product(v1.angular,inertia_tensor*v2.angular);}

    MATRIX<T,0> World_Space_Inertia_Tensor(const ROTATION<VECTOR<T,1> > orientation) const
    {STATIC_ASSERT((AND<NOT<world_space>::value,TV::m==1>::value));return MATRIX<T,0>();}

    MATRIX<T,0> World_Space_Inertia_Tensor_Inverse(const ROTATION<VECTOR<T,1> > orientation) const
    {STATIC_ASSERT((AND<NOT<world_space>::value,TV::m==1>::value));return MATRIX<T,0>();}

    MATRIX<T,1> World_Space_Inertia_Tensor(const ROTATION<VECTOR<T,2> >& orientation) const
    {STATIC_ASSERT((AND<NOT<world_space>::value,TV::m==2>::value));return inertia_tensor;}

    MATRIX<T,1> World_Space_Inertia_Tensor_Inverse(const ROTATION<VECTOR<T,2> >& orientation) const
    {STATIC_ASSERT((AND<NOT<world_space>::value,TV::m==2>::value));return inertia_tensor.Inverse();}

    SYMMETRIC_MATRIX<T,3> World_Space_Inertia_Tensor(const ROTATION<VECTOR<T,3> >& orientation) const
    {STATIC_ASSERT((AND<NOT<world_space>::value,TV::m==3>::value));
    MATRIX<T,3> object_to_world_transformation=orientation.Rotation_Matrix();
    return SYMMETRIC_MATRIX<T,3>::Conjugate(object_to_world_transformation,inertia_tensor);}

    SYMMETRIC_MATRIX<T,3> World_Space_Inertia_Tensor_Inverse(const ROTATION<VECTOR<T,3> >& orientation) const
    {STATIC_ASSERT((AND<NOT<world_space>::value,TV::m==3>::value));
    MATRIX<T,3> object_to_world_transformation=orientation.Rotation_Matrix();
    return SYMMETRIC_MATRIX<T,3>::Conjugate(object_to_world_transformation,inertia_tensor.Inverse());}

    MATRIX<T,0> World_Space_Inertia_Tensor(const FRAME<TV>& frame,const VECTOR<T,1>& reference_point) const // relative to a reference point
    {STATIC_ASSERT((AND<NOT<world_space>::value,TV::m==1>::value));return MATRIX<T,0>();}

    MATRIX<T,1> World_Space_Inertia_Tensor(const FRAME<TV>& frame,const VECTOR<T,2>& reference_point) const // relative to a reference point
    {STATIC_ASSERT((AND<NOT<world_space>::value,TV::m==2>::value));TV offset=reference_point-frame.t;
    return inertia_tensor+mass*offset.Magnitude_Squared();}

    SYMMETRIC_MATRIX<T,3> World_Space_Inertia_Tensor(const FRAME<TV>& frame,const VECTOR<T,3>& reference_point) const // relative to a reference point
    {STATIC_ASSERT((AND<NOT<world_space>::value,TV::m==3>::value));TV offset=reference_point-frame.t;
    return World_Space_Inertia_Tensor(frame.r)+mass*(offset.Magnitude_Squared()-SYMMETRIC_MATRIX<T,3>::Outer_Product(offset));}

    VECTOR<T,0> World_Space_Inertia_Tensor_Times(const ROTATION<VECTOR<T,1> > orientation,const VECTOR<T,0> angular_velocity) const
    {STATIC_ASSERT((AND<NOT<world_space>::value,TV::m==1>::value));return VECTOR<T,0>();}

    VECTOR<T,0> World_Space_Inertia_Tensor_Inverse_Times(const ROTATION<VECTOR<T,1> > orientation,const VECTOR<T,0> angular_velocity) const
    {STATIC_ASSERT((AND<NOT<world_space>::value,TV::m==1>::value));return inertia_tensor.Solve_Linear_System(angular_velocity);}

    VECTOR<T,1> World_Space_Inertia_Tensor_Times(const ROTATION<VECTOR<T,2> >& orientation,const VECTOR<T,1> angular_velocity) const
    {STATIC_ASSERT((AND<NOT<world_space>::value,TV::m==2>::value));return inertia_tensor*angular_velocity;}

    VECTOR<T,1> World_Space_Inertia_Tensor_Inverse_Times(const ROTATION<VECTOR<T,2> >& orientation,const VECTOR<T,1> angular_momentum) const
    {STATIC_ASSERT((AND<NOT<world_space>::value,TV::m==2>::value));return inertia_tensor.Solve_Linear_System(angular_momentum);}

    VECTOR<T,3> World_Space_Inertia_Tensor_Times(const ROTATION<VECTOR<T,3> >& orientation,const VECTOR<T,3>& angular_velocity) const
    {STATIC_ASSERT((AND<NOT<world_space>::value,TV::m==3>::value));return orientation.Rotate(inertia_tensor*orientation.Inverse_Rotate(angular_velocity));}

    VECTOR<T,3> World_Space_Inertia_Tensor_Inverse_Times(const ROTATION<VECTOR<T,3> >& orientation,const VECTOR<T,3>& angular_momentum) const
    {STATIC_ASSERT((AND<NOT<world_space>::value,TV::m==3>::value));return orientation.Rotate(inertia_tensor.Solve_Linear_System(orientation.Inverse_Rotate(angular_momentum)));}

    T Rotational_Kinetic_Energy(const ROTATION<VECTOR<T,1> > orientation,const VECTOR<T,0> angular_momentum) const
    {STATIC_ASSERT((AND<NOT<world_space>::value,TV::m==1>::value));return 0;}

    T Rotational_Kinetic_Energy(const ROTATION<VECTOR<T,2> >& orientation,const VECTOR<T,1> angular_momentum) const
    {STATIC_ASSERT((AND<NOT<world_space>::value,TV::m==2>::value));return (T).5*inertia_tensor.Inverse_Inner_Product(angular_momentum,angular_momentum);}

    T Rotational_Kinetic_Energy(const ROTATION<VECTOR<T,3> >& orientation,const VECTOR<T,3>& angular_momentum) const
    {STATIC_ASSERT((AND<NOT<world_space>::value,TV::m==3>::value));TV object_space_angular_momentum=orientation.Inverse_Rotate(angular_momentum);
    return (T).5*inertia_tensor.Inverse_Inner_Product(object_space_angular_momentum,object_space_angular_momentum);}

//#####################################################################
};

//#####################################################################
template<class TV,bool b> struct PRODUCT<RIGID_BODY_MASS<TV,b>,TWIST<TV> > {typedef TWIST<TV> TYPE;};
//#####################################################################
}
#include <PhysBAM_Solids/PhysBAM_Rigids/Read_Write/Rigid_Bodies/READ_WRITE_RIGID_BODY_MASS.h>
#endif
