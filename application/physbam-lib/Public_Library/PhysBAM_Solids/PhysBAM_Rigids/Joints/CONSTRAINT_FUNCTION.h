//#####################################################################
// Copyright 2007, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __CONSTRAINT_FUNCTION__
#define __CONSTRAINT_FUNCTION__

#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_ID.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_POLICY.h>
namespace PhysBAM{

template<class TV> class LINEAR_CONSTRAINT_FUNCTION;
template<class TV> class ANGULAR_CONSTRAINT_FUNCTION;
template<class TV> class LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION;
template<class TV> class ARTICULATED_RIGID_BODY;
template<class TV> class JOINT;
template<class TV> class RIGID_BODY;
//#####################################################################
// Class LINEAR_CONSTRAINT_FUNCTION
//#####################################################################
template<class T>
class LINEAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;typedef typename TV::SPIN T_SPIN;

    const JOINT<TV>& joint;
public:
    typedef TV T_IMPULSE;
    typedef TV T_CONSTRAINT_ERROR;

    const RIGID_BODY<TV>* rigid_body[2];
    const T dt,epsilon_scale;

    T_SPIN dt_angular_velocity[2];
    TV &location,p[2],c,r[2];
    T one_over_m;
    MATRIX<T,1,2> inverse_inertia_rhat_star[2];
    MATRIX<T,2> metric_tensor;

    LINEAR_CONSTRAINT_FUNCTION(const ARTICULATED_RIGID_BODY<TV>& arb,const JOINT_ID joint_id,const T dt_input,const T epsilon_scale_input,TV& location_input);

    T Inner_Product(const T_IMPULSE& j1,const T_IMPULSE& j2) const
    {return T_IMPULSE::Dot_Product(j1,metric_tensor*j2);}

    T Magnitude_Squared(const T_IMPULSE& j) const
    {return metric_tensor.Symmetric_Conjugate(j);}

    T Magnitude(const T_IMPULSE& j) const
    {return sqrt(Magnitude_Squared(j));}

    T Convergence_Norm_Squared(const T_CONSTRAINT_ERROR& f) const
    {return f.Magnitude_Squared();} // units of length^2

//#####################################################################
    T_CONSTRAINT_ERROR F(const T_IMPULSE& j) const;
    MATRIX<T,2> Jacobian(const T_IMPULSE& j) const;
    T_CONSTRAINT_ERROR F_Helper(const T_IMPULSE& j,const int i) const;
    MATRIX<T,2> Jacobian_Helper(const T_IMPULSE& j,const int i) const;
private:
    void Initialize();
//#####################################################################
};
//#####################################################################
// Class ANGULAR_CONSTRAINT_FUNCTION
//#####################################################################
template<class T>
class ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;typedef typename TV::SPIN T_SPIN;typedef typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR T_INERTIA_TENSOR;

    const JOINT<TV>& joint;
    const RIGID_BODY<TV>* rigid_body[2];

public:
    typedef T_SPIN T_IMPULSE;
    typedef VECTOR<T,1> T_CONSTRAINT_ERROR;

    const T dt,epsilon_scale;
    T_INERTIA_TENSOR inverse_inertia[2],metric_tensor;
    T_SPIN dt_angular_velocity[2],q_w_old[2]; // ={q_wp_old_q_pj_q_t,q_wc_old_q_cj}
    T length_scale_squared;

    ANGULAR_CONSTRAINT_FUNCTION(const ARTICULATED_RIGID_BODY<TV>& arb,const JOINT_ID joint_id,const T dt_input,const T epsilon_scale_input);

    MATRIX<T,1> Jacobian(const T_IMPULSE& j_tau) const
    {return Jacobian_Helper(0)+Jacobian_Helper(1);}

    MATRIX<T,1> Jacobian_Helper(const int i) const
    {return MATRIX<T,1>(inverse_inertia[i]);}

    T Inner_Product(const T_IMPULSE& j1,const T_IMPULSE& j2) const
    {return T_IMPULSE::Dot_Product(j1,metric_tensor*j2);}

    T Magnitude_Squared(const T_IMPULSE& j) const
    {return metric_tensor.Symmetric_Conjugate(j);}

    T Magnitude(const T_IMPULSE& j) const
    {return sqrt(Magnitude_Squared(j));}

    T Convergence_Norm_Squared(const T_CONSTRAINT_ERROR& f) const
    {return length_scale_squared*f.Magnitude_Squared();} // units of length^2

//#####################################################################
    T_CONSTRAINT_ERROR F(const T_IMPULSE& j_tau) const;
    T_SPIN F_Helper(const T_IMPULSE& j_tau,const int i) const;
private:
    void Initialize();
//#####################################################################
};
//#####################################################################
// Class LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION
//#####################################################################
template<class T>
class LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,2> >
{
    typedef VECTOR<T,2> TV;typedef typename TV::SPIN T_SPIN;

    LINEAR_CONSTRAINT_FUNCTION<TV> linear;
    ANGULAR_CONSTRAINT_FUNCTION<TV> angular;
    TV rhat[2];
    DIAGONAL_MATRIX<T,2> one_over_m_matrix;
    MATRIX<T,1,2> rhat_star[2];
    MATRIX<T,3> metric_tensor;

public:
    typedef VECTOR<T,3> T_IMPULSE;
    typedef VECTOR<T,3> T_CONSTRAINT_ERROR;

    const T dt,epsilon_scale;

    LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION(const ARTICULATED_RIGID_BODY<TV>& arb,const JOINT_ID joint_id,const T dt_input,const T epsilon_scale_input,TV& location);

    T Inner_Product(const T_IMPULSE& j1,const T_IMPULSE& j2) const
    {return T_IMPULSE::Dot_Product(j1,metric_tensor*j2);}

    T Magnitude_Squared(const T_IMPULSE& j) const
    {return metric_tensor.Symmetric_Conjugate(j);}

    T Magnitude(const T_IMPULSE& j) const
    {return sqrt(Magnitude_Squared(j));}

    T Convergence_Norm_Squared(const T_CONSTRAINT_ERROR& f) const
    {TV f_linear;T_SPIN f_angular;f.Split(f_linear,f_angular);
    return linear.Convergence_Norm_Squared(f_linear)+angular.Convergence_Norm_Squared(f_angular);} // units of length^2

//#####################################################################
    T_CONSTRAINT_ERROR F(const T_IMPULSE& j) const;
    MATRIX<T,3> Jacobian(const T_IMPULSE& j) const;
private:
    void Initialize();
    TV F_Linear_Helper(const T_SPIN& j,const int i) const;
    MATRIX<T,2,1> Jacobian_Linear_Helper(const T_SPIN& j,const int i) const;
//#####################################################################
};
//#####################################################################
// Class LINEAR_CONSTRAINT_FUNCTION
//#####################################################################
template<class T>
class LINEAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >
{
    typedef VECTOR<T,3> TV;typedef typename TV::SPIN T_SPIN;

    const JOINT<TV>& joint;
public:
    typedef TV T_IMPULSE;
    typedef TV T_CONSTRAINT_ERROR;

    const RIGID_BODY<TV>* rigid_body[2];
    const T dt,epsilon_scale;
    TV p[2]; // time n position
    ROTATION<TV> q[2]; // time n orientation
    T_SPIN dt_angular_velocity[2],c,r[2];
    TV& location;
    T one_over_m;
    MATRIX<T,3> inverse_inertia_rhat_star[2];
    MATRIX<T,3> metric_tensor;

    LINEAR_CONSTRAINT_FUNCTION(const ARTICULATED_RIGID_BODY<TV>& arb,const JOINT_ID joint_id,const T dt_input,const T epsilon_scale_input,TV& location_input);

    T Inner_Product(const T_IMPULSE& j1,const T_IMPULSE& j2) const
    {return T_IMPULSE::Dot_Product(j1,metric_tensor*j2);}

    T Magnitude_Squared(const T_IMPULSE& j) const
    {return metric_tensor.Symmetric_Conjugate(j);}

    T Magnitude(const T_IMPULSE& j) const
    {return sqrt(Magnitude_Squared(j));}

    T Convergence_Norm_Squared(const T_CONSTRAINT_ERROR& f) const
    {return f.Magnitude_Squared();} // units of length^2

//#####################################################################
    T_CONSTRAINT_ERROR F(const T_IMPULSE& j) const;
    MATRIX<T,3> Jacobian(const T_IMPULSE& j) const;
    T_CONSTRAINT_ERROR F_Helper(const T_IMPULSE& j,const int i) const;
    MATRIX<T,3> Jacobian_Helper(const T_IMPULSE& j,const int i) const;
private:
    void Initialize();
//#####################################################################
};
//#####################################################################
// Class ANGULAR_CONSTRAINT_FUNCTION
//#####################################################################
template<class T>
class ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >
{
    typedef VECTOR<T,3> TV;typedef typename TV::SPIN T_SPIN;

    const JOINT<TV>& joint;
    const RIGID_BODY<TV>* rigid_body[2];

public:
    typedef TV T_IMPULSE;
    typedef TV T_CONSTRAINT_ERROR;

    const T dt,epsilon_scale;
    SYMMETRIC_MATRIX<T,3> inverse_inertia[2];
    T_SPIN dt_angular_velocity[2];
    ROTATION<TV> q_w_old[2]; // ={q_wp_old_q_pj_q_t,q_wc_old_q_cj}
    MATRIX<T,3> modified_b_star[2];
    MATRIX<T,3> metric_tensor;
    T length_scale_squared;

    ANGULAR_CONSTRAINT_FUNCTION(const ARTICULATED_RIGID_BODY<TV>& arb,const JOINT_ID joint_id,const T dt_input,const T epsilon_scale_input);

    MATRIX<T,3> Jacobian(const T_IMPULSE& j_tau) const
    {return Jacobian_Helper(j_tau,j_tau,0)-Jacobian_Helper(j_tau,j_tau,1);}

    T Inner_Product(const T_IMPULSE& j1,const T_IMPULSE& j2) const
    {return T_IMPULSE::Dot_Product(j1,metric_tensor*j2);}

    T Magnitude_Squared(const T_IMPULSE& j) const
    {return metric_tensor.Symmetric_Conjugate(j);}

    T Magnitude(const T_IMPULSE& j) const
    {return sqrt(Magnitude_Squared(j));}

    T Convergence_Norm_Squared(const T_CONSTRAINT_ERROR& f) const
    {return length_scale_squared*f.Magnitude_Squared();} // units of length^2

//#####################################################################
    T_CONSTRAINT_ERROR F(const T_IMPULSE& j_tau) const;
    ROTATION<TV> F_Helper(const T_IMPULSE& j_tau,const int i) const;
    MATRIX<T,3> Jacobian_Helper(const T_IMPULSE& j_tau_i,const T_IMPULSE& j_tau_1_m_i,const int i) const;
    MATRIX<T,4,3> Jacobian_Old_Helper(const T_IMPULSE& j_tau,const int i) const;
private:
    void Initialize();
//#####################################################################
};
//#####################################################################
// Class LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION
//#####################################################################
template<class T>
class LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION<VECTOR<T,3> >
{
    typedef VECTOR<T,3> TV;typedef typename TV::SPIN T_SPIN;

    LINEAR_CONSTRAINT_FUNCTION<TV> linear;
    ANGULAR_CONSTRAINT_FUNCTION<TV> angular;
    TV rhat[2];
    DIAGONAL_MATRIX<T,3> one_over_m_matrix;
    MATRIX<T,6> metric_tensor;

public:
    typedef VECTOR<T,6> T_IMPULSE;
    typedef VECTOR<T,6> T_CONSTRAINT_ERROR;

    const T dt,epsilon_scale;

    LINEAR_AND_ANGULAR_CONSTRAINT_FUNCTION(const ARTICULATED_RIGID_BODY<TV>& arb,const JOINT_ID joint_id,const T dt_input,const T epsilon_scale_input,TV& location);

    T Inner_Product(const T_IMPULSE& j1,const T_IMPULSE& j2) const
    {return T_IMPULSE::Dot_Product(j1,metric_tensor*j2);}

    T Magnitude_Squared(const T_IMPULSE& j) const
    {return metric_tensor.Symmetric_Conjugate(j);}

    T Magnitude(const T_IMPULSE& j) const
    {return sqrt(Magnitude_Squared(j));}

    T Convergence_Norm_Squared(const T_CONSTRAINT_ERROR& f) const
    {TV f_linear;T_SPIN f_angular;f.Split(f_linear,f_angular);
    return linear.Convergence_Norm_Squared(f_linear)+angular.Convergence_Norm_Squared(f_angular);} // units of length^2

//#####################################################################
    T_CONSTRAINT_ERROR F(const T_IMPULSE& j) const;
    MATRIX<T,6> Jacobian(const T_IMPULSE& j) const;
private:
    void Initialize();
    TV F_Linear_Helper(const TV& j,const int i) const;
    MATRIX<T,3> Jacobian_Linear_Helper(const TV& j,const int i) const;
//#####################################################################
};
}
#endif
