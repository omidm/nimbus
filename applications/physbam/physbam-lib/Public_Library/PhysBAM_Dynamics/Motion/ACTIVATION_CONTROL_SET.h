//#####################################################################
// Copyright 2004-2005, Igor Neverov, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ACTIVATION_CONTROL_SET
//#####################################################################
#ifndef __ACTIVATION_CONTROL_SET__
#define __ACTIVATION_CONTROL_SET__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/SOLIDS_FORCES.h>
#include <PhysBAM_Dynamics/Motion/FACE_CONTROL_SET.h>
namespace PhysBAM{

template<class T>
class ACTIVATION_CONTROL_SET:public FACE_CONTROL_SET<T>
{
    typedef VECTOR<T,3> TV;
public:
    typedef typename FACE_CONTROL_SET<T>::TYPE T_TYPE;
    using FACE_CONTROL_SET<T>::ACTIVATION;

    ARRAY<T> activations,activations_save;
    ARRAY<std::string> activation_names;
    ARRAY<bool> activation_active;
    mutable int single_activation_used_for_force_derivative;
    SOLIDS_FORCES<TV>* muscle_force;
    T activation_cutoff,max_optimization_step_length,activation_penalty_coefficient;

    ACTIVATION_CONTROL_SET();
    virtual ~ACTIVATION_CONTROL_SET();

//#####################################################################
    void Set_Attachment_Positions(ARRAY<TV>& X) const PHYSBAM_OVERRIDE {}
    void Kinematically_Update_Positions(ARRAY<TV>& X) const PHYSBAM_OVERRIDE {}
    void Kinematically_Update_Jacobian(ARRAY<TV>& dX) const PHYSBAM_OVERRIDE {}
    int Size() const PHYSBAM_OVERRIDE;
    T_TYPE Type() const PHYSBAM_OVERRIDE;
    int Add_Activation(const std::string name="",const bool active=true);
    T operator()(const int control_id) const PHYSBAM_OVERRIDE;
    T& operator()(const int control_id) PHYSBAM_OVERRIDE;
    T Identity(const int control_id) const PHYSBAM_OVERRIDE;
    void Maximal_Controls() PHYSBAM_OVERRIDE;
    bool Control_Active(const int control_id) const PHYSBAM_OVERRIDE;
    bool& Control_Active(const int control_id) PHYSBAM_OVERRIDE;
    bool Positions_Determined_Kinematically(const int control_id) const PHYSBAM_OVERRIDE;
    void Force_Derivative(ARRAY<TV>& dFdl,ARRAY<TWIST<TV> >& dFrdl,const int control_id) const PHYSBAM_OVERRIDE;
    T Penalty() const PHYSBAM_OVERRIDE;
    T Penalty_Gradient(const int control_id) const PHYSBAM_OVERRIDE;
    T Penalty_Hessian(const int control_id1,const int control_id2) const PHYSBAM_OVERRIDE;
    void Save_Controls() PHYSBAM_OVERRIDE;
    void Project_Parameters_To_Allowable_Range(const bool active_controls_only) PHYSBAM_OVERRIDE;
    void Interpolate(const T interpolation_fraction) PHYSBAM_OVERRIDE;
    void Scale(const T scale) PHYSBAM_OVERRIDE;
    T Distance(const ARRAY<T>& weights,int base) PHYSBAM_OVERRIDE;
    void Print_Diagnostics(std::ostream& output) const PHYSBAM_OVERRIDE;
    static T Piecewise_Quadratic_Penalty(const T a,const T a_min,const T a_max);
    static T Piecewise_Quadratic_Penalty_Prime(const T a,const T a_min,const T a_max);
    static T Piecewise_Quadratic_Penalty_Double_Prime(const T a,const T a_min,const T a_max);
    void Read_Configuration(const STREAM_TYPE& stream_type,std::istream& input_stream);
    void Write_Configuration(const STREAM_TYPE& stream_type,std::ostream& output_stream) const;
//#####################################################################
};
}
#endif
