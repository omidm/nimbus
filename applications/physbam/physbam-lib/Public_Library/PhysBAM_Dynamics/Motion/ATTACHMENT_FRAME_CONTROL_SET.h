//#####################################################################
// Copyright 2004-2007, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ATTACHMENT_FRAME_CONTROL_SET
//#####################################################################
#ifndef __ATTACHMENT_FRAME_CONTROL_SET__
#define __ATTACHMENT_FRAME_CONTROL_SET__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/SOLIDS_FORCES.h>
#include <PhysBAM_Dynamics/Motion/FACE_CONTROL_SET.h>
#include <PhysBAM_Dynamics/Motion/QUASI_RIGID_TRANSFORM_3D.h>
namespace PhysBAM{

template<class T>
class ATTACHMENT_FRAME_CONTROL_SET:public FACE_CONTROL_SET<T>
{
    typedef VECTOR<T,3> TV;
public:
    typedef typename FACE_CONTROL_SET<T>::TYPE T_TYPE;
    using FACE_CONTROL_SET<T>::ATTACHMENT_FRAME;

    QUASI_RIGID_TRANSFORM_3D<T> cranium_transform,jaw_transform,cranium_transform_save,jaw_transform_save;
    ARRAY<bool> coefficient_active;
    ARRAY<TV>& X;
    ARRAY<TV> X_save;
    SOLIDS_FORCES<TV>* muscle_force;
    ARRAY<ARRAY<int> >& attached_nodes;
    int jaw_attachment_index;
    T rigidity_penalty_coefficient,jaw_constraint_penalty_coefficient;
    TV jaw_midpoint,jaw_normal,jaw_axis,jaw_front;
    T jaw_axis_length,jaw_sliding_length,max_opening_angle;
    FRAME<TV> cranium_frame_save,jaw_frame_save;

    ATTACHMENT_FRAME_CONTROL_SET(ARRAY<TV>& X_input,ARRAY<ARRAY<int> >& attached_nodes_input,const int jaw_attachment_index_input);
    ~ATTACHMENT_FRAME_CONTROL_SET();

    FRAME<TV> Cranium_Frame()
    {return cranium_transform.Frame()*cranium_frame_save;}

    FRAME<TV> Jaw_Frame()
    {return QUASI_RIGID_TRANSFORM_3D<T>::Composite_Transform(cranium_transform,jaw_transform).Frame()*jaw_frame_save;}

//#####################################################################
    int Size() const PHYSBAM_OVERRIDE;
    T_TYPE Type() const PHYSBAM_OVERRIDE;
    T operator()(const int control_id) const PHYSBAM_OVERRIDE;
    T& operator()(const int control_id) PHYSBAM_OVERRIDE;
    T Identity(const int control_id) const PHYSBAM_OVERRIDE;
    void Maximal_Controls() PHYSBAM_OVERRIDE;
    bool Control_Active(const int control_id) const PHYSBAM_OVERRIDE;
    bool& Control_Active(const int control_id) PHYSBAM_OVERRIDE;
    bool Positions_Determined_Kinematically(const int control_id) const PHYSBAM_OVERRIDE;
    void Force_Derivative(ARRAY<TV>& dFdl,ARRAY<TWIST<TV> >& dFrdl,const int control_id) const PHYSBAM_OVERRIDE;
    void Position_Derivative(ARRAY<TV>& dXdl,const int control_id) const PHYSBAM_OVERRIDE;
    void Set_Attachment_Positions(ARRAY<TV>&X) const PHYSBAM_OVERRIDE;
    void Save_Controls() PHYSBAM_OVERRIDE;
    void Kinematically_Update_Positions(ARRAY<TV>&X) const PHYSBAM_OVERRIDE;
    void Kinematically_Update_Jacobian(ARRAY<TV>&dX) const PHYSBAM_OVERRIDE;
    T Penalty() const PHYSBAM_OVERRIDE;
    T Penalty_Gradient(const int control_id) const PHYSBAM_OVERRIDE;
    T Penalty_Hessian(const int control_id1,const int control_id2) const PHYSBAM_OVERRIDE;
    T Jaw_Constraint_Penalty() const;
    T Jaw_Constraint_Penalty_Gradient(const int control_id) const;
    T Jaw_Constraint_Penalty_Hessian(const int control_id1,const int control_id2) const;
    void Read_Jaw_Joint_From_File(const STREAM_TYPE& stream_type,const std::string& filename);
    void Write_Jaw_Joint_To_File(const STREAM_TYPE& stream_type,const std::string& filename) const;
    void Read_Configuration(const STREAM_TYPE& stream_type,std::istream& input_stream);
    void Write_Configuration(const STREAM_TYPE& stream_type,std::ostream& output_stream) const;
    void Set_Original_Attachment_Configuration(const FRAME<TV>& cranium_frame_input,const FRAME<TV>& jaw_frame_input);
    void Project_Parameters_To_Allowable_Range(const bool active_controls_only) PHYSBAM_OVERRIDE;
    void Set_From_Generalized_Coordinates(const T left_sliding_parameter,const T right_sliding_parameter,const T opening_angle);
    T Opening_Angle() const;
    void Increment_Opening_Angle(const T angle);
    T Left_Condyle_Sliding_Parameter() const;
    T Right_Condyle_Sliding_Parameter() const;
    void Interpolate(const T interpolation_fraction) PHYSBAM_OVERRIDE;
    void Scale(const T scale) PHYSBAM_OVERRIDE;
    T Distance(const ARRAY<T>& weights,int base) PHYSBAM_OVERRIDE;
    void Print_Diagnostics(std::ostream& output) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
