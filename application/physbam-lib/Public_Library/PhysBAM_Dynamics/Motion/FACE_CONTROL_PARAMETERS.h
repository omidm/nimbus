//#####################################################################
// Copyright 2004-2007, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_CONTROL_PARAMETERS
//#####################################################################
#ifndef __FACE_CONTROL_PARAMETERS__    
#define __FACE_CONTROL_PARAMETERS__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
namespace PhysBAM{

template<class T> class FACE_CONTROL_SET;
template<class T>
class FACE_CONTROL_PARAMETERS
{
    typedef VECTOR<T,3> TV;
public:
    ARRAY<FACE_CONTROL_SET<T>*> list;

    FACE_CONTROL_PARAMETERS();
    ~FACE_CONTROL_PARAMETERS();
    int Size() const;
    int Active_Size() const;
    int Active_Kinematic_Size() const;
    int Active_Nonkinematic_Size() const;
    void Seek(const int control_index,int& set_subindex,int& control_subindex) const;
    T& operator()(const int control_index);
    T operator()(const int control_index) const;
    bool& Active(const int control_index);
    bool Active(const int control_index) const;
    bool Kinematic(const int control_index) const;
    bool Active_Kinematic(const int control_index) const;
    bool Active_Nonkinematic(const int control_index) const;
    ARRAY<int> Active_Subset() const;
    ARRAY<int> Active_Kinematic_Subset() const;
    ARRAY<int> Active_Nonkinematic_Subset() const;
    void Get(VECTOR_ND<T>& values) const;
    void Get(VECTOR_ND<T>& values,const ARRAY<int>& subset) const;
    void Get_Active(ARRAY<bool>& active) const;
    void Get_Active(ARRAY<bool>& active,const ARRAY<int>& subset) const;
    void Set(const VECTOR_ND<T>& values);
    void Set(const VECTOR_ND<T>& values,const ARRAY<int>& subset);
    void Set(const T value_input,const ARRAY<int>& subset);
    void Set_Active(const ARRAY<bool>& active);
    void Set_Active(const ARRAY<bool>& active,const ARRAY<int>& subset);
    void Set_Active(const bool active_input,const ARRAY<int>& subset);
    T Penalty() const;
    VECTOR_ND<T> Penalty_Gradient() const;
    MATRIX_MXN<T> Penalty_Hessian() const;
    void Save_Controls();
    void Print_Diagnostics(std::ostream& output) const;
    void Kinematically_Update_Positions(ARRAY<TV>&X) const;
    void Kinematically_Update_Jacobian(ARRAY<TV>&dX) const;
    void Force_Derivative(ARRAY<TV>& dFdl,ARRAY<TWIST<TV> >& dFrdl,const int control_index) const;
    void Position_Derivative(ARRAY<TV>& dXdl,const int control_index) const;
    void Project_Parameters_To_Allowable_Range(const bool active_controls_only=false);
    void Interpolate(const T interpolation_fraction);
    void Scale(const T scale);
    T Distance(const ARRAY<T>& weights);
    T Identity(const int control_index) const;
    void Identity_Controls(VECTOR_ND<T>& values) const;
    void Maximal_Controls();
    template<class RW> void Read(std::istream& input_stream);
    template<class RW> void Write(std::ostream& output_stream) const;
    void Read_Configuration_From_File(const STREAM_TYPE& stream_type,const std::string filename);
    void Write_Configuration_To_File(const STREAM_TYPE& stream_type,const std::string filename) const;
    void Write_Control_Set_Types(const STREAM_TYPE& stream_type,std::ostream& output_stream) const;
//#####################################################################
};
}
#endif
