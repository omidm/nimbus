//#####################################################################
// Copyright 2004-2007, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Dynamics/Motion/ACTIVATION_CONTROL_SET.h>
#include <PhysBAM_Dynamics/Motion/ATTACHMENT_FRAME_CONTROL_SET.h>
#include <PhysBAM_Dynamics/Motion/FACE_CONTROL_PARAMETERS.h>
#include <PhysBAM_Dynamics/Motion/FACE_CONTROL_SET.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> FACE_CONTROL_PARAMETERS<T>::
FACE_CONTROL_PARAMETERS()
{
}
template<class T> FACE_CONTROL_PARAMETERS<T>::
~FACE_CONTROL_PARAMETERS()
{
    list.Delete_Pointers_And_Clean_Memory();
}
template<class T> int FACE_CONTROL_PARAMETERS<T>::
Size() const
{
    int n=0;
    for(int s=1;s<=list.m;s++) n+=list(s)->Size();
    return n;
}
template<class T> int FACE_CONTROL_PARAMETERS<T>::
Active_Size() const
{
    int n=0;
    for(int s=1;
    s<=list.m;
    s++) for(int c=1;c<=list(s)->Size();c++) if(list(s)->Control_Active(c)) n++;
    return n;
}
template<class T> int FACE_CONTROL_PARAMETERS<T>::
Active_Kinematic_Size() const
{
    int n=0;
    for(int s=1;
    s<=list.m;
    s++) for(int c=1;c<=list(s)->Size();c++) if(list(s)->Control_Active(c)&&list(s)->Positions_Determined_Kinematically(c)) n++;
    return n;
}
template<class T> int FACE_CONTROL_PARAMETERS<T>::
Active_Nonkinematic_Size() const
{
    int n=0;
    for(int s=1;s<=list.m;s++) for(int c=1;c<=list(s)->Size();c++) if(list(s)->Control_Active(c)&&!list(s)->Positions_Determined_Kinematically(c)) n++;
    return n;
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Seek(const int control_index,int& set_subindex,int& control_subindex) const
{
    assert(1<=control_index&&control_index<=Size());
    control_subindex=control_index;
    for(set_subindex=1;control_subindex>list(set_subindex)->Size();set_subindex++) control_subindex-=list(set_subindex)->Size();
}
template<class T> T& FACE_CONTROL_PARAMETERS<T>::
operator()(const int control_index)
{
    int s,c;
    Seek(control_index,s,c);
    return (*list(s))(c);
}
template<class T> T FACE_CONTROL_PARAMETERS<T>::
operator()(const int control_index) const
{
    int s,c;
    Seek(control_index,s,c);
    return (*list(s))(c);
}
template<class T> bool& FACE_CONTROL_PARAMETERS<T>::
Active(const int control_index)
{
    int s,c;
    Seek(control_index,s,c);
    return list(s)->Control_Active(c);
}
template<class T> bool FACE_CONTROL_PARAMETERS<T>::
Active(const int control_index) const
{
    int s,c;
    Seek(control_index,s,c);
    return list(s)->Control_Active(c);
}
template<class T> bool FACE_CONTROL_PARAMETERS<T>::
Kinematic(const int control_index) const
{
    int s,c;
    Seek(control_index,s,c);
    return list(s)->Positions_Determined_Kinematically(c);
}
template<class T> bool FACE_CONTROL_PARAMETERS<T>::
Active_Kinematic(const int control_index) const
{
    int s,c;
    Seek(control_index,s,c);
    return list(s)->Control_Active(c)&&list(s)->Positions_Determined_Kinematically(c);
}
template<class T> bool FACE_CONTROL_PARAMETERS<T>::
Active_Nonkinematic(const int control_index) const
{
    int s,c;
    Seek(control_index,s,c);
    return list(s)->Control_Active(c)&&!list(s)->Positions_Determined_Kinematically(c);
}
template<class T> ARRAY<int> FACE_CONTROL_PARAMETERS<T>::
Active_Subset() const
{
    ARRAY<int> result(Active_Size());
    int n=0;
    for(int i=1;i<=Size();i++) if(Active(i)) result(++n)=i;
    return result;
}
template<class T> ARRAY<int> FACE_CONTROL_PARAMETERS<T>::
Active_Kinematic_Subset() const
{
    ARRAY<int> result(Active_Kinematic_Size());
    int n=0;
    for(int i=1;i<=Size();i++) if(Active_Kinematic(i)) result(++n)=i;
    return result;
}
template<class T> ARRAY<int> FACE_CONTROL_PARAMETERS<T>::
Active_Nonkinematic_Subset() const
{
    ARRAY<int> result(Active_Nonkinematic_Size());
    int n=0;
    for(int i=1;i<=Size();i++) if(Active_Nonkinematic(i)) result(++n)=i;
    return result;
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Get(VECTOR_ND<T>& values) const
{
    values=VECTOR_ND<T>(Size());
    for(int i=1;i<=Size();i++) values(i)=(*this)(i);
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Get(VECTOR_ND<T>& values,const ARRAY<int>& subset) const
{
    values=VECTOR_ND<T>(subset.m);
    for(int i=1;i<=subset.m;i++) values(i)=(*this)(subset(i));
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Get_Active(ARRAY<bool>& active) const
{
    active.Resize(Size());
    for(int i=1;i<=Size();i++) active(i)=Active(i);
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Get_Active(ARRAY<bool>& active,const ARRAY<int>& subset) const
{
    active.Resize(subset.m);
    for(int i=1;i<=subset.m;i++) active(i)=Active(subset(i));
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Set(const VECTOR_ND<T>& values)
{
    assert(values.n==Size());
    for(int i=1;i<=Size();i++) (*this)(i)=values(i);
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Set(const VECTOR_ND<T>& values,const ARRAY<int>& subset)
{
    assert(values.n==subset.m);
    for(int i=1;i<=subset.m;i++) (*this)(subset(i))=values(i);
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Set(const T value_input,const ARRAY<int>& subset)
{
    for(int i=1;i<=subset.m;i++) (*this)(subset(i))=value_input;
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Set_Active(const ARRAY<bool>& active)
{
    assert(active.m==Size());
    for(int i=1;i<=Size();i++) Active(i)=active(i);
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Set_Active(const ARRAY<bool>& active,const ARRAY<int>& subset)
{
    assert(active.m==subset.m);
    for(int i=1;i<=subset.m;i++) Active(subset(i))=active(i);
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Set_Active(const bool active_input,const ARRAY<int>& subset)
{
    for(int i=1;i<=subset.m;i++) Active(subset(i))=active_input;
}
template<class T> T FACE_CONTROL_PARAMETERS<T>::
Penalty() const
{
    T result=0;
    for(int s=1;s<=list.m;s++) result+=list(s)->Penalty();
    return result;
}
template<class T> VECTOR_ND<T> FACE_CONTROL_PARAMETERS<T>::
Penalty_Gradient() const
{
    ARRAY<int> subset(Active_Subset());
    VECTOR_ND<T> result(subset.m);
    for(int i=1;i<=subset.m;i++){
        int s,c;
        Seek(subset(i),s,c);
        result(i)=list(s)->Penalty_Gradient(c);}
    return result;
}
template<class T> MATRIX_MXN<T> FACE_CONTROL_PARAMETERS<T>::
Penalty_Hessian() const
{
    ARRAY<int> subset(Active_Subset());
    MATRIX_MXN<T> result(subset.m,subset.m);
    for(int i1=1;i1<=subset.m;i1++) for(int i2=1;i2<=subset.m;i2++){
        int s1,c1,s2,c2;
        Seek(subset(i1),s1,c1);
        Seek(subset(i2),s2,c2);
        if(s1==s2) result(i1,i2)=list(s1)->Penalty_Hessian(c1,c2);}
    return result;
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Save_Controls()
{
    for(int s=1;s<=list.m;s++) list(s)->Save_Controls();
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Print_Diagnostics(std::ostream& output) const
{
    for(int s=1;s<=list.m;s++) list(s)->Print_Diagnostics(output);
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Kinematically_Update_Positions(ARRAY<TV>&X) const
{
    for(int s=1;s<=list.m;s++) list(s)->Kinematically_Update_Positions(X);
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Kinematically_Update_Jacobian(ARRAY<TV>&dX) const
{
    for(int s=1;s<=list.m;s++) list(s)->Kinematically_Update_Jacobian(dX);
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Force_Derivative(ARRAY<TV>& dFdl,ARRAY<TWIST<TV> >& dFrdl,const int control_index) const
{
    int s,c;
    Seek(control_index,s,c);
    list(s)->Force_Derivative(dFdl,dFrdl,c);
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Position_Derivative(ARRAY<TV>& dXdl,const int control_index) const
{
    int s,c;
    Seek(control_index,s,c);
    list(s)->Position_Derivative(dXdl,c);
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Project_Parameters_To_Allowable_Range(const bool active_controls_only)
{
    for(int s=1;s<=list.m;s++) list(s)->Project_Parameters_To_Allowable_Range(active_controls_only);
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Interpolate(const T interpolation_fraction)
{
    for(int s=1;s<=list.m;s++) list(s)->Interpolate(interpolation_fraction);
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Scale(const T scale)
{
    if(scale==1)return;
    for(int s=1;s<=list.m;s++) list(s)->Scale(scale);
}
template<class T> T FACE_CONTROL_PARAMETERS<T>::
Distance(const ARRAY<T>& weights)
{
    T distance=0;
    int base=0;
    for(int s=1;s<=list.m;s++){
        distance+=list(s)->Distance(weights,base);
        base+=list(s)->Size();}
    return distance;
}
template<class T> T FACE_CONTROL_PARAMETERS<T>::
Identity(const int control_index) const
{
    int s,c;
    Seek(control_index,s,c);
    return (*list(s)).Identity(c);
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Identity_Controls(VECTOR_ND<T>& values) const
{
    values=VECTOR_ND<T>(Size());
    for(int i=1;i<=Size();i++) values(i)=Identity(i);
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Maximal_Controls()
{
    for(int s=1;
    s<=list.m;
    s++)list(s)->Maximal_Controls();
}
template<class T> template<class RW> void FACE_CONTROL_PARAMETERS<T>::
Read(std::istream& input_stream)
{
    VECTOR_ND<T> values;
    ARRAY<bool> active;
    Read_Binary<RW>(input_stream,values,active);
    assert(values.n==Size()&&active.m==Size());
    Set(values);
    Set_Active(active);
}
template<class T> template<class RW> void FACE_CONTROL_PARAMETERS<T>::
Write(std::ostream& output_stream) const
{
    VECTOR_ND<T> values;
    ARRAY<bool> active;
    Get(values);
    Get_Active(active);
    Write_Binary<RW>(output_stream,values,active);
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Read_Configuration_From_File(const STREAM_TYPE& stream_type,const std::string filename)
{
    list.Delete_Pointers_And_Clean_Memory();
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename);
    TYPED_ISTREAM typed_input(*input,stream_type);
    int number_of_control_sets;
    Read_Binary(typed_input,number_of_control_sets);
    ARRAY<int> control_set_types(number_of_control_sets);
    for(int i=1;i<=number_of_control_sets;i++) Read_Binary(typed_input,control_set_types(i));
    for(int i=1;i<=number_of_control_sets;i++){
        switch(control_set_types(i)){
            case FACE_CONTROL_SET<T>::ACTIVATION:{
                ACTIVATION_CONTROL_SET<T>* activation_control_set=new ACTIVATION_CONTROL_SET<T>();
                activation_control_set->Read_Configuration(stream_type,*input);
                list.Append(activation_control_set);
                break;}
            case FACE_CONTROL_SET<T>::ATTACHMENT_FRAME:{
                ATTACHMENT_FRAME_CONTROL_SET<T>* attachment_frame_control_set=new ATTACHMENT_FRAME_CONTROL_SET<T>(*(new ARRAY<TV>),*(new ARRAY<ARRAY<int> >),0);
                attachment_frame_control_set->Read_Configuration(stream_type,*input);
                list.Append(attachment_frame_control_set);
                break;}
            default: throw READ_ERROR("Unsupported FACE_CONTROL_SET type");}}
    delete input;
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Write_Configuration_To_File(const STREAM_TYPE& stream_type,const std::string filename) const
{
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(filename);
    Write_Control_Set_Types(stream_type,*output);
    for(int i=1;i<=list.m;i++)
        switch(list(i)->Type()){
            case FACE_CONTROL_SET<T>::ACTIVATION: ((ACTIVATION_CONTROL_SET<T>*)list(i))->Write_Configuration(stream_type,*output);break;
            case FACE_CONTROL_SET<T>::ATTACHMENT_FRAME: ((ATTACHMENT_FRAME_CONTROL_SET<T>*)list(i))->Write_Configuration(stream_type,*output);break;
            default: PHYSBAM_FATAL_ERROR("Unsupported FACE_CONTROL_SET type");}
    delete output;
}
template<class T> void FACE_CONTROL_PARAMETERS<T>::
Write_Control_Set_Types(const STREAM_TYPE& stream_type,std::ostream& output_stream) const
{
    TYPED_OSTREAM typed_output(output_stream,stream_type);
    Write_Binary(typed_output,list.m);
    for(int i=1;i<=list.m;i++)Write_Binary(typed_output,(int)list(i)->Type());
}
template class FACE_CONTROL_PARAMETERS<float>;
template void FACE_CONTROL_PARAMETERS<float>::Read<float>(std::basic_istream<char,std::char_traits<char> >&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FACE_CONTROL_PARAMETERS<double>;
template void FACE_CONTROL_PARAMETERS<double>::Read<double>(std::basic_istream<char,std::char_traits<char> >&);
#endif
