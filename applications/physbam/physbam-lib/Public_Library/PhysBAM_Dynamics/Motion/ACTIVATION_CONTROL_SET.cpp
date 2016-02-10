//#####################################################################
// Copyright 2004-2005, Igor Neverov, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Dynamics/Motion/ACTIVATION_CONTROL_SET.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> ACTIVATION_CONTROL_SET<T>::
ACTIVATION_CONTROL_SET()
    :single_activation_used_for_force_derivative(0),muscle_force(0),activation_cutoff((T)1),max_optimization_step_length((T).4),activation_penalty_coefficient((T)10)
{
    Save_Controls();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> ACTIVATION_CONTROL_SET<T>::
~ACTIVATION_CONTROL_SET()
{
}
template<class T> int ACTIVATION_CONTROL_SET<T>::
Add_Activation(const std::string name,const bool active)
{
    assert(activations.m==activation_active.m);
    activations.Append(0),activation_active.Append(active);
    activation_names.Append(name);
    Save_Controls();
    return activations.m;
}
template<class T> int ACTIVATION_CONTROL_SET<T>::
Size() const
{
    return activations.m;
}
template<class T> typename FACE_CONTROL_SET<T>::TYPE ACTIVATION_CONTROL_SET<T>::
Type() const
{
    return ACTIVATION;
}
template<class T> T ACTIVATION_CONTROL_SET<T>::
operator()(const int control_id) const
{
    return activations(control_id);
}
template<class T> T& ACTIVATION_CONTROL_SET<T>::
operator()(const int control_id)
{
    return activations(control_id);
}
template<class T> T ACTIVATION_CONTROL_SET<T>::
Identity(const int control_id) const
{
    return 0;
}
template<class T> void ACTIVATION_CONTROL_SET<T>::
Maximal_Controls()
{
    for(int i=1;i<=activations.m;i++) activations(i)=max(activations(i),activations_save(i));
}
template<class T> bool ACTIVATION_CONTROL_SET<T>::
Control_Active(const int control_id) const
{
    return activation_active(control_id);
}
template<class T> bool& ACTIVATION_CONTROL_SET<T>::
Control_Active(const int control_id)
{
    return activation_active(control_id);
}
template<class T> bool ACTIVATION_CONTROL_SET<T>::
Positions_Determined_Kinematically(const int control_id) const
{
    return false;
}
template<class T> void ACTIVATION_CONTROL_SET<T>::
Force_Derivative(ARRAY<TV>& dFdl,ARRAY<TWIST<TV> >& dFrdl,const int control_id) const
{
    assert(muscle_force);
    assert(Control_Active(control_id));
    single_activation_used_for_force_derivative=control_id;
    ARRAYS_COMPUTATIONS::Fill(dFdl,TV());
    ARRAYS_COMPUTATIONS::Fill(dFrdl,TWIST<TV>());
    muscle_force->Add_Velocity_Independent_Forces(dFdl,dFrdl,0);
    single_activation_used_for_force_derivative=0;
}
template<class T> T ACTIVATION_CONTROL_SET<T>::
Penalty() const
{
    assert(activations.m==activations_save.m);
    T penalty=0;
    for(int i=1;i<=activations.m;i++) penalty+=Piecewise_Quadratic_Penalty(activations(i),clamp<T>(activations_save(i)-max_optimization_step_length,0,activation_cutoff),
        clamp<T>(activations_save(i)+max_optimization_step_length,0,activation_cutoff));
    return activation_penalty_coefficient*penalty;
}
template<class T> T ACTIVATION_CONTROL_SET<T>::
Penalty_Gradient(const int control_id) const
{
    assert(1<=control_id && control_id<=Size());
    return activation_penalty_coefficient*Piecewise_Quadratic_Penalty_Prime(activations(control_id),clamp<T>(activations_save(control_id)-max_optimization_step_length,0,activation_cutoff),
        clamp<T>(activations_save(control_id)+max_optimization_step_length,0,activation_cutoff));
}
template<class T> T ACTIVATION_CONTROL_SET<T>::
Penalty_Hessian(const int control_id1,const int control_id2) const
{
    assert(1<=control_id1 && control_id1<=Size() && 1<=control_id2 && control_id2<=Size());
    if(control_id1!=control_id2) return 0;
    return activation_penalty_coefficient*Piecewise_Quadratic_Penalty_Double_Prime(activations(control_id1),
        clamp<T>(activations_save(control_id1)-max_optimization_step_length,0,activation_cutoff),
        clamp<T>(activations_save(control_id1)+max_optimization_step_length,0,activation_cutoff));
}
template<class T> void ACTIVATION_CONTROL_SET<T>::
Save_Controls()
{
    activations_save=activations;
}
template<class T> void ACTIVATION_CONTROL_SET<T>::
Project_Parameters_To_Allowable_Range(const bool active_controls_only)
{
    for(int i=1;i<=activations.m;i++) if(!active_controls_only || activation_active(i))
        activations(i)=clamp<T>(activations(i),max(activations_save(i)-max_optimization_step_length,(T)0),min(activations_save(i)+max_optimization_step_length,activation_cutoff));
}
template<class T> void ACTIVATION_CONTROL_SET<T>::
Interpolate(const T interpolation_fraction)
{
    for(int i=1;i<=activations.m;i++) activations(i)=LINEAR_INTERPOLATION<T,T>::Linear(activations_save(i),activations(i),interpolation_fraction);
}
template<class T> void ACTIVATION_CONTROL_SET<T>::
Scale(const T scale)
{
    for(int i=1;i<=activations.m;i++) activations(i)=min((T)1,activations(i)*scale);
}
template<class T> T ACTIVATION_CONTROL_SET<T>::
Distance(const ARRAY<T>& weights,int base)
{
    T distance=0;
    for(int i=1;i<=activations.m;i++){
        T mean=(T).5*abs((activations_save(i)+activations(i)));
        distance+=weights(base+i)*pow(activations_save(i)-activations(i),2)/(mean?mean:1);}
    return sqrt(distance);
}
template<class T> void ACTIVATION_CONTROL_SET<T>::Print_Diagnostics(std::ostream& output) const
{
    for(int i=1;i<=activations.m;i++) output<<STRING_UTILITIES::string_sprintf("%2d - %-35s [%s] : %g\n",i,activation_names(i).c_str(),(activation_active(i)?"active":"inactive"),activations(i));
    output<<"Penalty at current activation levels : "<<Penalty()<<std::endl;
}
template<class T> T ACTIVATION_CONTROL_SET<T>::
Piecewise_Quadratic_Penalty(const T a,const T a_min,const T a_max)
{
    if(a<a_min) return sqr(a-a_min);
    if(a<=a_max) return 0;
    return sqr(a-a_max);
}
template<class T> T ACTIVATION_CONTROL_SET<T>::
Piecewise_Quadratic_Penalty_Prime(const T a,const T a_min,const T a_max)
{
    if(a<a_min) return 2*(a-a_min);
    if(a<=a_max) return 0;
    return 2*(a-a_max);
}
template<class T> T ACTIVATION_CONTROL_SET<T>::
Piecewise_Quadratic_Penalty_Double_Prime(const T a,const T a_min,const T a_max)
{
    if(a<a_min) return 2;
    if(a<=a_max) return 0;
    return 2;
}
template<class T> void ACTIVATION_CONTROL_SET<T>::
Read_Configuration(const STREAM_TYPE& stream_type,std::istream& input_stream)
{
    TYPED_ISTREAM typed_input(input_stream,stream_type);
    Read_Binary(typed_input,activation_names);
    activations.Resize(activation_names.m);activations_save.Resize(activation_names.m),activation_active.Resize(activation_names.m);
}
template<class T> void ACTIVATION_CONTROL_SET<T>::
Write_Configuration(const STREAM_TYPE& stream_type,std::ostream& output_stream) const
{
    TYPED_OSTREAM typed_output(output_stream,stream_type);
    Write_Binary(typed_output,activation_names);
}
template class ACTIVATION_CONTROL_SET<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ACTIVATION_CONTROL_SET<double>;
#endif
