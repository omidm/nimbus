//#####################################################################
// Copyright 2005-2007, Kevin Der, Ron Fedkiw, Eran Guendelman, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MUSCLE
//#####################################################################
#ifndef __MUSCLE__
#define __MUSCLE__

#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/ANALYTIC_SURFACE_MUSCLE_SEGMENT.h>
namespace PhysBAM{

template<class TV> class MUSCLE_SEGMENT;
template<class TV> class MUSCLE_FORCE_CURVE;
template<class TV> class RIGID_BODY_COLLECTION;
template<class TV> class RIGID_BODY;
template<class TV> class ATTACHMENT_POINT;

template<class TV>
class MUSCLE
{
    typedef typename TV::SCALAR T;
    typedef typename MUSCLE_SEGMENT<TV>::MUSCLE_SEGMENT_TYPE T_MUSCLE_SEGMENT_TYPE;
    typedef typename ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::CURVE_TYPE T_MUSCLE_SEGMENT_CURVE_TYPE;
    typedef TRIPLE<T_MUSCLE_SEGMENT_TYPE,T_MUSCLE_SEGMENT_CURVE_TYPE,ARRAY<T> > T_MUSCLE_SEGMENT_DATA;
public:
    typedef int HAS_TYPED_READ_WRITE;

    const MUSCLE_FORCE_CURVE<T>& force_curve;
    ATTACHMENT_POINT<TV> *attachment_point_1,*attachment_point_2;
    T optimal_length,peak_force,pennation_angle,tendon_slack_length,max_shortening_velocity;
    ARRAY<ATTACHMENT_POINT<TV>*> via_points;
    ARRAY<MUSCLE_SEGMENT<TV>*> muscle_segments;
    int id;
    TV local_fiber_direction;
    std::string name;
protected:
    T peak_force_times_cos_pennation_angle,length_normalization,velocity_normalization,tendon_normalization,one_over_peak_force,cos_pennation_angle,one_over_cos_pennation_angle; // acceleration constants

public:
    MUSCLE(const MUSCLE_FORCE_CURVE<T>& force_curve_input);

    virtual ~MUSCLE(); // assume we own the constrained points!

    virtual void Clean_Memory();

    void Set_Attachment_Point_1(ATTACHMENT_POINT<TV>* attachment_point_1_input)
    {attachment_point_1=attachment_point_1_input;}

    void Set_Attachment_Point_2(ATTACHMENT_POINT<TV>* attachment_point_2_input)
    {attachment_point_2=attachment_point_2_input;}

    void Set_Optimal_Length(const T optimal_length_input)
    {optimal_length=optimal_length_input;length_normalization=(T)1/optimal_length;Update_Acceleration_Constants();}
    
    void Set_Peak_Force(const T peak_force_input)
    {peak_force=peak_force_input;one_over_peak_force=(T)1/peak_force;Update_Acceleration_Constants();}

    void Set_Pennation_Angle(const T pennation_angle_input)
    {pennation_angle=pennation_angle_input;cos_pennation_angle=cos(pennation_angle);one_over_cos_pennation_angle=(T)1/cos_pennation_angle;Update_Acceleration_Constants();}

    void Set_Tendon_Slack_Length(const T tendon_slack_length_input)
    {tendon_slack_length=tendon_slack_length_input;tendon_normalization=(T)1/tendon_slack_length;Update_Acceleration_Constants();}

    void Set_Max_Shortening_Velocity(const T max_shortening_velocity_input)
    {max_shortening_velocity=max_shortening_velocity_input;Update_Acceleration_Constants();}

    void Update_Acceleration_Constants() // Update constants that depend on more than one parameter
    {peak_force_times_cos_pennation_angle=peak_force*cos_pennation_angle;
    if(optimal_length!=0 && max_shortening_velocity!=0) velocity_normalization=(T)1/(max_shortening_velocity*optimal_length);}

    void Add_Via_Point(ATTACHMENT_POINT<TV>* via_point)
    {via_points.Append(via_point);} // at the moment these must be added in order from attachment point 1

    T Passive_Force(const T muscle_length) const
    {return force_curve.Passive_Force(muscle_length*length_normalization)*peak_force_times_cos_pennation_angle;}

    T Passive_Force_Slope(const T muscle_length) const
    {return force_curve.Passive_Force_Slope(muscle_length*length_normalization)*peak_force_times_cos_pennation_angle*length_normalization;}

    T Peak_Active_Force(const T muscle_length,const T muscle_velocity) const // scaled by activation when actually used
    {return force_curve.Active_Force(muscle_length*length_normalization)*force_curve.Velocity_Scale(muscle_velocity*velocity_normalization)*peak_force_times_cos_pennation_angle;}

    T Peak_Active_Force_Quasistatic(const T muscle_length) const
    {return force_curve.Active_Force(muscle_length*length_normalization)*peak_force_times_cos_pennation_angle;}

    T Peak_Active_Force_Quasistatic_Slope(const T muscle_length) const // does not include force velocity relationship
    {return force_curve.Active_Force_Slope(muscle_length*length_normalization)*peak_force_times_cos_pennation_angle*length_normalization;}

    T Velocity_Scale(const T muscle_velocity) const // multiplies Peak_Active_Force_Quasistatic to give Peak_Active_Force
    {return force_curve.Velocity_Scale(muscle_velocity*velocity_normalization);}

    T Tendon_Force(const T tendon_length) const
    {return force_curve.Tendon_Force((tendon_length-tendon_slack_length)*tendon_normalization)*peak_force;}

    T Tendon_Force_Slope(const T tendon_length) const
    {return force_curve.Tendon_Force_Slope((tendon_length-tendon_slack_length)*tendon_normalization)*peak_force*tendon_normalization;}

    T Total_Force(const T muscle_length,const T muscle_velocity,const T activation,const T tendon_length) const
    {return activation*Peak_Active_Force(muscle_length,muscle_velocity)+Passive_Force(muscle_length)+Tendon_Force(tendon_length);}

    virtual T Force(const T activation) const // TODO: this hack assumes that there are no tendons
    {T muscle_length=Total_Length(),muscle_velocity=Total_Velocity(),tendon_length=0; // TODO: should use nonlinear root finding to separate Total_Length into muscle and tendon length
    return Total_Force(muscle_length,muscle_velocity,activation,tendon_length);}

    void Set_Name(const std::string& name_input)
    {name=name_input;}

    void Write_Constrained_Point(TYPED_OSTREAM& output_stream,const ATTACHMENT_POINT<TV>* constrained_point) const
    {Write_Binary(output_stream,constrained_point->object_space_position,constrained_point->rigid_body.particle_index);}

//#####################################################################
    void Initialize();
    void Initialize(ARRAY<T_MUSCLE_SEGMENT_DATA>& segment_data);
    MUSCLE_SEGMENT<TV>* Create_Muscle_Segment(ATTACHMENT_POINT<TV>* point_1,ATTACHMENT_POINT<TV>* point_2,const T_MUSCLE_SEGMENT_TYPE& segment_type,
        const T_MUSCLE_SEGMENT_CURVE_TYPE& curve_type,const ARRAY<T>& parameters);
    T Total_Length() const;
    T Total_Velocity() const;
    T Calculate_Activation(const T desired_force) const;
    void Update_Segments();
    void Set_Segment_Activations(const T activation);
    void Apply_Fixed_Impulse_At_All_Points(const T impulse);
    ATTACHMENT_POINT<TV>* Read_Constrained_Point(TYPED_ISTREAM& input_stream,RIGID_BODY_COLLECTION<TV>& rigid_body_collection);
    virtual void Read(TYPED_ISTREAM& input_stream,RIGID_BODY_COLLECTION<TV>& rigid_body_collection); // Currently assumes attachment points are also constrained point in rigid body!!
    virtual void Write(TYPED_OSTREAM& output_stream) const;
//#####################################################################
};
}
#endif
