//#####################################################################
// Copyright 2006-2007, Kevin Der, Ranjitha Kumar, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_SURFACE_MUSCLE_SEGMENT
//#####################################################################
#ifndef __ANALYTIC_SURFACE_MUSCLE_SEGMENT__
#define __ANALYTIC_SURFACE_MUSCLE_SEGMENT__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_SEGMENT.h>
namespace PhysBAM{

template<class T> class TETRAHEDRALIZED_VOLUME;
template<class TV> class ATTACHMENT_POINT;
template<class T_input>
class ANALYTIC_SURFACE_MUSCLE_SEGMENT:public MUSCLE_SEGMENT<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef MUSCLE_SEGMENT<TV> BASE;
    using BASE::point_1;using BASE::point_2;using BASE::segment_type;using BASE::frame;using BASE::activations;using BASE::activation_memory;
    using BASE::Length;using BASE::Set_Current_Activation;

    enum CURVE_TYPE {CURVE_NONE=1,CURVE_COSINE,CURVE_ELLIPSE};
    CURVE_TYPE curve_type;
    T initial_volume,initial_length,initial_curve_thickness,initial_curve_offset_thickness,curve_thickness,curve_offset_thickness;
    T tendon_fraction_1,tendon_fraction_2,tendon_slack_length,activation_factor;
    static const T small_number;
    ARRAY<ARRAY<int> > muscle_tets;
    ARRAY<TV> inside_particle_rest_positions;
    T initial_curve_length;
    int num_segments_over_2;
    ARRAY<T> initial_segment_volumes, fractional_segment_endpoints;
    ARRAY<int> inside_particle_segments;

    ANALYTIC_SURFACE_MUSCLE_SEGMENT();
    ANALYTIC_SURFACE_MUSCLE_SEGMENT(CURVE_TYPE curve_type_input,ATTACHMENT_POINT<TV>* point_1_input,ATTACHMENT_POINT<TV>* point_2_input,
        const T curve_thickness_input,const T curve_offset_thickness_input,const T tendon_slack_length_input,const T tendon_fraction_1_input,const T tendon_fraction_2_input);
    virtual ~ANALYTIC_SURFACE_MUSCLE_SEGMENT();

//#####################################################################
    void Initialize() PHYSBAM_OVERRIDE;
    void Read_And_Set_Parameters(TYPED_ISTREAM& input) PHYSBAM_OVERRIDE;
    void Write_Parameters(TYPED_OSTREAM& output) const PHYSBAM_OVERRIDE;
    T Maximum_Radius() const PHYSBAM_OVERRIDE;
    static ANALYTIC_SURFACE_MUSCLE_SEGMENT* Create();
    virtual std::string Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    static std::string Static_Name();
    T Curve_Length() const;
    T Tendon_Length() const;
    T Get_Fractional_Curve_Value(const T fraction,const bool initial=false) const;
    T Compute_Volume() const;
    void Set_Current_Activation(const T activation) PHYSBAM_OVERRIDE;
    void Update_Parameters() PHYSBAM_OVERRIDE;
    bool Analytic_Inside_Test(const TV& local_position) const;
    T Get_Elongation_At_Local_Position(const TV& initial_local_position,const int segment_number) const;
    void Get_Local_Positions_For_Particles(const int m,const int n,GEOMETRY_PARTICLES<TV>& particles);
    void Initialize_Inside_Particles(const TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume);
    T Scaled_Volume(const T segment_endpoint) const;
    T Scaled_Volume_Derivative(const T segment_endpoint) const;
//#####################################################################
};
template<class T> const T ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::small_number=(T)1e-6;
}
#endif
