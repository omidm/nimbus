//#####################################################################
// Copyright 2006-2008, Ronald Fedkiw, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MULTILINEAR_SPRINGS
//#####################################################################
#ifndef __MULTILINEAR_SPRINGS__
#define __MULTILINEAR_SPRINGS__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
namespace PhysBAM{

template<class TV>
class MULTILINEAR_SPRINGS:public LINEAR_SPRINGS<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef LINEAR_SPRINGS<TV> BASE;
    using BASE::visual_restlength;using BASE::restlength;using BASE::constant_youngs_modulus;using BASE::youngs_modulus;using BASE::particles;using BASE::segment_mesh;
    using BASE::damping;using BASE::force_segments;using BASE::states;using BASE::current_lengths;
    using BASE::Invalidate_CFL;
    typedef typename BASE::SEGMENT_ITERATOR SEGMENT_ITERATOR;

    ARRAY<T,VECTOR<int,1> > intervals;
    ARRAY<T,VECTOR<int,1> > youngs_modulus_scaling;
    ARRAY<T,VECTOR<int,1> > correction_force_over_youngs_modulus; // correction force to ensure contuinuity
    ARRAY<ARRAY<T> ,VECTOR<int,1> > springs_damping;
    ARRAY<int,VECTOR<int,1> > spring_count;
    ARRAY<T> base_youngs_modulus;
    T constant_base_youngs_modulus;
    ARRAY<T> correction_force;

    MULTILINEAR_SPRINGS(PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh)
        :LINEAR_SPRINGS<TV>(particles,segment_mesh,false)
    {}

    int Find_Interval(const T relative_deformation)
    {if(relative_deformation<0){
        int i=-1;for(;i>=intervals.domain.min_corner.x && relative_deformation<intervals(i);i--);
        return i+1;}
    else if(relative_deformation>0){
        int i=1;for(;i<=intervals.domain.max_corner.x && relative_deformation>intervals(i);i++);
        return i-1;}
    else return 0;}

//#####################################################################
    void Set_Spring_Phases(const ARRAY<VECTOR<T,2> >& compression_intervals_input,const ARRAY<VECTOR<T,2> >& stretching_intervals_input);
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Set_Overdamping_Fraction(const T overdamping_fraction) PHYSBAM_OVERRIDE;
    void Set_Overdamping_Fraction(ARRAY_VIEW<const T> overdamping_fraction) PHYSBAM_OVERRIDE;
    void Set_Damping(const T damping_input) PHYSBAM_OVERRIDE;
    void Set_Damping(ARRAY_VIEW<const T> damping_input) PHYSBAM_OVERRIDE;
    void Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction);
private:
    void Set_All_Springs_To_Phase(const int phase_index);
//#####################################################################
};

template<class T,class TV> MULTILINEAR_SPRINGS<TV>*
Create_Multilinear_Springs(PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh,const ARRAY<VECTOR<T,2> >& compression_intervals,
    const ARRAY<VECTOR<T,2> >& stretching_intervals,const T stiffness=2e3,const T overdamping_fraction=1,const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,
    const bool use_rest_state_for_strain_rate=true,const T restlength_enlargement_fraction=0,const bool verbose=true)
{
    MULTILINEAR_SPRINGS<TV>* ls=new MULTILINEAR_SPRINGS<TV>(particles,segment_mesh);
    ls->Set_Restlength_From_Particles();
    ls->Set_Stiffness(stiffness);
    ls->Set_Spring_Phases(compression_intervals,stretching_intervals);
    ls->Set_Overdamping_Fraction(overdamping_fraction);
    ls->Limit_Time_Step_By_Strain_Rate(limit_time_step_by_strain_rate,max_strain_per_time_step);
    ls->Use_Rest_State_For_Strain_Rate(use_rest_state_for_strain_rate);
    if(restlength_enlargement_fraction) ls->Clamp_Restlength_With_Fraction_Of_Springs(restlength_enlargement_fraction);
    if(verbose) ls->Print_Restlength_Statistics();
    return ls;
}

template<class T,class T_OBJECT> MULTILINEAR_SPRINGS<typename T_OBJECT::VECTOR_T>*
Create_Multilinear_Springs(T_OBJECT& object,const ARRAY<VECTOR<T,2> >& compression_intervals,const ARRAY<VECTOR<T,2> >& stretching_intervals,const T stiffness=2e3,
    const T overdamping_fraction=1,const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,const bool use_rest_state_for_strain_rate=true,
    const T restlength_enlargement_fraction=0,const bool verbose=true)
{
    return Create_Multilinear_Springs(object.particles,object.Get_Segment_Mesh(),compression_intervals,stretching_intervals,stiffness,overdamping_fraction,
        limit_time_step_by_strain_rate,max_strain_per_time_step,use_rest_state_for_strain_rate,restlength_enlargement_fraction,verbose);
}

}
#endif
