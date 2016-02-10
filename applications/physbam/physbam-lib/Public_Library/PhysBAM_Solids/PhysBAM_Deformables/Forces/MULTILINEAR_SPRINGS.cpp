//#####################################################################
// Copyright 2006-2008, Ronald Fedkiw, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MULTILINEAR_SPRINGS
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/MULTILINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Set_Spring_Phases
//#####################################################################
template<class TV> void MULTILINEAR_SPRINGS<TV>::
Set_Spring_Phases(const ARRAY<VECTOR<T,2> >& compression_intervals_input,const ARRAY<VECTOR<T,2> >& stretching_intervals_input)
{
    intervals.Resize(-compression_intervals_input.m,stretching_intervals_input.m,false,false);
    youngs_modulus_scaling.Resize(-compression_intervals_input.m,stretching_intervals_input.m,false,false);
    correction_force_over_youngs_modulus.Resize(-compression_intervals_input.m,stretching_intervals_input.m,false,false);
    intervals(0)=0;youngs_modulus_scaling(0)=1;correction_force_over_youngs_modulus(0)=0;
    for(int i=1;i<=stretching_intervals_input.m;i++){ // stretching phases
        intervals(i)=stretching_intervals_input(i)[1];youngs_modulus_scaling(i)=stretching_intervals_input(i)[2];
        correction_force_over_youngs_modulus(i)=correction_force_over_youngs_modulus(i-1)+intervals(i)*(youngs_modulus_scaling(i-1)-youngs_modulus_scaling(i));}
    for(int i=1;i<=compression_intervals_input.m;i++){ // compression phases
        intervals(-i)=compression_intervals_input(i)[1];youngs_modulus_scaling(-i)=compression_intervals_input(i)[2];
        correction_force_over_youngs_modulus(-i)=correction_force_over_youngs_modulus(1-i)+intervals(-i)*(youngs_modulus_scaling(1-i)-youngs_modulus_scaling(-i));}

    if(constant_youngs_modulus){
        constant_base_youngs_modulus=constant_youngs_modulus;constant_youngs_modulus=0;
        youngs_modulus.Resize(segment_mesh.elements.m);}
    else base_youngs_modulus=youngs_modulus;
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void MULTILINEAR_SPRINGS<TV>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    states.Resize(segment_mesh.elements.m,false,false);current_lengths.Resize(segment_mesh.elements.m,false,false);
    correction_force.Resize(segment_mesh.elements.m,false,false);
    spring_count.Resize(intervals.domain.min_corner.x,intervals.domain.max_corner.x,false,false);
    ARRAY_VIEW<const TV> X(particles.X);
    Invalidate_CFL();
    spring_count.Fill(0);
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        typename BASE::STATE& state=states(s);
        const VECTOR<int,2>& nodes=segment_mesh.elements(s);
        state.nodes=nodes;
        state.direction=X(nodes[2])-X(nodes[1]);
        current_lengths(s)=state.direction.Normalize();
        T relative_deformation=(current_lengths(s)-visual_restlength(s))/restlength(s);
        int index=Find_Interval(relative_deformation);
        spring_count(index)++;
        T ym;if(constant_base_youngs_modulus) ym=constant_base_youngs_modulus;else ym=base_youngs_modulus(s);
        youngs_modulus(s)=youngs_modulus_scaling(index)*ym;
        correction_force(s)=correction_force_over_youngs_modulus(index)*youngs_modulus(s);
        damping(s)=springs_damping(index)(s);
        state.coefficient=damping(s)/restlength(s);}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void MULTILINEAR_SPRINGS<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    LINEAR_SPRINGS<TV>::Add_Velocity_Independent_Forces(F,time);
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const VECTOR<int,2>& nodes=segment_mesh.elements(s);
        TV force=correction_force(s)*states(s).direction;
        F(nodes[1])+=force;F(nodes[2])-=force;}
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class TV> void MULTILINEAR_SPRINGS<TV>::
Set_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    springs_damping.Resize(intervals.domain.min_corner.x,intervals.domain.max_corner.x,false,false);
    for(int i=springs_damping.domain.min_corner.x;i<=springs_damping.domain.max_corner.x;i++){
        Set_All_Springs_To_Phase(i);
        LINEAR_SPRINGS<TV>::Set_Overdamping_Fraction(overdamping_fraction);
        springs_damping(i)=damping;}
}
//#####################################################################
// Function Set_Damping
//#####################################################################
template<class TV> void MULTILINEAR_SPRINGS<TV>::
Set_Damping(const T constant_damping_input)
{
    springs_damping.Resize(intervals.domain.min_corner.x,intervals.domain.max_corner.x,false,false);
    for(int i=springs_damping.domain.min_corner.x;i<=springs_damping.domain.max_corner.x;i++){
        springs_damping(i).Resize(segment_mesh.elements.m,false,false);
        ARRAYS_COMPUTATIONS::Fill(springs_damping(i),constant_damping_input);}
}
//#####################################################################
// Function Set_Damping
//#####################################################################
template<class TV> void MULTILINEAR_SPRINGS<TV>::
Set_Damping(ARRAY_VIEW<const T> damping_input)
{
    springs_damping.Resize(intervals.domain.min_corner.x,intervals.domain.max_corner.x,false,false);
    for(int i=springs_damping.domain.min_corner.x;i<=springs_damping.domain.max_corner.x;i++) springs_damping(i)=damping_input;
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class TV> void MULTILINEAR_SPRINGS<TV>::
Set_Overdamping_Fraction(ARRAY_VIEW<const T> overdamping_fraction) // 1 is critically damped
{
    springs_damping.Resize(intervals.domain.min_corner.x,intervals.domain.max_corner.x,false,false);
    for(int i=springs_damping.domain.min_corner.x;i<=springs_damping.domain.max_corner.x;i++){
        Set_All_Springs_To_Phase(i);
        LINEAR_SPRINGS<TV>::Set_Overdamping_Fraction(overdamping_fraction);
        springs_damping(i)=damping;}
}
//#####################################################################
// Function Ensure_Minimum_Overdamping_Fraction
//#####################################################################
template<class TV> void MULTILINEAR_SPRINGS<TV>::
Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    for(int i=springs_damping.domain.min_corner.x;i<=springs_damping.domain.max_corner.x;i++){
        Set_All_Springs_To_Phase(i);
        LINEAR_SPRINGS<TV>::Set_Overdamping_Fraction(overdamping_fraction);
        for(int k=1;k<=damping.m;k++) springs_damping(i)(k)=max(springs_damping(i)(k),damping(k));}
}
//#####################################################################
// Function Set_All_Springs_To_Phase
//#####################################################################
template<class TV> void MULTILINEAR_SPRINGS<TV>::
Set_All_Springs_To_Phase(const int phase_index) // 1 is critically damped
{
    if(constant_base_youngs_modulus) for(int s=1;s<=segment_mesh.elements.m;s++) youngs_modulus(s)=youngs_modulus_scaling(phase_index)*constant_base_youngs_modulus;
    else for(int s=1;s<=segment_mesh.elements.m;s++) youngs_modulus(s)=youngs_modulus_scaling(phase_index)*base_youngs_modulus(s);
}
//#####################################################################
template class MULTILINEAR_SPRINGS<VECTOR<float,2> >;
template class MULTILINEAR_SPRINGS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MULTILINEAR_SPRINGS<VECTOR<double,2> >;
template class MULTILINEAR_SPRINGS<VECTOR<double,3> >;
#endif
