//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Andrew Selle, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using ::std::sqrt;
using namespace PhysBAM;
//#####################################################################
// Function Set_Stiffness_Based_On_Reduced_Mass
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_2D<T>::
Set_Stiffness_Based_On_Reduced_Mass(const T scaling_coefficient) // assumes mass and restlength are already defined
{
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        for(int s=1;s<=3;s++){int node1,node2,node3;
            switch(s){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}
            T one_over_triangle_mass=(T).25*(particles.one_over_effective_mass(node2)+particles.one_over_effective_mass(node3)),one_over_particle_mass=particles.one_over_effective_mass(node1),
               harmonic_mass=Pseudo_Inverse(one_over_particle_mass+one_over_triangle_mass);
            parameters(t)(s).youngs_modulus=scaling_coefficient*harmonic_mass/parameters(t)(s).restlength;}}
}
//#####################################################################
// Function Set_Restlength_From_Particles
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_2D<T>::
Set_Restlength_From_Particles()
{
    Set_Restlength_From_Material_Coordinates(particles.X);
}
//#####################################################################
// Function Set_Restlength_From_Material_Coordinates
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_2D<T>::
Set_Restlength_From_Material_Coordinates(ARRAY_VIEW<TV> material_coordinates)
{
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        parameters(t)(1).restlength=SEGMENT_2D<T>(material_coordinates(j),material_coordinates(k)).Distance_From_Point_To_Line(material_coordinates(i));
        parameters(t)(2).restlength=SEGMENT_2D<T>(material_coordinates(i),material_coordinates(k)).Distance_From_Point_To_Line(material_coordinates(j));
        parameters(t)(3).restlength=SEGMENT_2D<T>(material_coordinates(i),material_coordinates(j)).Distance_From_Point_To_Line(material_coordinates(k));
        for(int k=1;k<=3;k++) parameters(t)(k).visual_restlength=parameters(t)(k).restlength;}
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_2D<T>::
Set_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        T harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(i)+(T).25*(particles.one_over_effective_mass(j)+particles.one_over_effective_mass(k)));
        parameters(t)(1).damping=overdamping_fraction*2*sqrt(parameters(t)(1).youngs_modulus*parameters(t)(1).restlength*harmonic_mass);
        harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(j)+(T).25*(particles.one_over_effective_mass(i)+particles.one_over_effective_mass(k)));
        parameters(t)(2).damping=overdamping_fraction*2*sqrt(parameters(t)(2).youngs_modulus*parameters(t)(2).restlength*harmonic_mass);
        harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(k)+(T).25*(particles.one_over_effective_mass(j)+particles.one_over_effective_mass(i)));
        parameters(t)(3).damping=overdamping_fraction*2*sqrt(parameters(t)(3).youngs_modulus*parameters(t)(3).restlength*harmonic_mass);}
    Invalidate_CFL();
}
//#####################################################################
// Function Ensure_Minimum_Overdamping_Fraction
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_2D<T>::
Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    ARRAY<VECTOR<T,3> > save_damping(parameters.m);for(int i=1;i<=parameters.m;i++) for(int k=1;k<=3;k++) save_damping(i)(k)=parameters(i)(k).damping;
    Set_Overdamping_Fraction(overdamping_fraction);
    for(int k=1;k<=parameters.m;k++) for(int i=1;i<=3;i++) parameters(k)(i).damping=max(parameters(k)(i).damping,save_damping(k)(i));
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_2D<T>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    ARRAY_VIEW<const TV> X(particles.X);int node1,node2,node3; // node1 is the isolated vertex and nodes2,3 are the segment
    int used_springs=0,total_elements=0;
    spring_states.Resize(mesh.elements.m);
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int i,j,k;mesh.elements(t).Get(i,j,k);
        int hmin=0;T cross_length_max=-FLT_MAX;
        for(int h=1;h<=3;h++){
            switch(h){case 1:node2=j;node3=k;break;case 2:node2=k;node3=i;break;default:node2=i;node3=j;}
            T cross_length=(X(node3)-X(node2)).Magnitude_Squared();if(cross_length>cross_length_max){hmin=h;cross_length_max=cross_length;}}
        switch(hmin){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}
        TV direction=SEGMENT_2D<T>::Normal(X(node3),X(node2));
        T rl=parameters(t)(hmin).restlength,current_length=TV::Dot_Product(X(node1)-X(node2),direction);
        SPRING_STATE& spring_state=spring_states(t);
        total_elements++;
        if(use_springs_compressed_beyond_threshold && current_length-rl>-spring_compression_fraction_threshold*rl) spring_state.node=0;
        else{
            used_springs++;
            spring_state.current_length=current_length;spring_state.coefficient=parameters(t)(hmin).damping/rl;
            spring_state.barycentric=SEGMENT_2D<T>::Clamped_Barycentric_Coordinates(X(node1),X(node2),X(node3));
            if(use_plasticity) Compute_Plasticity(hmin,t,current_length);            
            spring_state.direction=direction;spring_state.node=hmin;}}
    if(print_number_used) {std::stringstream ss;ss<<"using "<<used_springs<<" of "<<total_elements<<" altitude springs"<<std::endl;LOG::filecout(ss.str());}
    if(compute_half_forces) for(int i=1;i<=spring_states.m;i++) spring_states(i).sqrt_coefficient=sqrt(spring_states(i).coefficient);
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_2D<T>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        const SPRING_STATE& state=spring_states(t);
        if(state.node){
            int i,j,k;mesh.elements(t).Get(i,j,k);
            int node1,node2,node3; // node1 is the isolated vertex and nodes2,3 is the segment
            const SPRING_PARAMETER& parameter=parameters(t)(state.node);
            switch(state.node){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}
            T rl=parameter.restlength,vrl=parameter.visual_restlength;
            TV force=parameter.youngs_modulus*(state.current_length-vrl)/rl*state.direction;
            F(node1)-=force;F(node2)+=state.barycentric.x*force;F(node3)+=state.barycentric.y*force;}}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_2D<T>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        const SPRING_STATE& spring_state=spring_states(t);
        if(spring_state.node){
            int i,j,k;mesh.elements(t).Get(i,j,k);
            int node1,node2,node3; // node1 is the isolated vertex and nodes2,3 are the segment
            switch(spring_state.node){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}
            TV force=(spring_state.coefficient*TV::Dot_Product(V(node1)-spring_state.barycentric.x*V(node2)-spring_state.barycentric.y*V(node3),spring_state.direction))*spring_state.direction;
            F(node1)-=force;F(node2)+=spring_state.barycentric.x*force;F(node3)+=spring_state.barycentric.y*force;}}
}
//#####################################################################
// Function Velocity_Dependent_Size
//#####################################################################
template<class T> int LINEAR_ALTITUDE_SPRINGS_2D<T>::
Velocity_Dependent_Forces_Size() const
{
    int aggregate_id=0;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){
        int t=iterator.Data();
        const SPRING_STATE& spring_state=spring_states(t);
        if(spring_state.node) aggregate_id++;}
    return aggregate_id;
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_2D<T>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    int aggregate_id=1;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        const SPRING_STATE& spring_state=spring_states(t);
        if(spring_state.node){
            int i,j,k;mesh.elements(t).Get(i,j,k);
            int node1,node2,node3; // node1 is the isolated vertex and nodes2,3 are the segment
            switch(spring_state.node){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}
            aggregate(aggregate_id++)+=spring_state.sqrt_coefficient*TV::Dot_Product(V(node1)-spring_state.barycentric.x*V(node2)-spring_state.barycentric.y*V(node3),spring_state.direction);}}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_2D<T>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    int aggregate_id=1;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        const SPRING_STATE& spring_state=spring_states(t);
        if(spring_state.node){
            int i,j,k;mesh.elements(t).Get(i,j,k);
            int node1,node2,node3; // node1 is the isolated vertex and nodes2,3 are the segment
            switch(spring_state.node){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}
            TV force=spring_state.sqrt_coefficient*aggregate(aggregate_id++)*spring_state.direction;
            F(node1)+=force;F(node2)-=spring_state.barycentric.x*force;F(node3)-=spring_state.barycentric.y*force;}}
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class T> T LINEAR_ALTITUDE_SPRINGS_2D<T>::
CFL_Strain_Rate() const
{
    T max_strain_rate=0,dx;TV direction,barycentric,v_interpolated;int node1,node2,node3; // 1 is vertex, 2,3 is segment
    ARRAY<VECTOR<int,3> >& elements=mesh.elements;ARRAY_VIEW<const TV> X(particles.X),V(particles.V);
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int i,j,k;elements(t).Get(i,j,k);
        int hmin=0;T cross_length_max=-FLT_MAX;
        for(int h=1;h<=3;h++){
            switch(h){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}
            T cross_length=(X(node3)-X(node2)).Magnitude_Squared();
            if(cross_length>cross_length_max){hmin=h;cross_length_max=cross_length;}}
        switch(hmin){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}
        direction=SEGMENT_2D<T>::Normal(X(node2),X(node3));
        T rl=parameters(t)(hmin).restlength;
        if(use_rest_state_for_strain_rate) dx=rl;else dx=TV::Dot_Product(direction,X(node1)-X(node2));
        if(use_springs_compressed_beyond_threshold){
            T length;if(use_rest_state_for_strain_rate) length=TV::Dot_Product(direction,X(node1)-X(node2));else length=dx;
            if(length-rl>-spring_compression_fraction_threshold*rl) continue;}
        barycentric=SEGMENT_2D<T>::Clamped_Barycentric_Coordinates(X(node1),X(node2),X(node3));
        v_interpolated=barycentric.x*V(node2)+barycentric.y*V(node3);
        max_strain_rate=max(max_strain_rate,abs(TV::Dot_Product(V(node1)-v_interpolated,direction)/dx));}
    return Robust_Divide(max_strain_per_time_step,max_strain_rate);
}
//#####################################################################
template class LINEAR_ALTITUDE_SPRINGS_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LINEAR_ALTITUDE_SPRINGS_2D<double>;
#endif
