//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_S3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using ::std::sqrt;
using namespace PhysBAM;
//#####################################################################
// Function Set_Stiffness_Based_On_Reduced_Mass
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
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
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Set_Restlength_From_Particles()
{
    Set_Restlength_From_Material_Coordinates(particles.X);
}
//#####################################################################
// Function Set_Restlength_From_Material_Coordinates
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Set_Restlength_From_Material_Coordinates(ARRAY_VIEW<const TV> X)
{
    Invalidate_CFL();
    for(int t=1;t<=mesh.elements.m;t++){
        int i=mesh.elements(t)(1),j=mesh.elements(t)(2),k=mesh.elements(t)(3);
        parameters(t)(1).restlength=SEGMENT_3D<T>(X(j),X(k)).Distance_From_Point_To_Line(X(i));
        parameters(t)(2).restlength=SEGMENT_3D<T>(X(i),X(k)).Distance_From_Point_To_Line(X(j));
        parameters(t)(3).restlength=SEGMENT_3D<T>(X(i),X(j)).Distance_From_Point_To_Line(X(k));
        for(int k=1;k<=3;k++) parameters(t)(k).visual_restlength=parameters(t)(k).restlength;}
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
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
// Function Set_Overdamping_Fraction
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Set_Overdamping_Fraction(const ARRAY<VECTOR<T,3> >& overdamping_fraction) // 1 is critically damped
{
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        T harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(i)+(T).25*(particles.one_over_effective_mass(j)+particles.one_over_effective_mass(k)));
        parameters(t)(1).damping=overdamping_fraction(t)(1)*2*sqrt(parameters(t)(1).youngs_modulus*parameters(t)(1).restlength*harmonic_mass);
        harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(j)+(T).25*(particles.one_over_effective_mass(i)+particles.one_over_effective_mass(k)));
        parameters(t)(2).damping=overdamping_fraction(t)(2)*2*sqrt(parameters(t)(2).youngs_modulus*parameters(t)(2).restlength*harmonic_mass);
        harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(k)+(T).25*(particles.one_over_effective_mass(j)+particles.one_over_effective_mass(i)));
        parameters(t)(3).damping=overdamping_fraction(t)(3)*2*sqrt(parameters(t)(3).youngs_modulus*parameters(t)(3).restlength*harmonic_mass);}
    Invalidate_CFL();
}
//#####################################################################
// Function Ensure_Minimum_Overdamping_Fraction
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    ARRAY<VECTOR<T,3> > save_damping(parameters.m);for(int i=1;i<=parameters.m;i++) for(int k=1;k<=3;k++) save_damping(i)(k)=parameters(i)(k).damping;
    Set_Overdamping_Fraction(overdamping_fraction);
    for(int k=1;k<=parameters.m;k++) for(int i=1;i<=3;i++) parameters(k)(i).damping=max(parameters(k)(i).damping,save_damping(k)(i));
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    ARRAY_VIEW<const TV> X(particles.X);int node1,node2,node3; // node1 is the isolated vertex and nodes2,3 are the segment
    int used_springs=0,total_elements=0;
    spring_states.Resize(mesh.elements.m);
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data(); // use shortest spring only
        int i,j,k;mesh.elements(t).Get(i,j,k);
        int hmin=0;T cross_length_max=-FLT_MAX;
        if(is_position_update){
            for(int h=1;h<=3;h++){
                switch(h){case 1:node2=j;node3=k;break;case 2:node2=k;node3=i;break;default:node2=i;node3=j;}
                T cross_length=(X(node3)-X(node2)).Magnitude_Squared();if(cross_length>cross_length_max){hmin=h;cross_length_max=cross_length;}}}
        else hmin=spring_states(t).node;
        switch(hmin){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}
        TV direction=SEGMENT_3D<T>::Normal(X(node2),X(node3),X(node1));
        if(triangle_inverted && (*triangle_inverted)(t)) direction=-direction;
        T rl=parameters(t)(hmin).restlength,current_length=TV::Dot_Product(X(node1)-X(node2),direction);
        SPRING_STATE& state=spring_states(t);
        total_elements++;
        if(use_springs_compressed_beyond_threshold && current_length-rl>-spring_compression_fraction_threshold*rl) state.node=0;
        else{
            used_springs++;
            state.current_length=current_length;state.coefficient=parameters(t)(hmin).damping/rl;
            state.barycentric=SEGMENT_3D<T>::Clamped_Barycentric_Coordinates(X(node1),X(node2),X(node3));
            if(use_plasticity) Compute_Plasticity(hmin,t,current_length);
            state.direction=direction;state.node=hmin;}}

    residual_PE.Resize(mesh.elements.m);
    if(!total_PE.m){
        total_PE.Resize(mesh.elements.m);
        for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            total_PE(t)=Potential_Energy(t,time);}}

    if(print_number_used) {std::stringstream ss;ss<<"using "<<used_springs<<" of "<<total_elements<<" altitude springs"<<std::endl;LOG::filecout(ss.str());}
    if(!mesh.incident_elements) mesh.Initialize_Incident_Elements();
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        const SPRING_STATE& state=spring_states(t);
        if(state.node){
            int i,j,k;mesh.elements(t).Get(i,j,k);
            int node1,node2,node3; // node1 is the isolated vertex and nodes2,3 is the segment
            switch(state.node){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}
            const SPRING_PARAMETER& parameter=parameters(t)(state.node);
            T rl=parameter.restlength,vrl=parameter.visual_restlength;
            TV force=parameter.youngs_modulus/rl*(state.current_length-vrl)*state.direction;
            F(node1)-=force;F(node2)+=state.barycentric.x*force;F(node3)+=state.barycentric.y*force;}}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        const SPRING_STATE& state=spring_states(t);
        if(state.node){
            int i,j,k;mesh.elements(t).Get(i,j,k);
            int node1,node2,node3; // node1 is the isolated vertex and nodes2,3 are the segment
            switch(state.node){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}
            TV force=(state.coefficient*TV::Dot_Product(V(node1)-state.barycentric.x*V(node2)-state.barycentric.y*V(node3),state.direction))*state.direction;
            F(node1)-=force;F(node2)+=state.barycentric.x*force;F(node3)+=state.barycentric.y*force;}}
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class T> T LINEAR_ALTITUDE_SPRINGS_S3D<T>::
CFL_Strain_Rate() const
{
    T max_strain_rate=0,dx;VECTOR<T,3> direction,v_interpolated; VECTOR<T,2> barycentric;int node1,node2,node3; // 1 is vertex, 2,3 is segment
    ARRAY<VECTOR<int,3> >& elements=mesh.elements;ARRAY_VIEW<const TV> X(particles.X),V(particles.V);
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data(); // use shortest spring only
        int i,j,k;elements(t).Get(i,j,k);
        int hmin=0;T cross_length_max=-FLT_MAX;
        for(int h=1;h<=3;h++){
            switch(h){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}
            T cross_length=(X(node3)-X(node2)).Magnitude_Squared();
            if(cross_length>cross_length_max){hmin=h;cross_length_max=cross_length;}}
        switch(hmin){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}
        direction=SEGMENT_3D<T>::Normal(X(node2),X(node3),X(node1));
        T rl=parameters(t)(hmin).restlength;
        if(use_rest_state_for_strain_rate) dx=rl;else dx=VECTOR<T,3>::Dot_Product(direction,X(node1)-X(node2));
        if(use_springs_compressed_beyond_threshold){
            T length;if(use_rest_state_for_strain_rate) length=VECTOR<T,3>::Dot_Product(direction,X(node1)-X(node2));else length=dx;
            if(length-rl>-spring_compression_fraction_threshold*rl) continue;}
        barycentric=SEGMENT_3D<T>::Clamped_Barycentric_Coordinates(X(node1),X(node2),X(node3));
        v_interpolated=barycentric.x*V(node2)+barycentric.y*V(node3);
        T strain_rate=VECTOR<T,3>::Dot_Product(V(node1)-v_interpolated,direction)/dx;max_strain_rate=max(max_strain_rate,abs(strain_rate));
        if(cache_strain){strains_of_spring(t)=VECTOR<T,2>(abs(strain_rate),abs((dx-parameters(t)(hmin).visual_restlength)/rl));}}
   return Robust_Divide(max_strain_per_time_step,max_strain_rate);
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class T> T LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Potential_Energy(const int t,const T time) const
{
    const SPRING_STATE& state=spring_states(t);
    if(state.node){
        int i,j,k;mesh.elements(t).Get(i,j,k);
        int node1,node2,node3; // node1 is the isolated vertex and nodes2,3 is the segment
        switch(state.node){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}
        const SPRING_PARAMETER& parameter=parameters(t)(state.node);
        T rl=parameter.restlength,vrl=parameter.visual_restlength;
        TV direction=SEGMENT_3D<T>::Normal(particles.X(node2),particles.X(node3),particles.X(node1));
        T current_length=TV::Dot_Product(particles.X(node1)-particles.X(node2),direction);
        return (T).5*parameter.youngs_modulus/rl*sqr(current_length-vrl);}
    else return 0;
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class T> T LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        potential_energy+=Potential_Energy(t,time);}
    return potential_energy;
}
//#####################################################################
// Function Residual_Energy
//#####################################################################
template<class T> T LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Residual_Energy(const T time) const
{
    T residual_energy=0;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        residual_energy+=residual_PE(t);}
    return residual_energy;
}
//#####################################################################
// Function Save_Potential_Energy
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Save_Potential_Energy(const T time)
{
    potential_energy_save.Resize(mesh.elements.m);
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        potential_energy_save(t)=Potential_Energy(t,time);}
}
//#####################################################################
// Function Setup_Set_Velocity_From_Positions
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Setup_Set_Velocity_From_Positions(const T time,const bool is_position_update,const bool reset_alphas)
{
    if(!is_position_update){
        // UPBS isn't called, so do most of the work, except keep spring direction the same
        ARRAY_VIEW<const TV> X(particles.X);int node1,node2,node3; // node1 is the isolated vertex and nodes2,3 are the segment
        spring_states.Resize(mesh.elements.m);
        for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data(); // use shortest spring only
            int i,j,k;mesh.elements(t).Get(i,j,k);
            SPRING_STATE& state=spring_states(t);
            int hmin=state.node;
            switch(hmin){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}
            TV direction=SEGMENT_3D<T>::Normal(X(node2),X(node3),X(node1));
            if(triangle_inverted && (*triangle_inverted)(t)) direction=-direction;
            T rl=parameters(t)(hmin).restlength,current_length=TV::Dot_Product(X(node1)-X(node2),direction);
            if(use_springs_compressed_beyond_threshold && current_length-rl>-spring_compression_fraction_threshold*rl) state.node=0;
            else{
                state.current_length=current_length;state.coefficient=parameters(t)(hmin).damping/rl;
                state.barycentric=SEGMENT_3D<T>::Clamped_Barycentric_Coordinates(X(node1),X(node2),X(node3));
                if(use_plasticity) Compute_Plasticity(hmin,t,current_length);
                state.direction=direction;state.node=hmin;}}}

    force_estimates.Resize(mesh.elements.m);
    incident_nodes.Resize(mesh.elements.m);
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int i,j,k;mesh.elements(t).Get(i,j,k);
        incident_nodes(t).Remove_All();
        incident_nodes(t).Append(i);incident_nodes(t).Append(j);incident_nodes(t).Append(k);}
    if(is_position_update){
        for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            if(reset_alphas) force_estimates(t)=0;
            residual_PE(t)=Potential_Energy(t,time)-total_PE(t);}
        Save_Potential_Energy(time);}
}
//#####################################################################
// Function Get_Element_Count
//#####################################################################
template<class T> int LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Get_Element_Count()
{
    return mesh.elements.m;
}
//#####################################################################
// Function Get_Force_Elements
//#####################################################################
template<class T> FORCE_ELEMENTS* LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Get_Force_Elements()
{
    return &force_elements;
}
//#####################################################################
// Function Incident_Nodes
//#####################################################################
template<class T> ARRAY<int>* LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Incident_Nodes(const int force_element)
{
    return &(incident_nodes(force_element));
}
//#####################################################################
// Function Get_Direction
//#####################################################################
template<class T> VECTOR<T,3> LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Get_Direction(const int force_element)
{
    const SPRING_STATE& state=spring_states(force_element);
    return state.direction;
}
//#####################################################################
// Function Get_Combined_One_Over_Mass
//#####################################################################
template<class T> T LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Get_Combined_One_Over_Mass(const int force_element)
{
    const SPRING_STATE& state=spring_states(force_element);
    if(state.node){
        int i,j,k;mesh.elements(force_element).Get(i,j,k);
        int node1,node2,node3; // node1 is the isolated vertex and nodes2,3 are the segment
        switch(state.node){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}
        return particles.one_over_mass(node1)+(sqr(state.barycentric.x)*particles.one_over_mass(node2)+sqr(state.barycentric.y)*particles.one_over_mass(node3));}
    return 0;
}
//#####################################################################
// Function Incident_Force_Elements
//#####################################################################
template<class T> ARRAY<int>* LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Incident_Force_Elements(const int particle)
{
    return &(*mesh.incident_elements)(particle);
}
//#####################################################################
// Function Choose_Solution
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Choose_Solution(const bool use_orig_force,const int force_element,const T dt,const T alpha1,const T alpha2,ARRAY<T>& v_n_hats)
{
    if(use_orig_force){
        const SPRING_STATE& state=spring_states(force_element);
        const SPRING_PARAMETER& parameter=parameters(force_element)(state.node);
        T rl=parameter.restlength,vrl=parameter.visual_restlength;
        T force=parameter.youngs_modulus/rl*(state.current_length-vrl);

        if(fabs(force-alpha1) < fabs(force-alpha2)) Set_Force(force_element,alpha1);
        else Set_Force(force_element,alpha2);}
    else{
        const SPRING_STATE& state=spring_states(force_element);
        if(state.node){
            int i,j,k;mesh.elements(force_element).Get(i,j,k);
            int node1,node2,node3; // node1 is the isolated vertex and nodes2,3 are the segment
            T v1_n_hat;
            T v2_n_hat;
            switch(state.node){
                case 1:{
                    node1=i;node2=j;node3=k;
                    v1_n_hat=state.barycentric.x*v_n_hats(2)+state.barycentric.y*v_n_hats(3);
                    v2_n_hat=v_n_hats(1);}
                    break;
                case 2:{
                    node1=j;node2=k;node3=i;
                    v1_n_hat=state.barycentric.x*v_n_hats(3)+state.barycentric.y*v_n_hats(1);
                    v2_n_hat=v_n_hats(2);}
                    break;
                default:{
                    node1=k;node2=i;node3=j;
                    v1_n_hat=state.barycentric.x*v_n_hats(1)+state.barycentric.y*v_n_hats(2);
                    v2_n_hat=v_n_hats(3);}}
            
            T v1_one_over_mass,v1_mass;
            if(particles.mass(node2)==FLT_MAX && particles.mass(node3)==FLT_MAX){
                v1_one_over_mass=0;
                v1_mass=FLT_MAX;}
            else{
                v1_one_over_mass=sqr(state.barycentric.x)*particles.one_over_mass(node2)+sqr(state.barycentric.y)*particles.one_over_mass(node3);
                v1_mass = (particles.mass(node2)*particles.mass(node3))/(sqr(state.barycentric.x)*particles.mass(node3)+sqr(state.barycentric.y)*particles.mass(node2));}
            T v1_alpha1_temp=v1_n_hat+dt*v1_one_over_mass*alpha1;
            T v2_alpha1_temp=v2_n_hat-dt*particles.one_over_mass(node1)*alpha1;
            
            T v1_alpha2_temp=v1_n_hat+dt*v1_one_over_mass*alpha2;
            T v2_alpha2_temp=v2_n_hat-dt*particles.one_over_mass(node1)*alpha2;
            
            T combined_mass=particles.mass(node1)+v1_mass;
            TV x1=state.barycentric.x*particles.X(node2)+state.barycentric.y*particles.X(node3);
            T x1u=TV::Dot_Product(x1,state.direction);
            T x2u=TV::Dot_Product(particles.X(node1),state.direction);
            
            const SPRING_PARAMETER& parameter=parameters(force_element)(state.node);
            T rl=parameter.restlength;
            T w=sqrt((parameter.youngs_modulus/rl)*Get_Combined_One_Over_Mass(force_element));
            T v_cm_0=(v1_mass*v1_n_hat+particles.mass(node1)*v2_n_hat)/combined_mass;
            T c1=(x2u-x1u)-rl;
            T c2=(v2_n_hat-v1_n_hat)/w;
            T v1_analytic=v_cm_0-(particles.mass(node1)/combined_mass)*(-w*c1*sin(w*dt)+w*c2*cos(w*dt));
            T v2_analytic=v_cm_0+(v1_mass/combined_mass)*(-w*c1*sin(w*dt)+w*c2*cos(w*dt));
            if((fabs(v1_alpha1_temp-v1_analytic)+fabs(v2_alpha1_temp-v2_analytic)) < (fabs(v1_alpha2_temp-v1_analytic)+fabs(v2_alpha2_temp-v2_analytic))) Set_Force(force_element,alpha1);
            else Set_Force(force_element,alpha2);}
        else Set_Force(force_element,0);}
}
//#####################################################################
// Function Get_Force
//#####################################################################
template<class T> VECTOR<T,3> LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Get_Force(const int force_element,const int particle,const bool use_original_force)
{
    const SPRING_STATE& state=spring_states(force_element);
    if(state.node){
        int i,j,k;mesh.elements(force_element).Get(i,j,k);
        int node1,node2,node3; // node1 is the isolated vertex and nodes2,3 are the segment
        switch(state.node){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}
        TV force;
        if(use_original_force){
            const SPRING_PARAMETER& parameter=parameters(force_element)(state.node);
            T rl=parameter.restlength,vrl=parameter.visual_restlength;
            force=parameter.youngs_modulus/rl*(state.current_length-vrl)*state.direction;}
        else force=force_estimates(force_element)*state.direction;
        if(particle==node1) return -force;
        else if(particle==node2) return state.barycentric.x*force;
        else if(particle==node3) return state.barycentric.y*force;}
    return TV();
}
//#####################################################################
// Function Get_Force
//#####################################################################
template<class T> T LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Get_Force(const int force_element)
{
    return force_estimates(force_element);
}
//#####################################################################
// Function Set_Force
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Set_Force(const int force_element,const T force)
{
    force_estimates(force_element)=force;
}
//#####################################################################
// Function Get_Damping_Force
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Get_Damping_Force(const int particle,TV& damping_force,const T dt,const bool use_coefficient)
{
    ARRAY<int>& node_incident_forces=(*mesh.incident_elements)(particle);
    for(int t=1;t<=node_incident_forces.m;t++){
        const SPRING_STATE& state=spring_states(node_incident_forces(t));
        if(state.node){
            int i,j,k;mesh.elements(node_incident_forces(t)).Get(i,j,k);
            int node1,node2,node3; // node1 is the isolated vertex and nodes2,3 are the segment
            switch(state.node){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}
            TV force=(TV::Dot_Product(particles.V(node1)-state.barycentric.x*particles.V(node2)-state.barycentric.y*particles.V(node3),state.direction))*state.direction;
            if(particle==node1) force*=(T)-1;
            else if(particle==node2) force*=state.barycentric.x;
            else force*=state.barycentric.y;
            if(use_coefficient) damping_force+=dt*particles.one_over_mass(particle)*state.coefficient*force;
            else damping_force+=dt*particles.one_over_mass(particle)*force;}}
}
//#####################################################################
// Function Update_Residual_Energy
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Update_Residual_Energy(const int force_element,const T residual_energy,const T time)
{
    residual_PE(force_element)=residual_energy;
    total_PE(force_element)=Potential_Energy(force_element,time)-residual_PE(force_element);
}
//#####################################################################
// Function Update_Residual_Energy
//#####################################################################
template<class T> T LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Get_Residual_Energy(const int force_element)
{
    return residual_PE(force_element);
}
//#####################################################################
// Function Compute_Quadratic_Contribution_For_Node
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Compute_Quadratic_Contribution_For_Force(T& A,T& a,T&c,const T dt,const int force_element,const T combined_one_over_mass,const bool ignore_PE_terms)
{
    const SPRING_STATE& state=spring_states(force_element);
    if(state.node){
        if(!ignore_PE_terms){
            int i,j,k;mesh.elements(force_element).Get(i,j,k);
            int node1,node2,node3; // node1 is the isolated vertex and nodes2,3 is the segment
            switch(state.node){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}

            T x1u=TV::Dot_Product(particles.X(node1),state.direction);
            T x2u=TV::Dot_Product(particles.X(node2),state.direction);
            T x3u=TV::Dot_Product(particles.X(node3),state.direction);
            
            a=x1u-(state.barycentric.x*x2u+state.barycentric.y*x3u);
            c=-(T).5*dt*dt*combined_one_over_mass;
            
            const SPRING_PARAMETER& parameter=parameters(force_element)(state.node);
            A+=(T).5*(parameter.youngs_modulus/parameter.restlength)*c*c;}
        
        A+=(T).5*dt*dt*combined_one_over_mass;}
}
//#####################################################################
// Function Compute_Quadratic_Contribution_For_Node
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Compute_Quadratic_Contribution_For_Node(T& B,T& C,T&b,const T dt,const int node,const int force_element,const T combined_one_over_mass,const T v_n_correction,const bool ignore_PE_terms)
{
    const SPRING_STATE& state=spring_states(force_element);
    if(state.node){
        int i,j,k;mesh.elements(force_element).Get(i,j,k);
        int node1,node2,node3; // node1 is the isolated vertex and nodes2,3 is the segment
        switch(state.node){case 1:node1=i;node2=j;node3=k;break;case 2:node1=j;node2=k;node3=i;break;default:node1=k;node2=i;node3=j;}

        T vu=TV::Dot_Product(particles.V(node),state.direction);
        
        if(node==node1){
            T term=dt*(vu+(T).5*v_n_correction);
            if(!ignore_PE_terms) b+=term;
            B+=-term;}
        else if(node==node2){
            T term=dt*state.barycentric.x*(vu+(T).5*v_n_correction);
            if(!ignore_PE_terms) b+=-term;
            B+=term;}
        else if(node==node3){
            T term=dt*state.barycentric.y*(vu+(T).5*v_n_correction);
            if(!ignore_PE_terms) b+=-term;
            B+=term;}}
}
//#####################################################################
// Function Compute_Quadratic_Contribution_For_Residual
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Compute_Quadratic_Contribution_For_Residual(T& B,T& C,T& a,T&b,T&c,const T dt,const T time,const int force_element,const bool ignore_PE_terms)
{
    const SPRING_STATE& state=spring_states(force_element);
    if(state.node){
        if(!ignore_PE_terms){
            const SPRING_PARAMETER& parameter=parameters(force_element)(state.node);
            B+=(parameter.youngs_modulus/parameter.restlength)*(a+b-parameter.visual_restlength)*c;
            C=(parameter.youngs_modulus/parameter.restlength)*(a+b/2-parameter.visual_restlength)*b;}
        else C=delta_PE(force_element);
        
        C+=residual_PE(force_element);}
}
//#####################################################################
// Function Store_Delta_PE
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Store_Delta_PE(const T time)
{
    delta_PE.Resize(potential_energy_save.m);
    T total_current_PE=0;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        delta_PE(t)=Potential_Energy(t,time)-potential_energy_save(t);
        total_current_PE+=Potential_Energy(t,time);}
    total_delta_PE=total_current_PE-ARRAYS_COMPUTATIONS::Sum(potential_energy_save);
}
//#####################################################################
// Function Get_Total_Delta_PE
//#####################################################################
template<class T> T LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Get_Total_Delta_PE()
{
    return total_delta_PE;
}
//#####################################################################
// Function Save_And_Reset_Elastic_Coefficient
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Save_And_Reset_Elastic_Coefficient()
{
    saved_youngs_modulus.Resize(parameters.m);
    for(int i=1;i<=parameters.m;i++)
        for(int j=1;j<=3;j++){
            saved_youngs_modulus(i)(j)=parameters(i)(j).youngs_modulus;
            parameters(i)(j).youngs_modulus=0;}
}
//#####################################################################
// Function Restore_Elastic_Coefficient
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_S3D<T>::
Restore_Elastic_Coefficient()
{
    for(int i=1;i<=parameters.m;i++) for(int j=1;j<=3;j++) parameters(i)(j).youngs_modulus=saved_youngs_modulus(i)(j);
}
//#####################################################################
template class LINEAR_ALTITUDE_SPRINGS_S3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LINEAR_ALTITUDE_SPRINGS_S3D<double>;
#endif
