//#####################################################################
// Copyright 2002-2010, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Neil Molino, Igor Neverov, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using ::std::sqrt;
using namespace PhysBAM;
template<class T> LINEAR_ALTITUDE_SPRINGS_3D<T>::
LINEAR_ALTITUDE_SPRINGS_3D(PARTICLES<TV>& particles,TETRAHEDRON_MESH& tetrahedron_mesh)
    :LINEAR_ALTITUDE_SPRINGS<VECTOR<T,3>,3>(particles,tetrahedron_mesh)
{
    use_shortest_spring_only=true;
}
template<class T> LINEAR_ALTITUDE_SPRINGS_3D<T>::
~LINEAR_ALTITUDE_SPRINGS_3D()
{
}
//#####################################################################
// Function Set_Stiffness_Based_On_Reduced_Mass
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Set_Stiffness_Based_On_Reduced_Mass(const T scaling_coefficient) // assumes mass and restlength are already defined
{
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        for(int s=1;s<=4;s++){int node1,node2,node3,node4;
            Fill_Node_Indices(i,j,k,l,s,node1,node2,node3,node4);
            T one_over_triangle_mass=(T)one_ninth*(particles.one_over_effective_mass(node2)+particles.one_over_effective_mass(node3)+particles.one_over_effective_mass(node4)),
                  one_over_particle_mass=particles.one_over_effective_mass(node1),harmonic_mass=Pseudo_Inverse(one_over_particle_mass+one_over_triangle_mass);
            parameters(t)(s).youngs_modulus=scaling_coefficient*harmonic_mass/parameters(t)(s).restlength;}}
}
//#####################################################################
// Function Set_Restlength_From_Particles
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Set_Restlength_From_Particles()
{
    Set_Restlength_From_Material_Coordinates(particles.X);
}
//#####################################################################
// Function Set_Restlength_From_Material_Coordinates
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Set_Restlength_From_Material_Coordinates(ARRAY_VIEW<const TV> material_coordinates)
{
    Invalidate_CFL();
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        // TODO: use a for loop instead
        TV barycentric=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(material_coordinates(i),material_coordinates(j),material_coordinates(l),material_coordinates(k));
        parameters(t)(1).restlength=(material_coordinates(i)-barycentric.x*material_coordinates(j)-barycentric.y*material_coordinates(l)-barycentric.z*material_coordinates(k)).Normalize();
        barycentric=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(material_coordinates(j),material_coordinates(i),material_coordinates(k),material_coordinates(l));
        parameters(t)(2).restlength=(material_coordinates(j)-barycentric.x*material_coordinates(i)-barycentric.y*material_coordinates(k)-barycentric.z*material_coordinates(l)).Normalize();
        barycentric=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(material_coordinates(k),material_coordinates(i),material_coordinates(l),material_coordinates(j));
        parameters(t)(3).restlength=(material_coordinates(k)-barycentric.x*material_coordinates(i)-barycentric.y*material_coordinates(l)-barycentric.z*material_coordinates(j)).Normalize();
        barycentric=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(material_coordinates(l),material_coordinates(i),material_coordinates(j),material_coordinates(k));
        parameters(t)(4).restlength=(material_coordinates(l)-barycentric.x*material_coordinates(i)-barycentric.y*material_coordinates(j)-barycentric.z*material_coordinates(k)).Normalize();
        assert(parameters(t)(1).restlength>0);assert(parameters(t)(2).restlength>0);assert(parameters(t)(3).restlength>0);assert(parameters(t)(4).restlength>0);
        for(int k=1;k<=4;k++) parameters(t)(k).visual_restlength=parameters(t)(k).restlength;}
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Set_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        // TODO: use a for loop instead
        T harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(i)+(T)one_ninth*(particles.one_over_effective_mass(j)+particles.one_over_effective_mass(k)+particles.one_over_effective_mass(l)));
        parameters(t)(1).damping=overdamping_fraction*2*sqrt(parameters(t)(1).youngs_modulus*parameters(t)(1).restlength*harmonic_mass);
        harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(j)+(T)one_ninth*(particles.one_over_effective_mass(i)+particles.one_over_effective_mass(k)+particles.one_over_effective_mass(l)));
        parameters(t)(2).damping=overdamping_fraction*2*sqrt(parameters(t)(2).youngs_modulus*parameters(t)(2).restlength*harmonic_mass);
        harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(k)+(T)one_ninth*(particles.one_over_effective_mass(j)+particles.one_over_effective_mass(i)+particles.one_over_effective_mass(l)));
        parameters(t)(3).damping=overdamping_fraction*2*sqrt(parameters(t)(3).youngs_modulus*parameters(t)(3).restlength*harmonic_mass);
        harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(l)+(T)one_ninth*(particles.one_over_effective_mass(j)+particles.one_over_effective_mass(k)+particles.one_over_effective_mass(i)));
        parameters(t)(4).damping=overdamping_fraction*2*sqrt(parameters(t)(4).youngs_modulus*parameters(t)(4).restlength*harmonic_mass);}
    Invalidate_CFL();
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Set_Overdamping_Fraction(const ARRAY<VECTOR<T,4> >& overdamping_fraction) // 1 is critically damped
{
    for(int t=1;t<=mesh.elements.m;t++){
        int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        // TODO: use a for loop instead
        T harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(i)+(T)one_ninth*(particles.one_over_effective_mass(j)+particles.one_over_effective_mass(k)+particles.one_over_effective_mass(l)));
        parameters(t)(1).damping=overdamping_fraction(t)(1)*2*sqrt(parameters(t)(1).youngs_modulus*parameters(t)(1).restlength*harmonic_mass);
        harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(j)+(T)one_ninth*(particles.one_over_effective_mass(i)+particles.one_over_effective_mass(k)+particles.one_over_effective_mass(l)));
        parameters(t)(2).damping=overdamping_fraction(t)(2)*2*sqrt(parameters(t)(2).youngs_modulus*parameters(t)(2).restlength*harmonic_mass);
        harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(k)+(T)one_ninth*(particles.one_over_effective_mass(j)+particles.one_over_effective_mass(i)+particles.one_over_effective_mass(l)));
        parameters(t)(3).damping=overdamping_fraction(t)(3)*2*sqrt(parameters(t)(3).youngs_modulus*parameters(t)(3).restlength*harmonic_mass);
        harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(l)+(T)one_ninth*(particles.one_over_effective_mass(j)+particles.one_over_effective_mass(k)+particles.one_over_effective_mass(i)));
        parameters(t)(4).damping=overdamping_fraction(t)(4)*2*sqrt(parameters(t)(4).youngs_modulus*parameters(t)(4).restlength*harmonic_mass);}
    Invalidate_CFL();
}
//#####################################################################
// Function Ensure_Minimum_Overdamping_Fraction
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    ARRAY<VECTOR<T,4> > save_damping(parameters.m);for(int i=1;i<=parameters.m;i++) for(int k=1;k<=4;k++) save_damping(i)(k)=parameters(i)(k).damping;
    Set_Overdamping_Fraction(overdamping_fraction);
    for(int k=1;k<=parameters.m;k++) for(int i=1;i<=4;i++) parameters(k)(i).damping=max(parameters(k)(i).damping,save_damping(k)(i));
}
//#####################################################################
// Function Fill_Node_Indices
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Fill_Node_Indices(int i,int j,int k,int l,int isolated_node_number,int& node1,int& node2,int& node3,int& node4) const
{
    // node1 is the isolated vertex and nodes2,3,4 are the triangle
    switch(isolated_node_number){
        case 1:node1=i;node2=j;node3=l;node4=k;break;
        case 2:node1=j;node2=i;node3=k;node4=l;break;
        case 3:node1=k,node2=i,node3=l,node4=j;break;
        default:node1=l;node2=i;node3=j;node4=k;}
}
//#####################################################################
// Function Fill_Spring_State
//#####################################################################
template<class T> bool LINEAR_ALTITUDE_SPRINGS_3D<T>::
Fill_Spring_State(int t,int isolated_node_number,int node1,int node2,int node3,int node4,SPRING_STATE& state)
{
    ARRAY_VIEW<const TV> X=particles.X;
    state.barycentric=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(X(node1),X(node2),X(node3),X(node4));
    state.direction=X(node1)-state.barycentric.x*X(node2)-state.barycentric.y*X(node3)-state.barycentric.z*X(node4);
    state.current_length=state.direction.Normalize();
    TV normal=TRIANGLE_3D<T>::Normal(X(node2),X(node3),X(node4));
    if(!state.current_length) state.direction=normal;
    else if(TV::Dot_Product(state.direction,normal)<0){
        state.direction=-state.direction;
        state.current_length=-state.current_length;}
    T rl=parameters(t)(isolated_node_number).restlength;
    if(use_springs_compressed_beyond_threshold && state.current_length-rl>-spring_compression_fraction_threshold*rl)
    {state.node=0;return false;}
    else{
        state.coefficient=parameters(t)(isolated_node_number).damping/rl;
        if(use_plasticity) Compute_Plasticity(isolated_node_number,t,state.current_length);
        state.node=isolated_node_number;
        return true;}
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    if(use_shortest_spring_only) spring_states.Resize(mesh.elements.m);
    else spring_states_all_springs.Resize(mesh.elements.m);
    int used_springs=0,total_elements=0;
    ARRAY_VIEW<const TV> X=particles.X;int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data(); // use shortest spring only
        total_elements++;
        int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        if(use_shortest_spring_only){
            int hmin=0;T cross_area_max=-(T)FLT_MAX;
            if(is_position_update){
                for(int h=1;h<=4;h++){
                    Fill_Node_Indices(i,j,k,l,h,node1,node2,node3,node4);
                    T cross_area=TV::Cross_Product(X(node3)-X(node2),X(node4)-X(node2)).Magnitude_Squared();
                    if(cross_area>cross_area_max){hmin=h;cross_area_max=cross_area;}}}
            else hmin=spring_states(t).node;
            Fill_Node_Indices(i,j,k,l,hmin,node1,node2,node3,node4);
            SPRING_STATE& state=spring_states(t);
            if(Fill_Spring_State(t,hmin,node1,node2,node3,node4,state)) used_springs++;}
        else{
            for(int h=1;h<=4;h++){
                Fill_Node_Indices(i,j,k,l,h,node1,node2,node3,node4);
                SPRING_STATE& state=spring_states_all_springs(t)[h];
                if(Fill_Spring_State(t,h,node1,node2,node3,node4,state)) used_springs++;}}}

    residual_PE.Resize(mesh.elements.m);
    if(!total_PE.m){
        total_PE.Resize(mesh.elements.m);
        for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            total_PE(t)=Potential_Energy(t,time);}}

    if(print_number_used) {std::stringstream ss;ss<<"using "<<used_springs<<" of "<<total_elements<<" altitude springs"<<std::endl;LOG::filecout(ss.str());}
    if(compute_half_forces){
        if(use_shortest_spring_only) for(int i=1;i<=spring_states.m;i++){
            spring_states(i).sqrt_coefficient=sqrt(spring_states(i).coefficient);}
        else for(int i=1;i<=spring_states_all_springs.m;i++) for(int h=1;h<=4;h++){
                spring_states_all_springs(i)(h).sqrt_coefficient=sqrt(spring_states_all_springs(i)(h).coefficient);}}
    if(!mesh.incident_elements) mesh.Initialize_Incident_Elements();
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int total_springs=use_shortest_spring_only?1:4;
        for(int spring_index=1;spring_index<=total_springs;spring_index++){
            const SPRING_STATE* state_ptr;
            if(use_shortest_spring_only) state_ptr=&spring_states(t);
            else state_ptr=&spring_states_all_springs(t)[spring_index];
            const SPRING_STATE& state=*state_ptr;
            if(state.node){
                int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
                int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
                const SPRING_PARAMETER& parameter=parameters(t)(state.node);
                Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);
                T rl=parameter.restlength,vrl=parameter.visual_restlength;
                TV force=parameter.youngs_modulus/rl*(state.current_length-vrl)*state.direction;
                F(node1)-=force;F(node2)+=state.barycentric.x*force;F(node3)+=state.barycentric.y*force;F(node4)+=state.barycentric.z*force;}}}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int total_springs=use_shortest_spring_only?1:4;
        for(int spring_index=1;spring_index<=total_springs;spring_index++){
            const SPRING_STATE* state_ptr;
            if(use_shortest_spring_only) state_ptr=&spring_states(t);
            else state_ptr=&spring_states_all_springs(t)[spring_index];
            const SPRING_STATE& state=*state_ptr;
            if(state.node){
                int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
                int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
                Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);
                TV force=(state.coefficient*TV::Dot_Product(V(node1)-state.barycentric.x*V(node2)-state.barycentric.y*V(node3)-state.barycentric.z*V(node4),
                        state.direction))*state.direction;
                F(node1)-=force;F(node2)+=state.barycentric.x*force;F(node3)+=state.barycentric.y*force;F(node4)+=state.barycentric.z*force;}}}
}
//#####################################################################
// Function Velocity_Dependent_Size
//#####################################################################
template<class T> int LINEAR_ALTITUDE_SPRINGS_3D<T>::
Velocity_Dependent_Forces_Size() const
{
    int aggregate_id=0;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int total_springs=use_shortest_spring_only?1:4;
        for(int spring_index=1;spring_index<=total_springs;spring_index++){
            const SPRING_STATE* state_ptr;
            if(use_shortest_spring_only) state_ptr=&spring_states(t);
            else state_ptr=&spring_states_all_springs(t)[spring_index];
            const SPRING_STATE& state=*state_ptr;
            if(state.node) aggregate_id++;}}
    return aggregate_id;
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    int aggregate_id=1;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int total_springs=use_shortest_spring_only?1:4;
        for(int spring_index=1;spring_index<=total_springs;spring_index++){
            const SPRING_STATE* state_ptr;
            if(use_shortest_spring_only) state_ptr=&spring_states(t);
            else state_ptr=&spring_states_all_springs(t)[spring_index];
            const SPRING_STATE& state=*state_ptr;
            if(state.node){
                int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
                int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
                Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);
                aggregate(aggregate_id++)+=state.sqrt_coefficient*TV::Dot_Product(V(node1)-state.barycentric.x*V(node2)-state.barycentric.y*V(node3)-state.barycentric.z*V(node4),state.direction);}}}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    int aggregate_id=1;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int total_springs=use_shortest_spring_only?1:4;
        for(int spring_index=1;spring_index<=total_springs;spring_index++){
            const SPRING_STATE* state_ptr;
            if(use_shortest_spring_only) state_ptr=&spring_states(t);
            else state_ptr=&spring_states_all_springs(t)[spring_index];
            const SPRING_STATE& state=*state_ptr;
            if(state.node){
                int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
                int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
                Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);
                TV force=state.sqrt_coefficient*aggregate(aggregate_id++)*state.direction;
                F(node1)+=force;F(node2)-=state.barycentric.x*force;F(node3)-=state.barycentric.y*force;F(node4)-=state.barycentric.z*force;}}}
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int total_springs=use_shortest_spring_only?1:4;
        for(int spring_index=1;spring_index<=total_springs;spring_index++){
            const SPRING_STATE* state_ptr;
            if(use_shortest_spring_only) state_ptr=&spring_states(t);
            else state_ptr=&spring_states_all_springs(t)[spring_index];
            const SPRING_STATE& state=*state_ptr;
            if(state.node){
                int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
                int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
                const SPRING_PARAMETER& parameter=parameters(t)(state.node);
                Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);
                T rl=parameter.restlength;
                TV force=parameter.youngs_modulus/rl*TV::Dot_Product(V(node1)-state.barycentric.x*V(node2)-state.barycentric.y*V(node3)-state.barycentric.z*V(node4),state.direction)*state.direction;
                F(node1)-=force;F(node2)+=state.barycentric.x*force;F(node3)+=state.barycentric.y*force;F(node4)+=state.barycentric.z*force;}}}
}
template<class T> bool LINEAR_ALTITUDE_SPRINGS_3D<T>::
Compute_Strain_Rate_And_Strain(int t,int isolated_node_number,int node1,int node2,int node3,int node4,T& strain_rate,T& strain) const
{
    ARRAY_VIEW<const TV> X(particles.X),V(particles.V);
    TV direction=TRIANGLE_3D<T>::Normal(X(node2),X(node3),X(node4));
    T rl=parameters(t)(isolated_node_number).restlength;
    T dx;if(use_rest_state_for_strain_rate) dx=rl;else dx=TV::Dot_Product(direction,X(node1)-X(node2));
    if(use_springs_compressed_beyond_threshold){
        T length;if(use_rest_state_for_strain_rate) length=TV::Dot_Product(direction,X(node1)-X(node2));else length=dx;
        if(length-rl>-spring_compression_fraction_threshold*rl) return false;}
    TV barycentric=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(X(node1),X(node2),X(node3),X(node4));
    TV v_interpolated=barycentric.x*V(node2)+barycentric.y*V(node3)+barycentric.z*V(node4);
    strain_rate=TV::Dot_Product(V(node1)-v_interpolated,direction)/dx;
    strain=(dx-parameters(t)(isolated_node_number).visual_restlength)/rl;
    return true;
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class T> T LINEAR_ALTITUDE_SPRINGS_3D<T>::
CFL_Strain_Rate() const
{
    T max_strain_rate=0;int node1,node2,node3,node4; // node1 is the vertex and nodes2,3,4 are the triangle
    T strain_rate,strain;
    ARRAY<VECTOR<int,4> >& elements=mesh.elements;ARRAY_VIEW<const TV> X(particles.X);
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int i,j,k,l;elements(t).Get(i,j,k,l);
        if(use_shortest_spring_only){
            int hmin=0;T cross_area_max=(T)-FLT_MAX;
            for(int h=1;h<=4;h++){
                Fill_Node_Indices(i,j,k,l,h,node1,node2,node3,node4);
                T cross_area=TV::Cross_Product(X(node3)-X(node2),X(node4)-X(node2)).Magnitude_Squared();
                if(cross_area>cross_area_max){hmin=h;cross_area_max=cross_area;}}
            Fill_Node_Indices(i,j,k,l,hmin,node1,node2,node3,node4);
            if(!Compute_Strain_Rate_And_Strain(t,hmin,node1,node2,node3,node4,strain_rate,strain)) continue;
            max_strain_rate=max(max_strain_rate,abs(strain_rate));
            if(cache_strain){strains_of_spring(t)=VECTOR<T,2>(abs(strain_rate),abs(strain));}}
        else{
            for(int h=1;h<=4;h++){
                Fill_Node_Indices(i,j,k,l,h,node1,node2,node3,node4);
                if(!Compute_Strain_Rate_And_Strain(t,h,node1,node2,node3,node4,strain_rate,strain)) continue;
                Compute_Strain_Rate_And_Strain(t,h,node1,node2,node3,node4,strain_rate,strain);
                max_strain_rate=max(max_strain_rate,abs(strain_rate));
                if(cache_strain){strains_of_spring_all_springs(t)[h]=VECTOR<T,2>(abs(strain_rate),abs(strain));}}}}
    return Robust_Divide(max_strain_per_time_step,max_strain_rate);
}
//#####################################################################
// Function Add_Force_Data
//#####################################################################
template<class TV> void LINEAR_ALTITUDE_SPRINGS_3D<TV>::
Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name) const
{
    if((use_shortest_spring_only && !spring_states.m) || (!use_shortest_spring_only && !spring_states_all_springs.m)) return;
    ARRAY_VIEW<const TV> X=particles.X;
    FORCE_DATA<TV> force_data;
    if(force_name.empty()) force_data.name="LINEAR_ALTITUDE_SPRINGS_3D";
    else force_data.name=force_name;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int total_springs=use_shortest_spring_only?1:4;
        for(int spring_index=1;spring_index<=total_springs;spring_index++){
            const SPRING_STATE* state_ptr;
            if(use_shortest_spring_only) state_ptr=&spring_states(t);
            else state_ptr=&spring_states_all_springs(t)[spring_index];
            const SPRING_STATE& state=*state_ptr;
            if(state.node){
                int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
                int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
                const SPRING_PARAMETER& parameter=parameters(t)(state.node);
                Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);

                force_data.state=(state.current_length-parameter.visual_restlength)/parameter.visual_restlength;
                force_data.first_action_point=X(node1);
                force_data.second_action_point=state.barycentric.x*X(node2)+state.barycentric.y*X(node3)+state.barycentric.z*X(node4);
                force_data_list.Append(force_data);}}}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class T> T LINEAR_ALTITUDE_SPRINGS_3D<T>::
Potential_Energy(const int t,const T time) const
{
    const SPRING_STATE& state=spring_states(t);
    if(state.node){
        int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
        const SPRING_PARAMETER& parameter=parameters(t)(state.node);
        Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);
        T rl=parameter.restlength,vrl=parameter.visual_restlength;

        ARRAY_VIEW<const TV> X=particles.X;
        VECTOR<T,3> barycentric=TRIANGLE_3D<T>::Clamped_Barycentric_Coordinates(X(node1),X(node2),X(node3),X(node4));
        TV direction=X(node1)-barycentric.x*X(node2)-barycentric.y*X(node3)-barycentric.z*X(node4);
        T current_length=direction.Normalize();
        TV normal=TRIANGLE_3D<T>::Normal(X(node2),X(node3),X(node4));
        if(!current_length) direction=normal;
        else if(TV::Dot_Product(direction,normal)<0){
            direction=-direction;
            current_length=-current_length;}
        return (T).5*parameter.youngs_modulus/rl*sqr(current_length-vrl);}
    else return 0;
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class T> T LINEAR_ALTITUDE_SPRINGS_3D<T>::
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
template<class T> T LINEAR_ALTITUDE_SPRINGS_3D<T>::
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
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Save_Potential_Energy(const T time)
{
    potential_energy_save.Resize(mesh.elements.m);
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        potential_energy_save(t)=Potential_Energy(t,time);}
}
//#####################################################################
// Function Setup_Set_Velocity_From_Positions
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Setup_Set_Velocity_From_Positions(const T time,const bool is_position_update,const bool reset_alphas)
{
    if(!is_position_update){
        // UPBS isn't called, so do most of the work, except keep spring direction the same
        if(use_shortest_spring_only) spring_states.Resize(mesh.elements.m);
        else spring_states_all_springs.Resize(mesh.elements.m);
        int used_springs=0,total_elements=0;
        ARRAY_VIEW<const TV> X=particles.X;int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
        for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data(); // use shortest spring only
            total_elements++;
            int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
            if(use_shortest_spring_only){
                SPRING_STATE& state=spring_states(t);
                int hmin=state.node;
                Fill_Node_Indices(i,j,k,l,hmin,node1,node2,node3,node4);
                if(Fill_Spring_State(t,hmin,node1,node2,node3,node4,state)) used_springs++;}
            else{
                for(int h=1;h<=4;h++){
                    Fill_Node_Indices(i,j,k,l,h,node1,node2,node3,node4);
                    SPRING_STATE& state=spring_states_all_springs(t)[h];
                    if(Fill_Spring_State(t,h,node1,node2,node3,node4,state)) used_springs++;}}}}

    force_estimates.Resize(mesh.elements.m);
    incident_nodes.Resize(mesh.elements.m);
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        int i,j,k,l;mesh.elements(t).Get(i,j,k,l);
        incident_nodes(t).Remove_All();
        incident_nodes(t).Append(i);incident_nodes(t).Append(j);incident_nodes(t).Append(k);incident_nodes(t).Append(l);}
    if(is_position_update){
        for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            if(reset_alphas) force_estimates(t)=0;
            residual_PE(t)=Potential_Energy(t,time)-total_PE(t);}
        Save_Potential_Energy(time);}
}
//#####################################################################
// Function Get_Element_Count
//#####################################################################
template<class T> int LINEAR_ALTITUDE_SPRINGS_3D<T>::
Get_Element_Count()
{
    return mesh.elements.m;
}
//#####################################################################
// Function Get_Force_Elements
//#####################################################################
template<class T> FORCE_ELEMENTS* LINEAR_ALTITUDE_SPRINGS_3D<T>::
Get_Force_Elements()
{
    return &force_elements;
}
//#####################################################################
// Function Incident_Nodes
//#####################################################################
template<class T> ARRAY<int>* LINEAR_ALTITUDE_SPRINGS_3D<T>::
Incident_Nodes(const int force_element)
{
    return &(incident_nodes(force_element));
}
//#####################################################################
// Function Get_Direction
//#####################################################################
template<class T> VECTOR<T,3> LINEAR_ALTITUDE_SPRINGS_3D<T>::
Get_Direction(const int force_element)
{
    const SPRING_STATE& state=spring_states(force_element);
    return state.direction;
}
//#####################################################################
// Function Get_Combined_One_Over_Mass
//#####################################################################
template<class T> T LINEAR_ALTITUDE_SPRINGS_3D<T>::
Get_Combined_One_Over_Mass(const int force_element)
{
    const SPRING_STATE& state=spring_states(force_element);
    if(state.node){
        int i,j,k,l;mesh.elements(force_element).Get(i,j,k,l);
        int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
        Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);
        return particles.one_over_mass(node1)+
            (sqr(state.barycentric.x)*particles.one_over_mass(node2)+sqr(state.barycentric.y)*particles.one_over_mass(node3)+sqr(state.barycentric.z)*particles.one_over_mass(node4));}
    return 0;
}
//#####################################################################
// Function Incident_Force_Elements
//#####################################################################
template<class T> ARRAY<int>* LINEAR_ALTITUDE_SPRINGS_3D<T>::
Incident_Force_Elements(const int particle)
{
    return &(*mesh.incident_elements)(particle);
}
//#####################################################################
// Function Get_Force
//#####################################################################
template<class T> VECTOR<T,3> LINEAR_ALTITUDE_SPRINGS_3D<T>::
Get_Force(const int force_element,const int particle,const bool use_original_force)
{
    const SPRING_STATE& state=spring_states(force_element);
    if(state.node){
        int i,j,k,l;mesh.elements(force_element).Get(i,j,k,l);
        int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
        Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);
        TV force;
        if(use_original_force){
            const SPRING_PARAMETER& parameter=parameters(force_element)(state.node);
            T rl=parameter.restlength,vrl=parameter.visual_restlength;
            force=parameter.youngs_modulus/rl*(state.current_length-vrl)*state.direction;}
        else force=force_estimates(force_element)*state.direction;
        if(particle==node1) return -force;
        else if(particle==node2) return state.barycentric.x*force;
        else if(particle==node3) return state.barycentric.y*force;
        else if(particle==node4) return state.barycentric.z*force;}
    return TV();
}
//#####################################################################
// Function Choose_Solution
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
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
            int i,j,k,l;mesh.elements(force_element).Get(i,j,k,l);
            int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
            T v1_n_hat;
            T v2_n_hat;

            switch(state.node){
                case 1:{
                    node1=i;node2=j;node3=l;node4=k;
                    v1_n_hat=state.barycentric.x*v_n_hats(2)+state.barycentric.y*v_n_hats(4)+state.barycentric.z*v_n_hats(3);
                    v2_n_hat=v_n_hats(1);}
                    break;
                case 2:{
                    node1=j;node2=i;node3=k;node4=l;
                    v1_n_hat=state.barycentric.x*v_n_hats(1)+state.barycentric.y*v_n_hats(3)+state.barycentric.z*v_n_hats(4);
                    v2_n_hat=v_n_hats(2);}
                    break;
                case 3:{
                    node1=k,node2=i,node3=l,node4=j;
                    v1_n_hat=state.barycentric.x*v_n_hats(1)+state.barycentric.y*v_n_hats(4)+state.barycentric.z*v_n_hats(2);
                    v2_n_hat=v_n_hats(3);}
                    break;
                default:{
                    node1=l;node2=i;node3=j;node4=k;
                    v1_n_hat=state.barycentric.x*v_n_hats(1)+state.barycentric.y*v_n_hats(2)+state.barycentric.z*v_n_hats(3);
                    v2_n_hat=v_n_hats(4);}}
            
            T v1_one_over_mass,v1_mass;
            if(particles.mass(node2)==FLT_MAX && particles.mass(node3)==FLT_MAX && particles.mass(node4)==FLT_MAX){
                v1_one_over_mass=0;
                v1_mass=FLT_MAX;}
            else{
                v1_one_over_mass=sqr(state.barycentric.x)*particles.one_over_mass(node2)+sqr(state.barycentric.y)*particles.one_over_mass(node3)+
                    sqr(state.barycentric.z)*particles.one_over_mass(node4);
                v1_mass = (particles.mass(node2)*particles.mass(node3)*particles.mass(node4))/(sqr(state.barycentric.x)*particles.mass(node3)*particles.mass(node4)+sqr(state.barycentric.y)*particles.mass(node2)*particles.mass(node4)+sqr(state.barycentric.z)*particles.mass(node2)*particles.mass(node3));}
            T v1_alpha1_temp=v1_n_hat+dt*v1_one_over_mass*alpha1;
            T v2_alpha1_temp=v2_n_hat-dt*particles.one_over_mass(node1)*alpha1;
            
            T v1_alpha2_temp=v1_n_hat+dt*v1_one_over_mass*alpha2;
            T v2_alpha2_temp=v2_n_hat-dt*particles.one_over_mass(node1)*alpha2;
            
            T combined_mass=particles.mass(node1)+v1_mass;
            TV x1=state.barycentric.x*particles.X(node2)+state.barycentric.y*particles.X(node3)+state.barycentric.z*particles.X(node4);
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
template<class T> T LINEAR_ALTITUDE_SPRINGS_3D<T>::
Get_Force(const int force_element)
{
    return force_estimates(force_element);
}
//#####################################################################
// Function Set_Force
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Set_Force(const int force_element,const T force)
{
    force_estimates(force_element)=force;
}
//#####################################################################
// Function Get_Damping_Force
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Get_Damping_Force(const int particle,TV& damping_force,const T dt,const bool use_coefficient)
{
    ARRAY<int>& node_incident_forces=(*mesh.incident_elements)(particle);
    for(int t=1;t<=node_incident_forces.m;t++){
        int total_springs=use_shortest_spring_only?1:4;
        for(int spring_index=1;spring_index<=total_springs;spring_index++){
            const SPRING_STATE* state_ptr;
            if(use_shortest_spring_only) state_ptr=&spring_states(node_incident_forces(t));
            else state_ptr=&spring_states_all_springs(node_incident_forces(t))[spring_index];
            const SPRING_STATE& state=*state_ptr;
            if(state.node){
                int i,j,k,l;mesh.elements(node_incident_forces(t)).Get(i,j,k,l);
                int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
                Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);
                TV force=(TV::Dot_Product(particles.V(node1)-state.barycentric.x*particles.V(node2)-state.barycentric.y*particles.V(node3)-state.barycentric.z*particles.V(node4),state.direction))*state.direction;
                if(particle==node1) force*=(T)-1;
                else if(particle==node2) force*=state.barycentric.x;
                else if(particle==node3) force*=state.barycentric.y;
                else force*=state.barycentric.z;
                if(use_coefficient) damping_force+=dt*particles.one_over_mass(particle)*state.coefficient*force;
                else damping_force+=dt*particles.one_over_mass(particle)*force;}}}
}
//#####################################################################
// Function Update_Residual_Energy
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Update_Residual_Energy(const int force_element,const T residual_energy,const T time)
{
    residual_PE(force_element)=residual_energy;
    total_PE(force_element)=Potential_Energy(force_element,time)-residual_PE(force_element);
}
//#####################################################################
// Function Update_Residual_Energy
//#####################################################################
template<class T> T LINEAR_ALTITUDE_SPRINGS_3D<T>::
Get_Residual_Energy(const int force_element)
{
    return residual_PE(force_element);
}
//#####################################################################
// Function Compute_Quadratic_Contribution_For_Node
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Compute_Quadratic_Contribution_For_Force(T& A,T& a,T&c,const T dt,const int force_element,const T combined_one_over_mass,const bool ignore_PE_terms)
{
    const SPRING_STATE& state=spring_states(force_element);
    if(state.node){
        if(!ignore_PE_terms){
            int i,j,k,l;mesh.elements(force_element).Get(i,j,k,l);
            int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
            Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);

            T x1u=TV::Dot_Product(particles.X(node1),state.direction);
            T x2u=TV::Dot_Product(particles.X(node2),state.direction);
            T x3u=TV::Dot_Product(particles.X(node3),state.direction);
            T x4u=TV::Dot_Product(particles.X(node4),state.direction);
            
            a=x1u-(state.barycentric.x*x2u+state.barycentric.y*x3u+state.barycentric.z*x4u);
            c=-(T).5*dt*dt*combined_one_over_mass;
            
            const SPRING_PARAMETER& parameter=parameters(force_element)(state.node);
            A+=(T).5*(parameter.youngs_modulus/parameter.restlength)*c*c;}
        
        A+=(T).5*dt*dt*combined_one_over_mass;}
}
//#####################################################################
// Function Compute_Quadratic_Contribution_For_Node
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Compute_Quadratic_Contribution_For_Node(T& B,T& C,T&b,const T dt,const int node,const int force_element,const T combined_one_over_mass,const T v_n_correction,const bool ignore_PE_terms)
{
    const SPRING_STATE& state=spring_states(force_element);
    if(state.node){
        int i,j,k,l;mesh.elements(force_element).Get(i,j,k,l);
        int node1,node2,node3,node4; // node1 is the isolated vertex and nodes2,3,4 are the triangle
        Fill_Node_Indices(i,j,k,l,state.node,node1,node2,node3,node4);

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
            B+=term;}
        else if(node==node4){
            T term=dt*state.barycentric.z*(vu+(T).5*v_n_correction);
            if(!ignore_PE_terms) b+=-term;
            B+=term;}}
}
//#####################################################################
// Function Compute_Quadratic_Contribution_For_Residual
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
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
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
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
template<class T> T LINEAR_ALTITUDE_SPRINGS_3D<T>::
Get_Total_Delta_PE()
{
    return total_delta_PE;
}
//#####################################################################
// Function Save_And_Reset_Elastic_Coefficient
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Save_And_Reset_Elastic_Coefficient()
{
    saved_youngs_modulus.Resize(parameters.m);
    for(int i=1;i<=parameters.m;i++)
        for(int j=1;j<=4;j++){
            saved_youngs_modulus(i)(j)=parameters(i)(j).youngs_modulus;
            parameters(i)(j).youngs_modulus=0;}
}
//#####################################################################
// Function Restore_Elastic_Coefficient
//#####################################################################
template<class T> void LINEAR_ALTITUDE_SPRINGS_3D<T>::
Restore_Elastic_Coefficient()
{
    for(int i=1;i<=parameters.m;i++) for(int j=1;j<=4;j++) parameters(i)(j).youngs_modulus=saved_youngs_modulus(i)(j);
}
//#####################################################################

//#####################################################################
// Function Create_Altitude_Springs
//#####################################################################
template<class T> LINEAR_ALTITUDE_SPRINGS_3D<T>* PhysBAM::
Create_Altitude_Springs(PARTICLES<VECTOR<T,3> >& particles,TETRAHEDRON_MESH& mesh,
    const T stiffness,const T overdamping_fraction,const bool use_compressed_by_threshold_only,const T fraction_compression,const bool limit_time_step_by_strain_rate,
    const T max_strain_per_time_step,const bool use_rest_state_for_strain_rate,const T restlength_enlargement_fraction,const bool verbose)
{
    return Create_Altitude_Springs_Base(particles,mesh,stiffness,overdamping_fraction,use_compressed_by_threshold_only,fraction_compression,
    limit_time_step_by_strain_rate,max_strain_per_time_step,use_rest_state_for_strain_rate,restlength_enlargement_fraction,verbose);
}
//#####################################################################
// Function Create_Altitude_Springs
//#####################################################################
template<class T> LINEAR_ALTITUDE_SPRINGS_3D<T>* PhysBAM::
Create_Altitude_Springs(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,
    const T stiffness,const T overdamping_fraction,const bool use_compressed_by_threshold_only,const T fraction_compression,
    const bool limit_time_step_by_strain_rate,const T max_strain_per_time_step,const bool use_rest_state_for_strain_rate,const T restlength_enlargement_fraction,
    const bool verbose)
{
    return Create_Altitude_Springs(dynamic_cast<PARTICLES<VECTOR<T,3> >&>(tetrahedralized_volume.particles),tetrahedralized_volume.mesh,stiffness,overdamping_fraction,use_compressed_by_threshold_only,
    fraction_compression,limit_time_step_by_strain_rate,max_strain_per_time_step,use_rest_state_for_strain_rate,restlength_enlargement_fraction,verbose);
}

template class LINEAR_ALTITUDE_SPRINGS_3D<float>;
template LINEAR_ALTITUDE_SPRINGS_3D<float>* PhysBAM::Create_Altitude_Springs<float>(TETRAHEDRALIZED_VOLUME<float>&,float,float,bool,float,bool,float,bool,float,bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LINEAR_ALTITUDE_SPRINGS_3D<double>;
template LINEAR_ALTITUDE_SPRINGS_3D<double>* PhysBAM::Create_Altitude_Springs<double>(TETRAHEDRALIZED_VOLUME<double>&,double,double,bool,double,bool,double,bool,double,bool);
#endif
