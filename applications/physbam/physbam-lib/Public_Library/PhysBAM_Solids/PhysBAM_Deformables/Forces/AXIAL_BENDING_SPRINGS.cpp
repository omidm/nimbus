//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Math_Tools/cyclic_shift.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_3D.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/AXIAL_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using ::std::sqrt;
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> AXIAL_BENDING_SPRINGS<T>::
AXIAL_BENDING_SPRINGS(PARTICLES<TV>& particles_input,TRIANGLE_MESH& triangle_mesh_input)
    :DEFORMABLES_FORCES<TV>(particles_input),triangle_mesh(triangle_mesh_input),use_kinetic_energy_fix(true),relaxation_fraction(1),
    use_gauss_seidel_in_energy_correction(false),allow_kd_direction_flip(false),verbose(false)
{
    Initialize();
    Set_Stiffness(0);
    ARRAYS_COMPUTATIONS::Fill(restlength,(T)0);
    ARRAYS_COMPUTATIONS::Fill(visual_restlength,(T)0);
    Set_Damping(0);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> AXIAL_BENDING_SPRINGS<T>::
~AXIAL_BENDING_SPRINGS()
{}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    triangle_mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    force_springs.Update(spring_particles,particle_is_simulated);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Initialize()
{
    bool adjacent_triangles_defined=(triangle_mesh.adjacent_elements!=0);if(!adjacent_triangles_defined) triangle_mesh.Initialize_Adjacent_Elements();

    int number_quadruples=0;
    for(int t=1;t<=triangle_mesh.elements.m;t++) for(int a=1;a<=(*triangle_mesh.adjacent_elements)(t).m;a++) if((*triangle_mesh.adjacent_elements)(t)(a)>t) number_quadruples++;
    spring_particles.Resize(number_quadruples,false);youngs_modulus.Resize(number_quadruples);restlength.Resize(number_quadruples);
    visual_restlength.Resize(number_quadruples);damping.Resize(number_quadruples);attached_edge_length.Resize(number_quadruples);attached_edge_restlength.Resize(number_quadruples);
    int index=0; // reset number
    for(int t=1;t<=triangle_mesh.elements.m;t++){
        int t1,t2,t3;triangle_mesh.elements(t).Get(t1,t2,t3);
        for(int a=1;a<=(*triangle_mesh.adjacent_elements)(t).m;a++){
            int s=(*triangle_mesh.adjacent_elements)(t)(a);
            if(s>t){
                int s1,s2,s3;triangle_mesh.elements(s).Get(s1,s2,s3);
                if(t1==s1 || t1==s2 || t1==s3){cyclic_shift(t1,t2,t3);if(t1==s1 || t1==s2 || t1==s3) cyclic_shift(t1,t2,t3);}
                spring_particles(++index).Set(t2,t3,t1,triangle_mesh.Other_Node(t2,t3,s));}}}

    if(!adjacent_triangles_defined){delete triangle_mesh.adjacent_elements;triangle_mesh.adjacent_elements=0;}
}
//#####################################################################
// Function Set_Restlength_From_Particles
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Set_Restlength_From_Particles()
{
    for(int s=1;s<=spring_particles.m;s++){const VECTOR<int,4>& nodes=spring_particles(s);
        TV axial_direction;VECTOR<T,2> weights;
        Axial_Vector(nodes,visual_restlength(s),axial_direction,weights,attached_edge_restlength(s));
        restlength(s)=visual_restlength(s);}
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Set_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    Invalidate_CFL();
    for(int s=1;s<=spring_particles.m;s++){
        T harmonic_mass=Pseudo_Divide((T)4,ARRAYS_COMPUTATIONS::Sum(particles.one_over_effective_mass.Subset(spring_particles(s))));
        damping(s)=overdamping_fraction*2*sqrt(youngs_modulus(s)*restlength(s)*harmonic_mass);}
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Set_Overdamping_Fraction(ARRAY_VIEW<const T> overdamping_fraction) // 1 is critically damped
{
    Invalidate_CFL();
    for(int s=1;s<=spring_particles.m;s++){
        T harmonic_mass=Pseudo_Divide((T)4,ARRAYS_COMPUTATIONS::Sum(particles.one_over_effective_mass.Subset(spring_particles(s))));
        damping(s)=overdamping_fraction(s)*2*sqrt(youngs_modulus(s)*restlength(s)*harmonic_mass);}
}
//#####################################################################
// Function Ensure_Minimum_Overdamping_Fraction
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    Invalidate_CFL();
    ARRAY<T> save_damping(damping);Set_Overdamping_Fraction(overdamping_fraction);
    for(int k=1;k<=damping.m;k++) damping(k)=max(damping(k),save_damping(k));
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    optimization_current_length.Resize(spring_particles.m,false,false);optimization_direction.Resize(spring_particles.m,false,false);
    optimization_weights.Resize(spring_particles.m,false,false);optimization_coefficient.Resize(spring_particles.m,false,false);
    extra_energy.Resize(spring_particles.m);
    force_correction.Resize(spring_particles.m);
    previously_applied_forces.Resize(spring_particles.m);

    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const VECTOR<int,4>& nodes=spring_particles(s);
        VECTOR<T,2> weights;
        Axial_Vector(nodes,optimization_current_length(s),optimization_direction(s),weights,attached_edge_length(s));
        optimization_weights(s).Set(1-weights.x,weights.x,1-weights.y,weights.y);
        optimization_coefficient(s)=damping(s)/restlength(s);}

    residual_PE.Resize(spring_particles.m);
    if(!incident_elements.m){
        incident_elements.Resize(particles.array_collection->Size());
        for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
            const VECTOR<int,4>& nodes=spring_particles(s);
            incident_elements(nodes(1)).Append(s);incident_elements(nodes(2)).Append(s);incident_elements(nodes(3)).Append(s);incident_elements(nodes(4)).Append(s);}}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        int node1,node2,node3,node4;spring_particles(s).Get(node1,node2,node3,node4);
        T w1,w2,w3,w4;optimization_weights(s).Get(w1,w2,w3,w4);
        TV force=youngs_modulus(s)/restlength(s)/*sqr(attached_edge_length(s)/attached_edge_restlength(s)-1)*/*(optimization_current_length(s)-visual_restlength(s))*optimization_direction(s);
        F(node1)-=w1*force;F(node2)-=w2*force;F(node3)+=w3*force;F(node4)+=w4*force;}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        int node1,node2,node3,node4;spring_particles(s).Get(node1,node2,node3,node4);
        T w1,w2,w3,w4;optimization_weights(s).Get(w1,w2,w3,w4);
        TV force=(optimization_coefficient(s)*TV::Dot_Product(w1*V(node1)+w2*V(node2)-w3*V(node3)-w4*V(node4),optimization_direction(s)))*optimization_direction(s);
        F(node1)-=w1*force;F(node2)-=w2*force;F(node3)+=w3*force;F(node4)+=w4*force;}
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        int node1,node2,node3,node4;spring_particles(s).Get(node1,node2,node3,node4);
        T w1,w2,w3,w4;optimization_weights(s).Get(w1,w2,w3,w4);
        TV dl=w1*V(node1)+w2*V(node2)-w3*V(node3)-w4*V(node4);
        TV force=youngs_modulus(s)/restlength(s)/*sqr(attached_edge_length(s)/attached_edge_restlength(s)-1)*/*dl.Projected_On_Unit_Direction(optimization_direction(s));
        F(node1)-=w1*force;F(node2)-=w2*force;F(node3)+=w3*force;F(node4)+=w4*force;}
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{
    T one_over_cfl_number=1/cfl_number,one_over_cfl_number_squared=sqr(one_over_cfl_number);
    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const VECTOR<int,4>& nodes=spring_particles(s);
        T one_over_mass_times_restlength=(T).25/restlength(s)*ARRAYS_COMPUTATIONS::Sum(particles.one_over_effective_mass.Subset(nodes));
        T elastic_hertz_squared=4*youngs_modulus(s)*one_over_mass_times_restlength*one_over_cfl_number_squared;
        T damping_hertz=2*damping(s)*one_over_mass_times_restlength*one_over_cfl_number;
        for(int k=1;k<=4;k++){frequency(nodes[k]).elastic_squared+=elastic_hertz_squared;frequency(nodes[k]).damping+=damping_hertz;}}
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
CFL_Strain_Rate() const
{
    ARRAY_VIEW<const TV> V(particles.V);
    T max_strain_rate=0;
    if(use_rest_state_for_strain_rate) for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const VECTOR<int,4>& nodes=spring_particles(s);
        TV dx;VECTOR<T,2> weights;T current_length;T current_edge_length;
        Axial_Vector(nodes,current_length,dx,weights,current_edge_length); // dx is normalized
        T strain_rate=TV::Dot_Product((1-weights.x)*V(nodes[1])+weights.x*V(nodes[2])-(1-weights.y)*V(nodes[3])-weights.y*V(nodes[4]),dx)/restlength(s);
        max_strain_rate=max(max_strain_rate,abs(strain_rate));}
    else for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const VECTOR<int,4>& nodes=spring_particles(s);
        TV dx;VECTOR<T,2> weights;T current_length;T current_edge_length;
        Axial_Vector(nodes,current_length,dx,weights,current_edge_length); // dx is normalized
        T strain_rate=TV::Dot_Product((1-weights.x)*V(nodes[1])+weights.x*V(nodes[2])-(1-weights.y)*V(nodes[3])-weights.y*V(nodes[4]),dx)/current_length;
        max_strain_rate=max(max_strain_rate,abs(strain_rate));}
    return Robust_Divide(max_strain_per_time_step,max_strain_rate);
}
//#####################################################################
// Function Set_Stiffness_Based_On_Reduced_Mass
//#####################################################################
// assumes mass and restlength are already defined
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Set_Stiffness_Based_On_Reduced_Mass(const T scaling_coefficient)
{
    for(int s=1;s<=spring_particles.m;s++){
        T harmonic_mass=Pseudo_Divide((T)4,ARRAYS_COMPUTATIONS::Sum(particles.one_over_effective_mass.Subset(spring_particles(s))));
        youngs_modulus(s)=scaling_coefficient*harmonic_mass/restlength(s);}
}
//#####################################################################
// Function Axial_Vector
//#####################################################################
// direction points from cross edge to common edge
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Axial_Vector(const VECTOR<int,4>& nodes,T& axial_length,TV& axial_direction,VECTOR<T,2>& weights,T& attached_edge_length) const
{
    axial_direction=SEGMENT_3D<T>(particles.X(nodes[1]),particles.X(nodes[2])).Shortest_Vector_Between_Segments(SEGMENT_3D<T>(particles.X(nodes[3]),particles.X(nodes[4])),weights);
    axial_length=axial_direction.Normalize();
    // TODO: How do we choose the sign of axial_direction?
    if(!axial_length) axial_direction=TV::Cross_Product(particles.X(nodes[2])-particles.X(nodes[1]),particles.X(nodes[4])-particles.X(nodes[3])).Normalized();
    attached_edge_length=(particles.X(nodes[1])-particles.X(nodes[3])).Magnitude()+
        (particles.X(nodes[1])-particles.X(nodes[4])).Magnitude()+
        (particles.X(nodes[2])-particles.X(nodes[3])).Magnitude()+
        (particles.X(nodes[2])-particles.X(nodes[4])).Magnitude();
}
//#####################################################################
// Function Prepare_Energy_Correction_Force
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Prepare_Energy_Correction_Force()
{
    energy_correction_forces.Resize(restlength.m);
    energy_correction_sign.Resize(restlength.m);

    // figure out damping directions
    if(verbose) for(int i=1;i<=particles.array_collection->Size();i++) {std::stringstream ss;ss<<"Initial body "<<i<<" velocity: "<<particles.V(i)<<std::endl;LOG::filecout(ss.str());}

    // velocities: v^n, v^{n+1} (unmodified), v^{n+1} (current)
    // at the beginning of jacobi step: v^{n+1} (current) = v^{n+1} (unmodified) + accumulated_impulses
    // during step: work is computed according to v^{n+1} (current) and v^n, along with previous_impulses and accumulated_impulses

    if(!allow_kd_direction_flip)
        for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
            TV new_V1=Endpoint_Velocity(s,1),new_V2=Endpoint_Velocity(s,2);
            if(TV::Dot_Product(new_V2-new_V1,optimization_direction(s))>0)
                energy_correction_sign(s)=-1;
            else
                energy_correction_sign(s)=1;}
}
//#####################################################################
// Function Compute_Energy_Correction_Force
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Compute_Energy_Correction_Force(ARRAY_VIEW<const TV> velocity_save,const int max_particle_degree,const T time,const T dt)
{
    T one_over_dt=1/dt;

    ARRAYS_COMPUTATIONS::Fill(energy_correction_forces,(T)0);

    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        TV old_V1=Endpoint_Velocity(velocity_save,s,1);
        TV old_V2=Endpoint_Velocity(velocity_save,s,2);
        TV new_V1=Endpoint_Velocity(s,1), new_V2=Endpoint_Velocity(s,2);
        
        TV V=new_V1-new_V2+old_V1-old_V2;
        T k=dt*Effective_Impulse_Factor(s);
        
        T current_potential_energy=Potential_Energy(s,time);
        T current_kinetic_energy=Endpoint_Kinetic_Energy(s);
        T delta_pe=current_potential_energy-potential_energy_save(s)+extra_energy(s);
        
        if(verbose) {std::stringstream ss;ss<<"AXIAL BENDING SPRING "<<s<<" ****************"<<std::endl;LOG::filecout(ss.str());}
        T current_work=Compute_Spring_Work(s,velocity_save,time,dt);
        if(verbose) {std::stringstream ss;ss<<"Direction "<<optimization_direction(s)<<" Total work done: "<<current_work<<std::endl;LOG::filecout(ss.str());}
        
        if(current_work+delta_pe<0)
            continue;
        
        if(verbose) {std::stringstream ss;ss<<"Saved potential energy: "<<potential_energy_save(s)<<" current potential energy: "<<current_potential_energy<<" delta pe: "<<delta_pe<<std::endl;LOG::filecout(ss.str());}
        if(verbose) {std::stringstream ss;ss<<"Current total energy: "<<current_potential_energy+current_kinetic_energy<<std::endl;LOG::filecout(ss.str());}
        
        T root1,root2;
        int roots=0;
        if(use_kinetic_energy_fix){
            T a_ke=k/2*dt;
            T b_ke=TV::Dot_Product(new_V1-new_V2,optimization_direction(s))*dt;
            T c_ke=(current_work+delta_pe)*relaxation_fraction;
            QUADRATIC<T> quadratic(a_ke,b_ke,c_ke);
            quadratic.Compute_Roots();
            roots=quadratic.roots;
            root1=quadratic.root1;
            root2=quadratic.root2;}
        else{
            TV previously_applied_force=previously_applied_forces(s)+force_correction(s)*optimization_direction(s);
            T a=k;
            T b=TV::Dot_Product(V+previously_applied_force*k,optimization_direction(s));
            T c=2*(current_work+delta_pe)*one_over_dt*relaxation_fraction;
            QUADRATIC<T> quadratic(a,b,c);
            quadratic.Compute_Roots();
            roots=quadratic.roots;
            root1=quadratic.root1;
            root2=quadratic.root2;}
        T force_magnitude;
        if(verbose) {std::stringstream ss;ss<<"ROOTS "<<roots<<std::endl;LOG::filecout(ss.str());}
        if(roots==2){
            if(abs(root1)<abs(root2))
                force_magnitude=root1;
            else
                force_magnitude=root2;
        }
        else if(roots==1)
            force_magnitude=root1;
        else{ // do the best job you can
            // it's always possible to add relative velocity, so that will never have no solution; this is only possible if we're trying to remove it
            assert(current_work+delta_pe>0);
            // maximize energy loss, i.e. kill relative velocity
            force_magnitude=TV::Dot_Product(new_V2-new_V1,optimization_direction(s))/k;}
        energy_correction_forces(s)=force_magnitude/max_particle_degree;
        if(verbose) {std::stringstream ss;ss<<"Force I want to apply: "<<energy_correction_forces(s)*optimization_direction(s)<<std::endl;LOG::filecout(ss.str());}
        if(!allow_kd_direction_flip){
            if(energy_correction_sign(s)==1){
                if(force_correction(s)+energy_correction_forces(s)<0)
                    energy_correction_forces(s)=-force_correction(s);}
            else{
                if(force_correction(s)+energy_correction_forces(s)>0)
                    energy_correction_forces(s)=-force_correction(s);}}
        if(verbose) {std::stringstream ss;ss<<"Correction force: "<<energy_correction_forces(s)*optimization_direction(s)<<std::endl;LOG::filecout(ss.str());}
        
        if(use_gauss_seidel_in_energy_correction) Apply_Energy_Correction_Impulse(s,energy_correction_forces(s),dt);}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Potential_Energy(int s,const T time) const
{
    const VECTOR<int,4>& nodes=spring_particles(s);
    VECTOR<T,2> weights;
    TV axial_direction=SEGMENT_3D<T>(particles.X(nodes[1]),particles.X(nodes[2])).Shortest_Vector_Between_Segments(SEGMENT_3D<T>(particles.X(nodes[3]),particles.X(nodes[4])),weights);
    T axial_length=axial_direction.Normalize();
    return (T).5*youngs_modulus(s)/restlength(s)*sqr(axial_length-visual_restlength(s));
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        potential_energy+=Potential_Energy(s,time);}
    return potential_energy;
}
//#####################################################################
// Function Apply_Energy_Correction
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Apply_Energy_Correction(const T time,const T dt)
{        
    if(!use_gauss_seidel_in_energy_correction)
        for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
            Apply_Energy_Correction_Impulse(s,energy_correction_forces(s),dt);}
    
    if(verbose) for(int i=1;i<=particles.array_collection->Size();i++) {std::stringstream ss;ss<<"Body "<<i<<" velocity: "<<particles.V(i)<<std::endl;LOG::filecout(ss.str());}

    if(verbose) 
        for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
            std::stringstream ss;ss<<"Force correction for spring "<<s<<": "<<force_correction(s)<<std::endl;
            ss<<"Final total energy: "<<Potential_Energy(s,time)+Endpoint_Kinetic_Energy(s)<<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
// Function Add_Connectivity
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Add_Connectivity(ARRAY<int>& particle_degree)
{
    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        for(int i=1;i<=4;i++) particle_degree(spring_particles(s)[i])++;}
}
//#####################################################################
// Function Endpoint_Velocity
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Endpoint_Mass(int s,int b) const
{
    switch(b){
        case 1:
            return particles.mass(spring_particles(s)(1))*optimization_weights(s)(1)+particles.mass(spring_particles(s)(2))*optimization_weights(s)(2);
        case 2:
            return particles.mass(spring_particles(s)(3))*optimization_weights(s)(3)+particles.mass(spring_particles(s)(4))*optimization_weights(s)(4);
        default:
            PHYSBAM_FATAL_ERROR("Invalid endpoint");}
}
//#####################################################################
// Function Endpoint_Velocity
//#####################################################################
template<class T> VECTOR<T,3> AXIAL_BENDING_SPRINGS<T>::
Endpoint_Velocity(int s,int b) const
{
    return Endpoint_Velocity(particles.V,s,b);
}
//#####################################################################
// Function Endpoint_Velocity
//#####################################################################
template<class T> VECTOR<T,3> AXIAL_BENDING_SPRINGS<T>::
Endpoint_Velocity(ARRAY_VIEW<const TV> velocity,int s,int b) const
{
    switch(b){
        case 1:{
            int node1=spring_particles(s)(1),node2=spring_particles(s)(2);
            return velocity(node1)*optimization_weights(s)(1)+velocity(node2)*optimization_weights(s)(2);}
        case 2:{
            int node3=spring_particles(s)(3),node4=spring_particles(s)(4);
            return velocity(node3)*optimization_weights(s)(3)+velocity(node4)*optimization_weights(s)(4);}
        default:
            PHYSBAM_FATAL_ERROR("Invalid endpoint");}
}
//#####################################################################
// Function Endpoint_Kinetic_Energy
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Endpoint_Kinetic_Energy(int s,int b) const
{
    return Endpoint_Kinetic_Energy(particles.V,s,b);
}
//#####################################################################
// Function Endpoint_Kinetic_Energy
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Endpoint_Kinetic_Energy(ARRAY_VIEW<const TV> velocity,int s,int b) const
{
    return (T).5*Endpoint_Mass(s,b)*Endpoint_Velocity(velocity,s,b).Magnitude_Squared();
}
//#####################################################################
// Function Endpoint_Kinetic_Energy
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Endpoint_Kinetic_Energy(int s) const
{
    return Endpoint_Kinetic_Energy(s,1)+Endpoint_Kinetic_Energy(s,2);
}
//#####################################################################
// Function Effective_Impulse_Factor
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Effective_Impulse_Factor(int s) const
{
    T mass1=Endpoint_Mass(s,1);
    T mass2=Endpoint_Mass(s,2);
    T one_over_denom=1/(mass1*mass2);
    return (mass1+mass2)*one_over_denom;
}
//#####################################################################
// Function Store_Potential_Energy
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Save_Potential_Energy(const T time)
{
    potential_energy_save.Resize(restlength.m);
    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        potential_energy_save(s)=Potential_Energy(s,time);}
}
//#####################################################################
// Function Compute_Spring_Work
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Compute_Spring_Work(int s,ARRAY_VIEW<const TV> velocity_save,const T time,const T dt) const
{
    TV new_V1=Endpoint_Velocity(s,1), new_V2=Endpoint_Velocity(s,2);
    TV force=previously_applied_forces(s);
    force+=force_correction(s)*optimization_direction(s);
    TV old_V1=Endpoint_Velocity(velocity_save,s,1);
    TV old_V2=Endpoint_Velocity(velocity_save,s,2);
    TV average_vel1=(new_V1+old_V1)*(T).5;
    TV average_vel2=(new_V2+old_V2)*(T).5;
    TV distance1=average_vel1*dt;
    TV distance2=average_vel2*dt;
    T first_node_work=-TV::Dot_Product(distance1,force);
    T second_node_work=TV::Dot_Product(distance2,force);
    if(verbose) {std::stringstream ss;ss<<"Work computation: force: "<<force<<" average v1: "<<average_vel1<<" average v2: "<<average_vel2<<std::endl;}
    return first_node_work+second_node_work;
}
//#####################################################################
// Function Compute_Previously_Applied_Forces
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Compute_Previously_Applied_Forces()
{
    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        previously_applied_forces(s)=optimization_coefficient(s)*TV::Dot_Product(Endpoint_Velocity(s,2)-Endpoint_Velocity(s,1),optimization_direction(s))*optimization_direction(s);
        previously_applied_forces(s)+=youngs_modulus(s)/restlength(s)*(optimization_current_length(s)-visual_restlength(s))*optimization_direction(s);
        force_correction(s)=0;}
}
//#####################################################################
// Function Apply_Energy_Correction_impulse
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Apply_Energy_Correction_Impulse(const int s,const T force,const T dt)
{
    force_correction(s)+=force;
    TV impulse=dt*force*optimization_direction(s);
    //LOG::cout<<"Force correction for spring "<<s<<": "<<energy_correction_forces(s)<<std::endl;
    
    int node1,node2,node3,node4;spring_particles(s).Get(node1,node2,node3,node4);
    T w1,w2,w3,w4;optimization_weights(s).Get(w1,w2,w3,w4);
    particles.V(node1)-=(impulse*particles.one_over_mass(node1)*w1);
    particles.V(node2)-=(impulse*particles.one_over_mass(node2)*w2);
    particles.V(node3)+=(impulse*particles.one_over_mass(node3)*w3);
    particles.V(node4)+=(impulse*particles.one_over_mass(node4)*w4);
    /*LOG::cout<<"Final total energy: "<<Potential_Energy(s,time)+body1.Kinetic_Energy()+body2.Kinetic_Energy()<<std::endl;*/
}
//#####################################################################
// Function Compute_Energy_Error
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Compute_Energy_Error(ARRAY_VIEW<const TV> velocity_save,const T time,const T dt)
{
    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        T current_potential_energy=Potential_Energy(s,time);
        T spring_extra_energy=(current_potential_energy-potential_energy_save(s))+extra_energy(s)+Compute_Spring_Work(s,velocity_save,time,dt);
        extra_energy(s)=max((T)0,spring_extra_energy);}
}
//#####################################################################
// Function Setup_Set_Velocity_From_Positions
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Setup_Set_Velocity_From_Positions(const T time,const bool is_position_update,const bool reset_alphas)
{
    if(!is_position_update) Update_Position_Based_State(time,is_position_update);

    force_estimates.Resize(spring_particles.m);
    incident_nodes.Resize(spring_particles.m);
    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const VECTOR<int,4>& nodes=spring_particles(s);
        incident_nodes(s).Remove_All();
        incident_nodes(s).Append(nodes(1));incident_nodes(s).Append(nodes(2));
        incident_nodes(s).Append(nodes(3));incident_nodes(s).Append(nodes(4));}
    if(is_position_update){
        if(reset_alphas) for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();force_estimates(s)=0;}
        Save_Potential_Energy(time);}
}
//#####################################################################
// Function Get_Element_Count
//#####################################################################
template<class T> int AXIAL_BENDING_SPRINGS<T>::
Get_Element_Count()
{
    return spring_particles.m;
}
//#####################################################################
// Function Get_Force_Elements
//#####################################################################
template<class T> FORCE_ELEMENTS* AXIAL_BENDING_SPRINGS<T>::
Get_Force_Elements()
{
    return &force_springs;
}
//#####################################################################
// Function Incident_Nodes
//#####################################################################
template<class T> ARRAY<int>* AXIAL_BENDING_SPRINGS<T>::
Incident_Nodes(const int force_element)
{
    return &(incident_nodes(force_element));
}
//#####################################################################
// Function Get_Direction
//#####################################################################
template<class T> VECTOR<T,3> AXIAL_BENDING_SPRINGS<T>::
Get_Direction(const int force_element)
{
    return optimization_direction(force_element);
}
//#####################################################################
// Function Get_Combined_One_Over_Mass
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Get_Combined_One_Over_Mass(const int force_element)
{
    int node1,node2,node3,node4;spring_particles(force_element).Get(node1,node2,node3,node4);
    T w1,w2,w3,w4;optimization_weights(force_element).Get(w1,w2,w3,w4);
    return (sqr(w1)*particles.one_over_mass(node1)+sqr(w2)*particles.one_over_mass(node2))+
        (sqr(w3)*particles.one_over_mass(node3)+sqr(w4)*particles.one_over_mass(node4));
}
//#####################################################################
// Function Incident_Force_Elements
//#####################################################################
template<class T> ARRAY<int>* AXIAL_BENDING_SPRINGS<T>::
Incident_Force_Elements(const int particle)
{
    return &(incident_elements)(particle);
}
//#####################################################################
// Function Choose_Solution
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Choose_Solution(const bool use_orig_force,const int force_element,const T dt,const T alpha1,const T alpha2,ARRAY<T>& v_n_hats)
{
    if(use_orig_force){
        T force=youngs_modulus(force_element)/restlength(force_element)*(optimization_current_length(force_element)-visual_restlength(force_element));
        
        if(fabs(force-alpha1) < fabs(force-alpha2)) Set_Force(force_element,alpha1);
        else Set_Force(force_element,alpha2);}
    else{
        int node1,node2,node3,node4;spring_particles(force_element).Get(node1,node2,node3,node4);
        T w1,w2,w3,w4;optimization_weights(force_element).Get(w1,w2,w3,w4);
        T v1_n_hat=w3*v_n_hats(3)+w4*v_n_hats(4);
        T v2_n_hat=w1*v_n_hats(1)+w2*v_n_hats(2);

        T v1_one_over_mass,v1_mass;
        if(particles.mass(node3)==FLT_MAX && particles.mass(node4)==FLT_MAX){
            v1_one_over_mass=0;
            v1_mass=FLT_MAX;}
        else{
            v1_one_over_mass=(sqr(w3)*particles.one_over_mass(node3)+sqr(w4)*particles.one_over_mass(node4));
            v1_mass = 1/v1_one_over_mass;}
        T v2_one_over_mass,v2_mass;
        if(particles.mass(node1)==FLT_MAX && particles.mass(node2)==FLT_MAX){
            v2_one_over_mass=0;
            v2_mass=FLT_MAX;}
        else{
            v2_one_over_mass=(sqr(w1)*particles.one_over_mass(node1)+sqr(w2)*particles.one_over_mass(node2));
            v2_mass = 1/v2_one_over_mass;}

        T v1_alpha1_temp=v1_n_hat+dt*v1_one_over_mass*alpha1;
        T v2_alpha1_temp=v2_n_hat-dt*v2_one_over_mass*alpha1;

        T v1_alpha2_temp=v1_n_hat+dt*v1_one_over_mass*alpha2;
        T v2_alpha2_temp=v2_n_hat-dt*v2_one_over_mass*alpha2;
        
        T combined_mass=v1_mass+v2_mass;
        TV x1=w3*particles.X(node3)+w4*particles.X(node4);
        TV x2=w1*particles.X(node1)+w2*particles.X(node2);
        T x1u=TV::Dot_Product(x1,optimization_direction(force_element));
        T x2u=TV::Dot_Product(x2,optimization_direction(force_element));
        
        T w=sqrt((youngs_modulus(force_element)/restlength(force_element))*Get_Combined_One_Over_Mass(force_element));
        T v_cm_0=(v1_mass*v1_n_hat+v2_mass*v2_n_hat)/combined_mass;
        T c1=(x2u-x1u)-restlength(force_element);PHYSBAM_FATAL_ERROR(); // this should probably be visual restlength
        T c2=(v2_n_hat-v1_n_hat)/w;
        T v1_analytic=v_cm_0-(v2_mass/combined_mass)*(-w*c1*sin(w*dt)+w*c2*cos(w*dt));
        T v2_analytic=v_cm_0+(v1_mass/combined_mass)*(-w*c1*sin(w*dt)+w*c2*cos(w*dt));
        if((fabs(v1_alpha1_temp-v1_analytic)+fabs(v2_alpha1_temp-v2_analytic)) < (fabs(v1_alpha2_temp-v1_analytic)+fabs(v2_alpha2_temp-v2_analytic))) Set_Force(force_element,alpha1);
        else Set_Force(force_element,alpha2);}
}
//#####################################################################
// Function Get_Force
//#####################################################################
template<class T> VECTOR<T,3> AXIAL_BENDING_SPRINGS<T>::
Get_Force(const int force_element,const int particle,const bool use_original_force)
{
    int node1,node2,node3,node4;spring_particles(force_element).Get(node1,node2,node3,node4);
    T w1,w2,w3,w4;optimization_weights(force_element).Get(w1,w2,w3,w4);
    TV force;
    if(use_original_force) force=youngs_modulus(force_element)/restlength(force_element)*(optimization_current_length(force_element)-visual_restlength(force_element))*optimization_direction(force_element);
    else force=force_estimates(force_element)*optimization_direction(force_element);

    if(particle==node1) return -w1*force;
    else if(particle==node2) return -w2*force;
    else if(particle==node3) return w3*force;
    else if(particle==node4) return w4*force;
    else PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Get_Force
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Get_Force(const int force_element)
{
    return force_estimates(force_element);
}
//#####################################################################
// Function Set_Force
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Set_Force(const int force_element,const T force)
{
    force_estimates(force_element)=force;
}
//#####################################################################
// Function Get_Damping_Force
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Get_Damping_Force(const int particle,TV& damping_force,const T dt,const bool use_coefficient)
{
    ARRAY<int>& node_incident_forces=(incident_elements)(particle);
    for(int s=1;s<=node_incident_forces.m;s++){
        int force_element=node_incident_forces(s);
        int node1,node2,node3,node4;spring_particles(force_element).Get(node1,node2,node3,node4);
        T w1,w2,w3,w4;optimization_weights(force_element).Get(w1,w2,w3,w4);
        TV force=TV::Dot_Product(w1*particles.V(node1)+w2*particles.V(node2)-w3*particles.V(node3)-w4*particles.V(node4),optimization_direction(force_element))*optimization_direction(force_element);
        if(particle==node1) force*=-w1;
        else if(particle==node2) force*=-w2;
        else if(particle==node3) force*=w3;
        else if(particle==node4) force*=w4;
        if(use_coefficient) damping_force+=dt*particles.one_over_mass(particle)*optimization_coefficient(force_element)*force;
        else damping_force+=dt*particles.one_over_mass(particle)*force;}
}
//#####################################################################
// Function Update_Residual_Energy
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Update_Residual_Energy(const int force_element,const T residual_energy,const T time)
{
    residual_PE(force_element)=residual_energy;
}
//#####################################################################
// Function Get_Residual_Energy
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Get_Residual_Energy(const int force_element)
{
    return residual_PE(force_element);
}
//#####################################################################
// Function Compute_Quadratic_Contribution_For_Node
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Compute_Quadratic_Contribution_For_Force(T& A,T& a,T&c,const T dt,const int force_element,const T combined_one_over_mass,const bool ignore_PE_terms)
{
    if(!ignore_PE_terms){
        int node1,node2,node3,node4;spring_particles(force_element).Get(node1,node2,node3,node4);
        T w1,w2,w3,w4;optimization_weights(force_element).Get(w1,w2,w3,w4);

        T x1u=TV::Dot_Product(particles.X(node1),optimization_direction(force_element));
        T x2u=TV::Dot_Product(particles.X(node2),optimization_direction(force_element));
        T x3u=TV::Dot_Product(particles.X(node3),optimization_direction(force_element));
        T x4u=TV::Dot_Product(particles.X(node4),optimization_direction(force_element));

        a=(w1*x1u+w2*x2u)-(w3*x3u+w4*x4u);
        c=-(T).5*dt*dt*combined_one_over_mass;

        A+=(T).5*(youngs_modulus(force_element)/restlength(force_element))*c*c;}

    A+=(T).5*dt*dt*combined_one_over_mass;
}
//#####################################################################
// Function Compute_Quadratic_Contribution_For_Node
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Compute_Quadratic_Contribution_For_Node(T& B,T& C,T&b,const T dt,const int node,const int force_element,const T combined_one_over_mass,const T v_n_correction,const bool ignore_PE_terms)
{
    int node1,node2,node3,node4;spring_particles(force_element).Get(node1,node2,node3,node4);
    T w1,w2,w3,w4;optimization_weights(force_element).Get(w1,w2,w3,w4);

    T vu=TV::Dot_Product(particles.V(node),optimization_direction(force_element));

    if(node==node1){
        T term=dt*w1*(vu+(T).5*v_n_correction);
        if(!ignore_PE_terms) b+=term;
        B+=-term;}
    if(node==node2){
        T term=dt*w2*(vu+(T).5*v_n_correction);
        if(!ignore_PE_terms) b+=term;
        B+=-term;}
    if(node==node3){
        T term=dt*w3*(vu+(T).5*v_n_correction);
        if(!ignore_PE_terms) b+=-term;
        B+=term;}
    if(node==node4){
        T term=dt*w4*(vu+(T).5*v_n_correction);
        if(!ignore_PE_terms) b+=-term;
        B+=term;}
}
//#####################################################################
// Function Compute_Quadratic_Contribution_For_Residual
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Compute_Quadratic_Contribution_For_Residual(T& B,T& C,T& a,T&b,T&c,const T dt,const T time,const int force_element,const bool ignore_PE_terms)
{
    if(!ignore_PE_terms){
        B+=(youngs_modulus(force_element)/restlength(force_element))*(a+b-visual_restlength(force_element))*c;
        C=(youngs_modulus(force_element)/restlength(force_element))*(a+b/2-visual_restlength(force_element))*b;}
    else C=delta_PE(force_element);

    C+=residual_PE(force_element);
}
//#####################################################################
// Function Store_Delta_PE
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Store_Delta_PE(const T time)
{
    delta_PE.Resize(potential_energy_save.m);
    T total_current_PE=0;
    for(SPRING_ITERATOR iterator(force_springs);iterator.Valid();iterator.Next()){int s=iterator.Data();
        delta_PE(s)=Potential_Energy(s,time)-potential_energy_save(s);
        total_current_PE+=Potential_Energy(s,time);}
    total_delta_PE=total_current_PE-ARRAYS_COMPUTATIONS::Sum(potential_energy_save);
}
//#####################################################################
// Function Get_Total_Delta_PE
//#####################################################################
template<class T> T AXIAL_BENDING_SPRINGS<T>::
Get_Total_Delta_PE()
{
    return total_delta_PE;
}
//#####################################################################
// Function Save_And_Reset_Elastic_Coefficients
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Save_And_Reset_Elastic_Coefficient()
{
    saved_youngs_modulus.Resize(youngs_modulus.m);
    for(int i=1;i<=saved_youngs_modulus.m;i++){
        saved_youngs_modulus(i)=youngs_modulus(i);
        youngs_modulus(i)=0;}
}
//#####################################################################
// Function Restore_Elastic_Coefficients
//#####################################################################
template<class T> void AXIAL_BENDING_SPRINGS<T>::
Restore_Elastic_Coefficient()
{
    for(int i=1;i<=saved_youngs_modulus.m;i++)
        youngs_modulus(i)=saved_youngs_modulus(i);
}
//#####################################################################
// Function Create_Axial_Bending_Springs
//#####################################################################
template<class T> AXIAL_BENDING_SPRINGS<T>* PhysBAM::
Create_Axial_Bending_Springs(PARTICLES<VECTOR<T,3> >& particles,TRIANGLE_MESH& triangle_mesh,const T clamped_restlength,const T stiffness,const T overdamping_fraction,
    const bool limit_time_step_by_strain_rate,const T max_strain_per_time_step,const bool use_rest_state_for_strain_rate,const bool verbose)
{
    AXIAL_BENDING_SPRINGS<T>* axial=new AXIAL_BENDING_SPRINGS<T>(particles,triangle_mesh);
    axial->verbose=verbose;
    axial->Set_Restlength_From_Particles();
    axial->Clamp_Restlength(clamped_restlength);
    axial->Set_Stiffness(stiffness);
    axial->Set_Overdamping_Fraction(overdamping_fraction);
    axial->Limit_Time_Step_By_Strain_Rate(limit_time_step_by_strain_rate,max_strain_per_time_step);
    axial->Use_Rest_State_For_Strain_Rate(use_rest_state_for_strain_rate);
    return axial;
}
//#####################################################################
// Function Create_Axial_Bending_Springs
//#####################################################################
template<class T> AXIAL_BENDING_SPRINGS<T>* PhysBAM::
Create_Axial_Bending_Springs(TRIANGULATED_SURFACE<T>& triangulated_surface,
    const T clamped_restlength,const T stiffness,const T overdamping_fraction,const bool limit_time_step_by_strain_rate,const T max_strain_per_time_step,
    const bool use_rest_state_for_strain_rate,const bool verbose)
{
    return Create_Axial_Bending_Springs(dynamic_cast<PARTICLES<VECTOR<T,3> >&>(triangulated_surface.particles),triangulated_surface.mesh,clamped_restlength,stiffness,overdamping_fraction,limit_time_step_by_strain_rate,
        max_strain_per_time_step,use_rest_state_for_strain_rate,verbose);
}
//#####################################################################
template class AXIAL_BENDING_SPRINGS<float>;
template AXIAL_BENDING_SPRINGS<float>* PhysBAM::Create_Axial_Bending_Springs<float>(TRIANGULATED_SURFACE<float>&,float,float,float,bool,float,bool,bool);
template AXIAL_BENDING_SPRINGS<float>* PhysBAM::Create_Axial_Bending_Springs<float>(PARTICLES<VECTOR<float,3> >&,TRIANGLE_MESH&,float,float,float,bool,float,bool,bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class AXIAL_BENDING_SPRINGS<double>;
template AXIAL_BENDING_SPRINGS<double>* PhysBAM::Create_Axial_Bending_Springs<double>(TRIANGULATED_SURFACE<double>&,double,double,double,bool,double,bool,bool);
template AXIAL_BENDING_SPRINGS<double>* PhysBAM::Create_Axial_Bending_Springs<double>(PARTICLES<VECTOR<double,3> >&,TRIANGLE_MESH&,double,double,double,bool,double,bool,bool);
#endif
