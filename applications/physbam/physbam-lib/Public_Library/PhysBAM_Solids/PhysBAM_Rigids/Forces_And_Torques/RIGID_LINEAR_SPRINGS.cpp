//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_LINEAR_SPRINGS
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGID_LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <cfloat>
using ::std::sqrt;
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_LINEAR_SPRINGS<TV>::
RIGID_LINEAR_SPRINGS(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
    :BASE(rigid_body_collection_input),use_kinetic_energy_fix(true),relaxation_fraction(1),use_gauss_seidel_in_energy_correction(false),allow_kd_direction_flip(false)
{
    Invalidate_CFL();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_LINEAR_SPRINGS<TV>::
~RIGID_LINEAR_SPRINGS()
{
}
//#####################################################################
// Function Set_Restlength_From_Material_Coordinates
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Set_Restlengths()
{
    restlength.Resize(attachment_radius.m);
    Invalidate_CFL();
    for(int i=1;i<=segment_mesh.elements.m;i++) restlength(i)=Spring_Length(i);
    visual_restlength=restlength;
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    segment_mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated)
{
    force_segments.Update(segment_mesh.elements,particle_is_simulated);
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Set_Overdamping_Fraction(int b,const T overdamping_fraction) // 1 is critically damped
{
    TV X1=Attachment_Location(b,1),X2=Attachment_Location(b,2);
    TV direction=(X2-X1).Normalized();
    RIGID_BODY<TV>& b1=Body(b,1);
    RIGID_BODY<TV>& b2=Body(b,2);

    T harmonic_mass=TV::Dot_Product((b1.Impulse_Factor(X1)+b2.Impulse_Factor(X2))*direction,direction);
    Set_Damping(b,overdamping_fraction*2*sqrt(youngs_modulus(b)*restlength(b)/harmonic_mass));
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Update_Position_Based_State(const T time)
{
    states.Resize(segment_mesh.elements.m);
    current_lengths.Resize(segment_mesh.elements.m);
    extra_energy.Resize(youngs_modulus.m);
    force_correction.Resize(youngs_modulus.m);
    previously_applied_forces.Resize(youngs_modulus.m);
    
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        TV X1=Attachment_Location(s,1),X2=Attachment_Location(s,2);
        STATE& state=states(s);
        state.nodes=segment_mesh.elements(s);
        state.direction=X2-X1;
        current_lengths(s)=state.direction.Normalize();
        state.coefficient=damping(s)/restlength(s);
        state.r(1)=X1-Body(s,1).X();
        state.r(2)=X2-Body(s,2).X();}
}
//#####################################################################
// Function Add_Force
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Add_Force(ARRAY_VIEW<TWIST<TV> > rigid_F,const STATE& state,const TV& force) const
{
    if(!rigid_body_collection.Rigid_Body(state.nodes(1)).Has_Infinite_Inertia()){
        rigid_F(state.nodes(1)).linear+=force;
        rigid_F(state.nodes(1)).angular+=TV::Cross_Product(state.r(1),force);}
    if(!rigid_body_collection.Rigid_Body(state.nodes(2)).Has_Infinite_Inertia()){
        rigid_F(state.nodes(2)).linear-=force;
        rigid_F(state.nodes(2)).angular-=TV::Cross_Product(state.r(2),force);}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        TV force=youngs_modulus(s)/restlength(s)*(current_lengths(s)-visual_restlength(s))*state.direction;
        Add_Force(rigid_F,state,force);}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    //LOG::Time("spring add velocity dependent");
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        VECTOR<TV,2> V;
        for(int i=1;i<=2;i++){
            const TWIST<TV>& twist=rigid_V(segment_mesh.elements(s)(i));
            V(i)=twist.linear+TV::Cross_Product(twist.angular,state.r(i));}
        TV force=(state.coefficient*TV::Dot_Product(V(2)-V(1),state.direction))*state.direction;
        Add_Force(rigid_F,state,force);}
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        VECTOR<TV,2> V;
        for(int i=1;i<=2;i++){
            const TWIST<TV>& twist=rigid_V(segment_mesh.elements(s)(i));
            V(i)=twist.linear+TV::Cross_Product(twist.angular,state.r(i));}
        TV dl=V(2)-V(1),dl_projected=dl.Projected_On_Unit_Direction(state.direction);
        TV force=youngs_modulus(s)/restlength(s)*dl_projected;
        Add_Force(rigid_F,state,force);}
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{
//     T one_over_cfl_number=1/cfl_number,one_over_cfl_number_squared=sqr(one_over_cfl_number);
//     for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
//         const VECTOR<int,2>& nodes=segment_mesh.elements(s);
//         T ym=youngs_modulus(s);
//         T d=damping(s);
//         for(int k=1;k<=2;k++){
//             frequency(nodes[k]).elastic_squared+=particles.one_over_effective_mass(nodes[k])/restlength(s)*4*ym*one_over_cfl_number_squared;
//             frequency(nodes[k]).damping+=particles.one_over_effective_mass(nodes[k])/restlength(s)*2*d*one_over_cfl_number;}}
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
CFL_Strain_Rate() const
{
//     T max_strain_rate=0,strain_rate;TV dx;
//     if(use_rest_state_for_strain_rate) for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
//         int i,j;segment_mesh.elements(s).Get(i,j);
//         dx=particles.X(j)-particles.X(i);T magnitude=dx.Magnitude();if(magnitude!=0) dx.Normalize();
//         strain_rate=TV::Dot_Product(particles.V(j)-particles.V(i),dx)/restlength(s);
//         max_strain_rate=max(max_strain_rate,abs(strain_rate));}
//     else for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
//         int i,j;segment_mesh.elements(s).Get(i,j);
//         dx=particles.X(j)-particles.X(i);
//         strain_rate=TV::Dot_Product(particles.V(j)-particles.V(i),dx)/TV::Dot_Product(dx,dx);
//         max_strain_rate=max(max_strain_rate,abs(strain_rate));}
//     return Robust_Divide(max_strain_per_time_step,max_strain_rate);
    return FLT_MAX;
}
//#####################################################################
// Function Average_Restlength
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
Average_Restlength() const
{
    return ARRAYS_COMPUTATIONS::Average(restlength);
}
//#####################################################################
// Function Print_Restlength_Statistics
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Print_Restlength_Statistics() const
{
    LOG::SCOPE scope("linear spring statistics","linear spring statistics");
    LOG::Stat("count",restlength.m);
    ARRAY<T> length(restlength);Sort(length);
    ARRAY<T> visual_length(visual_restlength);Sort(visual_length);
    if(length.m){
        LOG::Stat("smallest restlength",length(1));LOG::Stat("smallest visual restlength",visual_length(1));
        LOG::Stat("one percent restlength",length((int)(.01*length.m)+1));LOG::Stat("one percent visual restlength",visual_length((int)(.01*length.m)+1));
        LOG::Stat("ten percent restlength",length((int)(.1*length.m)+1));LOG::Stat("ten percent visual restlength",visual_length((int)(.1*length.m)+1));
        LOG::Stat("median restlength",length((int)(.5*length.m)+1));LOG::Stat("median visual restlength",visual_length((int)(.5*length.m)+1));}
}
//#####################################################################
// Function Print_Deformation_Statistics
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Print_Deformation_Statistics() const
{
    LOG::SCOPE scope("linear spring deformation","linear spring deformation");
    ARRAY<T> deformation(segment_mesh.elements.m,false);
    for(int s=1;s<=segment_mesh.elements.m;s++){
        int i,j;segment_mesh.elements(s).Get(i,j);
        T length=Spring_Length(i),rl=visual_restlength(s);
        deformation(s)=rl?abs(length-rl)/rl:length==0?0:FLT_MAX;}
    Sort(deformation);
    LOG::Stat("maximum deformation",deformation(deformation.m));
    LOG::Stat("one percent deformation",deformation((int)(.99*deformation.m)+1));
    LOG::Stat("ten percent deformation",deformation((int)(.9*deformation.m)+1));
    LOG::Stat("twenty percent deformation",deformation((int)(.8*deformation.m)+1));
    LOG::Stat("median deformation",deformation((int)(.5*deformation.m)+1));
}
//#####################################################################
// Function Maximum_Compression_Or_Expansion_Fraction
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
Maximum_Compression_Or_Expansion_Fraction(int* index) const
{
    T max_compression=0;int max_index=1;
    for(int s=1;s<=segment_mesh.elements.m;s++){
        int i,j;segment_mesh.elements(s).Get(i,j);
        T length=Spring_Length(s);
        T rl=visual_restlength(s);
        T compression=(rl)?(abs(length-rl)/rl):((length==0)?0:FLT_MAX);
        if(compression>max_compression){max_compression=compression;max_index=s;}}
    if(index) (*index)=max_index;
    return max_compression;
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
Potential_Energy(int s,const T time) const
{
    return (T).5*youngs_modulus(s)/restlength(s)*sqr(Spring_Length(s)-visual_restlength(s));
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        potential_energy+=Potential_Energy(s,time);}
    return potential_energy;
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
Compute_Total_Energy(const T time) const
{
    T total_energy=0;
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        total_energy+=Potential_Energy(s,time);}
    for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
        total_energy+=rigid_body_collection.Rigid_Body(i).Kinetic_Energy();
    }
    total_energy+=ARRAYS_COMPUTATIONS::Sum(extra_energy);
    return total_energy;
}
//#####################################################################
// Function Store_Potential_Energy
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Save_Potential_Energy(const T time)
{
    potential_energy_save.Resize(youngs_modulus.m);
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        potential_energy_save(s)=Potential_Energy(s,time);}
}
//#####################################################################
// Function Compute_Spring_Work
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
Compute_Spring_Work(int s,ARRAY_VIEW<const TWIST<TV> > rigid_velocity_save,const T time,const T dt) const
{
    const STATE& state=states(s);
    TV new_V1=Endpoint_Velocity(s,1), new_V2=Endpoint_Velocity(s,2);
    TV force=previously_applied_forces(s);
    force+=force_correction(s)*state.direction;
    const RIGID_BODY<TV> &body1=Body(s,1), &body2=Body(s,2);
    TV old_V1=body1.Pointwise_Object_Velocity(rigid_velocity_save(state.nodes(1)),body1.X(),Attachment_Location(s,1));
    TV old_V2=body2.Pointwise_Object_Velocity(rigid_velocity_save(state.nodes(2)),body2.X(),Attachment_Location(s,2));
    TV average_vel1=(new_V1+old_V1)*(T).5;
    TV average_vel2=(new_V2+old_V2)*(T).5;
    TV distance1=average_vel1*dt;
    TV distance2=average_vel2*dt;
    T first_node_work=TV::Dot_Product(distance1,force);
    T second_node_work=-TV::Dot_Product(distance2,force);
    {std::stringstream ss;ss<<"Work computation: force: "<<force<<" average v1: "<<average_vel1<<" average v2: "<<average_vel2<<std::endl;LOG::filecout(ss.str());}
    return first_node_work+second_node_work;
}
//#####################################################################
// Function Compute_Energy_Error
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Compute_Energy_Error(ARRAY_VIEW<const TWIST<TV> > rigid_velocity_save,const T time,const T dt)
{
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        T current_potential_energy=Potential_Energy(s,time);
        T spring_extra_energy=(current_potential_energy-potential_energy_save(s))+extra_energy(s)+Compute_Spring_Work(s,rigid_velocity_save,time,dt);
        extra_energy(s)=max((T)0,spring_extra_energy);}
}
//#####################################################################
// Function Compute_Previously_Applied_Forces
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Compute_Previously_Applied_Forces()
{
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        previously_applied_forces(s)=state.coefficient*TV::Dot_Product(Endpoint_Velocity(s,2)-Endpoint_Velocity(s,1),state.direction)*state.direction;
        previously_applied_forces(s)+=youngs_modulus(s)/restlength(s)*(current_lengths(s)-visual_restlength(s))*state.direction;
        force_correction(s)=0;}
}
//#####################################################################
// Function Print_All_Energy_Debts
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Print_All_Energy_Debts(ARRAY_VIEW<const TWIST<TV> > rigid_velocity_save,const T time,const T dt) const
{
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){
        int s=iterator.Data();
        T current_potential_energy=Potential_Energy(s,time);
        T delta_pe=current_potential_energy-potential_energy_save(s)+extra_energy(s);    
        T current_work=Compute_Spring_Work(s,rigid_velocity_save,time,dt);
        {std::stringstream ss;ss << "spring " << s << " debt " << current_work+delta_pe << std::endl;LOG::filecout(ss.str());}}
}
//#####################################################################
// Function Apply_Energy_Correction_impulse
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Apply_Energy_Correction_Impulse(const int s,const T force,const T dt)
{
    force_correction(s)+=force;
    TV impulse=dt*force*states(s).direction;
    //{std::stringstream ss;ss<<"Force correction for spring "<<s<<": "<<energy_correction_forces(s)<<std::endl;LOG::filecout(ss.str());}
    
    RIGID_BODY<TV> &body1=Body(s,1), &body2=Body(s,2);
    body1.Apply_Impulse_To_Body(Attachment_Location(s,1),impulse,typename TV::SPIN());
    body2.Apply_Impulse_To_Body(Attachment_Location(s,2),-impulse,typename TV::SPIN());
    /*{std::stringstream ss;ss<<"Final total energy: "<<Potential_Energy(s,time)+body1.Kinetic_Energy()+body2.Kinetic_Energy()<<std::endl;LOG::filecout(ss.str());}*/
}
//#####################################################################
// Function Add_Energy_Correction_Force
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Add_Energy_Correction_Force(ARRAY_VIEW<const TWIST<TV> > rigid_velocity_save,const T time,const T dt)
{
    energy_correction_forces.Resize(youngs_modulus.m);

    T one_over_dt=1/dt;
    int iterations=20;

    // figure out damping directions
    ARRAY<int> energy_correction_sign(youngs_modulus.m);
    for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
        {std::stringstream ss;ss<<"Initial body "<<i<<" velocity: "<<rigid_body_collection.Rigid_Body(i).Twist()<<std::endl;LOG::filecout(ss.str());}
    }

    // velocities: v^n, v^{n+1} (unmodified), v^{n+1} (current)
    // at the beginning of jacobi step: v^{n+1} (current) = v^{n+1} (unmodified) + accumulated_impulses
    // during step: work is computed according to v^{n+1} (current) and v^n, along with previous_impulses and accumulated_impulses

    if(!allow_kd_direction_flip)
        for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
            const STATE& state=states(s);
            TV new_V1=Endpoint_Velocity(s,1),new_V2=Endpoint_Velocity(s,2);
            if(TV::Dot_Product(new_V2-new_V1,state.direction)>0)
                energy_correction_sign(s)=1;
            else
                energy_correction_sign(s)=-1;}

    for(int i=1;i<=iterations;i++){
        ARRAYS_COMPUTATIONS::Fill(energy_correction_forces,(T)0);

        for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
            const STATE& state=states(s);
            int first_node=state.nodes(1);
            int second_node=state.nodes(2);
            RIGID_BODY<TV> &body1=Body(s,1), &body2=Body(s,2);
            TV old_V1=body1.Pointwise_Object_Velocity(rigid_velocity_save(first_node),body1.X(),Attachment_Location(s,1));
            TV old_V2=body2.Pointwise_Object_Velocity(rigid_velocity_save(second_node),body2.X(),Attachment_Location(s,2));
            TV new_V1=Endpoint_Velocity(s,1), new_V2=Endpoint_Velocity(s,2);

            TV V=new_V1-new_V2+old_V1-old_V2;
            T k=dt*Effective_Impulse_Factor(s);

            T current_potential_energy=Potential_Energy(s,time);
            T current_kinetic_energy=body1.Kinetic_Energy()+body2.Kinetic_Energy();
            T delta_pe=current_potential_energy-potential_energy_save(s)+extra_energy(s);

            {std::stringstream ss;ss<<"SPRING "<<s<<" ****************"<<std::endl;LOG::filecout(ss.str());}
            T current_work=Compute_Spring_Work(s,rigid_velocity_save,time,dt);
            {std::stringstream ss;ss<<"Total work done: "<<current_work<<std::endl;LOG::filecout(ss.str());}

            if(current_work+delta_pe<0)
                continue;

            {std::stringstream ss;ss<<"Saved potential energy: "<<potential_energy_save(s)<<" current potential energy: "<<current_potential_energy<<" delta pe: "<<delta_pe<<std::endl;LOG::filecout(ss.str());}
            {std::stringstream ss;ss<<"Current total energy: "<<current_potential_energy+current_kinetic_energy<<std::endl;LOG::filecout(ss.str());}

            T root1,root2;
            int roots=0;
            if(use_kinetic_energy_fix){
                T a_ke=k/2*dt;
                T b_ke=TV::Dot_Product(new_V1-new_V2,state.direction)*dt;
                T c_ke=(current_work+delta_pe)*relaxation_fraction;
                QUADRATIC<T> quadratic(a_ke,b_ke,c_ke);
                quadratic.Compute_Roots();
                roots=quadratic.roots;
                root1=quadratic.root1;
                root2=quadratic.root2;}
            else{
                TV previously_applied_force=previously_applied_forces(s)+force_correction(s)*state.direction;
                T a=k;
                T b=TV::Dot_Product(V+previously_applied_force*k,state.direction);
                T c=2*(current_work+delta_pe)*one_over_dt*relaxation_fraction;
                QUADRATIC<T> quadratic(a,b,c);
                quadratic.Compute_Roots();
                roots=quadratic.roots;
                root1=quadratic.root1;
                root2=quadratic.root2;}
            T force_magnitude;
            {std::stringstream ss;ss<<"ROOTS "<<roots<<std::endl;LOG::filecout(ss.str());}
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
                force_magnitude=TV::Dot_Product(new_V2-new_V1,state.direction)/k;}
            energy_correction_forces(s)=force_magnitude;
            {std::stringstream ss;ss<<"Force I want to apply: "<<force_magnitude*state.direction<<std::endl;LOG::filecout(ss.str());}
            if(!allow_kd_direction_flip){
                if(energy_correction_sign(s)==1){
                    if(force_correction(s)+energy_correction_forces(s)<0)
                        energy_correction_forces(s)=-force_correction(s);}
                else{
                    if(force_correction(s)+energy_correction_forces(s)>0)
                        energy_correction_forces(s)=-force_correction(s);}}
            {std::stringstream ss;ss<<"Correction force: "<<energy_correction_forces(s)*state.direction<<std::endl;LOG::filecout(ss.str());}

            if(use_gauss_seidel_in_energy_correction) Apply_Energy_Correction_Impulse(s,energy_correction_forces(s),dt);}
        
        if(!use_gauss_seidel_in_energy_correction)
            for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
                Apply_Energy_Correction_Impulse(s,energy_correction_forces(s),dt);}

        for(int i=1;i<=rigid_body_collection.rigid_body_particle.array_collection->Size();i++){
            {std::stringstream ss;ss<<"Body "<<i<<" velocity: "<<rigid_body_collection.Rigid_Body(i).Twist()<<std::endl;LOG::filecout(ss.str());}
        }
    }

    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        {std::stringstream ss;ss<<"Force correction for spring "<<s<<": "<<force_correction(s)<<std::endl;LOG::filecout(ss.str());}

        RIGID_BODY<TV> &body1=Body(s,1), &body2=Body(s,2);
        {std::stringstream ss;ss<<"Final total energy: "<<Potential_Energy(s,time)+body1.Kinetic_Energy()+body2.Kinetic_Energy()<<std::endl;LOG::filecout(ss.str());}}
}
//#####################################################################
// Function Attachment_Location
//#####################################################################
template<class TV> TV RIGID_LINEAR_SPRINGS<TV>::
Attachment_Location(int s,int b) const
{
    return Body(s,b).World_Space_Point(attachment_radius(s)(b));
}
template<class TV> TV RIGID_LINEAR_SPRINGS<TV>::
Endpoint_Velocity(int s,int b) const
{
    return Body(s,b).Pointwise_Object_Velocity(Attachment_Location(s,b));
}
//#####################################################################
// Function Spring_Length
//#####################################################################
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
Spring_Length(int s) const
{
    return (Attachment_Location(s,1)-Attachment_Location(s,2)).Magnitude();
}
//#####################################################################
// Function Set_Stiffness
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Set_Stiffness(int b,T stiffness)
{
    if(youngs_modulus.m<b) youngs_modulus.Resize(segment_mesh.elements.m);
    youngs_modulus(b)=stiffness;
    Invalidate_CFL();
}
//#####################################################################
// Function Set_Stiffness
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Set_Damping(int b,T damp)
{
    if(damping.m<b) damping.Resize(segment_mesh.elements.m);
    damping(b)=damp;
    Invalidate_CFL();
}
//#####################################################################
// Function Set_Restlength
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Set_Restlength(int b,T length,T visual)
{
    if(restlength.m<b) restlength.Resize(segment_mesh.elements.m);
    if(visual_restlength.m<b) visual_restlength.Resize(segment_mesh.elements.m);
    restlength(b)=length>=0?length:Spring_Length(b);
    visual_restlength(b)=visual>=0?visual:restlength(b);
    Invalidate_CFL();
}
//#####################################################################
// Function Add_Spring
//#####################################################################
template<class TV> void RIGID_LINEAR_SPRINGS<TV>::
Add_Spring(int body1,int body2,const TV& r1,const TV& r2)
{
    attachment_radius.Append(VECTOR<TV,2>(r1,r2));
    segment_mesh.elements.Append(VECTOR<int,2>(body1,body2));
    Invalidate_CFL();
}
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
Effective_Impulse_Factor(int s,int b) const
{
    return TV::Dot_Product(Body(s,b).Impulse_Factor(Attachment_Location(s,b))*states(s).direction,states(s).direction);
}
template<class TV> typename TV::SCALAR RIGID_LINEAR_SPRINGS<TV>::
Effective_Impulse_Factor(int s) const
{
    return Effective_Impulse_Factor(s,1)+Effective_Impulse_Factor(s,2);
}
template<class TV> const RIGID_BODY<TV>& RIGID_LINEAR_SPRINGS<TV>::
Body(int s,int b) const
{
    return rigid_body_collection.Rigid_Body(segment_mesh.elements(s)(b));
}
template<class TV> RIGID_BODY<TV>& RIGID_LINEAR_SPRINGS<TV>::
Body(int s,int b)
{
    return rigid_body_collection.Rigid_Body(segment_mesh.elements(s)(b));
}
//#####################################################################
template class RIGID_LINEAR_SPRINGS<VECTOR<float,1> >;
template class RIGID_LINEAR_SPRINGS<VECTOR<float,2> >;
template class RIGID_LINEAR_SPRINGS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_LINEAR_SPRINGS<VECTOR<double,1> >;
template class RIGID_LINEAR_SPRINGS<VECTOR<double,2> >;
template class RIGID_LINEAR_SPRINGS<VECTOR<double,3> >;
#endif
