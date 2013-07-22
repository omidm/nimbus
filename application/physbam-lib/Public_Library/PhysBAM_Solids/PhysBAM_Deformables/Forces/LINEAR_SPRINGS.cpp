//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_SPRINGS
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Nonlinear_Equations/LINE_SEARCH.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS_SYSTEM.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <cfloat>
using ::std::sqrt;
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> LINEAR_SPRINGS<TV>::
LINEAR_SPRINGS(PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh_input,const bool implicit)
    :DEFORMABLES_FORCES<TV>(particles),segment_mesh(segment_mesh_input),use_plasticity(false),cache_strain(false),
    use_kinetic_energy_fix(true),relaxation_fraction(1),use_gauss_seidel_in_energy_correction(false),allow_kd_direction_flip(false),verbose(false)
{
    Set_Stiffness(0);Set_Damping(0);
    use_implicit_velocity_independent_forces=implicit;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> LINEAR_SPRINGS<TV>::
~LINEAR_SPRINGS()
{
}
//#####################################################################
// Function Set_Restlength_From_Particles
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Set_Restlength_From_Particles()
{
    Set_Restlength_From_Material_Coordinates(particles.X);
}
//#####################################################################
// Function Set_Restlength_From_Material_Coordinates
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Set_Restlength_From_Material_Coordinates(ARRAY_VIEW<TV> material_coordinates)
{
    restlength.Resize(segment_mesh.elements.m,false);Invalidate_CFL();
    for(int i=1;i<=segment_mesh.elements.m;i++) restlength(i)=(material_coordinates(segment_mesh.elements(i)(1))-material_coordinates(segment_mesh.elements(i)(2))).Magnitude();
    visual_restlength=restlength;
}
//#####################################################################
// Function Clamp_Restlength
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Clamp_Restlength(const T clamped_restlength)
{
    Invalidate_CFL();
    for(int i=1;i<=restlength.m;i++) restlength(i)=max(visual_restlength(i),clamped_restlength);
}
//#####################################################################
// Function Enable_Plasticity
//#####################################################################
template<class TV> template<class T_FIELD> void LINEAR_SPRINGS<TV>::
Enable_Plasticity(const T_FIELD& plastic_yield_strain_input,const T_FIELD& plastic_hardening_input,const T plasticity_clamp_ratio_input)
{
    use_plasticity=true;plasticity_clamp_ratio=plasticity_clamp_ratio_input;
    plastic_yield_strain.Resize(segment_mesh.elements.m,false,false);ARRAYS_COMPUTATIONS::Fill(plastic_yield_strain,plastic_yield_strain_input);
    plastic_hardening.Resize(segment_mesh.elements.m,false,false);ARRAYS_COMPUTATIONS::Fill(plastic_hardening,plastic_hardening_input);
    plastic_visual_restlength=visual_restlength;
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    segment_mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    force_segments.Update(segment_mesh.elements,particle_is_simulated);
    if(cache_strain) strains_of_segment.Resize(segment_mesh.elements.m,false,false);
}
//#####################################################################
// Function Set_Stiffness_Based_On_Reduced_Mass
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Set_Stiffness_Based_On_Reduced_Mass(const T scaling_coefficient)
{
    constant_youngs_modulus=0;youngs_modulus.Resize(segment_mesh.elements.m,false);Invalidate_CFL();
    for(int i=1;i<=segment_mesh.elements.m;i++){
        int end1,end2;segment_mesh.elements(i).Get(end1,end2);T reduced_mass=Pseudo_Inverse(particles.one_over_effective_mass(end1)+particles.one_over_effective_mass(end2));
        youngs_modulus(i)=scaling_coefficient*reduced_mass/restlength(i);}
}
//#####################################################################
// Function Set_Stiffness_Based_On_Reduced_Mass
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Set_Stiffness_Based_On_Reduced_Mass(ARRAY_VIEW<const T> scaling_coefficient)
{
    PHYSBAM_ASSERT(scaling_coefficient.Size()==segment_mesh.elements.m);
    constant_youngs_modulus=0;youngs_modulus.Resize(segment_mesh.elements.m,false);Invalidate_CFL();
    for(int i=1;i<=segment_mesh.elements.m;i++){
        int end1,end2;segment_mesh.elements(i).Get(end1,end2);T reduced_mass=Pseudo_Inverse(particles.one_over_effective_mass(end1)+particles.one_over_effective_mass(end2));
        youngs_modulus(i)=scaling_coefficient(i)*reduced_mass/restlength(i);}
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Set_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    constant_damping=0;damping.Resize(segment_mesh.elements.m,false,false);Invalidate_CFL();
    for(int i=1;i<=segment_mesh.elements.m;i++){
        T harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(segment_mesh.elements(i)(1))+particles.one_over_effective_mass(segment_mesh.elements(i)(2)));
        T ym;if(!youngs_modulus.m) ym=constant_youngs_modulus;else ym=youngs_modulus(i);
        damping(i)=overdamping_fraction*2*sqrt(ym*restlength(i)*harmonic_mass);}
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Set_Overdamping_Fraction(ARRAY_VIEW<const T> overdamping_fraction) // 1 is critically damped
{
    constant_damping=0;damping.Resize(segment_mesh.elements.m,false,false);Invalidate_CFL();
    for(int i=1;i<=segment_mesh.elements.m;i++){
        T harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(segment_mesh.elements(i)(1))+particles.one_over_effective_mass(segment_mesh.elements(i)(2)));
        T ym;if(!youngs_modulus.m) ym=constant_youngs_modulus;else ym=youngs_modulus(i);
        damping(i)=overdamping_fraction(i)*2*sqrt(ym*restlength(i)*harmonic_mass);}
}
//#####################################################################
// Function Ensure_Minimum_Overdamping_Fraction
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    constant_damping=0;damping.Resize(segment_mesh.elements.m);Invalidate_CFL();
    ARRAY<T> save_damping(damping);Set_Overdamping_Fraction(overdamping_fraction);
    for(int k=1;k<=damping.m;k++) damping(k)=max(damping(k),save_damping(k));
}
//#####################################################################
// Function Clamp_Restlength_With_Fraction_Of_Springs
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Clamp_Restlength_With_Fraction_Of_Springs(const T fraction)
{
    ARRAY<T> length(restlength);Sort(length);
    T minimum_restlength=length(min((int)(fraction*length.m)+1,length.m));std::stringstream ss;ss<<"Enlarging the restlength of all linear springs below "<<minimum_restlength<<std::endl;LOG::filecout(ss.str());
    Clamp_Restlength(minimum_restlength);
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    states.Resize(segment_mesh.elements.m,false,false);
    current_lengths.Resize(segment_mesh.elements.m,false,false);
    extra_energy.Resize(segment_mesh.elements.m);
    force_correction.Resize(segment_mesh.elements.m);
    previously_applied_forces.Resize(segment_mesh.elements.m);
    residual_PE.Resize(segment_mesh.elements.m);
    
    ARRAY_VIEW<const TV> X(particles.X);
    if(!damping.m) for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        STATE& state=states(s);
        state.nodes=VECTOR<int,2>(node1,node2);
        state.direction=X(node2)-X(node1);
        current_lengths(s)=state.direction.Normalize();
        state.coefficient=constant_damping/restlength(s);}
    else for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        STATE& state=states(s);
        state.nodes=VECTOR<int,2>(node1,node2);
        state.direction=X(node2)-X(node1);
        current_lengths(s)=state.direction.Normalize();
        state.coefficient=damping(s)/restlength(s);}
    if(use_plasticity) for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        T strain=(current_lengths(s)-visual_restlength(s))/restlength(s);T strain_sign=sign(strain),strain_magnitude=abs(strain);
        if(strain_magnitude>plastic_yield_strain(s)){
            visual_restlength(s)=clamp(current_lengths(s)-strain_sign*plastic_yield_strain(s)*restlength(s),plastic_visual_restlength(s)/plasticity_clamp_ratio,
                plastic_visual_restlength(s)*plasticity_clamp_ratio);
            //visual_restlength(s)=clamp(optimization_current_length(s)-strain_sign*plastic_yield_strain(s)*restlength(s),restlength(s)/plasticity_clamp_ratio,restlength(s)*plasticity_clamp_ratio);
            plastic_yield_strain(s)+=plastic_hardening(s)*(strain_magnitude-plastic_yield_strain(s));}}
    if(compute_half_forces) for(int i=1;i<=states.m;i++) states(i).sqrt_coefficient=sqrt(states(i).coefficient);
    if(!segment_mesh.incident_elements) segment_mesh.Initialize_Incident_Elements();
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    if(!youngs_modulus.m) for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        TV force=constant_youngs_modulus/restlength(s)*(current_lengths(s)-visual_restlength(s))*state.direction;
        F(state.nodes[1])+=force;F(state.nodes[2])-=force;}
    else for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        TV force=youngs_modulus(s)/restlength(s)*(current_lengths(s)-visual_restlength(s))*state.direction;
        F(state.nodes[1])+=force;F(state.nodes[2])-=force;}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    //int springs_processed=0;
    //LOG::Time("spring add velocity dependent");
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        //springs_processed++;
        TV force=(state.coefficient*TV::Dot_Product(V(state.nodes[1])-V(state.nodes[2]),state.direction))*state.direction;
        F(state.nodes[1])-=force;F(state.nodes[2])+=force;}
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    if(!youngs_modulus.m) for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        TV dl=V(node2)-V(node1),dl_projected=dl.Projected_On_Unit_Direction(state.direction);
        TV dforce=constant_youngs_modulus/restlength(s)*dl_projected;
        F(node1)+=dforce;F(node2)-=dforce;}
    else for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        TV dl=V(node2)-V(node1),dl_projected=dl.Projected_On_Unit_Direction(state.direction);
        TV dforce=youngs_modulus(s)/restlength(s)*dl_projected;
        F(node1)+=dforce;F(node2)-=dforce;}
}
//#####################################################################
// Function Velocity_Dependent_Size
//#####################################################################
template<class TV> int LINEAR_SPRINGS<TV>::
Velocity_Dependent_Forces_Size() const
{
    int aggregate_id=0;
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()) aggregate_id++;
    return aggregate_id;
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    int aggregate_id=1;
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        aggregate(aggregate_id++)+=state.sqrt_coefficient*TV::Dot_Product(V(state.nodes[1])-V(state.nodes[2]),state.direction);}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    int aggregate_id=1;
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        TV force=state.sqrt_coefficient*aggregate(aggregate_id++)*state.direction;
        F(state.nodes[1])+=force;F(state.nodes[2])-=force;}
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const
{
    if(!youngs_modulus.m) for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        TV dl=dX(node2)-dX(node1),dl_projected=dl.Projected_On_Unit_Direction(state.direction);
        TV dforce=constant_youngs_modulus/restlength(s)*dl_projected; // Component of dF along direction of the spring
        if(current_lengths(s)>visual_restlength(s)) // Component of dF normal to the spring
            dforce+=constant_youngs_modulus/restlength(s)*(1-visual_restlength(s)/current_lengths(s))*(dl-dl_projected);
        dF(node1)+=dforce;dF(node2)-=dforce;}
    else for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        TV dl=dX(node2)-dX(node1),dl_projected=dl.Projected_On_Unit_Direction(state.direction),dforce;
        dforce+=youngs_modulus(s)/restlength(s)*dl_projected; // Component of dF along direction of the spring
        if(current_lengths(s)>visual_restlength(s)) // Component of dF normal to the spring
            dforce+=youngs_modulus(s)/restlength(s)*(1-visual_restlength(s)/current_lengths(s))*(dl-dl_projected); // Component of dF normal to the spring
        dF(node1)+=dforce;dF(node2)-=dforce;}
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{
    T one_over_cfl_number=1/cfl_number,one_over_cfl_number_squared=sqr(one_over_cfl_number);
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const VECTOR<int,2>& nodes=segment_mesh.elements(s);
        T ym=youngs_modulus.m?youngs_modulus(s):constant_youngs_modulus;
        T d=damping.m?damping(s):constant_damping;
        for(int k=1;k<=2;k++){
            frequency(nodes[k]).elastic_squared+=particles.one_over_effective_mass(nodes[k])/restlength(s)*4*ym*one_over_cfl_number_squared;
            frequency(nodes[k]).damping+=particles.one_over_effective_mass(nodes[k])/restlength(s)*2*d*one_over_cfl_number;}}
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
CFL_Strain_Rate() const
{
    T max_strain_rate=0,strain_rate;TV dx;
    if(use_rest_state_for_strain_rate) for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        int i,j;segment_mesh.elements(s).Get(i,j);
        dx=particles.X(j)-particles.X(i);T magnitude=dx.Magnitude();if(magnitude!=0) dx.Normalize();
        strain_rate=TV::Dot_Product(particles.V(j)-particles.V(i),dx)/restlength(s);
        max_strain_rate=max(max_strain_rate,abs(strain_rate));
        if(cache_strain){strains_of_segment(s)=VECTOR<T,2>(abs(strain_rate),abs((magnitude-visual_restlength(s))/restlength(s)));}}
    else for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        int i,j;segment_mesh.elements(s).Get(i,j);
        dx=particles.X(j)-particles.X(i);
        strain_rate=TV::Dot_Product(particles.V(j)-particles.V(i),dx)/TV::Dot_Product(dx,dx);
        max_strain_rate=max(max_strain_rate,abs(strain_rate));
        if(cache_strain){strains_of_segment(s)=VECTOR<T,2>(abs(strain_rate),abs((dx.Magnitude()-visual_restlength(s))/restlength(s)));}}
    return Robust_Divide(max_strain_per_time_step,max_strain_rate);
}
//#####################################################################
// Function Average_Restlength
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Average_Restlength() const
{
    return ARRAYS_COMPUTATIONS::Average(restlength);
}
//#####################################################################
// Function Print_Restlength_Statistics
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
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
template<class TV> void LINEAR_SPRINGS<TV>::
Print_Deformation_Statistics() const
{
    LOG::SCOPE scope("linear spring deformation","linear spring deformation");
    ARRAY<T> deformation(segment_mesh.elements.m,false);
    for(int s=1;s<=segment_mesh.elements.m;s++){
        int i,j;segment_mesh.elements(s).Get(i,j);
        T length=(particles.X(i)-particles.X(j)).Magnitude(),rl=visual_restlength(s);
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
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Maximum_Compression_Or_Expansion_Fraction(int* index) const
{
    T max_compression=0;int max_index=1;
    for(int s=1;s<=segment_mesh.elements.m;s++){
        int i,j;segment_mesh.elements(s).Get(i,j);
        T length=(particles.X(i)-particles.X(j)).Magnitude();
        T rl=visual_restlength(s);
        T compression=(rl)?(abs(length-rl)/rl):((length==0)?0:FLT_MAX);
        if(compression>max_compression){max_compression=compression;max_index=s;}}
    if(index) (*index)=max_index;
    return max_compression;
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Potential_Energy(const int s,const T time) const
{
    T current_length=(particles.X(segment_mesh.elements(s)(1))-particles.X(segment_mesh.elements(s)(2))).Magnitude();
    T spring_youngs_modulus=youngs_modulus.m?youngs_modulus(s):constant_youngs_modulus;
    return (T).5*spring_youngs_modulus/restlength(s)*sqr(current_length-visual_restlength(s));
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        potential_energy+=Potential_Energy(s,time);}
    return potential_energy;
}
//#####################################################################
// Function Residual_Energy
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Residual_Energy(const T time) const
{
    T residual_energy=0;
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        residual_energy+=residual_PE(s);}
    return residual_energy;
}
//#####################################################################
// Function Add_Force_Data
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Add_Force_Data(ARRAY<FORCE_DATA<TV> >& force_data_list,const std::string& force_name) const
{   
    ARRAY_VIEW<const TV> X=particles.X;
    FORCE_DATA<TV> force_data;
    if(force_name.empty()) force_data.name="LINEAR_SPRINGS";
    else force_data.name=force_name;
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        force_data.first_action_point=X(node1);
        force_data.second_action_point=X(node2);
        T current_length=(X(node2)-X(node1)).Magnitude();
        force_data.state=visual_restlength(s)?(current_length-visual_restlength(s))/visual_restlength(s):(T)0;
        force_data_list.Append(force_data);}
}
//#####################################################################
// Function Prepare_Energy_Correction_Force
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Prepare_Energy_Correction_Force()
{
    energy_correction_forces.Resize(segment_mesh.elements.m);
    energy_correction_sign.Resize(segment_mesh.elements.m);

    // figure out damping directions
    if(verbose) for(int i=1;i<=particles.array_collection->Size();i++) {std::stringstream ss;ss<<"Initial body "<<i<<" velocity: "<<particles.V(i)<<std::endl;LOG::filecout(ss.str());}

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
}
//#####################################################################
// Function Compute_Energy_Correction_Force
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Compute_Energy_Correction_Force(ARRAY_VIEW<const TV> velocity_save,const int max_particle_degree,const T time,const T dt)
{
    T one_over_dt=1/dt;

    ARRAYS_COMPUTATIONS::Fill(energy_correction_forces,(T)0);

    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        int first_node=state.nodes(1);
        int second_node=state.nodes(2);
        TV old_V1=velocity_save(first_node);
        TV old_V2=velocity_save(second_node);
        TV new_V1=Endpoint_Velocity(s,1), new_V2=Endpoint_Velocity(s,2);
        
        TV V=new_V1-new_V2+old_V1-old_V2;
        T k=dt*Effective_Impulse_Factor(s);
        
        T current_potential_energy=Potential_Energy(s,time);
        T current_kinetic_energy=Endpoint_Kinetic_Energy(s);
        T delta_pe=current_potential_energy-potential_energy_save(s)+extra_energy(s);
        
        if(verbose) {std::stringstream ss;ss<<"SPRING "<<s<<" ****************"<<std::endl;LOG::filecout(ss.str());}
        T current_work=Compute_Spring_Work(s,velocity_save,time,dt);
        if(verbose) {std::stringstream ss;ss<<"Total work done: "<<current_work<<std::endl;LOG::filecout(ss.str());}
        
        if(current_work+delta_pe<0)
            continue;
        
        if(verbose) {std::stringstream ss;ss<<"Saved potential energy: "<<potential_energy_save(s)<<" current potential energy: "<<current_potential_energy<<" delta pe: "<<delta_pe<<std::endl;LOG::filecout(ss.str());}
        if(verbose) {std::stringstream ss;ss<<"Current total energy: "<<current_potential_energy+current_kinetic_energy<<std::endl;LOG::filecout(ss.str());}
        
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
            force_magnitude=TV::Dot_Product(new_V2-new_V1,state.direction)/k;}
        energy_correction_forces(s)=force_magnitude/max_particle_degree;
        if(verbose) {std::stringstream ss;ss<<"Force I want to apply: "<<energy_correction_forces(s)*state.direction<<std::endl;LOG::filecout(ss.str());}
        if(!allow_kd_direction_flip){
            if(energy_correction_sign(s)==1){
                if(force_correction(s)+energy_correction_forces(s)<0)
                    energy_correction_forces(s)=-force_correction(s);}
            else{
                if(force_correction(s)+energy_correction_forces(s)>0)
                    energy_correction_forces(s)=-force_correction(s);}}
        if(verbose) {std::stringstream ss;ss<<"Correction force: "<<energy_correction_forces(s)*state.direction<<std::endl;LOG::filecout(ss.str());}
        
        if(use_gauss_seidel_in_energy_correction) Apply_Energy_Correction_Impulse(s,energy_correction_forces(s),dt);}
}
//#####################################################################
// Function Apply_Energy_Correction
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Apply_Energy_Correction(const T time,const T dt)
{        
    if(!use_gauss_seidel_in_energy_correction)
        for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
            Apply_Energy_Correction_Impulse(s,energy_correction_forces(s),dt);}
    
    if(verbose) for(int i=1;i<=particles.array_collection->Size();i++) {std::stringstream ss;ss<<"Body "<<i<<" velocity: "<<particles.V(i)<<std::endl;LOG::filecout(ss.str());}

    if(verbose) 
        for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
            std::stringstream ss;ss<<"Force correction for spring "<<s<<": "<<force_correction(s)<<std::endl;
            ss<<"Final total energy: "<<Potential_Energy(s,time)+Endpoint_Kinetic_Energy(s)<<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
// Function Add_Connectivity
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Add_Connectivity(ARRAY<int>& particle_degree)
{
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        particle_degree(state.nodes(1))++;
        particle_degree(state.nodes(2))++;}
}
//#####################################################################
// Function Endpoint_Velocity
//#####################################################################
template<class TV> TV LINEAR_SPRINGS<TV>::
Endpoint_Velocity(int s,int b) const
{
    return particles.V(segment_mesh.elements(s)(b));
}
//#####################################################################
// Function Endpoint_Kinetic_Energy
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Endpoint_Kinetic_Energy(int s,int b) const
{
    return (T).5*particles.mass(segment_mesh.elements(s)(b))*particles.V(segment_mesh.elements(s)(b)).Magnitude_Squared();
}
//#####################################################################
// Function Endpoint_Kinetic_Energy
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Endpoint_Kinetic_Energy(int s) const
{
    return Endpoint_Kinetic_Energy(s,1)+Endpoint_Kinetic_Energy(s,2);
}
//#####################################################################
// Function Effective_Impulse_Factor
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Effective_Impulse_Factor(int s) const
{
    T one_over_denom=particles.one_over_mass(segment_mesh.elements(s)(1))*particles.one_over_mass(segment_mesh.elements(s)(2));
    return (particles.mass(segment_mesh.elements(s)(1))+particles.mass(segment_mesh.elements(s)(1)))*one_over_denom;
}
//#####################################################################
// Function Store_Potential_Energy
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Save_Potential_Energy(const T time)
{
    potential_energy_save.Resize(segment_mesh.elements.m);
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        potential_energy_save(s)=Potential_Energy(s,time);}
}
//#####################################################################
// Function Compute_Spring_Work
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Compute_Spring_Work(int s,ARRAY_VIEW<const TV> velocity_save,const T time,const T dt) const
{
    const STATE& state=states(s);
    TV new_V1=Endpoint_Velocity(s,1), new_V2=Endpoint_Velocity(s,2);
    TV force=previously_applied_forces(s);
    force+=force_correction(s)*state.direction;
    TV old_V1=velocity_save(state.nodes(1));
    TV old_V2=velocity_save(state.nodes(2));
    TV average_vel1=(new_V1+old_V1)*(T).5;
    TV average_vel2=(new_V2+old_V2)*(T).5;
    TV distance1=average_vel1*dt;
    TV distance2=average_vel2*dt;
    T first_node_work=TV::Dot_Product(distance1,force);
    T second_node_work=-TV::Dot_Product(distance2,force);
    if(verbose) {std::stringstream ss;ss<<"Work computation: force: "<<force<<" average v1: "<<average_vel1<<" average v2: "<<average_vel2<<std::endl;LOG::filecout(ss.str());}
    return first_node_work+second_node_work;
}
//#####################################################################
// Function Compute_Previously_Applied_Forces
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Compute_Previously_Applied_Forces()
{
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        const STATE& state=states(s);
        const T ym=youngs_modulus.m?youngs_modulus(s):constant_youngs_modulus;
        previously_applied_forces(s)=state.coefficient*TV::Dot_Product(Endpoint_Velocity(s,2)-Endpoint_Velocity(s,1),state.direction)*state.direction;
        previously_applied_forces(s)+=ym/restlength(s)*(current_lengths(s)-visual_restlength(s))*state.direction;
        force_correction(s)=0;}
}
//#####################################################################
// Function Apply_Energy_Correction_impulse
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Apply_Energy_Correction_Impulse(const int s,const T force,const T dt)
{
    force_correction(s)+=force;
    TV impulse=dt*force*states(s).direction;
    //LOG::cout<<"Force correction for spring "<<s<<": "<<energy_correction_forces(s)<<std::endl;
    
    int node1=states(s).nodes(1),node2=states(s).nodes(2);
    particles.V(node1)+=(impulse*particles.one_over_mass(node1));
    particles.V(node2)-=(impulse*particles.one_over_mass(node2));
    /*LOG::cout<<"Final total energy: "<<Potential_Energy(s,time)+body1.Kinetic_Energy()+body2.Kinetic_Energy()<<std::endl;*/
}
//#####################################################################
// Function Compute_Energy_Error
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Compute_Energy_Error(ARRAY_VIEW<const TV> velocity_save,const T time,const T dt)
{
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        T current_potential_energy=Potential_Energy(s,time);
        T spring_extra_energy=(current_potential_energy-potential_energy_save(s))+extra_energy(s)+Compute_Spring_Work(s,velocity_save,time,dt);
        extra_energy(s)=max((T)0,spring_extra_energy);}
}
//#####################################################################
// Function Setup_Set_Velocity_From_Positions
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Setup_Set_Velocity_From_Positions(const T time,const bool is_position_update,const bool reset_alphas)
{
    if(!is_position_update) Update_Position_Based_State(time,is_position_update);

    force_estimates.Resize(segment_mesh.elements.m);
    incident_nodes.Resize(segment_mesh.elements.m);
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        incident_nodes(s).Remove_All();
        incident_nodes(s).Append(node1);incident_nodes(s).Append(node2);}
    if(is_position_update){
        if(reset_alphas) for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();force_estimates(s)=0;}
        Save_Potential_Energy(time);}
}
//#####################################################################
// Function Get_Element_Count
//#####################################################################
template<class TV> int LINEAR_SPRINGS<TV>::
Get_Element_Count()
{
    return segment_mesh.elements.m;
}
//#####################################################################
// Function Get_Force_Elements
//#####################################################################
template<class TV> FORCE_ELEMENTS* LINEAR_SPRINGS<TV>::
Get_Force_Elements()
{
    return &force_segments;
}
//#####################################################################
// Function Incident_Nodes
//#####################################################################
template<class TV> ARRAY<int>* LINEAR_SPRINGS<TV>::
Incident_Nodes(const int force_element)
{
    return &(incident_nodes(force_element));
}
//#####################################################################
// Function Get_Direction
//#####################################################################
template<class TV> TV LINEAR_SPRINGS<TV>::
Get_Direction(const int force_element)
{
    const STATE& state=states(force_element);
    return state.direction;
}
//#####################################################################
// Function Get_Combined_One_Over_Mass
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Get_Combined_One_Over_Mass(const int force_element)
{
    int node1,node2;segment_mesh.elements(force_element).Get(node1,node2);
    return particles.one_over_mass(node1)+particles.one_over_mass(node2);
}
//#####################################################################
// Function Incident_Force_Elements
//#####################################################################
template<class TV> ARRAY<int>* LINEAR_SPRINGS<TV>::
Incident_Force_Elements(const int particle)
{
    return &(*segment_mesh.incident_elements)(particle);
}
//#####################################################################
// Function Choose_Solution
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Choose_Solution(const bool use_orig_force,const int force_element,const T dt,const T alpha1,const T alpha2,ARRAY<T>& v_n_hats)
{
    if(use_orig_force){
        T force=constant_youngs_modulus/restlength(force_element)*(current_lengths(force_element)-restlength(force_element));
        
        if(fabs(force-alpha1) < fabs(force-alpha2)) Set_Force(force_element,alpha1);
        else Set_Force(force_element,alpha2);}
    else{
        const STATE& state=states(force_element);
        T v1_alpha1_temp=v_n_hats(1)+dt*particles.one_over_mass(state.nodes(1))*alpha1;
        T v2_alpha1_temp=v_n_hats(2)-dt*particles.one_over_mass(state.nodes(2))*alpha1;
        
        T v1_alpha2_temp=v_n_hats(1)+dt*particles.one_over_mass(state.nodes(1))*alpha2;
        T v2_alpha2_temp=v_n_hats(2)-dt*particles.one_over_mass(state.nodes(2))*alpha2;
        
        T combined_mass=particles.mass(state.nodes(1))+particles.mass(state.nodes(2));
        T x1u=TV::Dot_Product(particles.X(state.nodes(1)),state.direction);
        T x2u=TV::Dot_Product(particles.X(state.nodes(2)),state.direction);
        
        T w;
        if(!youngs_modulus.m) w=sqrt((constant_youngs_modulus/restlength(force_element))*Get_Combined_One_Over_Mass(force_element));
        else w=sqrt((youngs_modulus(force_element)/restlength(force_element))*Get_Combined_One_Over_Mass(force_element));
        T v_cm_0=(particles.mass(state.nodes(1))*v_n_hats(1)+particles.mass(state.nodes(2))*v_n_hats(2))/combined_mass;
        T c1=(x2u-x1u)-restlength(force_element);
        T c2=(v_n_hats(2)-v_n_hats(1))/w;
        T v1_analytic=v_cm_0-(particles.mass(state.nodes(2))/(combined_mass))*(-w*c1*sin(w*dt)+w*c2*cos(w*dt));
        T v2_analytic=v_cm_0+(particles.mass(state.nodes(1))/(combined_mass))*(-w*c1*sin(w*dt)+w*c2*cos(w*dt));
        if((fabs(v1_alpha1_temp-v1_analytic)+fabs(v2_alpha1_temp-v2_analytic)) < (fabs(v1_alpha2_temp-v1_analytic)+fabs(v2_alpha2_temp-v2_analytic))) Set_Force(force_element,alpha1);
        else Set_Force(force_element,alpha2);}
}
//#####################################################################
// Function Get_Force
//#####################################################################
template<class TV> TV LINEAR_SPRINGS<TV>::
Get_Force(const int force_element,const int particle,const bool use_original_force)
{
    const STATE& state=states(force_element);
    TV force;
    if(use_original_force){
        if(youngs_modulus.m)
            force=youngs_modulus(force_element)/restlength(force_element)*(current_lengths(force_element)-visual_restlength(force_element))*state.direction;
        else
            force=constant_youngs_modulus/restlength(force_element)*(current_lengths(force_element)-visual_restlength(force_element))*state.direction;}
    else force=force_estimates(force_element)*state.direction;
    return (particle==segment_mesh.elements(force_element)(1))?force:-force;
}
//#####################################################################
// Function Get_Force
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Get_Force(const int force_element)
{
    return force_estimates(force_element);
}
//#####################################################################
// Function Set_Force
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Set_Force(const int force_element,const T force)
{
    force_estimates(force_element)=force;
}
//#####################################################################
// Function Get_Damping_Force
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Get_Damping_Force(const int particle,TV& damping_force,const T dt,const bool use_coefficient)
{
    ARRAY<int>& node_incident_forces=(*segment_mesh.incident_elements)(particle);
    for(int s=1;s<=node_incident_forces.m;s++){
        const STATE& state=states(node_incident_forces(s));
        TV force=TV::Dot_Product(particles.V(state.nodes[1])-particles.V(state.nodes[2]),state.direction)*state.direction;
        if(particle==state.nodes[1]) force*=(T)-1;
        if(use_coefficient) damping_force+=dt*particles.one_over_mass(particle)*state.coefficient*force;
        else damping_force+=dt*particles.one_over_mass(particle)*force;}
}
//#####################################################################
// Function Update_Residual_Energy
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Update_Residual_Energy(const int force_element,const T residual_energy,const T time)
{
    residual_PE(force_element)=residual_energy;
}
//#####################################################################
// Function Get_Residual_Energy
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Get_Residual_Energy(const int force_element)
{
    return residual_PE(force_element);
}
//#####################################################################
// Function Compute_Quadratic_Contribution_For_Node
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Compute_Quadratic_Contribution_For_Force(T& A,T& a,T&c,const T dt,const int force_element,const T combined_one_over_mass,const bool ignore_PE_terms)
{
    if(!ignore_PE_terms){
        const STATE& state=states(force_element);

        T x1u=TV::Dot_Product(particles.X(state.nodes(1)),state.direction);
        T x2u=TV::Dot_Product(particles.X(state.nodes(2)),state.direction);

        a=x2u-x1u;
        c=-(T).5*dt*dt*combined_one_over_mass;

        if(!youngs_modulus.m) A+=(T).5*(constant_youngs_modulus/restlength(force_element))*c*c;
        else A+=(T).5*(youngs_modulus(force_element)/restlength(force_element))*c*c;}

    A+=(T).5*dt*dt*combined_one_over_mass;
}
//#####################################################################
// Function Compute_Quadratic_Contribution_For_Node
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Compute_Quadratic_Contribution_For_Node(T& B,T& C,T&b,const T dt,const int node,const int force_element,const T combined_one_over_mass,const T v_n_correction,const bool ignore_PE_terms)
{
    const STATE& state=states(force_element);

    T vu=TV::Dot_Product(particles.V(node),state.direction);

    if(!ignore_PE_terms) b+=((node==state.nodes(1))?-dt:dt)*(vu+(T).5*v_n_correction);
    B+=((node==state.nodes(1))?dt:-dt)*(vu+(T).5*v_n_correction);
}
//#####################################################################
// Function Compute_Quadratic_Contribution_For_Residual
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Compute_Quadratic_Contribution_For_Residual(T& B,T& C,T& a,T&b,T&c,const T dt,const T time,const int force_element,const bool ignore_PE_terms)
{
    if(!ignore_PE_terms){
        if(!youngs_modulus.m){
            B+=(constant_youngs_modulus/restlength(force_element))*(a+b-restlength(force_element))*c;
            C=(constant_youngs_modulus/restlength(force_element))*(a+b/2-restlength(force_element))*b;}
        else{
            B+=(youngs_modulus(force_element)/restlength(force_element))*(a+b-restlength(force_element))*c;
            C=(youngs_modulus(force_element)/restlength(force_element))*(a+b/2-restlength(force_element))*b;}}
    else C=delta_PE(force_element);

    C+=residual_PE(force_element);
}
//#####################################################################
// Function Store_Delta_PE
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Store_Delta_PE(const T time)
{
    delta_PE.Resize(potential_energy_save.m);
    T total_current_PE=0;
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        delta_PE(s)=Potential_Energy(s,time)-potential_energy_save(s);
        total_current_PE+=Potential_Energy(s,time);}
    total_delta_PE=total_current_PE-ARRAYS_COMPUTATIONS::Sum(potential_energy_save);
}
//#####################################################################
// Function Get_Total_Delta_PE
//#####################################################################
template<class TV> typename TV::SCALAR LINEAR_SPRINGS<TV>::
Get_Total_Delta_PE()
{
    return total_delta_PE;
}
//#####################################################################
// Function Save_And_Reset_Elastic_Coefficients
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Save_And_Reset_Elastic_Coefficient()
{
    saved_constant_youngs_modulus=constant_youngs_modulus;
    constant_youngs_modulus=0;
}
//#####################################################################
// Function Restore_Elastic_Coefficients
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Restore_Elastic_Coefficient()
{
    constant_youngs_modulus=saved_constant_youngs_modulus;
}
//#####################################################################
// Function Store_Velocities
//#####################################################################
template<class TV> void LINEAR_SPRINGS<TV>::
Store_Velocities()
{
    saved_V.Resize(particles.array_collection->Size());
    for(int i=1;i<=particles.array_collection->Size();i++) saved_V(i)=particles.V(i);
}
//#####################################################################
// Function Create_Edge_Springs
//#####################################################################
template<class TV> LINEAR_SPRINGS<TV>* PhysBAM::
Create_Edge_Springs(PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh,const typename TV::SCALAR stiffness,const typename TV::SCALAR overdamping_fraction,
    const bool limit_time_step_by_strain_rate,const typename TV::SCALAR max_strain_per_time_step,const bool use_rest_state_for_strain_rate,
    const typename TV::SCALAR restlength_enlargement_fraction,const bool verbose,const bool implicit)
{
    LINEAR_SPRINGS<TV>* ls=new LINEAR_SPRINGS<TV>(particles,segment_mesh,implicit);
    ls->Set_Restlength_From_Particles();
    if(restlength_enlargement_fraction) ls->Clamp_Restlength_With_Fraction_Of_Springs(restlength_enlargement_fraction);
    ls->Set_Stiffness(stiffness);
    ls->Set_Overdamping_Fraction(overdamping_fraction);
    ls->Limit_Time_Step_By_Strain_Rate(limit_time_step_by_strain_rate,max_strain_per_time_step);
    ls->Use_Rest_State_For_Strain_Rate(use_rest_state_for_strain_rate);
    if(verbose) ls->Print_Restlength_Statistics();
    ls->verbose=verbose;
    return ls;
}
//#####################################################################
// Function Create_Edge_Springs
//#####################################################################
template<class T_OBJECT> LINEAR_SPRINGS<typename T_OBJECT::VECTOR_T>* PhysBAM::
Create_Edge_Springs(T_OBJECT& object,
    const typename T_OBJECT::SCALAR stiffness,const typename T_OBJECT::SCALAR overdamping_fraction,const bool limit_time_step_by_strain_rate,
    const typename T_OBJECT::SCALAR max_strain_per_time_step,const bool use_rest_state_for_strain_rate,const typename T_OBJECT::SCALAR restlength_enlargement_fraction,
    const bool verbose,const bool implicit)
{
    return Create_Edge_Springs(dynamic_cast<PARTICLES<typename T_OBJECT::VECTOR_T>&>(object.particles),object.Get_Segment_Mesh(),stiffness,overdamping_fraction,limit_time_step_by_strain_rate,max_strain_per_time_step,
        use_rest_state_for_strain_rate,restlength_enlargement_fraction,verbose,implicit);
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template class LINEAR_SPRINGS<VECTOR<T,2> >; \
    template class LINEAR_SPRINGS<VECTOR<T,3> >; \
    template LINEAR_SPRINGS<SEGMENTED_CURVE<VECTOR<T,1> >::VECTOR_T>* PhysBAM::Create_Edge_Springs<SEGMENTED_CURVE<VECTOR<T,1> > >(SEGMENTED_CURVE<VECTOR<T,1> >&, \
        SEGMENTED_CURVE<VECTOR<T,1> >::SCALAR,SEGMENTED_CURVE<VECTOR<T,1> >::SCALAR,bool,SEGMENTED_CURVE<VECTOR<T,1> >::SCALAR,bool,SEGMENTED_CURVE<VECTOR<T,1> >::SCALAR,bool,bool); \
    template LINEAR_SPRINGS<SEGMENTED_CURVE<VECTOR<T,2> >::VECTOR_T>* PhysBAM::Create_Edge_Springs<SEGMENTED_CURVE<VECTOR<T,2> > >(SEGMENTED_CURVE<VECTOR<T,2> >&, \
        SEGMENTED_CURVE<VECTOR<T,2> >::SCALAR,SEGMENTED_CURVE<VECTOR<T,2> >::SCALAR,bool,SEGMENTED_CURVE<VECTOR<T,2> >::SCALAR,bool,SEGMENTED_CURVE<VECTOR<T,2> >::SCALAR,bool,bool); \
    template LINEAR_SPRINGS<SEGMENTED_CURVE<VECTOR<T,3> >::VECTOR_T>* PhysBAM::Create_Edge_Springs<SEGMENTED_CURVE<VECTOR<T,3> > >(SEGMENTED_CURVE<VECTOR<T,3> >&, \
        SEGMENTED_CURVE<VECTOR<T,3> >::SCALAR,SEGMENTED_CURVE<VECTOR<T,3> >::SCALAR,bool,SEGMENTED_CURVE<VECTOR<T,3> >::SCALAR,bool,SEGMENTED_CURVE<VECTOR<T,3> >::SCALAR,bool,bool); \
    template LINEAR_SPRINGS<TETRAHEDRALIZED_VOLUME<T>::VECTOR_T>* PhysBAM::Create_Edge_Springs<TETRAHEDRALIZED_VOLUME<T> >(TETRAHEDRALIZED_VOLUME<T>&,TETRAHEDRALIZED_VOLUME<T>::SCALAR, \
        TETRAHEDRALIZED_VOLUME<T>::SCALAR,bool,TETRAHEDRALIZED_VOLUME<T>::SCALAR,bool,TETRAHEDRALIZED_VOLUME<T>::SCALAR,bool,bool); \
    template LINEAR_SPRINGS<TRIANGULATED_SURFACE<T>::VECTOR_T>* PhysBAM::Create_Edge_Springs<TRIANGULATED_SURFACE<T> >(TRIANGULATED_SURFACE<T>&,TRIANGULATED_SURFACE<T>::SCALAR, \
        TRIANGULATED_SURFACE<T>::SCALAR,bool,TRIANGULATED_SURFACE<T>::SCALAR,bool,TRIANGULATED_SURFACE<T>::SCALAR,bool,bool); \
    template LINEAR_SPRINGS<VECTOR<T,3> >* PhysBAM::Create_Edge_Springs<VECTOR<T,3> >(PARTICLES<VECTOR<T,3> >&,SEGMENT_MESH&,VECTOR<T,3>::SCALAR,VECTOR<T,3>::SCALAR,bool, \
        VECTOR<T,3>::SCALAR,bool,VECTOR<T,3>::SCALAR,bool,bool); \
    template LINEAR_SPRINGS<SEGMENTED_CURVE_2D<T>::VECTOR_T>* PhysBAM::Create_Edge_Springs<SEGMENTED_CURVE_2D<T> >(SEGMENTED_CURVE_2D<T>&,SEGMENTED_CURVE_2D<T>::SCALAR, \
        SEGMENTED_CURVE_2D<T>::SCALAR,bool,SEGMENTED_CURVE_2D<T>::SCALAR,bool,SEGMENTED_CURVE_2D<T>::SCALAR,bool,bool); \
    template LINEAR_SPRINGS<TRIANGULATED_AREA<T>::VECTOR_T>* PhysBAM::Create_Edge_Springs<TRIANGULATED_AREA<T> >(TRIANGULATED_AREA<T>&,TRIANGULATED_AREA<T>::SCALAR, \
        TRIANGULATED_AREA<T>::SCALAR,bool,TRIANGULATED_AREA<T>::SCALAR,bool,TRIANGULATED_AREA<T>::SCALAR,bool,bool);

INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
