//#####################################################################
// Copyright 2006-2007, Craig Schroeder, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_ZERO_LENGTH_SPRINGS
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/IMPLICIT_ZERO_LENGTH_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using ::std::sqrt;
using namespace PhysBAM;
template<class TV> IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
IMPLICIT_ZERO_LENGTH_SPRINGS(PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh_input)
    :DEFORMABLES_FORCES<TV>(particles),segment_mesh(segment_mesh_input)
{
    Set_Stiffness(0);Set_Damping(0);
    use_implicit_velocity_independent_forces=true;
}
template<class TV> IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
~IMPLICIT_ZERO_LENGTH_SPRINGS()
{
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    segment_mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    force_segments.Update(segment_mesh.elements,particle_is_simulated);
}
//#####################################################################
// Function Set_Stiffness_Based_On_Reduced_Mass
//#####################################################################
template<class TV> void IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Set_Stiffness_Based_On_Reduced_Mass(const T scaling_coefficient) // assumes mass is already defined
{
    constant_stiffness=0;stiffness.Resize(segment_mesh.elements.m,false);
    for(int i=1;i<=segment_mesh.elements.m;i++){
        int end1,end2;segment_mesh.elements(i).Get(end1,end2);
        T reduced_mass=Pseudo_Inverse(particles.one_over_effective_mass(end1)+particles.one_over_effective_mass(end2));
        stiffness(i)=scaling_coefficient*reduced_mass;}
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class TV> void IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Set_Overdamping_Fraction(const T overdamping_fraction) // 1 is critically damped
{
    constant_damping=0;damping.Resize(segment_mesh.elements.m,false,false);
    for(int i=1;i<=segment_mesh.elements.m;i++){
        T harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(segment_mesh.elements(i)(1))+particles.one_over_effective_mass(segment_mesh.elements(i)(2)));
        T ym;if(!stiffness.m) ym=constant_stiffness;else ym=stiffness(i);
        damping(i)=overdamping_fraction*2*sqrt(ym*harmonic_mass);}
}
//#####################################################################
// Function Set_Overdamping_Fraction
//#####################################################################
template<class TV> void IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Set_Overdamping_Fraction(ARRAY_VIEW<const T> overdamping_fraction) // 1 is critically damped
{
    constant_damping=0;damping.Resize(segment_mesh.elements.m,false,false);
    for(int i=1;i<=segment_mesh.elements.m;i++){
        T harmonic_mass=Pseudo_Inverse(particles.one_over_effective_mass(segment_mesh.elements(i)(1))+particles.one_over_effective_mass(segment_mesh.elements(i)(2)));
        T ym;if(!stiffness.m) ym=constant_stiffness;else ym=stiffness(i);
        damping(i)=overdamping_fraction(i)*2*sqrt(ym*harmonic_mass);}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    Add_Force_Differential(particles.X,F,time);
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    if(!damping.m) for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        TV force=constant_damping*(V(node1)-V(node2));
        F(node1)-=force;F(node2)+=force;}
    else for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        TV force=damping(s)*(V(node1)-V(node2));
        F(node1)-=force;F(node2)+=force;}
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class TV> void IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    if(!stiffness.m) for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        TV force=constant_stiffness*(V(node2)-V(node1));
        F(node1)+=force;F(node2)-=force;}
    else for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        TV force=stiffness(s)*(V(node2)-V(node1));
        F(node1)+=force;F(node2)-=force;}
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class TV> void IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const
{
    Add_Implicit_Velocity_Independent_Forces(dX,dF,time);
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR IMPLICIT_ZERO_LENGTH_SPRINGS<TV>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    for(SEGMENT_ITERATOR iterator(force_segments);iterator.Valid();iterator.Next()){int s=iterator.Data();
        int node1,node2;segment_mesh.elements(s).Get(node1,node2);
        T stiff=stiffness.m?stiffness(s):constant_stiffness;
        potential_energy+=(T).5*stiff*(particles.X(node2)-particles.X(node1)).Magnitude_Squared();}
    return potential_energy;
}
//#####################################################################

//#####################################################################
// Function Create_Edge_Zero_Length_Springs
//#####################################################################
template<class T,class TV> IMPLICIT_ZERO_LENGTH_SPRINGS<TV>* PhysBAM::
Create_Edge_Zero_Length_Springs(PARTICLES<TV>& particles,SEGMENT_MESH& segment_mesh,const T stiffness,const T overdamping_fraction,const bool verbose)
{
    IMPLICIT_ZERO_LENGTH_SPRINGS<TV>* zls=new IMPLICIT_ZERO_LENGTH_SPRINGS<TV>(particles,segment_mesh);
    zls->Set_Stiffness(stiffness);
    zls->Set_Overdamping_Fraction(overdamping_fraction);
    return zls;
}

template class IMPLICIT_ZERO_LENGTH_SPRINGS<VECTOR<float,2> >;
template class IMPLICIT_ZERO_LENGTH_SPRINGS<VECTOR<float,3> >;
template IMPLICIT_ZERO_LENGTH_SPRINGS<VECTOR<float,3> >* PhysBAM::Create_Edge_Zero_Length_Springs<float,VECTOR<float,3> >(PARTICLES<VECTOR<float,3> >&,SEGMENT_MESH&,float,float,bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class IMPLICIT_ZERO_LENGTH_SPRINGS<VECTOR<double,2> >;
template class IMPLICIT_ZERO_LENGTH_SPRINGS<VECTOR<double,3> >;
template IMPLICIT_ZERO_LENGTH_SPRINGS<VECTOR<double,3> >* PhysBAM::Create_Edge_Zero_Length_Springs<double,VECTOR<double,3> >(PARTICLES<VECTOR<double,3> >&,SEGMENT_MESH&,double,double,bool);
#endif
