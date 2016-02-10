//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Craig Schroeder, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Data_Structures/UNION_FIND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/cyclic_shift.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class T> void SEGMENT_BENDING_ELEMENTS<T>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    for(int t=1;t<=bending_triples.m;t++)
        for(int i=1;i<=2;i++) for(int j=i+1;j<=3;j++) dependency_mesh.Add_Element_If_Not_Already_There(VECTOR<int,2>(bending_triples(t)[i],bending_triples(t)[j]));
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class T> void SEGMENT_BENDING_ELEMENTS<T>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    force_triples.Update(bending_triples,particle_is_simulated);
}
//#####################################################################
// Function Set_Triples_From_Segment_Mesh
//#####################################################################
template<class T> void SEGMENT_BENDING_ELEMENTS<T>::
Set_Triples_From_Segment_Mesh(SEGMENT_MESH& mesh)
{
    if(!mesh.adjacent_elements) mesh.Initialize_Adjacent_Elements();

    // allocate proper array sizes
    int number_triples=0;
    for(int t=1;t<=mesh.elements.m;t++) for(int a=1;a<=(*mesh.adjacent_elements)(t).m;a++) if((*mesh.adjacent_elements)(t)(a)>t) number_triples++;
    bending_triples.Resize(number_triples);length_scale.Resize(number_triples);stiffness.Resize(number_triples);
    sine_half_rest_angle.Resize(number_triples);damping.Resize(number_triples);

    int index=0; // reset number
    for(int t=1;t<=mesh.elements.m;t++){
        VECTOR<int,2> segment1=mesh.elements(t);
        for(int a=1;a<=(*mesh.adjacent_elements)(t).m;a++){
            int s=(*mesh.adjacent_elements)(t)(a);
            if(s>t){
                VECTOR<int,2> segment2=mesh.elements(s);
                if(segment2.Contains(segment1.x)) cyclic_shift(segment1);
                if(segment1.Contains(segment2.y)) cyclic_shift(segment2);
                bending_triples(++index).Set(segment1.x,segment1.y,segment2.y);}}}
}
//#####################################################################
// Function Set_Constants_From_Particles
//#####################################################################
template<class T> void SEGMENT_BENDING_ELEMENTS<T>::
Set_Constants_From_Particles(const T material_stiffness,const T material_damping)
{
    for(int q=1;q<=bending_triples.m;q++){
        int i,j,k;bending_triples(q).Get(i,j,k);
        TV n1=(particles.X(j)-particles.X(i)).Rotate_Clockwise_90(),n2=(particles.X(k)-particles.X(j)).Rotate_Clockwise_90();
        length_scale(q)=(T).5*(n1.Normalize()+n2.Normalize());
        sine_half_rest_angle(q)=Sine_Half_Angle_Between(n1,n2);}
    ARRAYS_COMPUTATIONS::Fill(stiffness,material_stiffness);
    ARRAYS_COMPUTATIONS::Fill(damping,material_damping);
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T> void SEGMENT_BENDING_ELEMENTS<T>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    elastic_s.Resize(bending_triples.m,false,false);damping_coefficient.Resize(bending_triples.m,false,false);force_directions.Resize(bending_triples.m,false,false);
    for(TRIPLE_ITERATOR iterator(force_triples);iterator.Valid();iterator.Next()){int q=iterator.Data();
        int i,j,k;bending_triples(q).Get(i,j,k);
        TV n1=(particles.X(j)-particles.X(i)).Rotate_Clockwise_90(),n2=(particles.X(k)-particles.X(j)).Rotate_Clockwise_90();
        T length1=n1.Normalize(),length2=n2.Normalize();
        T sine_half_psi=Sine_Half_Angle_Between(n1,n2);
        elastic_s(q)=stiffness(q)/sqr(length_scale(q))*(sine_half_psi-sine_half_rest_angle(q));
        damping_coefficient(q)=damping(q)*stiffness(q)/cube(length_scale(q));
        T scale=1/max((T)1e-8,(T).5*(length1+length2)); // scale by length1*length2/length_average so that forces degrade gracefully when edges collapse
        force_directions(q).Set(scale*length2*n1,scale*length1*n2);}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> void SEGMENT_BENDING_ELEMENTS<T>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(TRIPLE_ITERATOR iterator(force_triples);iterator.Valid();iterator.Next()){int q=iterator.Data();
        int i,j,k;bending_triples(q).Get(i,j,k);
        TV n1,n2;force_directions(q).Get(n1,n2);TV dj=-n1-n2;
        T s=elastic_s(q);
        F(i)+=s*n1;F(j)+=s*dj;F(k)+=s*n2;}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T> void SEGMENT_BENDING_ELEMENTS<T>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(TRIPLE_ITERATOR iterator(force_triples);iterator.Valid();iterator.Next()){int q=iterator.Data();
        int i,j,k;bending_triples(q).Get(i,j,k);
        TV n1,n2;force_directions(q).Get(n1,n2);n1=n2=TV(0,-1);TV dj=-n1-n2;
        T s=-damping_coefficient(q)*(TV::Dot_Product(n1,V(i))+TV::Dot_Product(dj,V(j))+TV::Dot_Product(n2,V(k)));
        F(i)+=s*n1;F(j)+=s*dj;F(k)+=s*n2;}
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class T> T SEGMENT_BENDING_ELEMENTS<T>::
CFL_Strain_Rate() const
{
    T max_dtheta_dt=0;
    ARRAY_VIEW<const TV> V(particles.V);
    for(TRIPLE_ITERATOR iterator(force_triples);iterator.Valid();iterator.Next()){int q=iterator.Data();
        int i,j,k;bending_triples(q).Get(i,j,k);
        TV n1,n2;force_directions(q).Get(n1,n2);TV dj=-n1-n2;
        T dtheta_dt=TV::Dot_Product(n1,V(i))+TV::Dot_Product(dj,V(j))+TV::Dot_Product(n2,V(k)); // TODO: check scaling
        max_dtheta_dt=max(max_dtheta_dt,abs(dtheta_dt));}
    return Robust_Divide((T)pi*max_strain_per_time_step,max_dtheta_dt);
}
//#####################################################################
template class SEGMENT_BENDING_ELEMENTS<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SEGMENT_BENDING_ELEMENTS<double>;
#endif
