//#####################################################################
// Copyright 2002-2007, Zhaosheng Bao, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/cyclic_shift.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class T> void TRIANGLE_BENDING_ELEMENTS<T>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    for(int q=1;q<=bending_quadruples.m;q++)
        for(int i=1;i<=3;i++) for(int j=i+1;j<=4;j++) dependency_mesh.Add_Element_If_Not_Already_There(VECTOR<int,2>(bending_quadruples(q)[i],bending_quadruples(q)[j]));
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class T> void TRIANGLE_BENDING_ELEMENTS<T>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    force_quadruples.Update(bending_quadruples,particle_is_simulated);
}
//#####################################################################
// Function Set_Area_Cutoff_With_Fraction_Of_Triangles
//#####################################################################
template<class T> void TRIANGLE_BENDING_ELEMENTS<T>::
Set_Area_Cutoff_With_Fraction_Of_Triangles(TRIANGULATED_SURFACE<T>& triangulated_surface,const T fraction)
{
    ARRAY<VECTOR<int,3> >& elements=triangulated_surface.mesh.elements;
    ARRAY<T> area(elements.m);for(int t=1;t<=elements.m;t++) area(t)=triangulated_surface.Area(t);Sort(area);
    area_cutoff=area(min((int)(fraction*elements.m)+1,area.m));
}
//#####################################################################
// Function Set_Quadruples_From_Triangle_Mesh
//#####################################################################
template<class T> void TRIANGLE_BENDING_ELEMENTS<T>::
Set_Quadruples_From_Triangle_Mesh(TRIANGLE_MESH& mesh)
{
    if(!mesh.adjacent_elements) mesh.Initialize_Adjacent_Elements();

    // allocate proper array sizes
    int number_quadruples=0;
    for(int t=1;t<=mesh.elements.m;t++) for(int a=1;a<=(*mesh.adjacent_elements)(t).m;a++) if((*mesh.adjacent_elements)(t)(a)>t) number_quadruples++;
    bending_quadruples.Resize(number_quadruples);bending_stiffness.Resize(number_quadruples);
    sine_half_rest_angle.Resize(number_quadruples);damping.Resize(number_quadruples);

    int index=0; // reset number
    for(int t=1;t<=mesh.elements.m;t++){
        int t1,t2,t3;mesh.elements(t).Get(t1,t2,t3);
        for(int a=1;a<=(*mesh.adjacent_elements)(t).m;a++){
            int s=(*mesh.adjacent_elements)(t)(a);
            if(s>t){
                int s1,s2,s3;mesh.elements(s).Get(s1,s2,s3);
                if(t1==s1 || t1==s2 || t1==s3){cyclic_shift(t1,t2,t3);if(t1==s1 || t1==s2 || t1==s3) cyclic_shift(t1,t2,t3);}
                bending_quadruples(++index).Set(t1,t2,t3,mesh.Other_Node(t2,t3,s));}}}
}
//#####################################################################
// Function Set_Quadruples_From_Reference_Triangle_Mesh
//#####################################################################
template<class T> void TRIANGLE_BENDING_ELEMENTS<T>::
Set_Quadruples_From_Reference_Triangle_Mesh(TRIANGLE_MESH& mesh,const ARRAY<int>& triangle_map_to_reference)
{
    if(!mesh.adjacent_elements) mesh.Initialize_Adjacent_Elements();

    // allocate proper array sizes
    int number_quadruples=0;
    for(int t=1;t<=mesh.elements.m;t++) for(int a=1;a<=(*mesh.adjacent_elements)(t).m;a++)
        if(triangle_map_to_reference((*mesh.adjacent_elements)(t)(a))>triangle_map_to_reference(t)) number_quadruples++;
    bending_quadruples.Resize(number_quadruples);bending_stiffness.Resize(number_quadruples);
    sine_half_rest_angle.Resize(number_quadruples);damping.Resize(number_quadruples);

    int index=0; // reset number
    for(int t=1;t<=mesh.elements.m;t++){
        int t1,t2,t3;mesh.elements(t).Get(t1,t2,t3);
        for(int a=1;a<=(*mesh.adjacent_elements)(t).m;a++){
            int s=(*mesh.adjacent_elements)(t)(a);
            if(triangle_map_to_reference(s)>triangle_map_to_reference(t)){
                int s1,s2,s3;mesh.elements(s).Get(s1,s2,s3);
                if(t1==s1 || t1==s2 || t1==s3){cyclic_shift(t1,t2,t3);if(t1==s1 || t1==s2 || t1==s3) cyclic_shift(t1,t2,t3);}
                bending_quadruples(++index).Set(t1,t2,t3,mesh.Other_Node(t2,t3,s));}}}
}
//#####################################################################
// Function Set_Constants_From_Particles
//#####################################################################
template<class T> void TRIANGLE_BENDING_ELEMENTS<T>::
Set_Constants_From_Particles(const T material_stiffness,const T material_damping)
{
    for(int q=1;q<=bending_quadruples.m;q++){
        int i,j,k,l;bending_quadruples(q).Get(i,j,k,l);
        TV n1=TV::Cross_Product(particles.X(i)-particles.X(j),particles.X(i)-particles.X(k)).Normalized(),
            n2=TV::Cross_Product(particles.X(l)-particles.X(k),particles.X(l)-particles.X(j)).Normalized(),e=particles.X(k)-particles.X(j);
        sine_half_rest_angle(q)=Sine_Half_Angle_Between(n1,n2,e);
        bending_stiffness(q)=material_stiffness;damping(q)=material_damping;}
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T> void TRIANGLE_BENDING_ELEMENTS<T>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    elastic_s.Resize(bending_quadruples.m,false,false);damping_coefficient.Resize(bending_quadruples.m,false,false);force_directions.Resize(bending_quadruples.m,false,false);

    int ignored_elements=0,total_elements=0;
    ARRAY_VIEW<const TV> X(particles.X);
    T twice_area_cutoff_squared=sqr(2*area_cutoff);
    for(QUADRUPLE_ITERATOR iterator(force_quadruples);iterator.Valid();iterator.Next()){int q=iterator.Data();
        int i,j,k,l;bending_quadruples(q).Get(i,j,k,l);
        TV xij=X(i)-X(j),xik=X(i)-X(k),xlk=X(l)-X(k),xlj=X(l)-X(j),n1=TV::Cross_Product(xij,xik),n2=TV::Cross_Product(xlk,xlj);
        T n1_magnitude_squared=n1.Magnitude_Squared(),n2_magnitude_squared=n2.Magnitude_Squared();
        total_elements++;
        if(area_cutoff && (n1_magnitude_squared<twice_area_cutoff_squared || n2_magnitude_squared<twice_area_cutoff_squared)){
            elastic_s(q)=damping_coefficient(q)=0;ignored_elements++;continue;}
        T n1_magnitude=sqrt(n1_magnitude_squared),n2_magnitude=sqrt(n2_magnitude_squared);
        TV e=X(k)-X(j);T e_magnitude_squared=e.Magnitude_Squared(),e_magnitude=sqrt(e_magnitude_squared),e_magnitude_cubed=e_magnitude_squared*e_magnitude;
        T sine_half_psi=Sine_Half_Angle_Between(n1/n1_magnitude,n2/n2_magnitude,e);
        T one_over_normal_magnitude=1/(n1_magnitude+n2_magnitude),strain;
        if(plastic_yield){
            strain=sine_half_psi-(*sine_half_elastic_angle)(q);
            T scale=e_magnitude*one_over_normal_magnitude,stress=scale*strain,yield=(*plastic_yield)(q);if(stress<0){stress=-stress;scale=-scale;}
            if(stress>yield){
                strain=yield/scale;(*sine_half_elastic_angle)(q)=sine_half_psi-strain;
                (*plastic_yield)(q)+=(*plastic_hardening)(q)*(stress-yield);}}
        else strain=sine_half_psi-sine_half_rest_angle(q);
        elastic_s(q)=bending_stiffness(q)*e_magnitude_cubed*one_over_normal_magnitude*strain;
        damping_coefficient(q)=-damping(q)*e_magnitude_cubed;
        n1/=n1_magnitude_squared;n2/=n2_magnitude_squared;e/=e_magnitude_squared;
        force_directions(q).Set(n1,TV::Dot_Product(xik,e)*n1+TV::Dot_Product(xlk,e)*n2,-TV::Dot_Product(xij,e)*n1-TV::Dot_Product(xlj,e)*n2,n2);}

    if(print_number_ignored) {std::stringstream ss;ss<<"ignored "<<ignored_elements<<" of "<<total_elements<<" bending elements"<<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> void TRIANGLE_BENDING_ELEMENTS<T>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(QUADRUPLE_ITERATOR iterator(force_quadruples);iterator.Valid();iterator.Next()){int q=iterator.Data();
        int i,j,k,l;bending_quadruples(q).Get(i,j,k,l);
        TV n1,dj,dk,n2;force_directions(q).Get(n1,dj,dk,n2);
        T s=elastic_s(q);
        F(i)+=s*n1;F(j)+=s*dj;F(k)+=s*dk;F(l)+=s*n2;}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T> void TRIANGLE_BENDING_ELEMENTS<T>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(QUADRUPLE_ITERATOR iterator(force_quadruples);iterator.Valid();iterator.Next()){int q=iterator.Data();
        int i,j,k,l;bending_quadruples(q).Get(i,j,k,l);
        TV n1,dj,dk,n2;force_directions(q).Get(n1,dj,dk,n2);
        T s=damping_coefficient(q)*(TV::Dot_Product(n1,V(i))+TV::Dot_Product(dj,V(j))+TV::Dot_Product(dk,V(k))+TV::Dot_Product(n2,V(l)));
        F(i)+=s*n1;F(j)+=s*dj;F(k)+=s*dk;F(l)+=s*n2;}
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class T> T TRIANGLE_BENDING_ELEMENTS<T>::
CFL_Strain_Rate() const
{
    T max_dtheta_dt=0,dtheta_dt;
    ARRAY_VIEW<const TV> X(particles.X),V(particles.V);
    for(QUADRUPLE_ITERATOR iterator(force_quadruples);iterator.Valid();iterator.Next()){int q=iterator.Data();
        int i,j,k,l;bending_quadruples(q).Get(i,j,k,l);
        TV n1,dj,dk,n2;force_directions(q).Get(n1,dj,dk,n2);
        TV e=X(k)-X(j);
        dtheta_dt=e.Magnitude()*(TV::Dot_Product(n1,V(i))+TV::Dot_Product(dj,V(j))+TV::Dot_Product(dk,V(k))+TV::Dot_Product(n2,V(l)));
        max_dtheta_dt=max(max_dtheta_dt,abs(dtheta_dt));}
    T dt=Robust_Divide((T)pi*max_strain_per_time_step,max_dtheta_dt);
    std::stringstream ss;ss<<"max dtheta_dt = "<<max_dtheta_dt<<std::endl;
    ss<<"bending CFL strain rate = "<<dt<<std::endl;LOG::filecout(ss.str());
    return dt;
}
//#####################################################################
// Function Compute_Discrete_Shell_Energy
//#####################################################################
// Use formula in Discrete Shells paper for integral of curvature squared
template<class T> T TRIANGLE_BENDING_ELEMENTS<T>::
Compute_Discrete_Shell_Energy()
{
    T energy=0;
    for(int q=1;q<=bending_quadruples.m;q++){
        int i,j,k,l;bending_quadruples(q).Get(i,j,k,l);
        TV n1=TV::Cross_Product(particles.X(i)-particles.X(j),particles.X(i)-particles.X(k)),
            n2=TV::Cross_Product(particles.X(l)-particles.X(k),particles.X(l)-particles.X(j)),
            e=particles.X(k)-particles.X(j);
        n1.Normalize();n2.Normalize();T e_magnitude=e.Magnitude();e/=e_magnitude;
        T theta=atan2(TV::Triple_Product(n1,n2,e),TV::Dot_Product(n1,n2));
        energy+=e_magnitude*sqr(theta);}
    return energy;
}
//#####################################################################
// Function Initialize_Reference_Quantities
//#####################################################################
template<class T> void TRIANGLE_BENDING_ELEMENTS<T>::
Initialize_Reference_Quantities(const int hash_multiple)
{
    delete reference_bending_quadruples_hashtable;
    reference_bending_quadruples_hashtable=new HASHTABLE<VECTOR<int,2>,int>(hash_multiple*bending_quadruples.m);
    for(int q=1;q<=bending_quadruples.m;q++) reference_bending_quadruples_hashtable->Insert(VECTOR<int,2>(bending_quadruples(q)(2),bending_quadruples(q)(3)),q);
    reference_sine_half_rest_angle=new ARRAY<T>(sine_half_rest_angle);
    reference_bending_stiffness=new ARRAY<T>(bending_stiffness);
    reference_damping=new ARRAY<T>(damping);
    if(plastic_hardening) reference_plastic_hardening=new ARRAY<T>(*plastic_hardening);
}
//#####################################################################
// Function Copy_Back_Reference_Quantities
//#####################################################################
template<class T> void TRIANGLE_BENDING_ELEMENTS<T>::
Copy_Back_Reference_Quantities(const ARRAY<int>& node_map_to_reference)
{
    sine_half_rest_angle.Resize(bending_quadruples.m,false,false);
    damping.Resize(bending_quadruples.m,false,false);
    if(plastic_hardening) plastic_hardening->Resize(bending_quadruples.m,false,false);
    for(int q=1;q<=bending_quadruples.m;q++){
        int reference_j=node_map_to_reference(bending_quadruples(q)(2));
        int reference_k=node_map_to_reference(bending_quadruples(q)(3));
        int q_reference=0;reference_bending_quadruples_hashtable->Get(VECTOR<int,2>(reference_j,reference_k),q_reference);assert(q_reference);
        sine_half_rest_angle(q)=(*reference_sine_half_rest_angle)(q_reference);
        bending_stiffness(q)=(*reference_bending_stiffness)(q_reference);
        damping(q)=(*reference_damping)(q_reference);
        if(plastic_hardening) (*plastic_hardening)(q)=(*reference_plastic_hardening)(q_reference);}
}
//#####################################################################
// Function Initialize_Save_Quantities
//#####################################################################
template<class T> void TRIANGLE_BENDING_ELEMENTS<T>::
Initialize_Save_Quantities()
{
    delete plastic_yield_save;plastic_yield_save=new ARRAY<T>(*plastic_yield);
    delete sine_half_elastic_angle_save;sine_half_elastic_angle_save=new ARRAY<T>(*sine_half_elastic_angle);
    delete bending_quadruples_save;bending_quadruples_save=new ARRAY<VECTOR<int,4> >(bending_quadruples);
}
//#####################################################################
// Function Copy_Back_Save_Quantities
//#####################################################################
template<class T> void TRIANGLE_BENDING_ELEMENTS<T>::
Copy_Back_Save_Quantities(const ARRAY<int>& node_map_to_saved)
{
    HASHTABLE<VECTOR<int,2>,int> save_bending_quadruples_hashtable(2*bending_quadruples_save->m);
    for(int q=1;q<=bending_quadruples_save->m;q++) save_bending_quadruples_hashtable.Insert(VECTOR<int,2>((*bending_quadruples_save)(q)(2),(*bending_quadruples_save)(q)(3)),q);
    plastic_yield->Resize(bending_quadruples.m,false,false);
    sine_half_elastic_angle->Resize(bending_quadruples.m,false,false);
    for(int q=1;q<=bending_quadruples.m;q++){
        int save_j=node_map_to_saved(bending_quadruples(q)(2)),save_k=node_map_to_saved(bending_quadruples(q)(3));
        int q_save=0;save_bending_quadruples_hashtable.Get(VECTOR<int,2>(save_j,save_k),q_save);assert(q_save);
        (*plastic_yield)(q)=(*plastic_yield_save)(q_save);
        (*sine_half_elastic_angle)(q)=(*sine_half_elastic_angle_save)(q_save);}
}
//#####################################################################
template class TRIANGLE_BENDING_ELEMENTS<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TRIANGLE_BENDING_ELEMENTS<double>;
#endif
