//#####################################################################
// Copyright 2002-2007, Zhaosheng Bao, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/cyclic_shift.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SIMPLE_TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> SIMPLE_TRIANGLE_BENDING_ELEMENTS<T>::
SIMPLE_TRIANGLE_BENDING_ELEMENTS(PARTICLES<TV>& particles,BINDING_LIST<TV>& binding_list_input,bool implicit)
    :LINEAR_SPRINGS<TV>(particles,spring_connectivity,implicit),
    bending_quadruples(bending_quadruples_default),binding_list(binding_list_input)
{
    Print_Number_Ignored();
    Limit_Time_Step_By_Strain_Rate(false);
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> SIMPLE_TRIANGLE_BENDING_ELEMENTS<T>::
SIMPLE_TRIANGLE_BENDING_ELEMENTS(PARTICLES<TV>& particles,ARRAY<VECTOR<int,4> >& bending_quadruples_input,BINDING_LIST<TV>& binding_list_input,bool implicit)
    :LINEAR_SPRINGS<TV>(particles,spring_connectivity,implicit),bending_quadruples(bending_quadruples_input),bending_stiffness(bending_quadruples_input.m),
    damping(bending_quadruples_input.m),binding_list(binding_list_input)
{
    Print_Number_Ignored();
    Limit_Time_Step_By_Strain_Rate(false);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> SIMPLE_TRIANGLE_BENDING_ELEMENTS<T>::
~SIMPLE_TRIANGLE_BENDING_ELEMENTS()
{
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class T> void SIMPLE_TRIANGLE_BENDING_ELEMENTS<T>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    for(int q=1;q<=bending_quadruples.m;q++)
        for(int i=1;i<=3;i++) for(int j=i+1;j<=4;j++) dependency_mesh.Add_Element_If_Not_Already_There(VECTOR<int,2>(bending_quadruples(q)[i],bending_quadruples(q)[j]));
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class T> void SIMPLE_TRIANGLE_BENDING_ELEMENTS<T>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    force_quadruples.Update(bending_quadruples,particle_is_simulated);
    force_segments.Update(spring_connectivity.elements);
}
//#####################################################################
// Function Set_Quadruples_From_Triangle_Mesh
//#####################################################################
template<class T> void SIMPLE_TRIANGLE_BENDING_ELEMENTS<T>::
Set_Quadruples_From_Triangle_Mesh(TRIANGLE_MESH& mesh)
{
    if(!mesh.adjacent_elements) mesh.Initialize_Adjacent_Elements();

    // allocate proper array sizes
    int number_quadruples=0;
    for(int t=1;t<=mesh.elements.m;t++) for(int a=1;a<=(*mesh.adjacent_elements)(t).m;a++) if((*mesh.adjacent_elements)(t)(a)>t) number_quadruples++;
    bending_quadruples.Resize(number_quadruples);bending_stiffness.Resize(number_quadruples);
    damping.Resize(number_quadruples);

    int index=0; // reset number
    for(int t=1;t<=mesh.elements.m;t++){
        int t1,t2,t3;mesh.elements(t).Get(t1,t2,t3);
        for(int a=1;a<=(*mesh.adjacent_elements)(t).m;a++){
            int s=(*mesh.adjacent_elements)(t)(a);
            if(s>t){
                int s1,s2,s3;mesh.elements(s).Get(s1,s2,s3);
                if(t1==s1 || t1==s2 || t1==s3){cyclic_shift(t1,t2,t3);if(t1==s1 || t1==s2 || t1==s3) cyclic_shift(t1,t2,t3);}
                bending_quadruples(++index).Set(t1,t2,t3,mesh.Other_Node(t2,t3,s));}}}

    linear_bindings.Resize(number_quadruples);
    ARRAY_VIEW<const TV> X(particles.X);
    for(int q=1;q<=bending_quadruples.m;q++){
        int i,j,k,l;bending_quadruples(q).Get(i,j,k,l);
        //TV xij=X(i)-X(j),xik=X(i)-X(k),xlk=X(l)-X(k),xlj=X(l)-X(j);
        //SEGMENT_3D<T> shared_edge(X(j),X(k));
        //SEGMENT_3D<T> fictitious_edge(X(i),X(l));
        VECTOR<T,2> weights;
        //TV direction=shared_edge.Shortest_Vector_Between_Segments(fictitious_edge,weights); // points from fictitious to shared
        // weights.x corresponds to weights on fictitious
        // shared edge
        int shared_edge_particle=particles.array_collection->Add_Element();
        //linear_bindings(q).x=binding_list.Add_Binding(new LINEAR_BINDING<TV,2>(particles,shared_edge_particle,VECTOR<int,2>(j,k),VECTOR<T,2>(weights.y,1-weights.y)));
        linear_bindings(q).x=binding_list.Add_Binding(new LINEAR_BINDING<TV,2>(particles,shared_edge_particle,VECTOR<int,2>(j,k),VECTOR<T,2>()));
        // fictitious edge
        int fictitious_edge_particle=particles.array_collection->Add_Element();
        //linear_bindings(q).y=binding_list.Add_Binding(new LINEAR_BINDING<TV,2>(particles,fictitious_edge_particle,VECTOR<int,2>(i,l),VECTOR<T,2>(weights.x,1-weights.x)));
        linear_bindings(q).y=binding_list.Add_Binding(new LINEAR_BINDING<TV,2>(particles,fictitious_edge_particle,VECTOR<int,2>(i,l),VECTOR<T,2>()));
        spring_connectivity.elements.Append(VECTOR<int,2>(shared_edge_particle,fictitious_edge_particle));
    }
    spring_connectivity.Set_Number_Nodes(particles.array_collection->Size());
}
//#####################################################################
// Function Set_Constants_From_Particles
//#####################################################################
template<class T> void SIMPLE_TRIANGLE_BENDING_ELEMENTS<T>::
Set_Constants_From_Particles(const T material_stiffness,const T material_damping)
{
    fictitious_edge_length.Resize(damping.m);
    for(int q=1;q<=bending_quadruples.m;q++){
        int i,j,k,l;bending_quadruples(q).Get(i,j,k,l);
        TV n1=TV::Cross_Product(particles.X(i)-particles.X(j),particles.X(i)-particles.X(k)).Normalized(),
            n2=TV::Cross_Product(particles.X(l)-particles.X(k),particles.X(l)-particles.X(j)).Normalized(),e=particles.X(k)-particles.X(j);
        bending_stiffness(q)=material_stiffness;damping(q)=material_damping;
        //fictitious_edge_length(q)=(particles.X(i)-particles.X(l)).Magnitude();
        fictitious_edge_length(q)=(particles.X(i)-particles.X(k)).Magnitude()+(particles.X(i)-particles.X(j)).Magnitude()+(particles.X(k)-particles.X(l)).Magnitude()+(particles.X(j)-particles.X(l)).Magnitude();
    }
    restlength.Resize(damping.m);visual_restlength.Resize(damping.m);
    ARRAYS_COMPUTATIONS::Fill(restlength,(T)1);
    ARRAYS_COMPUTATIONS::Fill(visual_restlength,(T)0);
    youngs_modulus.Resize(damping.m);base_youngs_modulus.Resize(damping.m);
    ARRAYS_COMPUTATIONS::Fill(youngs_modulus,material_stiffness);
    ARRAYS_COMPUTATIONS::Fill(base_youngs_modulus,material_stiffness);
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T> void SIMPLE_TRIANGLE_BENDING_ELEMENTS<T>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    damping_coefficient.Resize(bending_quadruples.m,false,false);
    states.Resize(spring_connectivity.elements.m,false,false);current_lengths.Resize(spring_connectivity.elements.m,false,false);
    //BASE::Set_Stiffness(material_stiffness);

    int ignored_elements=0,total_elements=0;
    ARRAY_VIEW<const TV> X(particles.X);
    for(QUADRUPLE_ITERATOR iterator(force_quadruples);iterator.Valid();iterator.Next()){int q=iterator.Data();
        VECTOR<T,2> weights;
        int shared_edge_binding_index,fictitious_edge_binding_index;linear_bindings(q).Get(shared_edge_binding_index,fictitious_edge_binding_index);
        LINEAR_BINDING<TV,2>* shared_edge_binding=dynamic_cast<LINEAR_BINDING<TV,2>*>(binding_list.bindings(shared_edge_binding_index));
        LINEAR_BINDING<TV,2>* fictitious_edge_binding=dynamic_cast<LINEAR_BINDING<TV,2>*>(binding_list.bindings(fictitious_edge_binding_index));
        int i=fictitious_edge_binding->parents.x,l=fictitious_edge_binding->parents.y;
        int j=shared_edge_binding->parents.x,k=shared_edge_binding->parents.y;
        total_elements++;
        damping_coefficient(q)=damping(q);

        SEGMENT_3D<T> shared_edge(X(j),X(k));
        SEGMENT_3D<T> fictitious_edge(X(i),X(l));
        
        TV direction=shared_edge.Shortest_Vector_Between_Segments(fictitious_edge,weights); // points from fictitious to shared
        shared_edge_binding->weights.x=1-weights.x;
        shared_edge_binding->weights.y=weights.x;
        fictitious_edge_binding->weights.x=1-weights.y;
        fictitious_edge_binding->weights.y=weights.y;

        //T current_fictitious_edge_length=fictitious_edge.Length();
        T current_fictitious_edge_length=(particles.X(i)-particles.X(k)).Magnitude()+(particles.X(i)-particles.X(j)).Magnitude()+(particles.X(k)-particles.X(l)).Magnitude()+(particles.X(j)-particles.X(l)).Magnitude();
        youngs_modulus(q)=base_youngs_modulus(q)*sqr(current_fictitious_edge_length/fictitious_edge_length(q)-1);

        typename BASE::STATE& state=states(q);
        state.nodes=VECTOR<int,2>(shared_edge_binding->particle_index,fictitious_edge_binding->particle_index);
        state.direction=direction;
        current_lengths(q)=state.direction.Normalize();
        state.coefficient=sqr(current_fictitious_edge_length/fictitious_edge_length(q)-1)*damping_coefficient(q)/restlength(q);}
    binding_list.Clamp_Particles_To_Embedded_Positions();
    binding_list.Clamp_Particles_To_Embedded_Velocities();

    if(print_number_ignored){ std::stringstream ss;ss<<"ignored "<<ignored_elements<<" of "<<total_elements<<" bending elements"<<std::endl;LOG::filecout(ss.str());}
}
template<class T> SIMPLE_TRIANGLE_BENDING_ELEMENTS<T>* PhysBAM::
Create_Simple_Bending_Elements(PARTICLES<VECTOR<T,3> >& particles,TRIANGLE_MESH& mesh,BINDING_LIST<VECTOR<T,3> >& binding_list,const T stiffness,const T damping,
    const bool limit_time_step_by_strain_rate,const T max_strain_per_time_step,const bool use_plasticity,const T plastic_yield,const T plastic_hardening,
    const T cutoff_fraction_of_minimum_area,const T cutoff_fraction_of_triangles,const bool verbose,const bool implicit) //see header for defaults
{
    SIMPLE_TRIANGLE_BENDING_ELEMENTS<T>* bend=new SIMPLE_TRIANGLE_BENDING_ELEMENTS<T>(particles,binding_list,implicit);
    bend->Set_Quadruples_From_Triangle_Mesh(mesh);
    bend->Set_Constants_From_Particles(stiffness,damping);
    bend->Limit_Time_Step_By_Strain_Rate(limit_time_step_by_strain_rate,max_strain_per_time_step);
    return bend;
}

template<class T> SIMPLE_TRIANGLE_BENDING_ELEMENTS<T>* PhysBAM::
Create_Simple_Bending_Elements(TRIANGULATED_SURFACE<T>& triangulated_surface,BINDING_LIST<VECTOR<T,3> >& binding_list,const T stiffness,
    const T damping,const bool limit_time_step_by_strain_rate,const T max_strain_per_time_step,const bool use_plasticity,const T plastic_yield,const T plastic_hardening,
    const T cutoff_fraction_of_minimum_area,const T cutoff_fraction_of_triangles,const bool verbose,const bool implicit) //see header for defaults
{
    return Create_Simple_Bending_Elements(dynamic_cast<PARTICLES<VECTOR<T,3> >&>(triangulated_surface.particles),triangulated_surface.mesh,binding_list,stiffness,damping,limit_time_step_by_strain_rate,max_strain_per_time_step,use_plasticity,
        plastic_yield,plastic_hardening,cutoff_fraction_of_minimum_area,cutoff_fraction_of_triangles,verbose,implicit);
}
//#####################################################################
template class SIMPLE_TRIANGLE_BENDING_ELEMENTS<float>;
template SIMPLE_TRIANGLE_BENDING_ELEMENTS<float>* PhysBAM::Create_Simple_Bending_Elements<float>(TRIANGULATED_SURFACE<float>&,BINDING_LIST<VECTOR<float,3> >&,float,float,bool,float,bool,
    float,float,float,float,bool,bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SIMPLE_TRIANGLE_BENDING_ELEMENTS<double>;
template SIMPLE_TRIANGLE_BENDING_ELEMENTS<double>* PhysBAM::Create_Simple_Bending_Elements<double>(TRIANGULATED_SURFACE<double>&,BINDING_LIST<VECTOR<double,3> >&,double,double,bool,double,bool,
    double,double,double,double,bool,bool);
#endif
