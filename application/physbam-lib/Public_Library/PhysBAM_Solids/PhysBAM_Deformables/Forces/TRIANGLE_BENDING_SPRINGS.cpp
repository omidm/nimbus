//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Unnur Gretarsdottir, Geoffrey Irving, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/cyclic_shift.h>
#include <PhysBAM_Tools/Math_Tools/Robust_Arithmetic.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> TRIANGLE_BENDING_SPRINGS<T>::
TRIANGLE_BENDING_SPRINGS(PARTICLES<TV>& particles,TRIANGLE_MESH& triangle_mesh,const bool implicit)
    :LINEAR_SPRINGS<TV>(particles,bending_segment_mesh,implicit)
{
    Initialize(triangle_mesh);
    Set_Stiffness(0);Set_Damping(0);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> TRIANGLE_BENDING_SPRINGS<T>::
~TRIANGLE_BENDING_SPRINGS()
{
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class T> void TRIANGLE_BENDING_SPRINGS<T>::
Initialize(TRIANGLE_MESH& triangle_mesh)
{
    if(!triangle_mesh.adjacent_elements) triangle_mesh.Initialize_Adjacent_Elements();
    int number_quadruples=0;
    for(int t=1;t<=triangle_mesh.elements.m;t++) for(int a=1;a<=(*triangle_mesh.adjacent_elements)(t).m;a++) if((*triangle_mesh.adjacent_elements)(t)(a) > t) number_quadruples++;
    ARRAY<VECTOR<int,2> > segment_list(number_quadruples);
    int index=0; // reset number
    for(int t=1;t<=triangle_mesh.elements.m;t++){
        int t1,t2,t3;triangle_mesh.elements(t).Get(t1,t2,t3);
        for(int a=1;a<=(*triangle_mesh.adjacent_elements)(t).m;a++){
            int s=(*triangle_mesh.adjacent_elements)(t)(a);
            if(s > t){
                int s1,s2,s3;triangle_mesh.elements(s).Get(s1,s2,s3);
                if(t1==s1||t1==s2||t1==s3){cyclic_shift(t1,t2,t3);if(t1==s1||t1==s2||t1==s3)cyclic_shift(t1,t2,t3);}
                segment_list(++index).Set(t1,triangle_mesh.Other_Node(t2,t3,s));}}}
    segment_mesh.Initialize_Mesh(triangle_mesh.number_nodes,segment_list);
}
//#####################################################################
// Function Create_Bending_Springs
//#####################################################################
template<class T> TRIANGLE_BENDING_SPRINGS<T>* PhysBAM::
Create_Bending_Springs(PARTICLES<VECTOR<T,3> >& particles,TRIANGLE_MESH& triangle_mesh,const T stiffness,const T overdamping_fraction,
    const bool limit_time_step_by_strain_rate,const T max_strain_per_time_step,const bool use_rest_state_for_strain_rate,const T restlength_enlargement_fraction,
    const bool verbose,const bool implicit)
{
    TRIANGLE_BENDING_SPRINGS<T>* bend=new TRIANGLE_BENDING_SPRINGS<T>(particles,triangle_mesh,implicit);
    bend->Set_Restlength_From_Particles();
    bend->Set_Stiffness(stiffness);
    bend->Set_Overdamping_Fraction(overdamping_fraction);
    bend->Limit_Time_Step_By_Strain_Rate(limit_time_step_by_strain_rate,max_strain_per_time_step);
    bend->Use_Rest_State_For_Strain_Rate(use_rest_state_for_strain_rate);
    if(restlength_enlargement_fraction) bend->Clamp_Restlength_With_Fraction_Of_Springs(restlength_enlargement_fraction);
    if(verbose) bend->Print_Restlength_Statistics();
    bend->verbose=verbose;
    return bend;
}
//#####################################################################
// Function Create_Bending_Springs
//#####################################################################
template<class T> TRIANGLE_BENDING_SPRINGS<T>* PhysBAM::
Create_Bending_Springs(TRIANGULATED_SURFACE<T>& triangulated_surface,const T stiffness,const T overdamping_fraction,const bool limit_time_step_by_strain_rate,
    const T max_strain_per_time_step,const bool use_rest_state_for_strain_rate,const T restlength_enlargement_fraction,const bool verbose,const bool implicit)
{
    return Create_Bending_Springs(dynamic_cast<PARTICLES<VECTOR<T,3> >&>(triangulated_surface.particles),triangulated_surface.mesh,stiffness,overdamping_fraction,limit_time_step_by_strain_rate,max_strain_per_time_step,
        use_rest_state_for_strain_rate,restlength_enlargement_fraction,verbose,implicit);
}
//#####################################################################
template class TRIANGLE_BENDING_SPRINGS<float>;
template TRIANGLE_BENDING_SPRINGS<float>* PhysBAM::Create_Bending_Springs<float>(TRIANGULATED_SURFACE<float>&,float,float,bool,float,bool,float,bool,bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TRIANGLE_BENDING_SPRINGS<double>;
template TRIANGLE_BENDING_SPRINGS<double>* PhysBAM::Create_Bending_Springs<double>(TRIANGULATED_SURFACE<double>&,double,double,bool,double,bool,double,bool,bool);
#endif
