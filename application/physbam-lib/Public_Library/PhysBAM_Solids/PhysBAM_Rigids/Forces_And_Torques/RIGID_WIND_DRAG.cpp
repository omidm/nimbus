//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/SUMMATIONS.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGID_WIND_DRAG.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
using namespace PhysBAM;
//#####################################################################
template<class TV> const LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM<GRID<TV>,typename TV::SCALAR> RIGID_WIND_DRAG<TV>::interpolation;
template<class TV> const LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM<GRID<TV>,TV> RIGID_WIND_DRAG<TV>::vector_interpolation;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_WIND_DRAG<TV>::
RIGID_WIND_DRAG(RIGID_BODY<TV>& rigid_body_input)
    :RIGIDS_FORCES<TV>(rigid_body_input.rigid_body_collection),rigid_body(&rigid_body_input),use_constant_wind(false),
    use_spatially_varying_wind(false),spatially_varying_wind(0),wind_density(0),spatially_varying_wind_density(0),spatially_varying_wind_pressure(0),linear_normal_viscosity(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_WIND_DRAG<TV>::
~RIGID_WIND_DRAG()
{}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void RIGID_WIND_DRAG<TV>::
Update_Position_Based_State(const T time)
{
    T_SIMPLICIAL_OBJECT& simplicial_object=*rigid_body->simplicial_object;
    T wind_viscosity=use_constant_wind?constant_wind_viscosity:spatially_varying_wind_viscosity;
    if(wind_viscosity || spatially_varying_wind_pressure || wind_density || spatially_varying_wind_density){
        optimization.Resize(simplicial_object.mesh.elements.m,false,false);
        for(int t=1;t<=rigid_body->simplicial_object->mesh.elements.m;t++){
            if(rigid_body->Has_Infinite_Inertia() || rigid_body->particle_index<0) return;
            T_SIMPLEX world_space_simplex=rigid_body->World_Space_Simplex(t);
            optimization(t).center=world_space_simplex.Center(),optimization(t).inward_normal=-world_space_simplex.Normal();
            optimization(t).area_over_m=world_space_simplex.Size()/TV::m;
            if(use_spatially_varying_wind) optimization(t).wind_velocity=Spatially_Varying_Wind_Velocity(optimization(t).center);}}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> TV RIGID_WIND_DRAG<TV>::
Add_Velocity_Independent_Forces_Helper(TV relative_velocity,int t) const
{
    // wind drag pressure - per unit area
    T wind_viscosity=use_constant_wind?constant_wind_viscosity:spatially_varying_wind_viscosity;
    TV force=wind_viscosity*(relative_velocity-TV::Dot_Product(relative_velocity,optimization(t).inward_normal)*optimization(t).inward_normal);
    // wind pressure 
    T pressure;
    if(spatially_varying_wind_pressure) pressure=Spatially_Varying_Wind_Pressure(optimization(t).center);
    else{
        T normal_velocity=TV::Dot_Product(relative_velocity,optimization(t).inward_normal);
        surface_area+=abs(TV::Dot_Product(constant_wind,optimization(t).inward_normal))*optimization(t).area_over_m;
        if(normal_velocity>0) pressure=sqr(normal_velocity);else pressure=-sqr(normal_velocity);
        if(wind_density) pressure*=wind_density;
        else pressure*=Spatially_Varying_Wind_Density(optimization(t).center);}
    force+=pressure*optimization(t).inward_normal;
    force*=optimization(t).area_over_m;
    return force;
}

template<class TV> void RIGID_WIND_DRAG<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    surface_area=0;
    if(use_constant_wind==use_spatially_varying_wind) PHYSBAM_FATAL_ERROR();
    T wind_viscosity=use_constant_wind?constant_wind_viscosity:spatially_varying_wind_viscosity;
    if(wind_viscosity || spatially_varying_wind_pressure || wind_density || spatially_varying_wind_density){
        if(!spatially_varying_wind_pressure && !wind_density && !spatially_varying_wind_density) PHYSBAM_FATAL_ERROR();
            if(rigid_body->Has_Infinite_Inertia() || rigid_body->particle_index<0) return;
            TWIST<TV>& wrench=rigid_F(rigid_body->particle_index);
            for(int t=1;t<=rigid_body->simplicial_object->mesh.elements.m;t++){
                TV wind_velocity=use_constant_wind?constant_wind:optimization(t).wind_velocity;
                TV relative_velocity=wind_velocity-rigid_body->Pointwise_Object_Velocity(optimization(t).center);
                TV simplex_force=Add_Velocity_Independent_Forces_Helper(relative_velocity,t);
                wrench.linear+=simplex_force;wrench.angular+=TV::Cross_Product(optimization(t).center-rigid_body->X(),simplex_force);}}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void RIGID_WIND_DRAG<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    if(use_constant_wind==use_spatially_varying_wind) PHYSBAM_FATAL_ERROR();
    T wind_viscosity=use_constant_wind?constant_wind_viscosity:spatially_varying_wind_viscosity;
    if(wind_viscosity){
        if(rigid_body->Has_Infinite_Inertia() || rigid_body->particle_index<0) return;
        TWIST<TV>& wrench=rigid_F(rigid_body->particle_index);
        for(int t=1;t<=rigid_body->simplicial_object->mesh.elements.m;t++){
            TV negative_V=-rigid_body->Pointwise_Object_Velocity(optimization(t).center);
            TV simplex_force=wind_viscosity*(negative_V-TV::Dot_Product(negative_V,optimization(t).inward_normal)*optimization(t).inward_normal);
            // total force
            simplex_force*=optimization(t).area_over_m; // one third of the triangle force is distriduted to each node
            wrench.linear+=simplex_force;wrench.angular+=TV::Cross_Product(optimization(t).center-rigid_body->X(),simplex_force);}}
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template class RIGID_WIND_DRAG<VECTOR<T,3> >; \
    template class RIGID_WIND_DRAG<VECTOR<T,2> >;
INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
