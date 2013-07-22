//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLE_WIND_DRAG_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
template<class T> const LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM<GRID<VECTOR<T,3> >,T> DEFORMABLE_WIND_DRAG_3D<T>::interpolation;
template<class T> const LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM<GRID<VECTOR<T,3> >,VECTOR<T,3> > DEFORMABLE_WIND_DRAG_3D<T>::vector_interpolation;
//#####################################################################
// Constructor
//#####################################################################
namespace{
template<class T> TRIANGULATED_SURFACE<T>& Triangulated_Surface(TRIANGULATED_SURFACE<T>& triangulated_surface)
{
    return triangulated_surface;
}
template<class T> TRIANGULATED_SURFACE<T>& Triangulated_Surface(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
{
    if(!tetrahedralized_volume.triangulated_surface) tetrahedralized_volume.Initialize_Triangulated_Surface();
    return *tetrahedralized_volume.triangulated_surface;
}}
template<class T> template<class T_OBJECT> DEFORMABLE_WIND_DRAG_3D<T>::
DEFORMABLE_WIND_DRAG_3D(T_OBJECT& object)
    :DEFORMABLES_FORCES<TV>(dynamic_cast<PARTICLES<TV>&>(object.particles)),triangulated_surface(&Triangulated_Surface(object)),use_constant_wind(false),
    use_spatially_varying_wind(false),spatially_varying_wind(0),wind_density(0),spatially_varying_wind_density(0),spatially_varying_wind_pressure(0),linear_normal_viscosity(0),mpi_solids(0)
{}
template<class T> DEFORMABLE_WIND_DRAG_3D<T>::
DEFORMABLE_WIND_DRAG_3D(PARTICLES<TV>& deformable_body_particles_input)
    :DEFORMABLES_FORCES<TV>(deformable_body_particles_input),triangulated_surface(0),use_constant_wind(false),
    use_spatially_varying_wind(false),spatially_varying_wind(0),wind_density(0),spatially_varying_wind_density(0),spatially_varying_wind_pressure(0),linear_normal_viscosity(0),mpi_solids(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> DEFORMABLE_WIND_DRAG_3D<T>::
~DEFORMABLE_WIND_DRAG_3D()
{}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class T> void DEFORMABLE_WIND_DRAG_3D<T>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    if(triangulated_surface) triangulated_surface->mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class T> void DEFORMABLE_WIND_DRAG_3D<T>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    force_elements.Update(triangulated_surface->mesh.elements);
    force_particles.Update(triangulated_surface->mesh.elements.Flattened());
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T> void DEFORMABLE_WIND_DRAG_3D<T>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    T_SIMPLICIAL_OBJECT& simplicial_object=*triangulated_surface;
    T wind_viscosity=use_constant_wind?constant_wind_viscosity:spatially_varying_wind_viscosity;
    if(wind_viscosity || spatially_varying_wind_pressure || wind_density || spatially_varying_wind_density){
        optimization.Resize(simplicial_object.mesh.elements.m,false,false);
        if(triangulated_surface) for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            int i,j,k;simplicial_object.mesh.elements(t).Get(i,j,k);
            const TV &Xi=particles.X(i),&Xj=particles.X(j),&Xk=particles.X(k);
            optimization(t).center=(T)one_third*(Xi+Xj+Xk);optimization(t).inward_normal=TRIANGLE_3D<T>::Normal_Direction(Xi,Xk,Xj);
            optimization(t).one_third_area=(T)one_sixth*optimization(t).inward_normal.Normalize();
            if(use_spatially_varying_wind) optimization(t).wind_velocity=Spatially_Varying_Wind_Velocity(optimization(t).center);}}
    if(linear_normal_viscosity && triangulated_surface){ // compute vertex normals for fragment
        vertex_normals.Resize(particles.array_collection->Size(),false,false);
        for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()) vertex_normals(iterator.Data())=TV();
        for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            const VECTOR<int,3>& nodes=simplicial_object.mesh.elements(t);
            vertex_normals.Subset(nodes)+=TRIANGLE_3D<T>::Normal_Direction(particles.X.Subset(nodes));}
        for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()) vertex_normals(iterator.Data()).Normalize();}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> VECTOR<T,3> DEFORMABLE_WIND_DRAG_3D<T>::
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
        if(normal_velocity>0) pressure=sqr(normal_velocity);else pressure=-sqr(normal_velocity);
        if(wind_density) pressure*=wind_density;
        else pressure*=Spatially_Varying_Wind_Density(optimization(t).center);}
    force+=pressure*optimization(t).inward_normal;
    force*=optimization(t).one_third_area;
    return force;
}

template<class T> void DEFORMABLE_WIND_DRAG_3D<T>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    if(use_constant_wind==use_spatially_varying_wind) PHYSBAM_FATAL_ERROR();
    T wind_viscosity=use_constant_wind?constant_wind_viscosity:spatially_varying_wind_viscosity;
    if(wind_viscosity || spatially_varying_wind_pressure || wind_density || spatially_varying_wind_density){
        if(!spatially_varying_wind_pressure && !wind_density && !spatially_varying_wind_density) PHYSBAM_FATAL_ERROR();
        if(triangulated_surface) for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            int i,j,k;triangulated_surface->mesh.elements(t).Get(i,j,k);
            TV wind_velocity=use_constant_wind?constant_wind:optimization(t).wind_velocity;
            TV relative_velocity=wind_velocity-(T)one_third*(particles.V(i)+particles.V(j)+particles.V(k));
            TV triangle_force=Add_Velocity_Independent_Forces_Helper(relative_velocity,t);
            F(i)+=triangle_force;F(j)+=triangle_force;F(k)+=triangle_force;}}
    if(linear_normal_viscosity && triangulated_surface){
        TV wind_acceleration=linear_normal_viscosity*constant_wind;
        for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
            if(use_constant_wind) F(p)+=particles.mass(p)*TV::Dot_Product(wind_acceleration,vertex_normals(p))*vertex_normals(p);
            else F(p)+=particles.mass(p)*linear_normal_viscosity*TV::Dot_Product(Spatially_Varying_Wind_Velocity(particles.X(p)),vertex_normals(p))*vertex_normals(p);}}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T> void DEFORMABLE_WIND_DRAG_3D<T>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    if(use_constant_wind==use_spatially_varying_wind) PHYSBAM_FATAL_ERROR();
    T wind_viscosity=use_constant_wind?constant_wind_viscosity:spatially_varying_wind_viscosity;
    if(wind_viscosity){
        if(triangulated_surface) for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            int i,j,k;triangulated_surface->mesh.elements(t).Get(i,j,k);
            // wind drag pressure - per unit area
            TV negative_V=-(T)one_third*(V(i)+V(j)+V(k));
            TV triangle_force=wind_viscosity*(negative_V-TV::Dot_Product(negative_V,optimization(t).inward_normal)*optimization(t).inward_normal);
            // total force
            triangle_force*=optimization(t).one_third_area; // one third of the triangle force is distriduted to each node
            F(i)+=triangle_force;F(j)+=triangle_force;F(k)+=triangle_force;}}
    if(linear_normal_viscosity && triangulated_surface){
        for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
            F(p)-=particles.mass(p)*linear_normal_viscosity*TV::Dot_Product(V(p),vertex_normals(p))*vertex_normals(p);}}
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template class DEFORMABLE_WIND_DRAG_3D<T>; \
    template DEFORMABLE_WIND_DRAG_3D<T>::DEFORMABLE_WIND_DRAG_3D(TRIANGULATED_SURFACE<T>&); \
    template DEFORMABLE_WIND_DRAG_3D<T>::DEFORMABLE_WIND_DRAG_3D(TETRAHEDRALIZED_VOLUME<T>&);
INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
