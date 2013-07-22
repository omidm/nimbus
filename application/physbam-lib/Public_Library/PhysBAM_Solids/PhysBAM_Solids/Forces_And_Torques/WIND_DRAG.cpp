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
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/WIND_DRAG.h>
using namespace PhysBAM;
//#####################################################################
template<class TV> const LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM<GRID<TV>,typename TV::SCALAR> WIND_DRAG<TV>::interpolation;
template<class TV> const LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM<GRID<TV>,TV> WIND_DRAG<TV>::vector_interpolation;
//#####################################################################
// Constructor
//#####################################################################
namespace{
template<class TV> typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::dimension-1>::OBJECT& Simplicial_Object(typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::dimension-1>::OBJECT& simplicial_object)
{
    return simplicial_object;
}
template<class TV> typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::dimension-1>::OBJECT& Simplicial_Object(TETRAHEDRALIZED_VOLUME<typename TV::SCALAR>& tetrahedralized_volume)
{
    if(!tetrahedralized_volume.triangulated_surface) tetrahedralized_volume.Initialize_Triangulated_Surface();
    return *tetrahedralized_volume.triangulated_surface;
}}
template<class TV> template<class T_OBJECT> WIND_DRAG<TV>::
WIND_DRAG(T_OBJECT& object,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input)
    :SOLIDS_FORCES<TV>(dynamic_cast<PARTICLES<TV>&>(object.particles),rigid_body_collection_input),deformable_simplicial_object(&Simplicial_Object<TV>(object)),rigid_body(0),use_constant_wind(false),
    use_spatially_varying_wind(false),spatially_varying_wind(0),wind_density(0),spatially_varying_wind_density(0),spatially_varying_wind_pressure(0),linear_normal_viscosity(0),mpi_solids(0)
{}
template<class TV> WIND_DRAG<TV>::
WIND_DRAG(RIGID_BODY<TV>& rigid_body_input,PARTICLES<TV>& deformable_body_particles_input)
    :SOLIDS_FORCES<TV>(deformable_body_particles_input,rigid_body_input.rigid_body_collection),deformable_simplicial_object(0),rigid_body(&rigid_body_input),use_constant_wind(false),
    use_spatially_varying_wind(false),spatially_varying_wind(0),wind_density(0),spatially_varying_wind_density(0),spatially_varying_wind_pressure(0),linear_normal_viscosity(0),mpi_solids(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> WIND_DRAG<TV>::
~WIND_DRAG()
{}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV> void WIND_DRAG<TV>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    if(deformable_simplicial_object) deformable_simplicial_object->mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void WIND_DRAG<TV>::
Update_Position_Based_State(const T time)
{
    T_SIMPLICIAL_OBJECT& simplicial_object=deformable_simplicial_object?*deformable_simplicial_object:*rigid_body->simplicial_object;
    T wind_viscosity=use_constant_wind?constant_wind_viscosity:spatially_varying_wind_viscosity;
    if(wind_viscosity || spatially_varying_wind_pressure || wind_density || spatially_varying_wind_density){
        optimization.Resize(simplicial_object.mesh.elements.m,false,false);
        if(deformable_simplicial_object) for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            const TV_INT& nodes=simplicial_object.mesh.elements(t);
            optimization(t).center=ARRAYS_COMPUTATIONS::Sum(particles.X.Subset(nodes))/TV::m;optimization(t).inward_normal=T_SIMPLEX::Normal(particles.X.Subset(nodes));
            optimization(t).area_over_m=T_SIMPLEX::Size(particles.X.Subset(nodes))/TV::m;
            if(use_spatially_varying_wind) optimization(t).wind_velocity=Spatially_Varying_Wind_Velocity(optimization(t).center);}
        else for(int t=1;t<=rigid_body->simplicial_object->mesh.elements.m;t++){
            if(rigid_body->Has_Infinite_Inertia() || rigid_body->particle_index<0) return;
            T_SIMPLEX world_space_simplex=rigid_body->World_Space_Simplex(t);
            optimization(t).center=world_space_simplex.Center(),optimization(t).inward_normal=-world_space_simplex.Normal();
            optimization(t).area_over_m=world_space_simplex.Size()/TV::m;
            if(use_spatially_varying_wind) optimization(t).wind_velocity=Spatially_Varying_Wind_Velocity(optimization(t).center);}}
    if(linear_normal_viscosity && deformable_simplicial_object){ // compute vertex normals for fragment
        vertex_normals.Resize(particles.array_collection->Size(),false,false);
        for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()) vertex_normals(iterator.Data())=TV();
        for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            const TV_INT& nodes=simplicial_object.mesh.elements(t);
            vertex_normals.Subset(nodes)+=T_SIMPLEX::Normal(particles.X.Subset(nodes));}
        for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()) vertex_normals(iterator.Data()).Normalize();}
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> TV WIND_DRAG<TV>::
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

template<class TV> void WIND_DRAG<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    surface_area=0;
    if(use_constant_wind==use_spatially_varying_wind) PHYSBAM_FATAL_ERROR();
    T wind_viscosity=use_constant_wind?constant_wind_viscosity:spatially_varying_wind_viscosity;
    if(wind_viscosity || spatially_varying_wind_pressure || wind_density || spatially_varying_wind_density){
        if(!spatially_varying_wind_pressure && !wind_density && !spatially_varying_wind_density) PHYSBAM_FATAL_ERROR();
        if(deformable_simplicial_object) for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            const TV_INT& nodes=deformable_simplicial_object->mesh.elements(t);
            TV wind_velocity=use_constant_wind?constant_wind:optimization(t).wind_velocity;
            TV relative_velocity=wind_velocity-ARRAYS_COMPUTATIONS::Sum(particles.V.Subset(nodes))/TV::m;
            TV triangle_force=Add_Velocity_Independent_Forces_Helper(relative_velocity,t);
            F(nodes(1))+=triangle_force;F(nodes(2))+=triangle_force;F(nodes(3))+=triangle_force;}
        else{
            if(rigid_body->Has_Infinite_Inertia() || rigid_body->particle_index<0) return;
            {std::stringstream ss;ss<<"M: Adding drag to "<<rigid_body->particle_index<<std::endl;LOG::filecout(ss.str());}
            TWIST<TV>& wrench=rigid_F(rigid_body->particle_index);
            for(int t=1;t<=rigid_body->simplicial_object->mesh.elements.m;t++){
                TV wind_velocity=use_constant_wind?constant_wind:optimization(t).wind_velocity;
                TV relative_velocity=wind_velocity-rigid_body->Pointwise_Object_Velocity(optimization(t).center);
                TV simplex_force=Add_Velocity_Independent_Forces_Helper(relative_velocity,t);
                wrench.linear+=simplex_force;wrench.angular+=TV::Cross_Product(optimization(t).center-rigid_body->X(),simplex_force);}
            {std::stringstream ss;ss<<"M: Added "<<wrench.linear<<std::endl;LOG::filecout(ss.str());}
            {std::stringstream ss;ss<<"M: Surface "<<surface_area<<std::endl;LOG::filecout(ss.str());}
        }}
    if(linear_normal_viscosity && deformable_simplicial_object){
        TV wind_acceleration=linear_normal_viscosity*constant_wind;
        for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
            if(use_constant_wind) F(p)+=particles.mass(p)*TV::Dot_Product(wind_acceleration,vertex_normals(p))*vertex_normals(p);
            else F(p)+=particles.mass(p)*linear_normal_viscosity*TV::Dot_Product(Spatially_Varying_Wind_Velocity(particles.X(p)),vertex_normals(p))*vertex_normals(p);}}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV> void WIND_DRAG<TV>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    if(use_constant_wind==use_spatially_varying_wind) PHYSBAM_FATAL_ERROR();
    T wind_viscosity=use_constant_wind?constant_wind_viscosity:spatially_varying_wind_viscosity;
    if(wind_viscosity){
        if(deformable_simplicial_object) for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
            const TV_INT& nodes=deformable_simplicial_object->mesh.elements(t);
            // wind drag pressure - per unit area
            TV negative_V=-ARRAYS_COMPUTATIONS::Sum(V.Subset(nodes))/TV::m;
            TV triangle_force=wind_viscosity*(negative_V-TV::Dot_Product(negative_V,optimization(t).inward_normal)*optimization(t).inward_normal);
            // total force
            triangle_force*=optimization(t).area_over_m; // one third of the triangle force is distriduted to each node
            F(nodes(1))+=triangle_force;F(nodes(2))+=triangle_force;F(nodes(3))+=triangle_force;}
        else{
            if(rigid_body->Has_Infinite_Inertia() || rigid_body->particle_index<0) return;
            TWIST<TV>& wrench=rigid_F(rigid_body->particle_index);
            for(int t=1;t<=rigid_body->simplicial_object->mesh.elements.m;t++){
                TV negative_V=-rigid_body->Pointwise_Object_Velocity(optimization(t).center);
                TV simplex_force=wind_viscosity*(negative_V-TV::Dot_Product(negative_V,optimization(t).inward_normal)*optimization(t).inward_normal);
                // total force
                simplex_force*=optimization(t).area_over_m; // one third of the triangle force is distriduted to each node
                wrench.linear+=simplex_force;wrench.angular+=TV::Cross_Product(optimization(t).center-rigid_body->X(),simplex_force);}}}
    if(linear_normal_viscosity && deformable_simplicial_object){
        for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
            F(p)-=particles.mass(p)*linear_normal_viscosity*TV::Dot_Product(V(p),vertex_normals(p))*vertex_normals(p);}}
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV> void WIND_DRAG<TV>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,const ARRAY<bool>& rigid_particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    force_elements.Update(deformable_simplicial_object->mesh.elements,particle_is_simulated);
    force_particles.Update(deformable_simplicial_object->mesh.elements.Flattened(),particle_is_simulated);
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template class WIND_DRAG<VECTOR<T,3> >; \
    template WIND_DRAG<VECTOR<T,3> >::WIND_DRAG(TRIANGULATED_SURFACE<T>&,RIGID_BODY_COLLECTION<VECTOR<T,3> >&); \
    template WIND_DRAG<VECTOR<T,3> >::WIND_DRAG(TETRAHEDRALIZED_VOLUME<T>&,RIGID_BODY_COLLECTION<VECTOR<T,3> >&); \
    template class WIND_DRAG<VECTOR<T,2> >; \
    template WIND_DRAG<VECTOR<T,2> >::WIND_DRAG(SEGMENTED_CURVE_2D<T>&,RIGID_BODY_COLLECTION<VECTOR<T,2> >&);
INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
