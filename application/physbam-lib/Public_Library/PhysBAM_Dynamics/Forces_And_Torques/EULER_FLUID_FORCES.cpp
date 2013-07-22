//#####################################################################
// Copyright 2007-2008, Nipun Kwatra, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_FLUID_FORCES
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Dynamics/Forces_And_Torques/EULER_FLUID_FORCES.h>
using ::std::sqrt;
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> EULER_FLUID_FORCES<T_GRID>::
EULER_FLUID_FORCES(const T_GRID& grid_input,const T_FACE_ARRAYS_SCALAR& pressure_at_faces_input,
    const T_FACE_ARRAYS_BOOL& solid_fluid_face_input,const T_ARRAYS_BOOL& cells_inside_fluid_input,
    const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<T_GRID>* collision_bodies_affecting_fluid_input,PARTICLES<TV>& particles_input,
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input):SOLIDS_FORCES<TV>(particles_input,rigid_body_collection_input),
    grid(grid_input),pressure_at_faces(pressure_at_faces_input),solid_fluid_face(solid_fluid_face_input),
    cells_inside_fluid(cells_inside_fluid_input),collision_bodies_affecting_fluid(collision_bodies_affecting_fluid_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> EULER_FLUID_FORCES<T_GRID>::
~EULER_FLUID_FORCES()
{}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T_GRID> void EULER_FLUID_FORCES<T_GRID>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const
{
    COLLISION_GEOMETRY_ID body_id;int simplex_id;
    T distance,max_distance=grid.dX.Magnitude();
    TV one_over_dx=grid.one_over_dX;T cell_size=grid.Cell_Size();
    TV_INT face_index,first_cell_index,second_cell_index;TV first_cell_location,second_cell_location;int axis;
    bool first_cell_inside_fluid,second_cell_inside_fluid;
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){axis=iterator.Axis();face_index=iterator.Face_Index();
        if(!solid_fluid_face.Component(axis)(face_index)) continue;
        first_cell_index=iterator.First_Cell_Index();second_cell_index=iterator.Second_Cell_Index();
        first_cell_inside_fluid=cells_inside_fluid.Valid_Index(first_cell_index) && cells_inside_fluid(first_cell_index);
        second_cell_inside_fluid=cells_inside_fluid.Valid_Index(second_cell_index) && cells_inside_fluid(second_cell_index);
        int direction=0;
        if(first_cell_inside_fluid && !second_cell_inside_fluid) direction=1; 
        else if(!first_cell_inside_fluid && second_cell_inside_fluid) direction=-1;
        if(direction!=0){ // face on solid-fluid boundary
            TV location=iterator.Location();
            collision_bodies_affecting_fluid->collision_geometry_collection.Closest_Boundary_Point(location,max_distance,distance,body_id,simplex_id);
            const COLLISION_GEOMETRY<TV>& collision_body=collision_bodies_affecting_fluid->collision_geometry_collection(body_id);
            if(const RIGID_COLLISION_GEOMETRY<TV>* rigid_collision_geometry=dynamic_cast<const RIGID_COLLISION_GEOMETRY<TV>*>(&collision_body)){
                // apply force and torque to this body from the pressure flux times area
                const RIGID_GEOMETRY<TV>& rigid_geometry=rigid_collision_geometry->rigid_geometry;
                int rigid_body_index=rigid_geometry.particle_index;
                T face_area=cell_size*one_over_dx[axis],face_pressure=pressure_at_faces.Component(axis)(face_index);
                TV center_of_mass=rigid_geometry.X(),force=face_area*face_pressure*direction*TV::Axis_Vector(axis);
                rigid_F(rigid_body_index).linear+=force;
                rigid_F(rigid_body_index).angular+=TV::Cross_Product(location-center_of_mass,force);}
            else PHYSBAM_FATAL_ERROR("deformable part not implemented");}}
}
//#####################################################################
template class EULER_FLUID_FORCES<GRID<VECTOR<float,1> > >;
template class EULER_FLUID_FORCES<GRID<VECTOR<float,2> > >;
template class EULER_FLUID_FORCES<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EULER_FLUID_FORCES<GRID<VECTOR<double,1> > >;
template class EULER_FLUID_FORCES<GRID<VECTOR<double,2> > >;
template class EULER_FLUID_FORCES<GRID<VECTOR<double,3> > >;
#endif
