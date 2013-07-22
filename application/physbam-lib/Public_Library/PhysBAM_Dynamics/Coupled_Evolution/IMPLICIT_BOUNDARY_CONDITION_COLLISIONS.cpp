//#####################################################################
// Copyright 2009, Nipun Kwatra, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_BOUNDARY_CONDITION_COLLISIONS
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION_COLLISIONS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<TV>::
IMPLICIT_BOUNDARY_CONDITION_COLLISIONS(COLLISION_GEOMETRY_COLLECTION<TV>& collision_geometry_collection_input,
    const bool use_implicit_geometry_input,FLUIDS_PARAMETERS_CALLBACKS<T_GRID>& callbacks)
    :collision_geometry_collection(collision_geometry_collection_input),use_implicit_geometry(use_implicit_geometry_input),callbacks(callbacks)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<TV>::
~IMPLICIT_BOUNDARY_CONDITION_COLLISIONS()
{}
//#####################################################################
// Function Update_Boundary_Conditions
//#####################################################################
template<class TV> void IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<TV>::
Update_Boundary_Conditions(const GRID<TV>& grid,ARRAY<bool,TV_INT>& psi_D,ARRAY<bool,FACE_INDEX<TV::dimension> >& psi_N,ARRAY<T,TV_INT>& p,
    ARRAY<T,FACE_INDEX<TV::dimension> >& face_velocities,const T time)
{
    T p_inside_solid=0; // TODO: set to something nasty and make sure uncovered cells get good values for initial guess.

    // TODO(jontg): Need this to be aware of time n vs. time np1 solid location...
    // if(callbacks.Get_Psi_D_Inside_Solids(psi_D)) return;

    if(use_implicit_geometry){
        COLLISION_GEOMETRY_ID body_id;
        for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,1);iterator.Valid();iterator.Next()){TV location=iterator.Location();
            TV_INT cell_index=iterator.Cell_Index();
            if(collision_geometry_collection.Implicit_Geometry_Lazy_Inside_Any_Body(location,body_id)){
                psi_D(cell_index)=true;p(cell_index)=p_inside_solid;}}
        return;}

    for(COLLISION_GEOMETRY_ID i(1);i<=collision_geometry_collection.Size();i++) if(collision_geometry_collection.Is_Active(i)){
        T collision_thickness_over_two=collision_geometry_collection.collision_body_thickness*(T).5;
        if(DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>* object=dynamic_cast<DEFORMABLE_OBJECT_FLUID_COLLISIONS<TV>*>(&collision_geometry_collection(i))){
            if(!object->volume_object){
                BOX<TV>& box(*(object->object.bounding_box));
                TV_INT min_index=grid.Clamp_To_Cell(box.Minimum_Corner(),1);
                TV_INT max_index=grid.Clamp_To_Cell(box.Maximum_Corner(),1);
                int dummy_index;
                for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,RANGE<TV_INT>(min_index,max_index));iterator.Valid();iterator.Next())
                    if(object->object.Inside_Any_Simplex(iterator.Location(),dummy_index,collision_thickness_over_two)){
                        psi_D(iterator.Cell_Index())=true;p(iterator.Cell_Index())=p_inside_solid;}}
            else{
                for(int i=1;i<=object->volume_object->mesh.elements.m;i++){const T_SIMPLEX simplex=object->volume_object->Get_Element(i);
                    RANGE<TV> box(simplex.Bounding_Box());
                    TV_INT min_index=grid.Clamp_To_Cell(box.Minimum_Corner(),1);
                    TV_INT max_index=grid.Clamp_To_Cell(box.Maximum_Corner(),1);
                    for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,RANGE<TV_INT>(min_index,max_index));iterator.Valid();iterator.Next())
                        if(simplex.Inside(iterator.Location(),collision_thickness_over_two)){
                            psi_D(iterator.Cell_Index())=true;p(iterator.Cell_Index())=p_inside_solid;}}}}
        else if(RIGID_COLLISION_GEOMETRY_BASE<TV>* object=dynamic_cast<RIGID_COLLISION_GEOMETRY_BASE<TV>*>(&collision_geometry_collection(i))){
            RIGID_BODY<TV>& rigid_body=dynamic_cast<RIGID_BODY<TV>&>(object->rigid_geometry);
            if(!rigid_body.simplicial_object->mesh.incident_elements) rigid_body.simplicial_object->mesh.Initialize_Incident_Elements();
            if(!rigid_body.simplicial_object->mesh.adjacent_elements) rigid_body.simplicial_object->mesh.Initialize_Adjacent_Elements();
            object->Update_Bounding_Box();
            RANGE<TV> box(object->Axis_Aligned_Bounding_Box());
            TV_INT min_index=grid.Clamp_To_Cell(box.Minimum_Corner(),1);
            TV_INT max_index=grid.Clamp_To_Cell(box.Maximum_Corner(),1);
            if(rigid_body.thin_shell){
                int dummy_index;
                for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,RANGE<TV_INT>(min_index,max_index));iterator.Valid();iterator.Next()){
                    if(collision_geometry_collection(i).Inside_Any_Simplex(iterator.Location(),dummy_index)){
                        psi_D(iterator.Cell_Index())=true;p(iterator.Cell_Index())=p_inside_solid;}}}
            else{for(UNIFORM_GRID_ITERATOR_CELL<TV> iterator(grid,RANGE<TV_INT>(min_index,max_index));iterator.Valid();iterator.Next()){
                if(collision_geometry_collection(i).Inside(iterator.Location(),collision_thickness_over_two)){
                    psi_D(iterator.Cell_Index())=true;p(iterator.Cell_Index())=p_inside_solid;}}}}
        else PHYSBAM_FATAL_ERROR("Unrecognized collision body type");}
}
//#####################################################################
template class IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<VECTOR<float,1> >;
template class IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<VECTOR<float,2> >;
template class IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<VECTOR<double,1> >;
template class IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<VECTOR<double,2> >;
template class IMPLICIT_BOUNDARY_CONDITION_COLLISIONS<VECTOR<double,3> >;
#endif
