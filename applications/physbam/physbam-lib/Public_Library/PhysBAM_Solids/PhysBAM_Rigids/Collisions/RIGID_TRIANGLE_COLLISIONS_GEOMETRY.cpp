//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_TRIANGLE_COLLISIONS_GEOMETRY
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_STRUCTURE_INTERACTION_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_TRIANGLE_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Parallel_Computation/MPI_RIGIDS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_TRIANGLE_COLLISIONS_GEOMETRY<TV>::
RIGID_TRIANGLE_COLLISIONS_GEOMETRY()
    :mpi_solids(0),mass_modifier(0)
{
    // set parameters
    Allow_Intersections(false);Set_Allow_Intersections_Tolerance();
    // output
    Output_Number_Checked(false);Set_Small_Number();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_TRIANGLE_COLLISIONS_GEOMETRY<TV>::
~RIGID_TRIANGLE_COLLISIONS_GEOMETRY()
{
    structure_geometries.Delete_Pointers_And_Clean_Memory();
}
//#####################################################################
// Function Build_Collision_Geometry
//#####################################################################
template<class TV> void RIGID_TRIANGLE_COLLISIONS_GEOMETRY<TV>::
Build_Collision_Geometry()
{
    structure_geometries.Delete_Pointers_And_Clean_Memory();
    structure_geometries.Resize(structures.m);
    interacting_structure_pairs.Remove_All();
    for(int k=1;k<=structures.m;k++){
        if(MESH_OBJECT<TV,TRIANGLE_MESH>* mesh=dynamic_cast<MESH_OBJECT<TV,TRIANGLE_MESH>*>(structures(k))) structure_geometries(k)=new RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>(mesh->particles);
        else PHYSBAM_FATAL_ERROR("Geometry type is not supported");
        structure_geometries(k)->Build_Collision_Geometry(*structures(k));}
    for(int i=1;i<=structures.m;i++) for(int j=i+1;j<=structures.m;j++) interacting_structure_pairs.Append(VECTOR<int,2>(i,j));
}
//#####################################################################
// Function Build_Topological_Structure_Of_Hierarchies
//#####################################################################
template<class TV> void RIGID_TRIANGLE_COLLISIONS_GEOMETRY<TV>::
Build_Topological_Structure_Of_Hierarchies()
{
    for(int k=1;k<=structure_geometries.m;k++){
        structure_geometries(k)->Build_Topological_Structure_Of_Hierarchies();}
        //if(mpi_solids) structure_geometries(k)->Update_Processor_Masks(mpi_solids->Partition(),
        //    mpi_solids->partition_id_from_particle_index);}
}
//#####################################################################
// Function Allow_Intersections
//#####################################################################
template<class TV> void RIGID_TRIANGLE_COLLISIONS_GEOMETRY<TV>::
Allow_Intersections(const bool allow_intersections_input)
{
    allow_intersections=allow_intersections_input;
    if(allow_intersections)
        for(int k=1;k<=structure_geometries.m;k++)
            if(!structure_geometries(k)->triangulated_surface->mesh.element_edges) structure_geometries(k)->triangulated_surface->mesh.Initialize_Element_Edges();
}
//#####################################################################
// Function Save_Self_Collision_Free_State
//#####################################################################
template<class TV> void RIGID_TRIANGLE_COLLISIONS_GEOMETRY<TV>::
Save_Current_State(const ARRAY_VIEW<TV>& X,const ARRAY_VIEW<ROTATION<TV> >& rotation,const ARRAY_VIEW<TV>& V) // assumes mass does not change
{
    for(int i=1;i<=structure_geometries.m;i++) structure_geometries(i)->Save_Current_State(X(i),rotation(i),V(i));
}
//#####################################################################
// Function Save_Self_Collision_Free_State
//#####################################################################
template<class TV> void RIGID_TRIANGLE_COLLISIONS_GEOMETRY<TV>::
Save_Self_Collision_Free_State(const ARRAY_VIEW<TV>& X,const ARRAY_VIEW<ROTATION<TV> >& rotation,const ARRAY_VIEW<TV>& V) // assumes mass does not change
{
    for(int i=1;i<=structure_geometries.m;i++) structure_geometries(i)->Save_Self_Collision_Free_State(X(i),rotation(i),V(i));
}
//#####################################################################
// Function Save_Self_Collision_Free_State
//#####################################################################
template<class TV> void RIGID_TRIANGLE_COLLISIONS_GEOMETRY<TV>::
Save_Self_Collision_Free_State() // assumes mass does not change
{
    for(int i=1;i<=structure_geometries.m;i++) structure_geometries(i)->Save_Self_Collision_Free_State();
}
//#####################################################################
// Function Restore_Self_Collision_Free_State
//#####################################################################
template<class TV> void RIGID_TRIANGLE_COLLISIONS_GEOMETRY<TV>::
Restore_Self_Collision_Free_State()
{
    for(int i=1;i<=structure_geometries.m;i++) structure_geometries(i)->Restore_Self_Collision_Free_State();
}
//#####################################################################
// Function Compute_Intersecting_Segment_Face_Pairs
//#####################################################################
template<class TV> void RIGID_TRIANGLE_COLLISIONS_GEOMETRY<TV>::
Compute_Intersecting_Segment_Face_Pairs()
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//####################################################################
#define INSTANTIATION_HELPER(T,d) \
    template class RIGID_TRIANGLE_COLLISIONS_GEOMETRY<VECTOR<T,d> >;
INSTANTIATION_HELPER(float,1);
INSTANTIATION_HELPER(float,2);
INSTANTIATION_HELPER(float,3);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double,1);
INSTANTIATION_HELPER(double,2);
INSTANTIATION_HELPER(double,3);
#endif
