//#####################################################################
// Copyright 2004-2008, Eran Guendelman, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class THIN_SHELLS_FLUID_COUPLING_UTILITIES
//#####################################################################
#ifndef __THIN_SHELLS_FLUID_COUPLING_UTILITIES__
#define __THIN_SHELLS_FLUID_COUPLING_UTILITIES__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Read_Write/Utilities/TYPED_STREAM.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_OBJECT_FORWARD.h>
#include <string>
namespace PhysBAM{

class PARAMETER_LIST;
template<class T_GRID> class SOLIDS_FLUIDS_EXAMPLE_UNIFORM;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template<class T_GRID> class SOLIDS_FLUIDS_EXAMPLE_DYADIC;
#endif
template<class TV> class GRID;
template<class T> class SEGMENTED_CURVE_2D;
template<class T> class TRIANGULATED_SURFACE;
template<class T> class PLANE;

template<class T>
class THIN_SHELLS_FLUID_COUPLING_UTILITIES
{
public:
//#####################################################################
    static int Add_Deformable_Object(DEFORMABLE_BODY_COLLECTION<VECTOR<T,2> >& deformable_body_collection_list,const int number_of_vertices,const VECTOR<T,2>& start_position,
        const VECTOR<T,2>& end_position);
    static int Add_Circle_Deformable_Object(DEFORMABLE_BODY_COLLECTION<VECTOR<T,2> >& deformable_body_collection_list,const int number_of_vertices,const VECTOR<T,2>& center,const T radius);
    static int Add_Grid_Deformable_Object(DEFORMABLE_BODY_COLLECTION<VECTOR<T,2> >& deformable_body_collection_list,const GRID<VECTOR<T,2> >& grid,const int edge_subdivision);
    static int Add_Deformable_Object(DEFORMABLE_BODY_COLLECTION<VECTOR<T,3> >& deformable_body_collection_list,ARRAY<int>& deformable_body_collection_enslaved_nodes,
        const GRID<VECTOR<T,2> >& cloth_grid,const MATRIX<T,4>& transform,const int constraint_mode=1);
    static int Add_Deformable_Object_From_File(const STREAM_TYPE stream_type,DEFORMABLE_BODY_COLLECTION<VECTOR<T,3> >& deformable_body_collection_list,
        ARRAY<int>& deformable_body_collection_enslaved_nodes,const std::string& filename,const MATRIX<T,4>& transform,PLANE<T>* enslaved_halfplane=0);
    static void Add_Rigid_Body_Walls(SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T,2> > >& example,const T coefficient_of_restitution=(T).5,const T coefficient_of_friction=(T).5,
        ARRAY<int>* walls_added=0);
    static void Add_Rigid_Body_Walls(SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T,3> > >& example,const T coefficient_of_restitution=(T).5,const T coefficient_of_friction=(T).5,
        ARRAY<int>* walls_added=0);
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    static void Add_Rigid_Body_Walls(SOLIDS_FLUIDS_EXAMPLE_DYADIC<QUADTREE_GRID<T> >& example,const T coefficient_of_restitution=(T).5,const T coefficient_of_friction=(T).5,
        ARRAY<int>* walls_added=0);
    static void Add_Rigid_Body_Walls(SOLIDS_FLUIDS_EXAMPLE_DYADIC<OCTREE_GRID<T> >& example,const T coefficient_of_restitution=(T).5,const T coefficient_of_friction=(T).5,
        ARRAY<int>* walls_added=0);
#endif
    static void Set_Deformable_Object_Parameters_2D(const int id,T& density,PARAMETER_LIST& parameter_list);
    static void Set_Deformable_Object_Parameters_3D(const int id,T& edge_stiffness_scaling,T& altitude_stiffness_scaling,T& density,PARAMETER_LIST& parameter_list);
    static void Set_Rigid_Body_Parameters_2D(const int id,T& density,PARAMETER_LIST& parameter_list);
    static void Set_Rigid_Body_Parameters_3D(const int id,T& density,PARAMETER_LIST& parameter_list);
    static void Set_Mass(SEGMENTED_CURVE_2D<T>& segmented_curve,const T mass,const bool use_constant_mass=false);
    static void Set_Density(SEGMENTED_CURVE_2D<T>& segmented_curve,const T density,const bool use_constant_mass=false);
    static void Set_Mass(TRIANGULATED_SURFACE<T>& triangulated_surface,const T mass,const bool use_constant_mass=false);
    static void Set_Density(TRIANGULATED_SURFACE<T>& triangulated_surface,const T density,const bool use_constant_mass=false);
    template<class T_GRID> static void Set_Parameters_From_Parameter_List(SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID>& example,PARAMETER_LIST& parameter_list);
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    template<class T_GRID> static void Set_Parameters_From_Parameter_List(SOLIDS_FLUIDS_EXAMPLE_DYADIC<T_GRID>& example,PARAMETER_LIST& parameter_list);
#endif
//#####################################################################
};
}
#endif
