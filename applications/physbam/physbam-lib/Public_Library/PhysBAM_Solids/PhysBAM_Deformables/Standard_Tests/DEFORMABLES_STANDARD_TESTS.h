//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLES_STANDARD_TESTS
//#####################################################################
#ifndef __DEFORMABLES_STANDARD_TESTS__
#define __DEFORMABLES_STANDARD_TESTS__

#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_STATE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_POLICY.h>
namespace PhysBAM{

template<class TV> class SOFT_BINDINGS;
template<class TV> class PARTICLES;
template<class TV> class GEOMETRY_PARTICLES;
template<class TV> class LEVELSET_IMPLICIT_OBJECT;
template<class TV> class EXAMPLE;
template<class T> class EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE;
template<class TV> class BINDING_LIST;
template<class TV,class T_PARTICLES> class POINT_CLOUD_SUBSET;
template<class TV,int d> class EMBEDDED_MATERIAL_SURFACE;
template<class TV> class TRIANGLE_COLLISION_PARAMETERS;
template<class TV> class GRID;
template<class TV> class DEFORMABLE_BODY_COLLECTION;

template<class T>
class TRIANGULATED_SURFACE_CLIPPING_HELPER
{
public:
    virtual ~TRIANGULATED_SURFACE_CLIPPING_HELPER() {};
    virtual void operator()(TRIANGULATED_SURFACE<T>& surface) const=0;
};

template<class TV>
class DEFORMABLES_STANDARD_TESTS
{
protected:
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::TRIANGULATED_OBJECT T_TRIANGULATED_OBJECT;
    typedef typename TOPOLOGY_BASED_GEOMETRY_POLICY<TV>::SEGMENTED_CURVE T_SEGMENTED_CURVE;
public:
    EXAMPLE<TV>& example;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection;

    DEFORMABLES_STANDARD_TESTS(EXAMPLE<TV>& example_input,DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection_input);
    virtual ~DEFORMABLES_STANDARD_TESTS(){}

    TRIANGULATED_AREA<T>& Create_Mattress(const GRID<TV>& mattress_grid,const bool use_constant_mass,const RIGID_GEOMETRY_STATE<TV>& initial_state)
    {return Create_Mattress(mattress_grid,use_constant_mass,&initial_state);}

    template<class T_OBJECT> void
    Substitute_Soft_Bindings_For_Embedded_Nodes(T_OBJECT& object,SOFT_BINDINGS<TV>& soft_bindings,HASHTABLE<int,int>* persistent_soft_bindings=0)
    {Substitute_Soft_Bindings_For_Nodes(object,soft_bindings,persistent_soft_bindings,true);}

    TRIANGULATED_SURFACE<T>& Create_Cloth_Panel(const int number_side_panels,const T side_length,const T aspect_ratio,const RIGID_GEOMETRY_STATE<TV>& initial_state,
        TRIANGULATED_SURFACE_CLIPPING_HELPER<T> *clipping_function,ARRAY<int>* particle_indices=0)
    {return Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,&initial_state,clipping_function,particle_indices);}

    TRIANGULATED_SURFACE<T>& Create_Cloth_Panel(const int number_side_panels,const T side_length,const T aspect_ratio,const RIGID_GEOMETRY_STATE<TV>* initial_state,ARRAY<int>* particle_indices=0)
    {return Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,initial_state,0,particle_indices);}

    TRIANGULATED_SURFACE<T>& Create_Cloth_Panel(const int number_side_panels,const T side_length,const T aspect_ratio,const RIGID_GEOMETRY_STATE<TV>& initial_state,ARRAY<int>* particle_indices=0)
    {return Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,&initial_state,particle_indices);}

//#####################################################################
    virtual void Add_Gravity();
    template<class T_STRUCTURE>
    T_STRUCTURE& Copy_And_Add_Structure(T_STRUCTURE& structure,ARRAY<int>* particle_indices=0);
    void Set_Initial_Particle_Configuration(GEOMETRY_PARTICLES<TV>& particles,const RIGID_GEOMETRY_STATE<TV>& state,const bool relative_to_box_center);
    T_SEGMENTED_CURVE& Create_Segmented_Curve(const GRID<VECTOR<T,1> >& square_grid,const RIGID_GEOMETRY_STATE<TV>& initial_state,const T density);
    T_SEGMENTED_CURVE& Create_Segmented_Curve(const int m,const RIGID_GEOMETRY_STATE<TV>& initial_state,const T initial_radius=(T)1,const T density=(T)100);
    TETRAHEDRALIZED_VOLUME<T>& Create_Tetrahedralized_Volume(const std::string& filename,const RIGID_GEOMETRY_STATE<TV>& initial_state,const bool relative_to_box_center,
        const bool use_constant_mass,const T density,const T scale=1);
    T_TRIANGULATED_OBJECT& Create_Triangulated_Object(const GRID<TV>& square_grid,const RIGID_GEOMETRY_STATE<TV>& initial_state,const T density);
    T_TRIANGULATED_OBJECT& Create_Triangulated_Object(const std::string& filename,const RIGID_GEOMETRY_STATE<TV>& initial_state,const bool relative_to_box_center,
        const bool use_constant_mass,const T scale=1);
    T_SEGMENTED_CURVE& Create_Segmented_Curve(const std::string& filename,const RIGID_GEOMETRY_STATE<TV>& initial_state,const bool relative_to_box_center,
        const bool use_constant_mass);
    TRIANGULATED_AREA<T>& Create_Mattress(const GRID<VECTOR<T,2> >& mattress_grid,const bool use_constant_mass,const RIGID_GEOMETRY_STATE<TV>* initial_state=0,const T density=(T)1000,const bool reverse_triangles=false);
    TETRAHEDRALIZED_VOLUME<T>& Create_Mattress(const GRID<VECTOR<T,3> >& mattress_grid,const bool use_constant_mass=true,const RIGID_GEOMETRY_STATE<TV>* initial_state=0,const T density=(T)1000);
    template<class T_SHAPE>
    EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T>& Create_Embedded_Tetrahedralized_Volume(const T_SHAPE& shape,const RIGID_GEOMETRY_STATE<TV>& initial_state,const bool relative_to_box_center);
    template<class T_OBJECT> void
    Substitute_Soft_Bindings_For_Nodes(T_OBJECT& object,SOFT_BINDINGS<TV>& soft_bindings,HASHTABLE<int,int>* persistent_soft_bindings=0,const bool embedded_only=false,
        const bool use_impulses_for_collisions=true);
    LEVELSET_IMPLICIT_OBJECT<TV>* Read_Or_Initialize_Implicit_Surface(const std::string& levelset_filename,TRIANGULATED_SURFACE<T>& undeformed_triangulated_surface) const;
    void Initialize_Tetrahedron_Collisions(const int id_number,TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,TRIANGLE_COLLISION_PARAMETERS<TV>& triangle_collision_parameters,
        TRIANGULATED_SURFACE<T>* triangulated_surface=0);
    TRIANGULATED_SURFACE<T>& Create_Drifted_Surface(const TRIANGULATED_SURFACE<T>& triangulated_surface,SOFT_BINDINGS<TV>& soft_bindings,const bool use_impulses_for_collisions=false) const;
    template <class T_OBJECT> static void Set_Mass_Of_Particles(const T_OBJECT& volume,const T density,const bool use_constant_mass=false);
    void PD_Curl(const T scale,const TV shift,const ROTATION<TV> orient,const T k_p,const int number_of_joints,const bool parent_static=true,const T friction=.5);
    TRIANGULATED_SURFACE<T>& Create_Cloth_Panel(const int number_side_panels,const T side_length,const T aspect_ratio,const RIGID_GEOMETRY_STATE<TV>* initial_state,
        TRIANGULATED_SURFACE_CLIPPING_HELPER<T> *clipping_function,ARRAY<int>* particle_indices);
    void Embed_Particles_In_Tetrahedralized_Volume(BINDING_LIST<VECTOR<T,3> >& binding_list,const POINT_CLOUD_SUBSET<VECTOR<T,3>,PARTICLES<VECTOR<T,3> > >& particles_to_embed,
        TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,const T thickness_over_two);
    void Mark_Hard_Bindings_With_Free_Particles();
//#####################################################################
};
}
#endif
