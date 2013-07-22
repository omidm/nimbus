//#####################################################################
// Copyright 2006-2007, Geoffrey Irving, Nipun Kwatra, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STRUCTURE_INTERACTION_GEOMETRY
//#####################################################################
#ifndef __STRUCTURE_INTERACTION_GEOMETRY__
#define __STRUCTURE_INTERACTION_GEOMETRY__

#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/choice.h>
#include <PhysBAM_Tools/Parallel_Computation/PARTITION_ID.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD_SUBSET.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{

template<class TV>
class STRUCTURE_INTERACTION_GEOMETRY
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    struct UNUSABLE{};
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,1>::OBJECT T_SEGMENTED_CURVE;
public:
    GEOMETRY_PARTICLES<TV>& full_particles;
    POINT_CLOUD_SUBSET<TV,GEOMETRY_PARTICLES<TV> > collision_particles;
    TRIANGULATED_SURFACE<T>* triangulated_surface;
    T_SEGMENTED_CURVE* segmented_curve;
    POINT_SIMPLICES_1D<T>* point_simplices;
    INDIRECT_ARRAY<ARRAY_VIEW<TV> > subset;
    PARTICLE_HIERARCHY<TV,INDIRECT_ARRAY<ARRAY_VIEW<TV> > > particle_hierarchy;
    ARRAY<char> triangulated_surface_processor_masks,segmented_curve_processor_masks,point_processor_masks;
    ARRAY<bool> triangulated_surface_modified,segmented_curve_modified,point_modified; // whether boxes below you have been modified
private:
    bool need_destroy_segmented_curve;
    UNUSABLE unusable;
public:

    STRUCTURE_INTERACTION_GEOMETRY(GEOMETRY_PARTICLES<TV>& full_particles_input)
        :full_particles(full_particles_input),collision_particles(full_particles),triangulated_surface(0),segmented_curve(0),point_simplices(0),
        subset(collision_particles.point_cloud.X,collision_particles.active_indices),particle_hierarchy(subset,false,0),need_destroy_segmented_curve(false)
    {}

    ~STRUCTURE_INTERACTION_GEOMETRY()
    {Clean_Memory();}

    void Clean_Memory()
    {if(need_destroy_segmented_curve) delete segmented_curve;need_destroy_segmented_curve=false;
    collision_particles.Clean_Memory();triangulated_surface=0;segmented_curve=0;point_simplices=0;particle_hierarchy.Clean_Memory();}

    void Build_Topological_Structure_Of_Hierarchies()
    {if(triangulated_surface) triangulated_surface->Initialize_Hierarchy(false);
    if(segmented_curve) segmented_curve->Initialize_Hierarchy(false);
    ARRAY_VIEW<TV> tmp_view(collision_particles.point_cloud.X.Size(),collision_particles.point_cloud.X.Get_Array_Pointer());
    subset.array.Exchange(tmp_view);
    particle_hierarchy.Initialize_Hierarchy_Using_KD_Tree();}

    const typename IF<d==2,ARRAY<int>,typename IF<d==3,ARRAY<VECTOR<int,2> >,UNUSABLE>::TYPE>::TYPE& Edges() const
    {return choice<d>(unusable,collision_particles.active_indices,segmented_curve->mesh.elements);}

    const typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d-1>::OBJECT* Face_Mesh_Object() const
    {return choice<d>(point_simplices,segmented_curve,triangulated_surface);}

    const ARRAY<bool>& Edge_Modified() const
    {return d==2?point_modified:segmented_curve_modified;}

    const ARRAY<bool>& Face_Modified() const
    {return d==2?segmented_curve_modified:triangulated_surface_modified;}

    typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d-1>::HIERARCHY& Face_Hierarchy() const
    {return *Face_Mesh_Object()->hierarchy;}

    bool Has_Edges() const
    {return d==2 || segmented_curve!=0;}

    const typename IF<(d<2),UNUSABLE,typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d<2?0:d-2>::HIERARCHY>::TYPE& Edge_Hierarchy() const
    {return choice<d>(unusable,particle_hierarchy,*segmented_curve->hierarchy);}

    const ARRAY<char>& Edge_Processor_Masks() const
    {return d==2?point_processor_masks:segmented_curve_processor_masks;}

    const ARRAY<char>& Face_Processor_Masks() const
    {return d==2?segmented_curve_processor_masks:triangulated_surface_processor_masks;}

//#####################################################################
    static TRIANGULATED_SURFACE<T>* Triangulated_Surface(STRUCTURE<TV>* structure);
    static T_SEGMENTED_CURVE* Segmented_Curve(STRUCTURE<TV>* structure);
    void Build_Collision_Geometry(STRUCTURE<TV>& structure);
    void Update_Faces_And_Hierarchies_With_Collision_Free_Positions(ARRAY_VIEW<const T> node_thickness,const T node_thickness_multiplier,ARRAY_VIEW<const TV> X_old_full);
    void Update_Processor_Masks(const PARTITION_ID processor,const ARRAY<PARTITION_ID>& partition_id_from_particle_index);
private:
    template<class T_OBJECT,class T_HIERARCHY>
    void Update_Processor_Masks_Helper(T_OBJECT& object,T_HIERARCHY& hierarchy,const PARTITION_ID processor,const ARRAY<PARTITION_ID>& partition_id_from_particle_index,ARRAY<char>& mask);
//#####################################################################
};
}
#endif
