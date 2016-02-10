#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2005, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OCTREE_REMOVED_PARTICLES_PROCESSING
//#####################################################################
#ifndef __OCTREE_REMOVED_PARTICLES_PROCESSING__
#define __OCTREE_REMOVED_PARTICLES_PROCESSING__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

template<class T> class OCTREE_GRID;
template<class TV> class PARTICLE_LEVELSET_REMOVED_PARTICLES;
template<class T,int d> class VECTOR;
template<class TV> class RANGE;
template<class T> class ORIENTED_BOX;
template<class T> class REMOVED_PARTICLES_BLENDER_3D;
template<class T_GRID> class GRID_BASED_COLLISION_GEOMETRY_DYADIC;

template<class T>
class OCTREE_REMOVED_PARTICLES_PROCESSING:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;
public:
    OCTREE_GRID<T>& octree_grid; // NOTE: this grid will be refined
    ARRAY<T>& water_phi; // NOTE: this array will be refined (in sync with octree_grid)
    ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>& particles_array;
    ARRAY<T> particle_octree_phi;
    GRID_BASED_COLLISION_GEOMETRY_DYADIC<OCTREE_GRID<T> >* collision_body_list;
    T blending_parameter;
    T scale;
    int octree_maximum_depth;
    T particle_power;
    bool use_velocity_scaling;
    T dt; // for velocity scaling
    bool preserve_volume;

    OCTREE_REMOVED_PARTICLES_PROCESSING(OCTREE_GRID<T>& octree_grid_input,ARRAY<T>& water_phi_input,
        ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*>& particles_array_input)
        :octree_grid(octree_grid_input),water_phi(water_phi_input),particles_array(particles_array_input),
        collision_body_list(0),
        blending_parameter((T).8),scale((T).25),octree_maximum_depth(2),particle_power(T(1)),
        use_velocity_scaling(true),dt((T)1/24),preserve_volume(true)
    {}

//#####################################################################
    void Set_Collision_Aware(GRID_BASED_COLLISION_GEOMETRY_DYADIC<OCTREE_GRID<T> >* collision_body_list_input);
    void Refine_And_Create_Particle_Phi();
    void Merge_Phi(ARRAY<T>& result);
    void Union_Phi(ARRAY<T>& result);
    void Blend_Phi(ARRAY<T>& result,const T blend_cells);
private:
    void Get_Ellipsoid(const PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>* particles,const int p,T& radius_x,T& radius_yz,TV& major_axis);
    void Get_Particle_Bounding_Boxes(ARRAY<ARRAY<ORIENTED_BOX<TV> > >& particles_bounding_box,REMOVED_PARTICLES_BLENDER_3D<T>& particle_blender,int& number_of_particles);
//#####################################################################
};
}
#endif
#endif
