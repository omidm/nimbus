#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_REMOVED_PARTICLES_PROCESSING
//#####################################################################
#ifndef __RLE_REMOVED_PARTICLES_PROCESSING__
#define __RLE_REMOVED_PARTICLES_PROCESSING__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
namespace PhysBAM{

template<class T> class OCTREE_GRID;
template<class TV> class PARTICLE_LEVELSET_REMOVED_PARTICLES;
template<class TV> class RANGE;
template<class T> class ORIENTED_BOX;
template<class T> class REMOVED_PARTICLES_BLENDER_3D;
template<class T_GRID> class GRID_BASED_COLLISION_GEOMETRY_DYADIC;

template<class T>
class RLE_REMOVED_PARTICLES_PROCESSING:public NONCOPYABLE
{
public:
    typedef RLE_GRID_3D<T> T_GRID;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename T_GRID::VECTOR_HORIZONTAL TV_HORIZONTAL;typedef typename T_GRID::BOX_HORIZONTAL_INT T_BOX_HORIZONTAL_INT;
    typedef typename T_GRID::RUN T_RUN;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::HORIZONTAL_GRID::CELL_ITERATOR HORIZONTAL_CELL_ITERATOR;

    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles;
    T_GRID grid;
    ARRAY<T> phi;
    T blending_parameter;
    T scale;
    int octree_maximum_depth;
    T particle_power;
    bool use_velocity_scaling;
    T dt; // for velocity scaling
    bool preserve_volume;

    RLE_REMOVED_PARTICLES_PROCESSING(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles_input)
        :particles(particles_input),blending_parameter((T).8),scale((T).25),octree_maximum_depth(2),particle_power(T(1)),use_velocity_scaling(true),dt((T)1/24),
        preserve_volume(true)
    {}

//#####################################################################
    void Build_Grid_And_Rasterize_Particles();
private:
    void Get_Ellipsoid(const int p,T& radius_x,T& radius_yz,TV& major_axis);
    void Get_Particle_Bounding_Boxes(ARRAY<RANGE<TV> >& particles_bounding_box,REMOVED_PARTICLES_BLENDER_3D<T>& particle_blender);
//#####################################################################
};
}
#endif
#endif
