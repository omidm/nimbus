//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_REMOVED_PARTICLES_PROCESSING
//#####################################################################
#ifndef __UNIFORM_REMOVED_PARTICLES_PROCESSING__
#define __UNIFORM_REMOVED_PARTICLES_PROCESSING__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

template<class TV> class PARTICLE_LEVELSET_REMOVED_PARTICLES;
template<class TV> class LEVELSET_IMPLICIT_OBJECT;

template<class T>
class UNIFORM_REMOVED_PARTICLES_PROCESSING:public NONCOPYABLE
{
public:
    typedef VECTOR<T,3> TV;typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;typedef typename GRID<TV>::NODE_ITERATOR NODE_ITERATOR;typedef VECTOR<int,3> TV_INT;

    GRID<TV>& grid;
    GRID<TV>* sim_grid;
    ARRAY<T,VECTOR<int,3> >& water_phi;  
    ARRAY<T,VECTOR<int,3> > particle_phi;
    ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> *,VECTOR<int,3> >& particle_array;
    T blending_parameter;
    T scale;
    T scale_factor;
    T particle_power;
    bool use_velocity_scaling;
    T dt; // for velocity scaling
    bool preserve_volume;

    UNIFORM_REMOVED_PARTICLES_PROCESSING(GRID<TV>& grid_input,ARRAY<T,VECTOR<int,3> >& water_phi_input,
        ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> *,VECTOR<int,3> >& particle_array_input)
        :grid(grid_input),sim_grid(&grid_input),water_phi(water_phi_input),particle_array(particle_array_input),blending_parameter((T).8),
        scale((T).25),scale_factor((T)1),particle_power(T(1)),use_velocity_scaling(true),dt((T)1/24),preserve_volume(true)
    {
        Initialize();
    }

    void Initialize()
    {
        particle_phi.Resize(grid.Domain_Indices(3),false,false);particle_phi.Fill(0);
    }

//#####################################################################
    void Refine_Grid_To_Particle_Size(const LEVELSET_IMPLICIT_OBJECT<TV>* water_levelset);
    void Incorporate_Removed_Negative_Particles();
    void Merge_Phi(ARRAY<T,VECTOR<int,3> >& result) const;
    void Union_Phi(ARRAY<T,VECTOR<int,3> >& result) const;
    void Blend_Phi(ARRAY<T,VECTOR<int,3> >& result,const T blend_cells) const;
private:
    void Get_Ellipsoid(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles,const int p,T& radius_x,T& radius_yz,TV& major_axis) const;
//#####################################################################
};
}
#endif
