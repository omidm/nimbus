//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REMOVED_PARTICLES_PROCESSING
//#####################################################################
#ifndef __REMOVED_PARTICLES_PROCESSING__
#define __REMOVED_PARTICLES_PROCESSING__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/KD_TREE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Basic_Geometry/ELLIPSOID.h>
namespace PhysBAM{

template<class TV> class LEVELSET_IMPLICIT_OBJECT;
template<class T> class REMOVED_PARTICLES_BLENDER_3D;
template<class TV> class PARTICLE_LEVELSET_REMOVED_PARTICLES;

template<class T>
class REMOVED_PARTICLES_PROCESSING:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;

    PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> particles;
    ARRAY<ELLIPSOID<T> > ellipsoids;
    ARRAY<SYMMETRIC_MATRIX<T,3> > metrics;
    REMOVED_PARTICLES_BLENDER_3D<T>* particle_blender;
    T blending_parameter;
    T scale;
    T relative_tolerance,tolerance;
    RANGE<TV> particle_domain;
    int grid_divisions;
    GRID<TV> particle_grid;
    ARRAY<ARRAY<int> ,VECTOR<int,3> > particle_array;
    KD_TREE<TV> particle_tree;
    
    REMOVED_PARTICLES_PROCESSING(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles_input)
        :blending_parameter((T).8),scale((T)1),relative_tolerance((T)0.01),tolerance((T)0),grid_divisions(150)
    {
        Initialize(particles_input);
    }

    REMOVED_PARTICLES_PROCESSING(GRID<TV>& grid_input,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV> *,VECTOR<int,3> >& particles_array_input)
        :blending_parameter((T).8),scale((T)1),relative_tolerance((T)0.01),tolerance((T)0),grid_divisions(150)
    {
        Initialize(grid_input,particles_array_input);
    }

private:
    void Initialize(PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles_input)
    {particles.array_collection->Initialize(*particles_input.array_collection);
    ellipsoids.Resize(particles.array_collection->Size());metrics.Resize(particles.array_collection->Size());}

    void Initialize(const GRID<TV>& grid,ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,VECTOR<int,3> >& particles_array)
    {int number_of_particles=0;
    for(CELL_ITERATOR it(grid,3);it.Valid();it.Next()) if(particles_array(it.Cell_Index())) number_of_particles+=particles_array(it.Cell_Index())->array_collection->Size();
    std::stringstream ss;ss<<"Processing "<<number_of_particles<<" removed particles"<<std::endl;LOG::filecout(ss.str());
    particles.array_collection->Preallocate(number_of_particles);ellipsoids.Resize(number_of_particles);metrics.Resize(number_of_particles);
    for(CELL_ITERATOR it(grid,3);it.Valid();it.Next())if(particles_array(it.Cell_Index())) particles.array_collection->Take(*particles_array(it.Cell_Index())->array_collection);}

//#####################################################################
public:
    void Setup_Processing();
    T Phi(const TV& position) const;
    TV Normal(const TV& position) const;
private:
    ELLIPSOID<T> Get_Ellipsoid(const int p) const;
//#####################################################################
};
}
#endif
