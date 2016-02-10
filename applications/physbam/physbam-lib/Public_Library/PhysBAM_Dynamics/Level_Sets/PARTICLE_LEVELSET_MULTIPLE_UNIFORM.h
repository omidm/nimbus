//#####################################################################
// Copyright 2005, Ron Fedkiw, Eran Guendelman, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_LEVELSET_MULTIPLE_UNIFORM  
//##################################################################### 
#ifndef __PARTICLE_LEVELSET_MULTIPLE_UNIFORM__
#define __PARTICLE_LEVELSET_MULTIPLE_UNIFORM__

#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_MULTIPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/PARTICLE_LEVELSET_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID>
class PARTICLE_LEVELSET_MULTIPLE_UNIFORM:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<PARTICLE_LEVELSET_PARTICLES<TV>*>::TYPE T_ARRAYS_PARTICLE_LEVELSET_PARTICLES;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
public:
    T min_collision_distance_factor,max_collision_distance_factor,max_minus_min_collision_distance_factor_over_max_short;

    LEVELSET_MULTIPLE_UNIFORM<T_GRID> levelset_multiple;
    ARRAY<PARTICLE_LEVELSET_UNIFORM<T_GRID>*> particle_levelsets;
    int number_of_ghost_cells;

    PARTICLE_LEVELSET_MULTIPLE_UNIFORM(T_GRID& grid_input,ARRAY<T_ARRAYS_SCALAR>& phis_input,const int number_of_ghost_cells_input)
        :levelset_multiple(grid_input,phis_input,true),number_of_ghost_cells(number_of_ghost_cells_input)
    {
        Set_Collision_Distance_Factors(); // TODO: use this from a normal particle levelset
    }

    void Initialize_Particle_Levelsets_And_Grid_Values(T_GRID& grid,ARRAY<T_ARRAYS_SCALAR>& phis,const int number_of_regions)
    {if(particle_levelsets.m!=number_of_regions){
        for(int i=1;i<=particle_levelsets.m;i++)delete particle_levelsets(i);
        particle_levelsets.Resize(number_of_regions);levelset_multiple.levelsets.Resize(particle_levelsets.m);
        for(int i=1;i<=particle_levelsets.m;i++){
            particle_levelsets(i)=new PARTICLE_LEVELSET_UNIFORM<T_GRID>(grid,phis(i),number_of_ghost_cells);
            particle_levelsets(i)->only_use_negative_particles=true;
            particle_levelsets(i)->Initialize_Particle_Levelset_Grid_Values();
            levelset_multiple.levelsets(i)=&particle_levelsets(i)->levelset;}}
    else for(int i=1;i<=particle_levelsets.m;i++) particle_levelsets(i)->Initialize_Particle_Levelset_Grid_Values();}

    void Store_Unique_Particle_Id()
    {for(int i=1;i<=particle_levelsets.m;i++) particle_levelsets(i)->Store_Unique_Particle_Id();}

    void Set_Band_Width(const T number_of_cells=6) 
    {for(int i=1;i<=particle_levelsets.m;i++) particle_levelsets(i)->Set_Band_Width(number_of_cells);}

    void Use_Removed_Negative_Particles(const bool use_removed_negative_particles_input=true)
    {for(int i=1;i<=particle_levelsets.m;i++) particle_levelsets(i)->Use_Removed_Negative_Particles(use_removed_negative_particles_input);}

    void Use_Removed_Positive_Particles(const bool use_removed_positive_particles_input=true)
    {for(int i=1;i<=particle_levelsets.m;i++) particle_levelsets(i)->Use_Removed_Positive_Particles(use_removed_positive_particles_input);}

    void Set_Collision_Distance_Factors(const T min_collision_distance_factor_input=(T).1,const T max_collision_distance_factor_input=(T)1)
    {min_collision_distance_factor=min_collision_distance_factor_input;max_collision_distance_factor=max_collision_distance_factor_input;
    max_minus_min_collision_distance_factor_over_max_short=(max_collision_distance_factor-min_collision_distance_factor)/USHRT_MAX;}

    T Particle_Collision_Distance(const unsigned short quantized_collision_distance)
    {return (min_collision_distance_factor+(T)quantized_collision_distance*max_minus_min_collision_distance_factor_over_max_short)*levelset_multiple.grid.Minimum_Edge_Length();}

    void Seed_Particles(const T time,const bool verbose=true)
    {for(int i=1;i<=particle_levelsets.m;i++) particle_levelsets(i)->Seed_Particles(time,verbose);}

    void Adjust_Particle_Radii()
    {for(int i=1;i<=particle_levelsets.m;i++) particle_levelsets(i)->Adjust_Particle_Radii();}

    void Modify_Levelset_Using_Escaped_Particles(T_FACE_ARRAYS_SCALAR* face_velocities)
    {for(int i=1;i<=particle_levelsets.m;i++){
        ARRAY<T_ARRAYS_PARTICLE_LEVELSET_PARTICLES*> other_positive_particles(particle_levelsets.m-1);
        int index=0;for(int j=1;j<=particle_levelsets.m;j++)if(i!=j)other_positive_particles(++index)=&particle_levelsets(j)->negative_particles;
        particle_levelsets(i)->Modify_Levelset_Using_Escaped_Particles(face_velocities,&other_positive_particles);}}

    void Euler_Step_Particles(const T_FACE_ARRAYS_SCALAR& V,const T dt,const T time=0,const bool use_second_order_for_nonremoved_particles=false,const bool update_particle_cells_after_euler_step=true,const bool verbose=true)
    {for(int i=1;i<=particle_levelsets.m;i++) particle_levelsets(i)->Euler_Step_Particles(V,dt,time,use_second_order_for_nonremoved_particles,update_particle_cells_after_euler_step,verbose);}

    void Euler_Step_Removed_Particles(const T dt,const T time,const bool update_particle_cells_after_euler_step=true,const bool verbose=true)
    {for(int i=1;i<=particle_levelsets.m;i++) particle_levelsets(i)->Euler_Step_Removed_Particles(dt,time,update_particle_cells_after_euler_step,verbose);}

    int Reseed_Particles(const T time,T_ARRAYS_BOOL* cell_centered_mask=0)
    {int new_particles=0;for(int i=1;i<=particle_levelsets.m;i++) new_particles+=particle_levelsets(i)->Reseed_Particles(time,cell_centered_mask);
    return new_particles;}

    void Delete_Particles_Outside_Grid()
    {for(int i=1;i<=particle_levelsets.m;i++) particle_levelsets(i)->Delete_Particles_Outside_Grid();}

    void Identify_And_Remove_Escaped_Particles(const T_ARRAYS_VECTOR& V,const T radius_fraction=1.5,const bool verbose=true)
    {for(int i=1;i<=particle_levelsets.m;i++) particle_levelsets(i)->Identify_And_Remove_Escaped_Particles(V,radius_fraction,verbose);}

    void Reincorporate_Removed_Particles(const T radius_fraction)
    {for(int i=1;i<=particle_levelsets.m;i++) particle_levelsets(i)->Reincorporate_Removed_Particles(radius_fraction);}

//#####################################################################
};
}
#endif
