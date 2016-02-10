//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Dynamics/Level_Sets/REMOVED_PARTICLES_BLENDER_3D.h>
#include <PhysBAM_Dynamics/Level_Sets/UNIFORM_REMOVED_PARTICLES_PROCESSING.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
using namespace PhysBAM;

//#####################################################################
// Function Refine_Grid
//#####################################################################
template<class T> void UNIFORM_REMOVED_PARTICLES_PROCESSING<T>::
Refine_Grid_To_Particle_Size(const LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* water_levelset)
{
    sim_grid=new GRID<TV>(grid);
    GRID<TV> temp_grid(TV_INT(scale_factor*TV(grid.counts)),grid.domain,true);
    ARRAY<T,VECTOR<int,3> > temp_phi(temp_grid.Domain_Indices(3));
    for(CELL_ITERATOR iterator(temp_grid,3);iterator.Valid();iterator.Next()) temp_phi(iterator.Cell_Index())=(*water_levelset)(iterator.Location());
    grid=temp_grid;water_phi.Resize(grid.Domain_Indices(3),false,false);ARRAY<T,VECTOR<int,3> >::Copy(temp_phi,water_phi);
    particle_phi.Resize(grid.Domain_Indices(3),false,false);particle_phi.Fill(0);
}
//#####################################################################
// Function Get_Ellipsoid
//#####################################################################
template<class T> void UNIFORM_REMOVED_PARTICLES_PROCESSING<T>::
Get_Ellipsoid(PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<T,3> >& particles,int p,T& radius_x,T& radius_yz,VECTOR<T,3>& major_axis) const
{
    T radius=scale*particles.radius(p);
    T velocity_magnitude_squared=particles.V(p).Magnitude_Squared();
    if(velocity_magnitude_squared>1e-8){ // ellipsoid
        T speed=sqrt(velocity_magnitude_squared);
        major_axis=particles.V(p)/speed;
        if(use_velocity_scaling){
            radius_x=radius+(T).5*dt*speed;
            if(preserve_volume){radius_yz=sqrt(cube(radius)/radius_x);}
            else{radius_yz=radius;}}
        else{radius_x=3*radius;radius_yz=radius;}}
    else{ // sphere
        major_axis=VECTOR<T,3>(1,0,0); // arbitrary axis
        radius_x=radius;radius_yz=radius;}
}
//#####################################################################
// Function Incorporate_Removed_Negative_Particles
//#####################################################################
template<class T> void UNIFORM_REMOVED_PARTICLES_PROCESSING<T>::
Incorporate_Removed_Negative_Particles()
{
    REMOVED_PARTICLES_BLENDER_3D<T> particle_blender(blending_parameter);
    T max_dX_times_particle_power=grid.dX.Max()*particle_power;
    for(NODE_ITERATOR it(*sim_grid);it.Valid();it.Next())if(particle_array(it.Node_Index())){
        PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<T,3> >& particles=*particle_array(it.Node_Index());
        for(int p=1;p<=particles.array_collection->Size();p++){
            T radius_x,radius_yz;VECTOR<T,3> major_axis;Get_Ellipsoid(particles,p,radius_x,radius_yz,major_axis);
            T one_over_radius_x_squared=(T)1/sqr(radius_x),one_over_radius_yz_squared=(T)1/sqr(radius_yz);
            RANGE<TV> box=particle_blender.Get_Bounding_Box(radius_x,radius_yz,particles.X(p),major_axis);
            VECTOR<int,3> min_index=grid.Clamped_Index_End_Minus_One(box.Minimum_Corner())+VECTOR<int,3>(1,1,1),max_index=grid.Clamped_Index(box.Maximum_Corner());
            for(CELL_ITERATOR jt(grid,RANGE<VECTOR<int,3> >(min_index,max_index));jt.Valid();jt.Next()){
                T distance=particle_blender.Get_Distance(one_over_radius_x_squared,one_over_radius_yz_squared,particles.X(p),major_axis,grid.X(jt.Cell_Index()));
                particle_phi(jt.Cell_Index())-=max_dX_times_particle_power*particle_blender.C(distance);}}}
}
//#####################################################################
// Function Merge_Phi
//#####################################################################
template<class T> void UNIFORM_REMOVED_PARTICLES_PROCESSING<T>::
Merge_Phi(ARRAY<T,VECTOR<int,3> >& result) const
{
    assert((ARRAY<T,VECTOR<int,3> >::Equal_Dimensions(water_phi,particle_phi)));
    result.Resize(grid.Domain_Indices(3),false);
    result.array=water_phi.array+particle_phi.array;
}
//#####################################################################
// Function Union_Phi
//#####################################################################
template<class T> void UNIFORM_REMOVED_PARTICLES_PROCESSING<T>::
Union_Phi(ARRAY<T,VECTOR<int,3> >& result) const
{
    assert((ARRAY<T,VECTOR<int,3> >::Equal_Dimensions(water_phi,particle_phi)));
    result.Resize(grid.Domain_Indices(3),false);
    T offset=grid.dX.Max()*blending_parameter*particle_power;
    for(int i=1;i<=water_phi.array.Size();i++) result.array(i)=min(particle_phi.array(i)+offset,water_phi.array(i));
}
//#####################################################################
// Function Blend_Phi
//#####################################################################
template<class T> void UNIFORM_REMOVED_PARTICLES_PROCESSING<T>::
Blend_Phi(ARRAY<T,VECTOR<int,3> >& result,const T blend_cells) const
{
    assert((ARRAY<T,VECTOR<int,3> >::Equal_Dimensions(water_phi,particle_phi)));
    result.Resize(grid.Domain_Indices(3),false);
    T offset=grid.dX.Max()*blending_parameter*particle_power;
    T scale=1/(blend_cells*grid.dX.Max());
    for(int i=1;i<=water_phi.array.Size();i++){
        T alpha=clamp(scale*water_phi.array(i),(T)0,(T)1);
        result.array(i)=(1-alpha)*(water_phi.array(i)+particle_phi.array(i))+alpha*(particle_phi.array(i)+offset);}
}
//#####################################################################
template class UNIFORM_REMOVED_PARTICLES_PROCESSING<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class UNIFORM_REMOVED_PARTICLES_PROCESSING<double>;
#endif
