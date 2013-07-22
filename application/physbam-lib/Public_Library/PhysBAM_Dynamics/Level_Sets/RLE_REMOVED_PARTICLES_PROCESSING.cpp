#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_3D.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/AVERAGING_RLE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Log/PROGRESS_INDICATOR.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_RLE.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <PhysBAM_Dynamics/Level_Sets/REMOVED_PARTICLES_BLENDER_3D.h>
#include <PhysBAM_Dynamics/Level_Sets/RLE_REMOVED_PARTICLES_PROCESSING.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Get_Ellipsoid
//#####################################################################
template<class T> void RLE_REMOVED_PARTICLES_PROCESSING<T>::
Get_Ellipsoid(const int p,T& radius_x,T& radius_yz,TV& major_axis)
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
// Function Get_Particle_Bounding_Boxes
//#####################################################################
template<class T> void RLE_REMOVED_PARTICLES_PROCESSING<T>::
Get_Particle_Bounding_Boxes(ARRAY<RANGE<TV> >& particle_boxes,REMOVED_PARTICLES_BLENDER_3D<T>& particle_blender)
{
    for(int p=1;p<=particles.array_collection->Size();p++){
        T radius_x,radius_yz;VECTOR<T,3> major_axis;Get_Ellipsoid(p,radius_x,radius_yz,major_axis);
        ORIENTED_BOX<TV> oriented_box=particle_blender.Get_Oriented_Bounding_Box(radius_x,radius_yz,particles.X(p),major_axis);
        particle_boxes(p)=oriented_box.Axis_Aligned_Bounding_Box();}
}
//#####################################################################
// Struct BOX_HIERARCHY_PHI
//#####################################################################
template<class T_input> struct BOX_HIERARCHY_PHI:public NONCOPYABLE, public RLE_INITIAL_PHI_HELPER<T_input,3>
{
    typedef T_input T;
    BOX_HIERARCHY<VECTOR<T,3> > box_hierarchy;
    const T thickness_over_two;
private:
    mutable ARRAY<int> intersection_list;
public:

    BOX_HIERARCHY_PHI(const ARRAY<RANGE<VECTOR<T,3> > >& boxes,const T thickness_over_two_input)
        :thickness_over_two(thickness_over_two_input)
    {
        LOG::Time("initializing box hierarchy");
        box_hierarchy.box_hierarchy=boxes;
        box_hierarchy.Initialize_Hierarchy_Using_KD_Tree();
        box_hierarchy.Update_Nonleaf_Boxes();
    }

    T operator()(const VECTOR<T,3>& X) const
    {intersection_list.Remove_All();
    box_hierarchy.Intersection_List(X,intersection_list,thickness_over_two);
    T phi=thickness_over_two;
    for(int i=1;i<=intersection_list.m;i++)phi=min(phi,(X-box_hierarchy.box_hierarchy(intersection_list(i)).Clamp(X)).Magnitude());
    return phi;}
};
//#####################################################################
// Function Build_Grid_And_Rasterize_Particles
//#####################################################################
template<class T> void RLE_REMOVED_PARTICLES_PROCESSING<T>::
Build_Grid_And_Rasterize_Particles()
{
    LOG::SCOPE scope("PARTICLES PROCESSING","particle processing");
    REMOVED_PARTICLES_BLENDER_3D<T> particle_blender(blending_parameter);
    LOG::cout<<"PARAMETERS:\n  scale = "<<scale<<"\n  blending = "<<blending_parameter<<"\n  octree depth = "<<octree_maximum_depth<<"\n  R = "<<particle_blender.R<<std::endl;
    LOG::cout<<"  particle_power = "<<particle_power<<"\n  use_velocity_scaling = "<<use_velocity_scaling<<"\n  dt = "<<dt<<"\n  preserve_volume = "<<preserve_volume<<std::endl;

    LOG::Time("computing particle boxes");
    ARRAY<RANGE<TV> > boxes(particles.array_collection->Size());
    Get_Particle_Bounding_Boxes(boxes,particle_blender);

    LOG::Time("building refined grid");
    {BOX_HIERARCHY_PHI<T> box_hierarchy_implicit_surface(boxes,3*grid.positive_bandwidth);
    grid.Initialize(box_hierarchy_implicit_surface,0,true);}
    LOG::cout<<"number of cells = "<<grid.number_of_cells<<", number of blocks = "<<grid.number_of_blocks<<std::endl;

    LOG::Time("precomputing short cell locations");
    {ARRAY<T> short_cell_y(grid.number_of_cells);
    for(CELL_ITERATOR cell(grid,grid.number_of_ghost_cells);cell;cell++)if(cell.Short()) short_cell_y(cell.Cell())=cell.y();

    LOG::Time("rasterizing particles");
    PROGRESS_INDICATOR progress(particles.array_collection->Size());
    phi.Resize(grid.number_of_cells);
    for(int p=1;p<=particles.array_collection->Size();p++){
        progress.Progress();
        // construct ellipsoid
        RANGE<TV>& box=boxes(p);
        T radius_x,radius_yz;VECTOR<T,3> major_axis;Get_Ellipsoid(p,radius_x,radius_yz,major_axis);
        T one_over_radius_x_squared=(T)1/sqr(radius_x),one_over_radius_yz_squared=(T)1/sqr(radius_yz);
        // find cells in box
        TV_INT min_index=grid.uniform_grid.Clamp_To_Cell(box.Minimum_Corner(),grid.number_of_ghost_cells); // TODO: this is too conservative
        TV_INT max_index=grid.uniform_grid.Clamp_To_Cell(box.Maximum_Corner(),grid.number_of_ghost_cells);
        T_BOX_HORIZONTAL_INT horizontal_box=RANGE<TV_INT>(min_index,max_index).Get_Horizontal_Box();
        for(HORIZONTAL_CELL_ITERATOR iterator(grid.horizontal_grid,horizontal_box);iterator.Valid();iterator.Next()){
            const ARRAY<T_RUN>& column=grid.columns(iterator.Cell_Index());
            const T_RUN* start_run=grid.Clamped_Run_In_Column(column,min_index.y);int start_cell=start_run->is_long?start_run->cell+2:start_run->cell+min_index.y-start_run->jmin;
            const T_RUN* end_run=grid.Clamped_Run_In_Column(column,max_index.y);int end_cell=end_run->is_long?end_run->cell-1:end_run->cell+max_index.y-end_run->jmin;
            // rasterize distance
            TV_HORIZONTAL X_horizontal=iterator.Location();
            for(int c=start_cell;c<=end_cell;c++){TV X(X_horizontal.x,short_cell_y(c),X_horizontal.y);if(box.Lazy_Inside(X)){
                T distance=particle_blender.Get_Distance(one_over_radius_x_squared,one_over_radius_yz_squared,particles.X(p),major_axis,X);
                phi(c)-=grid.Minimum_Edge_Length()*particle_blender.C(distance)*particle_power;}}}}}

    LOG::Time("rebuilding grid");
    {T_GRID new_grid;
    {ARRAY<bool> cell_should_be_long(grid.number_of_cells);
    for(int c=1;c<=grid.number_of_cells;c++)cell_should_be_long(c)=!phi(c);
    T_GRID::Rebuild(grid,new_grid,cell_should_be_long,1,0);}
    LOG::cout<<"number of cells = "<<new_grid.number_of_cells<<", number of blocks = "<<new_grid.number_of_blocks<<std::endl;
    LOG::Time("transferring phi");
    LEVELSET_RLE<T_GRID> levelset(grid,phi);
    levelset.Transfer_Phi(new_grid);
    T_GRID::Transfer(new_grid,grid);}
    for(CELL_ITERATOR cell(grid,grid.number_of_ghost_cells);cell;cell++)if(!cell.Short()){int c=cell.Cell();phi(c)=phi(c+1)=grid.positive_bandwidth;}
}
//#####################################################################
template class RLE_REMOVED_PARTICLES_PROCESSING<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RLE_REMOVED_PARTICLES_PROCESSING<double>;
#endif
#endif
