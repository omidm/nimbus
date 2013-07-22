//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Data_Structures/KD_TREE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
#include <PhysBAM_Dynamics/Level_Sets/REMOVED_PARTICLES_BLENDER_3D.h>
#include <PhysBAM_Dynamics/Level_Sets/REMOVED_PARTICLES_PROCESSING.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Setup_Processing
//#####################################################################
template<class T> void REMOVED_PARTICLES_PROCESSING<T>::
Setup_Processing()
{
    particle_blender=new REMOVED_PARTICLES_BLENDER_3D<T>(blending_parameter);
    particle_tree.Create_Left_Balanced_KD_Tree(particles.X);
    tolerance=relative_tolerance*scale*ARRAYS_COMPUTATIONS::Min(particles.radius);
    particle_domain=RANGE<TV>(particles.X(1));
    ARRAY<RANGE<TV> > particle_boxes(particles.array_collection->Size());
    for(int p=1;p<=particles.array_collection->Size();p++){
        ellipsoids(p)=Get_Ellipsoid(p);metrics(p)=ellipsoids(p).Metric_Tensor();
        particle_boxes(p)=particle_blender->Get_Bounding_Box(ellipsoids(p));
        particle_domain.Enlarge_To_Include_Box(particle_boxes(p));}
    particle_grid=GRID<TV>(grid_divisions,grid_divisions,grid_divisions,particle_domain);
    ARRAY<ARRAY<int> ,VECTOR<int,3> > conservative_array(particle_grid.Domain_Indices(1));
    LOG::filecout("Rasterizing particles to grid...\n");
    for(int p=1;p<=particles.array_collection->Size();p++){
        RANGE<VECTOR<int,3> > particle_box=particle_grid.Clamp_To_Cell(particle_boxes(p),1);
        for(CELL_ITERATOR it(particle_grid,particle_box);it.Valid();it.Next())for(int n=1;n<=8;n++){
            T distance=REMOVED_PARTICLES_BLENDER_3D<T>::Get_Distance(ellipsoids(p).center,metrics(p),particle_grid.Node(it.Cell_Node_Index(n)));
            if(particle_blender->C(distance)>0){conservative_array(it.Cell_Index()).Append(p);break;}}}
    particle_array.Resize(particle_grid.Domain_Indices());
    for(CELL_ITERATOR it(particle_grid);it.Valid();it.Next()){particle_array(it.Cell_Index())=conservative_array(it.Cell_Index());
    for(int i=1;i<=26;i++)particle_array(it.Cell_Index()).Append_Unique_Elements(conservative_array(particle_grid.One_Ring_Neighbor(it.Cell_Index(),i)));}
}
//#####################################################################
// Function Phi
//#####################################################################
template<class T> T REMOVED_PARTICLES_PROCESSING<T>::
Phi(const TV& position) const
{
    T phi=0,lipschitz=0;
    TV_INT cell=particle_grid.Clamp_To_Cell(position);
    if(particle_array(cell).m==0) return particle_grid.min_dX;
    for(int p=1;p<=particle_array(cell).m;p++){
        int index=particle_array(cell)(p);
        T distance=REMOVED_PARTICLES_BLENDER_3D<T>::Get_Distance(ellipsoids(index).center,metrics(index),position);
        phi-=particle_blender->C(distance);lipschitz+=ellipsoids(index).Lipschitz_Constant();} 
    T distance=(T)2/3*(blending_parameter+phi)*particle_blender->R/lipschitz;
    return sign(distance)*min(abs(distance),particle_grid.min_dX);
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> REMOVED_PARTICLES_PROCESSING<T>::
Normal(const TV& position) const
{
    TV normal;
    TV_INT cell=particle_grid.Clamp_To_Cell(position);
    for(int p=1;p<=particle_array(cell).m;p++){
        int index=particle_array(cell)(p);
        T distance=REMOVED_PARTICLES_BLENDER_3D<T>::Get_Distance(ellipsoids(index).center,metrics(index),position);
        normal+=particle_blender->C_Prime(distance)*ellipsoids(index).Normal(position);}
    return normal.Normalized();
}
//#####################################################################
// Function Get_Ellipsoid
//#####################################################################
template<class T> ELLIPSOID<T> REMOVED_PARTICLES_PROCESSING<T>::
Get_Ellipsoid(const int p) const
{
    // TODO: generalize, optimize, and add options -- still largely a placeholder function
    T radius=scale*particles.radius(p);
    int number_of_points_in_estimate=100;
    T density_threshold=50/(radius);

    ARRAY<int> points_found(number_of_points_in_estimate);ARRAY<T> squared_distance_to_points_found(number_of_points_in_estimate);
    int number_of_points_found;T max_squared_distance_to_points_found;
    particle_tree.Locate_Nearest_Neighbors(particles.X(p),FLT_MAX,points_found,squared_distance_to_points_found,number_of_points_found,max_squared_distance_to_points_found,particles.X);
    
    ARRAY<TV> positions(number_of_points_in_estimate);
    for(int i=1;i<=number_of_points_in_estimate;i++)positions(i)=particles.X(points_found(i));

    ELLIPSOID<T> covariance=ELLIPSOID<T>::Covariance_Ellipsoid(positions);
    T density=number_of_points_in_estimate/((T)four_thirds_pi*cube(sqrt(max_squared_distance_to_points_found)));

    if(density<density_threshold){
        TV x_axis=particles.V(p).Normalized(),y_axis=x_axis.Unit_Orthogonal_Vector(),z_axis=TV::Cross_Product(x_axis,y_axis).Normalized();
        return ELLIPSOID<T>(particles.X(p),DIAGONAL_MATRIX<T,3>(3*radius,radius,radius),x_axis,y_axis,z_axis);}
    else{covariance.center=particles.X(p);T clamped_radius=min(3*radius,sqrt(max_squared_distance_to_points_found));
    covariance.radii=DIAGONAL_MATRIX<T,3>(clamped_radius,clamped_radius,(T).5*radius);return covariance;}
}
//#####################################################################
template class REMOVED_PARTICLES_PROCESSING<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class REMOVED_PARTICLES_PROCESSING<double>;
#endif
