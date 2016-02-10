//#####################################################################
// Copyright 2011-
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_CHIMERA_IMPLICIT_SURFACE
//#####################################################################
#ifndef __RENDERING_CHIMERA_IMPLICIT_SURFACE__
#define __RENDERING_CHIMERA_IMPLICIT_SURFACE__

#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <PhysBAM_Geometry/Tessellation/IMPLICIT_OBJECT_TESSELLATION.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/CHIMERA_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{
#define GRID(I) implicit_surfaces(I)->levelset.grid

template<class T_GRID> struct GRID_ARRAYS_POLICY;



template<class T>
class RENDERING_CHIMERA_IMPLICIT_SURFACE:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;using RENDERING_OBJECT<T>::inverse_transform;
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection; // silence -Woverloaded-virtual

    ARRAY<CHIMERA_LEVELSET_IMPLICIT_OBJECT<TV>*> implicit_surfaces;
    ARRAY<FRAME<TV>*> grid_frames;

    RENDERING_CHIMERA_IMPLICIT_SURFACE(ARRAY<CHIMERA_LEVELSET_IMPLICIT_OBJECT<TV>*>& implicit_surfaces_input,ARRAY<FRAME<TV>*>& grid_frames_input)
        :implicit_surfaces(implicit_surfaces_input),grid_frames(grid_frames_input)
    {
        for(int grid_index=2;grid_index<=implicit_surfaces.m;++grid_index){
            implicit_surfaces(grid_index)->coarse_levelset=&implicit_surfaces(grid_index-1)->levelset;
            implicit_surfaces(grid_index)->coarse_grid_frame=grid_frames(grid_index-1)->Inverse_Times(*grid_frames(grid_index));
        }
        implicit_surfaces(1)->intersect_with_box=true; //count intersection with box for the largest grid
    }

    virtual ~RENDERING_CHIMERA_IMPLICIT_SURFACE()
    {
        for(int grid_index=1;grid_index<=implicit_surfaces.m;++grid_index) delete implicit_surfaces(grid_index);
    }

    bool Intersection(RAY<TV>& ray) const PHYSBAM_OVERRIDE
    {
        bool has_any_intersection=false;
        ARRAY<RAY<TV> > rays;ARRAY<bool> have_intersect;
        const RAY<TV> temp_ray=Object_Space_Ray(ray);
        for(int grid_index=1;grid_index<=implicit_surfaces.m;++grid_index){
            RAY<TV> object_space_ray(grid_frames(grid_index)->Inverse_Times(temp_ray.endpoint),grid_frames(grid_index)->Inverse().r.Rotate(temp_ray.direction));
            object_space_ray.semi_infinite=temp_ray.semi_infinite;object_space_ray.t_max=temp_ray.t_max;object_space_ray.aggregate_id=temp_ray.aggregate_id;
            rays.Append(object_space_ray);have_intersect.Append(false);}
        for(int grid_index=1;grid_index<=implicit_surfaces.m;++grid_index){
            const RAY<TV> ray_save=rays(grid_index);
            if(implicit_surfaces(grid_index)->Intersection(rays(grid_index),small_number)){
                if(grid_index<implicit_surfaces.m && GRID(grid_index+1).domain.Inside(rays(grid_index+1).Point(rays(grid_index).t_max),small_number))
                    rays(grid_index)=ray_save;
                else have_intersect(grid_index)=has_any_intersection=true;}}
        int closest_intersection_grid_index=1;
        T smallest_t_max=FLT_MAX;
        for(int grid_index=implicit_surfaces.m;grid_index>=1;--grid_index)
            if(have_intersect(grid_index) && rays(grid_index).t_max<smallest_t_max-small_number){
                closest_intersection_grid_index=grid_index;
                smallest_t_max=rays(grid_index).t_max;}
        if(has_any_intersection){
            ray.semi_infinite=false;
            ray.t_max=smallest_t_max;
            ray.aggregate_id=rays(closest_intersection_grid_index).aggregate_id;}
        return has_any_intersection;
    }

    TV Normal(const TV& location,const int aggregate=0) const PHYSBAM_OVERRIDE
    {if(aggregate!=-1)return implicit_surfaces(1)->box.Normal(aggregate);
    for(int grid_index=implicit_surfaces.m;grid_index>=2;--grid_index){
        TV object_space_location=grid_frames(grid_index)->Inverse_Times(Object_Space_Point(location));
        const T& blend_thickness=implicit_surfaces(grid_index)->blend_thickness;
        if(GRID(grid_index).domain.Inside(object_space_location,blend_thickness)){
            TV normal=implicit_surfaces(grid_index)->Normal(object_space_location,aggregate);
            return (World_Space_Vector(grid_frames(grid_index)->r.Rotate(normal))).Normalized();}
        else if(GRID(grid_index).domain.Lazy_Inside(object_space_location)){
            T alpha=min((object_space_location-GRID(grid_index).domain.min_corner).Min(),(GRID(grid_index).domain.max_corner-object_space_location).Min())/blend_thickness;
            return (World_Space_Vector(alpha*grid_frames(grid_index)->r.Rotate(implicit_surfaces(grid_index)->Normal(object_space_location,aggregate))
              +(1-alpha)*grid_frames(grid_index-1)->r.Rotate(implicit_surfaces(grid_index-1)->Normal(implicit_surfaces(grid_index)->coarse_grid_frame*object_space_location,aggregate)))).Normalized();
        }
    }
    TV object_space_location=grid_frames(1)->Inverse_Times(Object_Space_Point(location));
    return World_Space_Vector(grid_frames(1)->r.Rotate(implicit_surfaces(1)->Normal(object_space_location,aggregate)));}

    bool Inside(const TV& location) const PHYSBAM_OVERRIDE
    {
        for(int grid_index=implicit_surfaces.m;grid_index>=1;--grid_index){
            TV object_space_location=grid_frames(grid_index)->Inverse_Times(Object_Space_Point(location));
            if(GRID(grid_index).domain.Inside(object_space_location,0.5*GRID(grid_index).dX.Max()))
                return implicit_surfaces(grid_index)->Inside(object_space_location,small_number);
        }
        return false;
    }

    bool Outside(const TV& location) const PHYSBAM_OVERRIDE
    {
        for(int grid_index=implicit_surfaces.m;grid_index>=1;--grid_index){
            TV object_space_location=grid_frames(grid_index)->Inverse_Times(Object_Space_Point(location));
            if(GRID(grid_index).domain.Inside(object_space_location,0.5*GRID(grid_index).dX.Max()))
                return implicit_surfaces(grid_index)->Outside(object_space_location,small_number);
        }
        return true;
    }

    bool Boundary(const TV& location) const PHYSBAM_OVERRIDE
    {
        for(int grid_index=implicit_surfaces.m;grid_index>=1;--grid_index){
            TV object_space_location=grid_frames(grid_index)->Inverse_Times(Object_Space_Point(location));
            if(GRID(grid_index).domain.Inside(object_space_location,0.5*GRID(grid_index).dX.Max()))
                return !implicit_surfaces(grid_index)->Inside(object_space_location,small_number) &&
                       !implicit_surfaces(grid_index)->Outside(object_space_location,small_number);
        }
        return false;
    }

    T Signed_Distance(const TV& location) const PHYSBAM_OVERRIDE
    {
        for(int grid_index=implicit_surfaces.m;grid_index>=2;--grid_index){
            TV object_space_location=grid_frames(grid_index)->Inverse_Times(Object_Space_Point(location));
            if(GRID(grid_index).domain.Inside(object_space_location,0.5*GRID(grid_index).dX.Max()))
                return (*implicit_surfaces(grid_index))(object_space_location);
        }
        TV object_space_location=grid_frames(1)->Inverse_Times(Object_Space_Point(location));
        return (*implicit_surfaces(1))(object_space_location);
    }

    TRIANGULATED_SURFACE<T>* Generate_Triangles() const PHYSBAM_OVERRIDE //TODO (yuey)
    {
    TRIANGULATED_SURFACE<T>* surface=TESSELLATION::Generate_Triangles(*implicit_surfaces(1));surface->Update_Triangle_List();return surface;}

    RANGE<TV> Object_Space_Bounding_Box() const PHYSBAM_OVERRIDE
    {
    return implicit_surfaces(1)->box;}

//#####################################################################
};
}
#endif
