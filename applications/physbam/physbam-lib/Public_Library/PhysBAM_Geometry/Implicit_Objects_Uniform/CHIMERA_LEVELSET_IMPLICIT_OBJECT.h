//#####################################################################
// Copyright 2002-2007, Doug Enright, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Neil Molino, Andrew Selle, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CHIMERA_LEVELSET_IMPLICIT_OBJECT
//#####################################################################
#ifndef __CHIMERA_LEVELSET_IMPLICIT_OBJECT__
#define __CHIMERA_LEVELSET_IMPLICIT_OBJECT__

#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
#include <PhysBAM_Geometry/Level_Sets/IMPLICIT_OBJECT_ON_A_RAY.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
namespace PhysBAM{

template<class TV> class GRID;

//#####################################################################
// Function Iterative_Solver_Tolerance
//#####################################################################
namespace{
template<class T> inline T Iterative_Solver_Tolerance(){STATIC_ASSERT((T)false);}
template<> inline float Iterative_Solver_Tolerance(){return (float).01;}
template<> inline double Iterative_Solver_Tolerance(){return .001;}
}

template<class TV>
class CHIMERA_LEVELSET_IMPLICIT_OBJECT:public LEVELSET_IMPLICIT_OBJECT<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;
    typedef typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX T_SYMMETRIC_MATRIX;
    enum WORKAROUND {d=TV::m};
    typedef VECTOR<T,d-1> T_PRINCIPAL_CURVATURES;
public:
    typedef LEVELSET_IMPLICIT_OBJECT<TV> BASE;
    using BASE::box;using BASE::levelset;
    using IMPLICIT_OBJECT<TV>::use_secondary_interpolation;

    T_LEVELSET* coarse_levelset;
    FRAME<TV> coarse_grid_frame;
    T blend_thickness;
    bool intersect_with_box;
public:

    CHIMERA_LEVELSET_IMPLICIT_OBJECT(GRID<TV>& grid_input,T_ARRAYS_SCALAR& phi_input)
    :LEVELSET_IMPLICIT_OBJECT<TV>(grid_input,phi_input),coarse_levelset(NULL),blend_thickness(12*grid_input.dX.Max()),intersect_with_box(false)
    {
    }

    static CHIMERA_LEVELSET_IMPLICIT_OBJECT<TV>* Create(){
        CHIMERA_LEVELSET_IMPLICIT_OBJECT* levelset_implicit_object=new CHIMERA_LEVELSET_IMPLICIT_OBJECT(*(new GRID<TV>),*(new T_ARRAYS_SCALAR));
        levelset_implicit_object->need_destroy_data=true;return levelset_implicit_object;
    }

    T operator()(const TV& location) const PHYSBAM_OVERRIDE{
        if(coarse_levelset && levelset.grid.domain.Outside(location,-blend_thickness)){
           T alpha=max((T)0,min((location-levelset.grid.domain.min_corner).Min(),(levelset.grid.domain.max_corner-location).Min())/blend_thickness);
           return alpha*levelset.Phi(location)+(1-alpha)*coarse_levelset->Phi(coarse_grid_frame*location);
        }
        return levelset.Phi(location);
    }


//#####################################################################
// Function Intersection
//#####################################################################
bool Intersection(RAY<TV>& ray,const T thickness) const
{
    bool save_semi_infinite=ray.semi_infinite;T save_t_max=ray.t_max;int save_aggregate=ray.aggregate_id;

    T t_start,t_end;
    bool exit_intersection=false;int exit_aggregate=0;T exit_t_max=0;
    int intersect_box=INTERSECTION::Intersects(ray,box,thickness);
    int outside_box=box.Outside(ray.endpoint,thickness);

    if(outside_box && !intersect_box) return false; // missed the box
    else if(outside_box){ // intersected the box from the outside
        TV point=ray.Point(ray.t_max); // intersection point with the box
        point=box.Thickened(-4*thickness).Clamp(point); // moves the point inside the box
        if(intersect_with_box && (*this)(point) <= 0) return true; // level set is on the edge of the box
        else{ // level set is not on the edge of the box
            t_start=ray.t_max; // intersection with the box
            ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate; // box intersection doesn't count
            RAY<TV> new_ray(point,ray.direction);INTERSECTION::Intersects(new_ray,box,thickness);
            if(ray.semi_infinite) t_end=t_start+new_ray.t_max;
            else t_end=min(t_start+new_ray.t_max,ray.t_max);}}
    else if(!intersect_box){t_start=0;t_end=ray.t_max;} // intersected some object inside the box
    else{ // intersects the box from inside
        t_start=0;t_end=ray.t_max;
        if(intersect_with_box) {exit_intersection=true;exit_t_max=ray.t_max;exit_aggregate=ray.aggregate_id;} // save for exiting rays
        ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate;} // box intersection doesn't count

    if(!use_secondary_interpolation){
        // set up marching
        IMPLICIT_OBJECT_ON_A_RAY<IMPLICIT_OBJECT<TV> > implicit_surface_on_a_ray(*this,ray);
        ITERATIVE_SOLVER<T> iterative_solver;iterative_solver.tolerance=Iterative_Solver_Tolerance<T>()*thickness;
        // start marching
        T t1=t_start+thickness,phi1=(*this)(ray.Point(t1)),t2=t1+Integration_Step(phi1);
        // march through the line segment
        while(t2 <= t_end){
            T phi2=(*this)(ray.Point(t2));
            if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){ray.semi_infinite=false;ray.t_max=iterative_solver.Bisection_Secant_Root(implicit_surface_on_a_ray,t1,t2);ray.aggregate_id=-1;return true;}
            else{t1=t2;phi1=phi2;t2=t1+Integration_Step(phi1);}}
        // check the last piece of the line segment
        t2=t_end;T phi2=(*this)(ray.Point(t_end));
        if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){ray.semi_infinite=false;ray.t_max=iterative_solver.Bisection_Secant_Root(implicit_surface_on_a_ray,t1,t2);ray.aggregate_id=-1;return true;}
        else if(exit_intersection && phi2 <= 0){ray.semi_infinite=false;ray.t_max=exit_t_max;ray.aggregate_id=exit_aggregate;return true;} // exiting ray
        else return false;}
    else{ // use_secondary_interpolation
        // set up marching
        //IMPLICIT_OBJECT_ON_A_RAY_SECONDARY_INTERPOLATION<T> implicit_surface_on_a_ray(*this,ray); // TODO: we should probably be using this instead
        IMPLICIT_OBJECT_ON_A_RAY<IMPLICIT_OBJECT<TV> > implicit_surface_on_a_ray(*this,ray);
        ITERATIVE_SOLVER<T> iterative_solver;iterative_solver.tolerance=Iterative_Solver_Tolerance<T>()*thickness;
        // start marching
        T t1=t_start+thickness,phi1=(*this)(ray.Point(t1)),t2=t1+Integration_Step(phi1);
        // march through the line segment

        while(t2 <= t_end){
            T phi2=(*this)(ray.Point(t2));
            if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){
                phi1=this->Phi_Secondary(ray.Point(t1));
                phi2=this->Phi_Secondary(ray.Point(t2));
                if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){
                    ray.semi_infinite=false;
                    ray.t_max=iterative_solver.Bisection_Secant_Root(implicit_surface_on_a_ray,t1,t2);
                    ray.aggregate_id=-1;return true;}
                else{t1=t2;phi1=phi2;t2=t1+Integration_Step(phi1);}}
            else{t1=t2;phi1=phi2;t2=t1+Integration_Step(phi1);}}

        // check the last piece of the line segment
        t2=t_end;T phi2=(*this)(ray.Point(t2));
        if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){
            phi1=this->Phi_Secondary(ray.Point(t1));
            phi2=this->Phi_Secondary(ray.Point(t2));
            if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){
                ray.semi_infinite=false;
                ray.t_max=iterative_solver.Bisection_Secant_Root(implicit_surface_on_a_ray,t1,t2);
                ray.aggregate_id=-1;return true;}}
        if(exit_intersection && phi2 <= 0){ray.semi_infinite=false;ray.t_max=exit_t_max;ray.aggregate_id=exit_aggregate;return true;} // exiting ray
        else return false;}
}
};
}
#endif
