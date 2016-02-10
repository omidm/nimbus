//#####################################################################
// Copyright 2011, Mridul Aanjaneya, Linhai Qiu. 
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/PROJECTED_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_FRAME.h>
#endif
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GRID.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/SEGMENT_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_SEGMENTED_CURVE_2D_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_TRIANGULATED_SURFACE_INTERSECTION.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> RIGID_GRID<T_GRID>::
RIGID_GRID(RIGID_GRID_COLLECTION<T_GRID>& rigid_grid_collection_input,int input_index,RIGID_GEOMETRY_STATE<TV> previous_state_input)
    :rigid_grid_collection(rigid_grid_collection_input),is_static(false),bounding_box_up_to_date(false),previous_state(previous_state_input)
{
    if(input_index) index=input_index;
    else index=rigid_grid_collection.particles.array_collection->Add_Element();
    assert(!rigid_grid_collection.particles.rigid_grid(index));
    rigid_grid_collection.particles.rigid_grid(index)=dynamic_cast<RIGID_GRID<T_GRID>*>(this);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> RIGID_GRID<T_GRID>::
~RIGID_GRID()
{
    rigid_grid_collection.particles.rigid_grid(index)=0;
}
//#####################################################################
// Function Update_Grid_From_Twist
//#####################################################################
template<class T_GRID> void RIGID_GRID<T_GRID>::
Update_Frame_From_Twist(const T dt, const T time)
{
    previous_state.frame=Frame();
    X()+=Twist().linear*dt;Rotation()=ROTATION<TV>::From_Rotation_Vector(Twist().angular*dt)*Rotation();Rotation().Normalize();
    Update_Bounding_Box();
}
//#####################################################################
// Function Update_Bounding_Box
//#####################################################################
template<class T_GRID> void RIGID_GRID<T_GRID>::
Update_Bounding_Box()
{
    bounding_box_up_to_date=true;
    const RANGE<TV>& box=grid.domain;
    oriented_box=T_ORIENTED_BOX(box,Frame());
    axis_aligned_bounding_box=oriented_box.Axis_Aligned_Bounding_Box();
}
//#####################################################################
// Function Pointwise_Object_Velocity
//#####################################################################
template<class T_GRID> typename T_GRID::VECTOR_T RIGID_GRID<T_GRID>::
Pointwise_Object_Velocity(const TV& X) const
{
    return Pointwise_Object_Velocity(Twist(),this->X(),X);
}
//#####################################################################
// Function Pointwise_Object_Velocity_At_Particle
//#####################################################################
template<class T_GRID> typename T_GRID::VECTOR_T RIGID_GRID<T_GRID>::
Pointwise_Object_Velocity_At_Particle(const TV& X,const int index) const
{
    return Pointwise_Object_Velocity(X);
}
//#####################################################################
// Function Compute_Velocity_Between_States
//#####################################################################
template<class T_GRID> void RIGID_GRID<T_GRID>::
Compute_Velocity_Between_States(const RIGID_GEOMETRY_STATE<TV>& state1,const RIGID_GEOMETRY_STATE<TV>& state2,RIGID_GEOMETRY_STATE<TV>& result_state)
{
    RIGID_GEOMETRY_STATE<TV>::Compute_Velocity_Between_States(state1,state2,result_state);
}
//#####################################################################
// Function Interpolate_Between_States
//#####################################################################
template<class T_GRID> void RIGID_GRID<T_GRID>::
Interpolate_Between_States(const RIGID_GEOMETRY_STATE<TV>& state1,const RIGID_GEOMETRY_STATE<TV>& state2,const T time,RIGID_GEOMETRY_STATE<TV>& interpolated_state)
{
    PHYSBAM_ASSERT(((time>=state1.time && time<=state2.time) || (time<=state1.time && time>=state2.time)) && state1.time<state2.time);
    T alpha=(time-state1.time)/(state2.time-state1.time);alpha=clamp(alpha,(T)0,(T)1);
    interpolated_state.frame.t=TV::Interpolate(state1.frame.t,state2.frame.t,alpha);
    interpolated_state.frame.r=ROTATION<TV>::Spherical_Linear_Interpolation(state1.frame.r,state2.frame.r,alpha);
    interpolated_state.twist.linear=(1-alpha)*state1.twist.linear+alpha*state2.twist.linear;
    interpolated_state.time=time;
}
//#####################################################################
template class RIGID_GRID<GRID<VECTOR<float,1> > >;
template class RIGID_GRID<GRID<VECTOR<float,2> > >;
template class RIGID_GRID<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_GRID<GRID<VECTOR<double,1> > >;
template class RIGID_GRID<GRID<VECTOR<double,2> > >;
template class RIGID_GRID<GRID<VECTOR<double,3> > >;
#endif
}
