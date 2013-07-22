//#####################################################################
// Copyright 2006-2007, Jeong-Mo Hong, Nipun Kwatra, Tamar Shinar, Jonatthan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_MULTIPLE_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/DETONATION_SHOCK_DYNAMICS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> DETONATION_SHOCK_DYNAMICS<T_GRID>::
DETONATION_SHOCK_DYNAMICS(T_GRID& grid_input,const T_LEVELSET& levelset_input,const int order_input)
    :grid(grid_input),levelset(levelset_input),Dn(grid),Dn_dot(grid),curvature(grid),curvature_old(grid),order(order_input),boundary(&boundary_default),
    boundary_vector(&boundary_vector_default)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> DETONATION_SHOCK_DYNAMICS<T_GRID>::
~DETONATION_SHOCK_DYNAMICS()
{
}
//#####################################################################
// Function Initialize_Grid
//#####################################################################
template<class T_GRID> void DETONATION_SHOCK_DYNAMICS<T_GRID>::
Initialize_Grid()
{
    assert(order>=1&&order<=3);
    Dn.Initialize_Array(3); // TODO
    if(order>=2) Dn_dot.Initialize_Array();
    if(order==3){curvature.Initialize_Array();curvature_old.Initialize_Array();}
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class T_GRID> void DETONATION_SHOCK_DYNAMICS<T_GRID>::
Advance_One_Time_Step(const T_FACE_ARRAYS_SCALAR& V,const T dt,const T time,const int number_of_ghost_cells)
{
    LOG::SCOPE scope("DSD Advance","DSD Advance (order=%d, dt=%f)",3,dt);

    LOG::Time("Advection of dsd variables");
    Dn.Set_Velocity(&V);Dn.Euler_Step(dt,time,number_of_ghost_cells);
    Dn_dot.Set_Velocity(&V);Dn_dot.Euler_Step(dt,time,number_of_ghost_cells);
    curvature_old.Set_Velocity(&V);curvature_old.Euler_Step(dt,time,number_of_ghost_cells);

    LOG::Time("Fill ghost phis and normals");
    T_ARRAYS_SCALAR& phi=levelset.phi;
    T_ARRAYS_VECTOR& normals=*levelset.normals;
    T_ARRAYS_SCALAR phi_ghost(grid.Cell_Indices(number_of_ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,number_of_ghost_cells);
    T_ARRAYS_VECTOR normals_ghost(grid.Cell_Indices(number_of_ghost_cells));boundary_vector->Fill_Ghost_Cells(grid,normals,normals_ghost,dt,time,number_of_ghost_cells);

    LOG::Time("Update state variables");
    T small_number=(T)0.01*grid.dX.x;int max_iteration=10;
    T_ARRAYS_BOOL interfacial(grid.Cell_Indices());
    Make_NB_Indices(grid,phi,indices_interface,dt,time,number_of_ghost_cells);
    T mean_curvature=(T)0;
    T mean_A=(T)0,mean_B=(T)0,mean_C=(T)0,mean_D=(T)0;
    T mean_Dn=(T)0,mean_Dn_dot=(T)0;
    for(int i=1;i<=indices_interface.m;i++){TV_INT index=indices_interface(i);
        interfacial(index)=true;
        // Determine Control Coefficients
        curvature.array(index)=(*levelset.curvature)(index)*grid.dX.x;
        mean_curvature+=curvature.array(index);// change to container?
        T delta=(Dn.array(index)-Dcj);// -1 < delta < 1
        T curvature_dot=(curvature.array(index)-curvature_old.array(index))/dt;
        curvature_old.array(index)=curvature.array(index);// Store old curvature
        // Determine control coefficients
        T A=A_coeff*exp((T)2*mutheta*delta);
        T B=B_coeff*exp(mutheta*delta);
        T C=C_coeff*exp((T)2*mutheta*delta);
        T D=D_coeff;
        T Lcj=(T)0;
        if(!use_log_Lcj)Lcj=curvature.array(index);
        else{
            T inside_log=dtheta*curvature.array(index)*exp(-mutheta*delta);
            inside_log=max((T)-0.999999,inside_log);
            Lcj=log((T)1+inside_log);}//nonlinear curvature relationship //Lcj(cell)=-(T)0.1*curvature(cell);//linear curvature relationship
        // Euler step of Dn and Dn_dot by Dn_ddot
        T Dn_ddot=-A*delta-B*Dn_dot.array(index)-C*Lcj-D*curvature_dot;        
        Dn_dot.array(index)+=Dn_ddot*dt;
        Dn.array(index)+=Dn_dot.array(index)*dt;
        // debugging data
        mean_A+=A*delta;mean_B+=B*Dn_dot.array(index);mean_C+=C*Lcj;mean_D+=D*curvature_dot;
        mean_Dn+=Dn.array(index);mean_Dn_dot+=Dn_dot.array(index);
        // Clamping
        if(Dn.array(index)<Dcj_min_clamp){Dn_dot.array(index)=max((T)0,Dn_dot.array(index));Dn.array(index)=max(Dn.array(index),Dcj_min_clamp);}
        else if(Dn.array(index)>Dcj_max_clamp){Dn_dot.array(index)=min((T)0,Dn_dot.array(index));Dn.array(index)=min(Dn.array(index),Dcj_max_clamp);}}

    LOG::Time("Extrapolation");
    T_ARRAYS_SCALAR Dn_ghost(grid.Cell_Indices(number_of_ghost_cells));boundary->Fill_Ghost_Cells(grid,Dn.array,Dn_ghost,dt,time,number_of_ghost_cells);
    T_ARRAYS_SCALAR Dn_dot_ghost(grid.Cell_Indices(number_of_ghost_cells));boundary->Fill_Ghost_Cells(grid,Dn_dot.array,Dn_dot_ghost,dt,time,number_of_ghost_cells);
    T_ARRAYS_SCALAR curvature_old_ghost(grid.Cell_Indices(number_of_ghost_cells));boundary->Fill_Ghost_Cells(grid,curvature_old.array,curvature_old_ghost,dt,time,number_of_ghost_cells);
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T> interpolation;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT index=iterator.Cell_Index();
        if(abs(phi(index))>(T)nb_width*grid.dX.x){Dn.array(index)=Dcj;continue;}//to make 'reaction_speed' files compact
        if(interfacial(index))continue;
        TV closest_point;
        if(Closest_Point_On_Boundary(phi_ghost,normals_ghost,grid.Center(index),closest_point,small_number,max_iteration)){
            Dn.array(index)=interpolation.Clamped_To_Array_Cell(grid,Dn_ghost,closest_point);
            Dn_dot.array(index)=interpolation.Clamped_To_Array_Cell(grid,Dn_dot_ghost,closest_point);
            curvature_old.array(index)=interpolation.Clamped_To_Array_Cell(grid,curvature_old_ghost,closest_point);}
        else{Dn.array(index)=Dcj;Dn_dot.array(index)=(T)0;curvature_old.array(index)=(T)0;}}

    Dn.boundary->Fill_Ghost_Cells(grid,Dn.array,Dn.array,dt,time,number_of_ghost_cells); // TODO

    // Debug data
    std::stringstream ss;
    if(indices_interface.m>0){
        T Dcj_min=Dn.array(indices_interface(1)),Dcj_max=Dcj_min;
        T Ddot_min=Dn_dot.array(indices_interface(1)),Ddot_max=Ddot_min;
        for(int i=1;i<=indices_interface.m;i++){TV_INT index=indices_interface(i);
            Dcj_min=min(Dn.array(index),Dcj_min);Dcj_max=max(Dn.array(index),Dcj_max);
            Ddot_min=min(Dn_dot.array(index),Ddot_min);Ddot_max=max(Dn_dot.array(index),Ddot_max);}
        ss<<"Dn ("<<Dcj_min<<","<<Dcj_max<<")"<<std::endl;
        ss<<"Dn_dot ("<<Ddot_min<<","<<Ddot_max<<")"<<std::endl;}
    if(indices_interface.m>0){
        mean_Dn/=(T)indices_interface.m;mean_Dn_dot/=(T)indices_interface.m;
        mean_curvature/=(T)indices_interface.m;
        mean_A/=(T)indices_interface.m;mean_B/=(T)indices_interface.m;mean_C/=(T)indices_interface.m;mean_D/=(T)indices_interface.m;
        ss<<"Mean Dn "<<mean_Dn<<" "<<mean_Dn_dot<<std::endl;
        ss<<"Mean curvature "<<mean_curvature<<std::endl;
        ss<<"Mean coefficients "<<mean_A<<" "<<mean_B<<" "<<mean_C<<" "<<mean_D<<std::endl;LOG::filecout(ss.str());}
}
//#####################################################################
// Function Make_NB_Indices
//#####################################################################
template<class T_GRID> void DETONATION_SHOCK_DYNAMICS<T_GRID>::
Make_NB_Indices(T_GRID &grid,T_ARRAYS_SCALAR &phi,ARRAY<TV_INT>& indices_interface,const T dt,const T time,int number_of_ghost_cells)
{
    T_ARRAYS_SCALAR phi_ghost;
    if(number_of_ghost_cells>0){   
        phi_ghost.Resize(grid.Cell_Indices(number_of_ghost_cells));
        boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,number_of_ghost_cells);}
    else phi_ghost=phi;
    indices_interface.Clean_Memory();// Clean memory of lists
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        T phi_1=phi_ghost(iterator.Cell_Index());
        for(int i=1;i<=T_GRID::number_of_faces_per_cell;i++){// painful symmetric debugging, T_GRID::dimension was a big mistake
            TV_INT neighbor_index=iterator.Cell_Neighbor(i);
            if(!phi_ghost.Valid_Index(neighbor_index))continue;
            T phi_2=phi_ghost(neighbor_index);
            if(LEVELSET_UTILITIES<T>::Interface(phi_1,phi_2)){indices_interface.Append(iterator.Cell_Index());break;}}}
}
//#####################################################################
// Function Closest_Point_On_Boundary
//#####################################################################
template<class T_GRID> bool DETONATION_SHOCK_DYNAMICS<T_GRID>::
Closest_Point_On_Boundary(T_ARRAYS_SCALAR &phi_ghost,T_ARRAYS_VECTOR &normals_ghost,const TV& location,TV& new_location,const T tolerance,const int max_iterations) const
{
    RANGE<TV> box=grid.Ghost_Domain(3);
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,T> interpolation;
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,TV> interpolation_vector;
    if(!tolerance){
        T phi=interpolation.Clamped_To_Array_Cell(grid,phi_ghost,location);
        TV normal=interpolation_vector.Clamped_To_Array_Cell(grid,normals_ghost,location);
        new_location=location-phi*normal;
        return box.Lazy_Inside(new_location);}
    else{
        int iterations=1;
        T phi=interpolation.Clamped_To_Array_Cell(grid,phi_ghost,location);
        TV normal=interpolation_vector.Clamped_To_Array_Cell(grid,normals_ghost,location);
        new_location=location-phi*normal;
        if(!box.Lazy_Inside(new_location))return false;
        while(iterations<max_iterations && abs(phi)>tolerance){
            phi=interpolation.Clamped_To_Array_Cell(grid,phi_ghost,new_location);
            TV normal=interpolation_vector.Clamped_To_Array_Cell(grid,normals_ghost,new_location);
            iterations++;new_location-=phi*normal;
            if(!box.Lazy_Inside(new_location))return false;}
        return true;}
}
//#####################################################################
// Function Normal_Flame_Speed
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR DETONATION_SHOCK_DYNAMICS<T_GRID>::
Normal_Flame_Speed(const int axis,const TV_INT& face_index) const
{
    TV_INT cell1,cell2;grid.Cells_Touching_Face(axis,face_index,cell1,cell2);
    if(abs(levelset.phi(cell1))>(T)nb_width*levelset.grid.dX.x || abs(levelset.phi(cell2))>(T)nb_width*levelset.grid.dX.x) return (T)0;
    return (T).5*(Dn.array(cell1)+Dn.array(cell2));
}
//#####################################################################
template class DETONATION_SHOCK_DYNAMICS<GRID<VECTOR<float,1> > >;
template class DETONATION_SHOCK_DYNAMICS<GRID<VECTOR<float,2> > >;
template class DETONATION_SHOCK_DYNAMICS<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DETONATION_SHOCK_DYNAMICS<GRID<VECTOR<double,1> > >;
template class DETONATION_SHOCK_DYNAMICS<GRID<VECTOR<double,2> > >;
template class DETONATION_SHOCK_DYNAMICS<GRID<VECTOR<double,3> > >;
#endif

