#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Frank Losasso, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Advection/ADVECTION.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/FACE_LOOKUP_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/EXTRAPOLATION_DYADIC.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_DYADIC.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> INCOMPRESSIBLE_DYADIC<T_GRID>::
INCOMPRESSIBLE_DYADIC(T_GRID& grid_input,const bool flame)
    :grid(grid_input),projection(grid,&flame_speed),boundary_default(*new BOUNDARY_DYADIC<T_GRID,T>)
{
    boundary=&boundary_default;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> INCOMPRESSIBLE_DYADIC<T_GRID>::
~INCOMPRESSIBLE_DYADIC()
{
    delete &boundary_default;
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_DYADIC<T_GRID>::
Euler_Step(const T dt,const T time)
{
    Advance_One_Time_Step(dt,time);
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_DYADIC<T_GRID>::
Advance_One_Time_Step(const T dt,const T time,const bool implicit_viscosity,const bool crank_nicolson,const bool first_step,const ARRAY<T>* phi_ghost)
{
    assert(!implicit_viscosity&&!crank_nicolson&&!first_step);
    LOG::Time("Incompressible Explicit Part");
    Advance_One_Time_Step_Explicit_Part(dt,time,implicit_viscosity,crank_nicolson,first_step,phi_ghost);
    LOG::Time("Incompressible Implicit Part");
    //Average_Delta_Node_Velocities_To_Faces();
    Advance_One_Time_Step_Implicit_Part(dt,time);
    //Average_Face_Velocities_To_Nodes();
    boundary->Apply_Boundary_Condition(grid,projection.face_velocities,time+dt);
}
//#####################################################################
// Function Advance_One_Time_Step_Explicit_Part
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_DYADIC<T_GRID>::
Advance_One_Time_Step_Explicit_Part(const T dt,const T time,const bool implicit_viscosity,const bool crank_nicolson,const bool first_step,const ARRAY<T>* phi_ghost)
{
    assert(!implicit_viscosity&&!crank_nicolson);
    // find ghost cells
    ARRAY<T>& face_velocities=projection.face_velocities;
    ARRAY<T> face_velocities_ghost(grid.number_of_faces);
    boundary->Fill_Ghost_Cells_Face(grid,face_velocities,face_velocities_ghost,time);
    // update convection
    advection->Update_Advection_Equation_Face(grid,face_velocities,face_velocities_ghost,face_velocities_ghost,*boundary,dt,time);
    
    // update gravity
    if(gravity){
        TV dt_times_gravity_vector=dt*gravity*downward_direction;
        for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next())
            face_velocities(iterator.Face_Index())+=dt_times_gravity_vector[iterator.Axis()+1];}

    // update body force
    if(use_force)
        for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next())
            face_velocities(iterator.Face_Index())+=dt*force(iterator.Face_Index());

    // TODO: CAN WE MAKE CONFINEMENT LOOK MORE LIKE THE UNIFORM CASE?
    if(vorticity_confinement || use_variable_vorticity_confinement){
        ARRAY<TV> F(grid.number_of_nodes,false);
        // TODO: make this cell based
        ARRAY<TV> V_ghost(grid.number_of_nodes,false);LINEAR_INTERPOLATION_DYADIC_HELPER<T_GRID>::Interpolate_From_Faces_To_Nodes(grid,face_velocities_ghost,V_ghost);
        Compute_Vorticity_Confinement_Force(grid,V_ghost,F);
        if(use_variable_vorticity_confinement){F*=dt;F*=variable_vorticity_confinement;}else F*=dt*vorticity_confinement;
        for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
            T face_force=0;for(int i=0;i<T_GRID::number_of_nodes_per_face;i++)face_force+=F(iterator.Face_Node(i))[iterator.Axis()+1];face_force*=T_GRID::one_over_number_of_nodes_per_face;
            face_velocities(iterator.Face_Index())+=face_force;}}

    boundary->Apply_Boundary_Condition(grid,face_velocities,time+dt);

    // Set ghost velocities to zero to avoid accumulating with no bound
    for(FACE_ITERATOR iterator(grid,grid.Map_Ghost_Faces());iterator.Valid();iterator.Next())face_velocities(iterator.Face_Index())=0;
}
//#####################################################################
// Function Advance_One_Time_Step_Implicit_Part
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_DYADIC<T_GRID>::
Advance_One_Time_Step_Implicit_Part(const T dt,const T time,const bool implicit_viscosity,const bool crank_nicolson)
{
    assert(!implicit_viscosity&&!crank_nicolson);
    // make divergence free
    if(projection.flame){
        PHYSBAM_NOT_IMPLEMENTED(); //projection.poisson->Set_Constant_beta(T(1)/projection.density_fuel,T(1)/projection.density);
        Calculate_Flame_Speed();
        Calculate_Pressure_Jump(dt,time);}
    projection.Make_Divergence_Free(dt,time);
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR INCOMPRESSIBLE_DYADIC<T_GRID>::
CFL(const bool inviscid,const bool viscous_only) const
{
    const ARRAY<T>& face_velocities=projection.face_velocities;
    T dx=grid.Minimum_Edge_Length();
    TV max_abs_V;
    PHYSBAM_NOT_IMPLEMENTED(); // need to update convection CFL to look like levelset CFL (do this after cleaning up the ghost values stuff)
    if(!viscous_only){
        if(projection.flame) for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
            T face_velocity=projection.Face_Velocity_With_Ghost_Value(projection.face_velocities,iterator.Face_Index(),1);
            T face_velocity_fuel=projection.Face_Velocity_With_Ghost_Value(projection.face_velocities,iterator.Face_Index(),-1);
            max_abs_V[iterator.Axis()+1]=max(max_abs_V[iterator.Axis()+1],abs(face_velocity),abs(face_velocity_fuel));}
        else for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next())
            max_abs_V[iterator.Axis()+1]=max(max_abs_V[iterator.Axis()+1],abs(face_velocities(iterator.Face_Index())));}
    T dt_convection=max_abs_V.L1_Norm()/dx;
    T dt_viscosity=0;
    if(!inviscid){ // TODO: CFL for variable viscosity (once we add variable_viscosity support to class
        if(viscosity) dt_viscosity=viscosity/projection.density*T_GRID::number_of_faces_per_cell/sqr(dx);
        if(projection.flame) PHYSBAM_NOT_IMPLEMENTED();/*dt_viscosity=max(dt_viscosity,viscosity_fuel/projection.density_fuel*T_GRID::number_of_faces_per_cell/sqr(dx));*/}
    if(viscous_only) return 1/max(dt_viscosity,1/max_time_step);
    TV max_force;
    if(use_force)for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next())
            max_force[iterator.Axis()+1]=max(max_force[iterator.Axis()+1],abs(force(iterator.Face_Index())));
    T dt_force=0;
    if(use_force) dt_force=max_force.L1_Norm()/dx;
    if(gravity) dt_force+=abs(gravity)*(downward_direction.L1_Norm())/dx;
    //if(strain) dt_force+=1/strain->CFL(); // TODO: uncomment once strain is implemented
    T dt_overall=(dt_convection+dt_viscosity+sqrt(sqr(dt_convection+dt_viscosity)+4*dt_force))/2; 
    return 1/max(dt_overall,1/max_time_step);
}
//#####################################################################
// Function Calculate_Pressure_Jump
//#####################################################################
// flame_speed must be up to date
template<class T_GRID> void INCOMPRESSIBLE_DYADIC<T_GRID>::
Calculate_Pressure_Jump(const T dt,const T time)
{
    PHYSBAM_NOT_IMPLEMENTED();
    //T dt_one_over_density_jump_density_fuel_squared=dt*(1/projection.density_fuel-1/projection.density)*sqr(projection.density_fuel); // [p]=dt*[1/density]*sqr(M) with M=-density_fuel*flame_speed
    //for(int i=1;i<=grid.number_of_cells;i++) projection.poisson->u_jump(i)=dt_one_over_density_jump_density_fuel_squared*sqr(flame_speed(i)); 
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_DYADIC<T_GRID>::
Set_Dirichlet_Boundary_Conditions(const ARRAY<T>& phi,const T pressure)
{
    for(DYADIC_GRID_ITERATOR_CELL<T_GRID> iterator(grid,0);iterator.Valid();iterator.Next())
        if(phi(iterator.Cell_Index())>0){projection.elliptic_solver->psi_D(iterator.Cell_Index())=true;projection.p(iterator.Cell_Index())=pressure;}
}
//#####################################################################
// Function Calculate_Flame_Speed
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_DYADIC<T_GRID>::
Calculate_Flame_Speed()
{
    PHYSBAM_NOT_IMPLEMENTED();
/*    if(!projection.flame) return;
    if(projection.normal_flame_speed_only)ARRAY<T>::copy(projection.normal_flame_speed,flame_speed);
    else{
        projection.poisson->levelset->Compute_Curvature();
        flame_speed=*projection.poisson->levelset->curvature;
        flame_speed*=projection.curvature_flame_speed;flame_speed+=projection.normal_flame_speed;}*/
}
//#####################################################################
// Function Extrapolate_Velocity_Across_Interface
//#####################################################################
template<class T_GRID> void INCOMPRESSIBLE_DYADIC<T_GRID>::
Extrapolate_Velocity_Across_Interface(ARRAY<T>& phi_ghost,const bool enforce_divergence_free,const T band_width,const T damping,const TV& air_speed,const T_FACE_NEIGHBORS* face_neighbors_visible)
{
    T delta=(T)band_width*grid.Minimum_Edge_Length(); 

    ARRAY<T> phi_face(grid.number_of_faces);ARRAY<bool> fixed_face(grid.number_of_faces);
    // this averaging will only give good values when the cells are the same size, which is fine because of the bandwidth
    for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
        T phi1=phi_ghost(iterator.Deepest_Cell_Index()),phi2=phi_ghost(iterator.Other_Cell_Index());
        phi_face(iterator.Face_Index())=(T).5*(phi1+phi2);
        if(phi1<=0||phi2<=0)fixed_face(iterator.Face_Index())=true;
        if(phi_face(iterator.Face_Index())>=delta)projection.face_velocities(iterator.Face_Index())=(T)0;}
    EXTRAPOLATION_DYADIC<T_GRID,TV> extrapolate(grid);extrapolate.Set_Band_Width(band_width);
    if(face_neighbors_visible) extrapolate.Set_Collision_Aware_Extrapolation(*face_neighbors_visible); // set for collision aware extrapolation
    extrapolate.Set_Custom_Seed_Done(&fixed_face);extrapolate.Extrapolate_Faces(phi_face,projection.face_velocities);
    //if(damping)for(i=1;i<=grid.number_of_nodes;i++)if(phi_nodes_ghost(i)>0&&phi_nodes_ghost(i)<delta)V(i)=(1-damping)*V(i)+damping*air_speed; // TODO

    // make extrapolated velocity divergence free
    if(enforce_divergence_free){
        ARRAY<T> p_new(grid.number_of_cells);ARRAY<bool> psi_D_new(grid.number_of_cells);ARRAY<bool> psi_N_new(grid.number_of_faces);
        ARRAY<T>::Exchange_Arrays(p_new,projection.p);ARRAY<bool>::Exchange_Arrays(psi_D_new,projection.elliptic_solver->psi_D);ARRAY<bool>::Exchange_Arrays(psi_N_new,projection.elliptic_solver->psi_N);
        projection.elliptic_solver->Set_Dirichlet_Outer_Boundaries();
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int index=iterator.Cell_Index();
            if(phi_ghost(index) >= delta) projection.elliptic_solver->psi_D(index)=true;
            else if(phi_ghost(index) <= 0){for(int i=0;i<T_GRID::number_of_faces_per_cell;i++)projection.elliptic_solver->psi_N(iterator.Cell_Pointer()->Face(i))=true;}
            else{
                bool local_maximum=true;
                for(int i=1;i<=T_GRID::number_of_neighbors_per_cell;i++) if(phi_ghost(index)<phi_ghost(iterator.Cell_Neighbor(i)->Cell())){local_maximum=false;break;}
                if(local_maximum)projection.elliptic_solver->psi_D(index)=true;}}
        projection.Make_Divergence_Free(0,0); // TODO: use real dt/time for divergence free extrapolation
        ARRAY<T>::Exchange_Arrays(p_new,projection.p);ARRAY<bool>::Exchange_Arrays(psi_D_new,projection.elliptic_solver->psi_D);ARRAY<bool>::Exchange_Arrays(psi_N_new,projection.elliptic_solver->psi_N);} // restore pressure for use as initial guess for incompressible projection
}
//#####################################################################
template class INCOMPRESSIBLE_DYADIC<OCTREE_GRID<float> >;
template class INCOMPRESSIBLE_DYADIC<QUADTREE_GRID<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class INCOMPRESSIBLE_DYADIC<OCTREE_GRID<double> >;
template class INCOMPRESSIBLE_DYADIC<QUADTREE_GRID<double> >;
#endif
#endif
