//#####################################################################
// Copyright 2002-2007, Doug Enright, Ronald Fedkiw, Jon Gretarsson, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_UNIFORM
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UTILITIES.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_1D_EIGENSYSTEM_F.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_1D_EIGENSYSTEM_F_ADVECTION_ONLY.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_2D_EIGENSYSTEM_F.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_2D_EIGENSYSTEM_F_ADVECTION_ONLY.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_2D_EIGENSYSTEM_G.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_2D_EIGENSYSTEM_G_ADVECTION_ONLY.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_3D_EIGENSYSTEM_F.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_3D_EIGENSYSTEM_F_ADVECTION_ONLY.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_3D_EIGENSYSTEM_G.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_3D_EIGENSYSTEM_G_ADVECTION_ONLY.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_3D_EIGENSYSTEM_H.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_3D_EIGENSYSTEM_H_ADVECTION_ONLY.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> EULER_UNIFORM<T_GRID>::
EULER_UNIFORM(const T_GRID& grid_input) :grid(grid_input),mpi_grid(0),U_ghost(U_ghost_private),psi_pointer(0),timesplit(false),
    use_sound_speed_for_cfl(false),perform_rungekutta_for_implicit_part(false),compute_pressure_fluxes(false),
    thinshell(true), use_sound_speed_based_dt_multiple_for_cfl(false),multiplication_factor_for_sound_speed_based_dt(1.),
    euler_projection(this),apply_cavitation_correction(false),euler_cavitation_density(*this,true,(T).0001),
    euler_cavitation_internal_energy(*this,false,(T).0001),e_min((T)1e-15),last_dt(0),initial_total_conserved_quantity(),
    accumulated_boundary_flux(),ghost_cells_valid(false),ghost_cells_valid_ring(0),need_to_remove_added_internal_energy(false)
{
    assert(!(use_sound_speed_for_cfl && use_sound_speed_based_dt_multiple_for_cfl));
    for(int i=1;i<=T_GRID::dimension;i++){pressure_flux_dimension_indices(i)=1+i;eigensystems[i]=0;eigensystems_default[i]=0;eigensystems_pressureonly[i]=0;}
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> EULER_UNIFORM<T_GRID>::
~EULER_UNIFORM()
{
    for(int i=1;i<=T_GRID::dimension;i++){
        delete eigensystems[i];
        delete eigensystems_default[i];
        delete eigensystems_pressureonly[i];}
    if(psi_pointer) delete psi_pointer;
}
//#####################################################################
// Function Set_Up_Cut_Out_Grid
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Set_Up_Cut_Out_Grid(T_ARRAYS_BOOL& psi_input)
{
    psi_pointer=&psi_input;cut_out_grid=true;
    psi=*psi_pointer;
}
//#####################################################################
// Function Set_Custom_Equation_Of_State
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Set_Custom_Equation_Of_State(EOS<T>& eos_input)
{
    for(int i=1;i<=T_GRID::dimension;i++){
        (dynamic_cast<EULER_EIGENSYSTEM<T_GRID>*>(eigensystems[i]))->Set_Custom_Equation_Of_State(eos_input);
        (dynamic_cast<EULER_EIGENSYSTEM<T_GRID>*>(eigensystems_default[i]))->Set_Custom_Equation_Of_State(eos_input);
        (dynamic_cast<EULER_EIGENSYSTEM<T_GRID>*>(eigensystems_pressureonly[i]))->Set_Custom_Equation_Of_State(eos_input);}
    BASE::Set_Custom_Equation_Of_State(eos_input);
}
//#####################################################################
// Function Set_Custom_Boundary
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Set_Custom_Boundary(T_BOUNDARY& boundary_input)
{
    BASE::Set_Custom_Boundary(boundary_input);
    euler_projection.Set_Constant_Extrapolation_Pressure_Boundary();
}
//#####################################################################
// Function Set_Body_Force
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Set_Body_Force(const bool use_force_input)
{
    use_force=use_force_input;
    if(use_force) force.Resize(grid,1);
    else force.Clean_Memory();
}
//#####################################################################
// Function Initialize_Domain
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Initialize_Domain(const T_GRID& grid_input)
{
    if(grid_input.Is_MAC_Grid()) grid=grid_input;
    else grid=grid_input.Get_MAC_Grid_At_Regular_Positions();
    U.Resize(grid.Domain_Indices()); U_save.Resize(grid.Domain_Indices());
    if(compute_pressure_fluxes && !timesplit) fluxes_pressure.Resize(grid,0);
    if(cut_out_grid) psi=*psi_pointer;
    else{psi.Resize(grid.Domain_Indices());psi.Fill(true);}
    added_internal_energy.Resize(grid.Domain_Indices());
    Set_Eigensystems(timesplit);
    euler_projection.Initialize_Grid();
}
//#####################################################################
// Function Save_State 
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Save_State(T_ARRAYS_DIMENSION_SCALAR& U_s,T_FACE_ARRAYS_SCALAR& face_velocities_s,bool& need_to_remove_added_internal_energy_s)
{
    T_ARRAYS_DIMENSION_SCALAR::Copy(U,U_s);
    euler_projection.Save_State(face_velocities_s);
    need_to_remove_added_internal_energy_s=need_to_remove_added_internal_energy;
}
//#####################################################################
// Function Restore_State 
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Restore_State(T_ARRAYS_DIMENSION_SCALAR& U_s,T_FACE_ARRAYS_SCALAR& face_velocities_s,bool& need_to_remove_added_internal_energy_s)
{
    T_ARRAYS_DIMENSION_SCALAR::Copy(U_s,U);
    euler_projection.Restore_State(face_velocities_s);
    need_to_remove_added_internal_energy=need_to_remove_added_internal_energy_s;
    Invalidate_Ghost_Cells(); // TODO(jontg): Cheaper to save and restore ghost values, too
}
//#####################################################################
// Function Get_Cell_Velocities
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Get_Cell_Velocities(const T dt,const T time,const int ghost_cells,T_ARRAYS_VECTOR& centered_velocities)
{
    Fill_Ghost_Cells(dt,time,ghost_cells);
    for(CELL_ITERATOR iterator(grid,ghost_cells);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        centered_velocities(cell_index)=EULER<T_GRID>::Get_Velocity(U_ghost(cell_index));}
}
//#####################################################################
// Compute_Total_Conserved_Quantity
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Compute_Total_Conserved_Quantity(const bool update_boundary_flux,const T dt,TV_DIMENSION& total_conserved_quantity)
{
    assert(conservation->save_fluxes);
    T cell_size=grid.Cell_Size();
    total_conserved_quantity=TV_DIMENSION();
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(psi(cell_index)) total_conserved_quantity+=U(cell_index)*cell_size;}
    if(update_boundary_flux) for(FACE_ITERATOR iterator(grid,0,T_GRID::BOUNDARY_REGION);iterator.Valid();iterator.Next()){
        const int axis=iterator.Axis();const TV_INT face_index=iterator.Face_Index();
        const T direction=iterator.First_Boundary()?(T)1:(T)-1;
        const TV_INT inside_cell_index=iterator.First_Boundary()?iterator.Second_Cell_Index():iterator.First_Cell_Index();
        if(psi(inside_cell_index)){
                accumulated_boundary_flux+=dt*direction*conservation->fluxes.Component(axis)(face_index)*iterator.Face_Size();
                if(timesplit) accumulated_boundary_flux+=dt*direction*euler_projection.fluxes->Component(axis)(face_index)*iterator.Face_Size();}}
    total_conserved_quantity-=accumulated_boundary_flux;
}
//#####################################################################
// Function Warn_For_Low_Internal_Energy
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Warn_For_Low_Internal_Energy() const
{
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        TV_DIMENSION U_cell=U(iterator.Cell_Index());
        T e=EULER<T_GRID>::e(U_cell);
        if(e<e_min){
            std::stringstream ss;ss<<"WARNING: low internal energy at cell "<<iterator.Cell_Index()<<" with e="<<e<<", state="<<U_cell<<std::endl;LOG::filecout(ss.str());}}
}
//#####################################################################
// Function Invalidate_Ghost_Cells
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Invalidate_Ghost_Cells()
{
    ghost_cells_valid=false;
}
//#####################################################################
// Function Equal_Real_Data
//#####################################################################
template<class T_GRID> bool EULER_UNIFORM<T_GRID>::
Equal_Real_Data(const T_ARRAYS_DIMENSION_SCALAR& U1,const T_ARRAYS_DIMENSION_SCALAR& U2) const
{
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(U1(cell_index)!=U2(cell_index)) return false;}
    return true;
}
//#####################################################################
// Function Apply_Boundary_Conditions
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Fill_Ghost_Cells(const T dt,const T time,const int ghost_cells) const
{
    if(ghost_cells_valid_ring >= ghost_cells && ghost_cells_valid){
        assert(Equal_Real_Data(U,U_ghost_private));
        return;}

    U_ghost_private.Resize(grid.Domain_Indices(ghost_cells));
    boundary->Fill_Ghost_Cells(grid,U,U_ghost_private,dt,time,ghost_cells);
    Clamp_Internal_Energy_Ghost(U_ghost_private,ghost_cells); // TODO(jontg): Is this always / ever necessary?
    ghost_cells_valid_ring=ghost_cells;ghost_cells_valid=true;
}
//#####################################################################
// Function Get_Dirichlet_Boundary_Conditions
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Get_Dirichlet_Boundary_Conditions(const T dt,const T time)
{
    if(timesplit){
        Fill_Ghost_Cells(dt,time,3);
        euler_projection.Get_Dirichlet_Boundary_Conditions(U_ghost);}
}
//#####################################################################
// Function Advance_One_Time_Step_Forces
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Advance_One_Time_Step_Forces(const T dt,const T time)
{
    // update gravity
    if(gravity) for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        for(int axis=1;axis<=T_GRID::dimension;axis++){
            if(downward_direction[axis]) U(cell_index)(axis+1)+=dt*U(cell_index)(1)*gravity*downward_direction[axis];}}

    // update body force
    if(use_force) for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        for(int axis=1;axis<=T_GRID::dimension;axis++) U(cell_index)(axis+1)+=dt*U(cell_index)(1)*force.Component(axis)(cell_index);}

    //TODO: Is this necessary?
    boundary->Apply_Boundary_Condition(grid,U,time+dt);
    Invalidate_Ghost_Cells();
}
//#####################################################################
// Function Compute_Cavitation_Velocity
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Compute_Cavitation_Velocity(T_ARRAYS_SCALAR& rho_n, T_FACE_ARRAYS_SCALAR& face_velocities_n, T_ARRAYS_DIMENSION_SCALAR& momentum_n)
{
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
        if(!euler_projection.elliptic_solver->psi_N.Component(axis)(face_index)){
            TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
            T rho_face_n=(T).5*(rho_n(first_cell_index)+rho_n(second_cell_index));T one_over_rho_face_n=1/rho_face_n;
            T rho_star = (T).5*(U_ghost(first_cell_index)(1)+U_ghost(second_cell_index)(1));
            T momentum_star=(T).5*(U_ghost(first_cell_index)(axis+1)+U_ghost(second_cell_index)(axis+1));
            T momentum_face_n=(T).5*(momentum_n(first_cell_index)(axis)+momentum_n(second_cell_index)(axis));
            face_velocities_n.Component(axis)(face_index)-=one_over_rho_face_n*(rho_star*face_velocities_n.Component(axis)(face_index)+momentum_star-2*momentum_face_n);}}
}
//#####################################################################
// Function Advance_One_Time_Step_Explicit_Part
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Advance_One_Time_Step_Explicit_Part(const T dt,const T time,const int rk_substep,const int rk_order)
{
    last_dt=dt;

    // Store time n density values
    T_ARRAYS_SCALAR rho_n(grid.Domain_Indices(0));T_ARRAYS_DIMENSION_SCALAR momentum_n(grid.Domain_Indices(0));T_FACE_ARRAYS_SCALAR face_velocities_n(grid);
    if(apply_cavitation_correction){
        for(CELL_ITERATOR iterator(grid,0);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            rho_n(cell_index)=U_ghost(cell_index)(1);}

        // Store time n velocities
        euler_projection.Compute_Density_Weighted_Face_Velocities(dt,time,euler_projection.elliptic_solver->psi_N);
        face_velocities_n=euler_projection.face_velocities;

        // Store time n momentum
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            for(int axis=1;axis<=T_GRID::dimension;axis++) momentum_n(cell_index)(axis)=U_ghost(cell_index)(axis+1);}}

    if(compute_pressure_fluxes && !timesplit){
        conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystems,eigensystems_default,euler_projection.elliptic_solver->psi_N,euler_projection.face_velocities,thinshell,
            open_boundaries,&eigensystems_pressureonly,&fluxes_pressure);}
    else{
        if(timesplit && thinshell && rk_substep==1) conservation->Save_Fluxes();
        conservation->Update_Conservation_Law(grid,U,U_ghost,psi,dt,eigensystems,eigensystems_default,euler_projection.elliptic_solver->psi_N,euler_projection.face_velocities,false,
            open_boundaries);}

    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(U_ghost(cell_index)(1)<=0){std::stringstream ss;ss<<"Density at cell "<<cell_index<<" is: "<<U_ghost(cell_index)(1)<<std::endl;LOG::filecout(ss.str());
            PHYSBAM_FATAL_ERROR("Density should not be zero or negative");}}

    boundary->Apply_Boundary_Condition(grid,U,time+dt);
    Invalidate_Ghost_Cells();

    if(apply_cavitation_correction){
        PHYSBAM_DEBUG_WRITE_SUBSTEP("Before cavitation correction",0,1);
        Compute_Cavitation_Velocity(rho_n,face_velocities_n,momentum_n);
        euler_cavitation_density.Apply_Cavitation_Correction(dt,time,face_velocities_n);
        PHYSBAM_DEBUG_WRITE_SUBSTEP("After cavitation density correction",0,1);
        //euler_cavitation_internal_energy.Apply_Cavitation_Correction(dt,time,face_velocities_n);
        //PHYSBAM_DEBUG_WRITE_SUBSTEP("After cavitation energy correction",0,1);
    }
}
//#####################################################################
// Function Advance_One_Time_Step_Implicit_Part
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Advance_One_Time_Step_Implicit_Part(const T dt,const T time)
{
    euler_projection.Compute_Density_Weighted_Face_Velocities(dt,time,euler_projection.elliptic_solver->psi_N);
    euler_projection.Project(euler_projection.face_velocities,dt,time);
    Invalidate_Ghost_Cells();
}
//#####################################################################
// Function Clamp_Internal_Energy
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Clamp_Internal_Energy(const T dt,const T time)
{
    if(apply_cavitation_correction||conservation->adaptive_time_step) return;
    //assert(!need_to_remove_added_internal_energy);
    need_to_remove_added_internal_energy=true;
    added_internal_energy.Fill(0);
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(!psi(cell_index)) continue;
        TV_DIMENSION U_cell=U(cell_index);
        T e=EULER<T_GRID>::e(U_cell);
        if(e<e_min){
            PHYSBAM_FATAL_ERROR("Do not clamp internal energy");
            std::stringstream ss;ss<<"Warning: Clamping internal energy at cell "<<cell_index<<" with e="<<e<<", state="<<U_cell<<std::endl;LOG::filecout(ss.str());
            EULER<T_GRID>::Set_Euler_State_From_rho_velocity_And_internal_energy(U,cell_index,U_cell(1),EULER<T_GRID>::Get_Velocity(U_cell),e_min);
            added_internal_energy(cell_index)=e_min-e;}}
    Invalidate_Ghost_Cells();
}
//#####################################################################
// Function Clamp_Internal_Energy_Ghost
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Clamp_Internal_Energy_Ghost(T_ARRAYS_DIMENSION_SCALAR& U_ghost,const int number_of_ghost_cells) const
{
    if(apply_cavitation_correction||conservation->adaptive_time_step) return;
    for(int axis=1;axis<=TV::dimension;++axis) for(int axis_side=1;axis_side<=2;++axis_side) if(!mpi_grid || !mpi_grid->Neighbor(axis,axis_side)){
        for(CELL_ITERATOR iterator(grid,number_of_ghost_cells,T_GRID::GHOST_REGION,2*axis-(2-axis_side));iterator.Valid();iterator.Next()){
            TV_DIMENSION U_cell=U_ghost(iterator.Cell_Index());
            T e=EULER<T_GRID>::e(U_cell);
            if(e<e_min){
                PHYSBAM_FATAL_ERROR("Do not clamp internal energy e2");
                std::stringstream ss;ss<<"Warning: Clamping internal energy at ghost cell "<<iterator.Cell_Index()<<" with e="<<e<<", state="<<U_cell<<std::endl;LOG::filecout(ss.str());
                EULER<T_GRID>::Set_Euler_State_From_rho_velocity_And_internal_energy(U_ghost,iterator.Cell_Index(),U_cell(1),EULER<T_GRID>::Get_Velocity(U_cell),e_min);}}}
}
//#####################################################################
// Function Remove_Added_Internal_Energy
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Remove_Added_Internal_Energy(const T dt,const T time)
{
    if(apply_cavitation_correction||conservation->adaptive_time_step) return;
    //assert(need_to_remove_added_internal_energy);
    need_to_remove_added_internal_energy=false;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(!psi(cell_index)) continue;
        T added_e=added_internal_energy(cell_index);
        if(added_e){
            PHYSBAM_FATAL_ERROR("Did not clamp internal energy");
            TV_DIMENSION U_cell=U(cell_index);
            T e=EULER<T_GRID>::e(U_cell);
            EULER<T_GRID>::Set_Euler_State_From_rho_velocity_And_internal_energy(U,cell_index,U_cell(1),EULER<T_GRID>::Get_Velocity(U_cell),e-added_e);}}
    Invalidate_Ghost_Cells();
}
//#####################################################################
// Function Euler_Step
//#####################################################################
// providing time is optional
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Euler_Step(const T dt,const T time,const bool is_time_n)
{
    Advance_One_Time_Step_Forces(dt,time);
    Advance_One_Time_Step_Explicit_Part(dt,time,1,1);
    if(timesplit) Advance_One_Time_Step_Implicit_Part(dt,time);
}
//#####################################################################
// Function CFL_Using_Sound_Speed
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR EULER_UNIFORM<T_GRID>::
CFL_Using_Sound_Speed() const
{
    T_ARRAYS_VECTOR velocity(grid.Domain_Indices()),velocity_minus_c(grid.Domain_Indices()),velocity_plus_c(grid.Domain_Indices());
    T max_sound_speed=0;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(psi(cell_index)){
            TV velocity_cell=Get_Velocity(U,cell_index);
            T sound_speed=eos->c(U(cell_index)(1),e(U,cell_index));
            max_sound_speed=max(max_sound_speed,sound_speed);
            velocity(cell_index)=velocity_cell;
            velocity_minus_c(cell_index)=velocity_cell-sound_speed*TV::All_Ones_Vector();velocity_plus_c(cell_index)=velocity_cell+sound_speed*TV::All_Ones_Vector();}}

    TV max_velocity_minus_c=velocity_minus_c.Componentwise_Maxabs();TV max_velocity_plus_c=velocity_plus_c.Componentwise_Maxabs();
    T dt_convect=0;for(int axis=1;axis<=T_GRID::dimension;axis++) dt_convect+=max(max_velocity_minus_c(axis),max_velocity_plus_c(axis))/grid.dX[axis];
    dt_convect=max(dt_convect,1/max_time_step);
    std::stringstream ss;ss<<"max sound speed="<<max_sound_speed<<std::endl;
    ss<<"max velocity="<<velocity.Componentwise_Maxabs()<<std::endl;LOG::filecout(ss.str());
    return 1/dt_convect;
}
//#####################################################################
// Function CFL
//#####################################################################
template<class T_GRID> typename T_GRID::SCALAR EULER_UNIFORM<T_GRID>::
CFL(const T time) const
{
    // TODO: add gravity and body force contribution
    static T sum_ratios=0;
    static int count_ratios=0;
    TV one_over_dx=grid.one_over_dX;
    if(!use_sound_speed_for_cfl){
        TV max_lambdas;TV max_grad_p_over_rho,max_grad_p;

        Fill_Ghost_Cells(last_dt,time,3);
        T_ARRAYS_SCALAR p_approx(grid.Domain_Indices(1));
        for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next())
            p_approx(iterator.Cell_Index())=eos->p(U_ghost(iterator.Cell_Index())(1),e(U_ghost,iterator.Cell_Index()));
        T_FACE_ARRAYS_SCALAR p_approx_face(grid);
        euler_projection.Compute_Face_Pressure_From_Cell_Pressures(grid,U_ghost,psi,p_approx_face,p_approx);
        T_ARRAYS_VECTOR grad_p_approx(grid.Domain_Indices());ARRAYS_UTILITIES<T_GRID,T>::Compute_Gradient_At_Cells_From_Face_Data(grid,grad_p_approx,p_approx_face);

        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){const TV_INT cell_index=iterator.Cell_Index();
            if(psi(cell_index)) for(int axis=1;axis<=T_GRID::dimension;axis++){
                max_lambdas[axis]=max(max_lambdas[axis],eigensystems[axis]->Maximum_Magnitude_Eigenvalue(U(cell_index)));
                max_grad_p[axis]=maxabs(max_grad_p[axis],grad_p_approx(cell_index)[axis]);
                max_grad_p_over_rho[axis]=maxabs(max_grad_p_over_rho[axis],grad_p_approx(cell_index)[axis]/U(cell_index)(1));}}

        T dt_convect=0;for(int axis=1;axis<=T_GRID::dimension;axis++){
            T max_u_over_dx=max_lambdas[axis]*one_over_dx[axis];
            dt_convect+=max_u_over_dx+sqrt(max_u_over_dx*max_u_over_dx+((T)4*max_grad_p_over_rho[axis])*one_over_dx[axis]);}
        dt_convect=(T).5*dt_convect;
        T one_over_dt=max(dt_convect,1/max_time_step);

        T dt_old=CFL_Using_Sound_Speed();
        T ratio_new_dt_to_old_dt=1/(one_over_dt*dt_old);
        sum_ratios+=ratio_new_dt_to_old_dt;count_ratios++;
        std::stringstream ss;ss<<"dt after CFL="<<1/one_over_dt<<std::endl;
        ss<<"Ratio of new dt to old dt="<<ratio_new_dt_to_old_dt<<", average ratio till now="<<(sum_ratios/(T)count_ratios)<<std::endl;
        ss<<"max_lambdas= ";for(int axis=1;axis<=T_GRID::dimension;axis++) ss<<max_lambdas[axis]<<", ";ss<<std::endl;
        ss<<"max_grad_p_over_rho= ";for(int axis=1;axis<=T_GRID::dimension;axis++) ss<<max_grad_p_over_rho[axis]<<", ";ss<<std::endl;
        ss<<"max_grad_p= ";for(int axis=1;axis<=T_GRID::dimension;axis++) ss<<max_grad_p[axis]<<", ";ss<<std::endl;LOG::filecout(ss.str());

        // Try going to a higher cfl number on time step computed using sound speed
        if(use_sound_speed_based_dt_multiple_for_cfl && (ratio_new_dt_to_old_dt>multiplication_factor_for_sound_speed_based_dt)){
            T dt_new=multiplication_factor_for_sound_speed_based_dt*dt_old;
            ss.clear();ss<<"clamping ratio to "<<multiplication_factor_for_sound_speed_based_dt<<". New dt after CFL="<<dt_new<<"ratio ="<<dt_new/dt_old<<std::endl;LOG::filecout(ss.str());
            return dt_new;}
        else return 1/one_over_dt;}
    else return CFL_Using_Sound_Speed();
}
//#####################################################################
// Function Set_Eigensystems
//#####################################################################
template<class T_GRID,class T> void Set_Eigensystems_Helper(VECTOR<EIGENSYSTEM<T,VECTOR<T,3> >*,1>& eigensystems_default,VECTOR<EIGENSYSTEM<T,VECTOR<T,3> >*,1>& eigensystems,VECTOR<EIGENSYSTEM<T,VECTOR<T,3> >*,1>& eigensystems_pressureonly,
    const EULER_PROJECTION_UNIFORM<T_GRID>& euler_projection,const bool advection_only)
{
    if(eigensystems_default[1]) delete eigensystems_default[1];eigensystems_default[1]=new EULER_1D_EIGENSYSTEM_F<T>();
    if(eigensystems_pressureonly[1]) delete eigensystems_pressureonly[1];eigensystems_pressureonly[1]=new EULER_1D_EIGENSYSTEM_F<T>(true);
    if(eigensystems[1]) delete eigensystems[1];
    if(advection_only) eigensystems[1]=new EULER_1D_EIGENSYSTEM_F_ADVECTION_ONLY<T>();
    else eigensystems[1]=new EULER_1D_EIGENSYSTEM_F<T>();
}
template<class T_GRID,class T> void Set_Eigensystems_Helper(VECTOR<EIGENSYSTEM<T,VECTOR<T,4> >*,2>& eigensystems_default,VECTOR<EIGENSYSTEM<T,VECTOR<T,4> >*,2>& eigensystems,VECTOR<EIGENSYSTEM<T,VECTOR<T,4> >*,2>& eigensystems_pressureonly,
    const EULER_PROJECTION_UNIFORM<T_GRID>& euler_projection,const bool advection_only)
{
    if(eigensystems_default[1]) delete eigensystems_default[1];eigensystems_default[1]=new EULER_2D_EIGENSYSTEM_F<T>();
    if(eigensystems_default[2]) delete eigensystems_default[2];eigensystems_default[2]=new EULER_2D_EIGENSYSTEM_G<T>();
    if(eigensystems_pressureonly[1]) delete eigensystems_pressureonly[1];eigensystems_pressureonly[1]=new EULER_2D_EIGENSYSTEM_F<T>(true);
    if(eigensystems_pressureonly[2]) delete eigensystems_pressureonly[2];eigensystems_pressureonly[2]=new EULER_2D_EIGENSYSTEM_G<T>(true);
    if(eigensystems[1]) delete eigensystems[1];
    if(eigensystems[2]) delete eigensystems[2];
    if(advection_only){
        eigensystems[1]=new EULER_2D_EIGENSYSTEM_F_ADVECTION_ONLY<T>();
        eigensystems[2]=new EULER_2D_EIGENSYSTEM_G_ADVECTION_ONLY<T>();}
    else{
        eigensystems[1]=new EULER_2D_EIGENSYSTEM_F<T>();
        eigensystems[2]=new EULER_2D_EIGENSYSTEM_G<T>();}
}
template<class T_GRID,class T> void Set_Eigensystems_Helper(VECTOR<EIGENSYSTEM<T,VECTOR<T,5> >*,3>& eigensystems_default,VECTOR<EIGENSYSTEM<T,VECTOR<T,5> >*,3>& eigensystems,VECTOR<EIGENSYSTEM<T,VECTOR<T,5> >*,3>& eigensystems_pressureonly,
    const EULER_PROJECTION_UNIFORM<T_GRID>& euler_projection,const bool advection_only)
{
    if(eigensystems_default[1]) delete eigensystems_default[1];eigensystems_default[1]=new EULER_3D_EIGENSYSTEM_F<T>();
    if(eigensystems_default[2]) delete eigensystems_default[2];eigensystems_default[2]=new EULER_3D_EIGENSYSTEM_G<T>();
    if(eigensystems_default[3]) delete eigensystems_default[3];eigensystems_default[3]=new EULER_3D_EIGENSYSTEM_H<T>();
    if(eigensystems_pressureonly[1]) delete eigensystems_pressureonly[1];eigensystems_pressureonly[1]=new EULER_3D_EIGENSYSTEM_F<T>(true);
    if(eigensystems_pressureonly[2]) delete eigensystems_pressureonly[2];eigensystems_pressureonly[2]=new EULER_3D_EIGENSYSTEM_G<T>(true);
    if(eigensystems_pressureonly[3]) delete eigensystems_pressureonly[3];eigensystems_pressureonly[3]=new EULER_3D_EIGENSYSTEM_H<T>(true);
    if(eigensystems[1]) delete eigensystems[1];
    if(eigensystems[2]) delete eigensystems[2];
    if(eigensystems[3]) delete eigensystems[3];
    if(advection_only){
        eigensystems[1]=new EULER_3D_EIGENSYSTEM_F_ADVECTION_ONLY<T>();
        eigensystems[2]=new EULER_3D_EIGENSYSTEM_G_ADVECTION_ONLY<T>();
        eigensystems[3]=new EULER_3D_EIGENSYSTEM_H_ADVECTION_ONLY<T>();}
    else{
        eigensystems[1]=new EULER_3D_EIGENSYSTEM_F<T>();
        eigensystems[2]=new EULER_3D_EIGENSYSTEM_G<T>();
        eigensystems[3]=new EULER_3D_EIGENSYSTEM_H<T>();}
}

template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Set_Eigensystems(const bool advection_only)
{
    Set_Eigensystems_Helper(eigensystems_default,eigensystems,eigensystems_pressureonly,euler_projection,advection_only);
    Set_Custom_Equation_Of_State(*eos);
}
//#####################################################################
// Function Log_Parameters
//#####################################################################
template<class T_GRID> void EULER_UNIFORM<T_GRID>::
Log_Parameters() const
{
    LOG::SCOPE scope("EULER_UNIFORM parameters");
    BASE::Log_Parameters();
    std::stringstream ss;ss<<"timesplit="<<timesplit<<std::endl;
    ss<<"use_sound_speed_for_cfl="<<use_sound_speed_for_cfl<<std::endl;
    ss<<"perform_rungekutta_for_implicit_part="<<perform_rungekutta_for_implicit_part<<std::endl;
    ss<<"use_sound_speed_based_dt_multiple_for_cfl="<<use_sound_speed_based_dt_multiple_for_cfl<<std::endl;
    ss<<"multiplication_factor_for_sound_speed_based_dt="<<multiplication_factor_for_sound_speed_based_dt<<std::endl;
    ss<<"e_min="<<e_min<<std::endl;
    ss<<"apply_cavitation_correction="<<apply_cavitation_correction<<std::endl;LOG::filecout(ss.str());
    if(timesplit) euler_projection.Log_Parameters();
}
//#####################################################################
template class EULER_UNIFORM<GRID<VECTOR<float,1> > >;
template class EULER_UNIFORM<GRID<VECTOR<float,2> > >;
template class EULER_UNIFORM<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EULER_UNIFORM<GRID<VECTOR<double,1> > >;
template class EULER_UNIFORM<GRID<VECTOR<double,2> > >;
template class EULER_UNIFORM<GRID<VECTOR<double,3> > >;
#endif
