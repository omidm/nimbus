//#####################################################################
// Copyright 2011, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Fluids/PhysBAM_Compressible/COMPRESSIBLE_DRIVER.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/COMPRESSIBLE_EXAMPLE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMPRESSIBLE_DRIVER<TV>::
COMPRESSIBLE_DRIVER(COMPRESSIBLE_EXAMPLE<TV>& example_input)
    : BASE(example_input),example(example_input),restart_dt(0),reset_with_restart(false)
{}
//#####################################################################
// Compute_Dt
//#####################################################################
template<class TV> typename TV::SCALAR COMPRESSIBLE_DRIVER<TV>::
Compute_Dt(const T time,const T target_time,bool& done)
{
    T dt=target_time-time;
    example.Limit_Dt(dt,time);
    if((restart_dt!=0) && reset_with_restart && example.adaptive_time_step){example.Restore_State();
        dt=restart_dt;restart_dt=0;reset_with_restart=false;
        if(dt<example.global_min_dt && example.clamp_fluxes){dt=example.global_min_dt;example.conservation_law_solver->clamp_fluxes=true;}
        std::stringstream ss;ss<<"Corrected dt to "<<dt<<std::endl;LOG::filecout(ss.str());}
    return dt;
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void COMPRESSIBLE_DRIVER<TV>::
Initialize()
{
    BASE::Initialize();
    example.Parse_Late_Options();
    example.Log_Parameters();

    if(example.restart) example.Read_Output_Files(example.restart_frame);
    else example.Initialize_Fluid_State();
    example.Initialize_Boundaries();

    // mpi
    if(example.mpi_grid){example.mpi_grid->Initialize(example.domain_walls);
        example.boundary=new BOUNDARY_MPI<T_GRID,TV_DIMENSION>(example.mpi_grid,*example.boundary_scalar);}
    else example.boundary=example.boundary_scalar;

    // collision objects
    example.Initialize_Bodies();
    example.conservation_law_solver->Set_Callbacks(&example);

    if(example.solid_affects_fluid){
        example.kinematic_coupling_utilities->Initialize_Solid_Fluid_Coupling(example.collision_bodies_affecting_fluid);
        example.kinematic_coupling_utilities->Set_Boundary(example.boundary);
        example.kinematic_coupling_utilities->Fill_Solid_Cells((T)0.,(T)0.);
        example.kinematic_coupling_utilities->Update_Cut_Out_Grid();}
}
//#####################################################################
// Update_Bodies
//#####################################################################
template<class TV> void COMPRESSIBLE_DRIVER<TV>::
Update_Bodies(const T dt,const T time)
{
    for(int i=1;i<=example.solid_body_collection.rigid_body_collection.kinematic_rigid_bodies.m;i++){
        RIGID_BODY<TV>& rigid_body=example.solid_body_collection.rigid_body_collection.Rigid_Body(i);
        rigid_body.X()+=rigid_body.V()*dt;
        rigid_body.Rotation()=ROTATION<TV>::From_Rotation_Vector(rigid_body.Angular_Velocity()*dt)*rigid_body.Rotation();
        rigid_body.Rotation().Normalize();
        rigid_body.Update_Bounding_Box();}
}
//#####################################################################
// Function Calculate_Maximum_Allowable_dt
//#####################################################################
template<class TV> void COMPRESSIBLE_DRIVER<TV>::
Calculate_Maximum_Allowable_dt(const T dt,T& min_dt,const int substep,RUNGEKUTTA<T_ARRAYS_DIMENSION_SCALAR>& rungekutta_u)
{
    T val=2;
    if((substep==2)&&(rungekutta_u.order==3)) val=4;
    else if((substep==3)&&(rungekutta_u.order==3)) val=1.5;

    ARRAY_VIEW<TV_DIMENSION,TV_INT> U_n(example.grid.Domain_Indices(),reinterpret_cast<TV_DIMENSION*>(rungekutta_u.u_copy.Get_Array_Pointer()));
    for(CELL_ITERATOR iterator(example.grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        T clamp_rho_cell=example.conservation_law_solver->clamp_rho*U_n(cell_index)(1);
        T clamp_e_cell=example.conservation_law_solver->clamp_e*EULER<T_GRID>::e(U_n,cell_index);
        if(example.U(cell_index)(1)<U_n(cell_index)(1) && abs(example.U(cell_index)(1)-U_n(cell_index)(1))>1e-5)
            min_dt=min(min_dt,((T)val*dt*(clamp_rho_cell-U_n(cell_index)(1)))/(example.U(cell_index)(1)-U_n(cell_index)(1)));
        assert(min_dt>0);
        if(EULER<T_GRID>::e(example.U,cell_index)<EULER<T_GRID>::e(U_n,cell_index) && abs(EULER<T_GRID>::e(example.U,cell_index)-EULER<T_GRID>::e(U_n,cell_index))>1e-5)
            min_dt=min(min_dt,((T)val*dt*(clamp_e_cell-EULER<T_GRID>::e(U_n,cell_index)))/(EULER<T_GRID>::e(example.U,cell_index)-EULER<T_GRID>::e(U_n,cell_index)));
        assert(min_dt>0);}
    if(example.mpi_grid) example.mpi_grid->Synchronize_Dt(min_dt);
    if(min_dt!=dt){rungekutta_u.RUNGEKUTTA_CORE<T>::u-=rungekutta_u.RUNGEKUTTA_CORE<T>::u_copy;
        rungekutta_u.RUNGEKUTTA_CORE<T>::u=min_dt*(rungekutta_u.RUNGEKUTTA_CORE<T>::u/(val*dt));
        rungekutta_u.RUNGEKUTTA_CORE<T>::u+=rungekutta_u.RUNGEKUTTA_CORE<T>::u_copy;}
}
//#####################################################################
// Advance_To_Target_Time
//#####################################################################
template<class TV> void COMPRESSIBLE_DRIVER<TV>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;
    for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP","substep %d",substep,1);
        Preprocess_Substep(current_frame,substep);
        T dt=Compute_Dt(time,target_time,done);
        if(example.mpi_grid) example.mpi_grid->Synchronize_Dt(dt);
        EXAMPLE<TV>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done);

        Update_Bodies(dt,time);

        RUNGEKUTTA<T_ARRAYS_DIMENSION_SCALAR> rungekutta_U(example.U);
        rungekutta_U.Set_Grid_And_Boundary_Condition(example.grid,*example.boundary);
        rungekutta_U.Set_Order(example.rk_order);rungekutta_U.Set_Time(time);rungekutta_U.Start(dt);T rk_time=time;
        example.Save_State();

        for(int rk_substep=1;rk_substep<=rungekutta_U.order;++rk_substep){
            example.boundary->Fill_Ghost_Cells(example.grid,example.U,example.U_ghost,dt,time,3);
            // example.conservation_law_solver->Update_Conservation_Law(example.grid,example.U,example.U_ghost,example.psi,dt,example.cfl_eigensystem);
            example.conservation_law_solver->Update_Conservation_Law(example.grid,example.U,example.U_ghost,example.psi,dt,example.cfl_eigensystem,example.cfl_eigensystem,example.psi_N,example.face_velocities);

            for(CELL_ITERATOR iterator(example.grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
                PHYSBAM_ASSERT(example.U(cell)(1)>0 && EULER<T_GRID>::e(example.U,cell)>0);}

            if(example.mpi_grid) example.mpi_grid->Synchronize_Dt(example.conservation_law_solver->min_dt);
            
            if((example.conservation_law_solver->min_dt<dt) && example.adaptive_time_step){
                if(rk_substep==1){restart_dt=example.conservation_law_solver->min_dt;reset_with_restart=true;break;}
                else if((rk_substep==2)&&(rungekutta_U.order==2)){T min_dt=dt;Calculate_Maximum_Allowable_dt(dt,min_dt,rk_substep,rungekutta_U);restart_dt=min_dt;break;}
                else if((rk_substep==2)&&(rungekutta_U.order==3)){T min_dt=dt;Calculate_Maximum_Allowable_dt(dt,min_dt,rk_substep,rungekutta_U);*const_cast<T*>(&dt)=min_dt;rk_time-=dt/2;continue;}
                else if((rk_substep==3)&&(rungekutta_U.order==3)){T min_dt=dt;Calculate_Maximum_Allowable_dt(dt,min_dt,rk_substep,rungekutta_U);restart_dt=min_dt;break;}}

            example.boundary->Apply_Boundary_Condition(example.grid,example.U,rk_time+dt);
            rk_time=rungekutta_U.Main();}
        if((restart_dt!=0) && reset_with_restart && example.adaptive_time_step) continue;
        
        if(example.solid_affects_fluid){
            example.kinematic_coupling_utilities->Update_Cut_Out_Grid();
            example.kinematic_coupling_utilities->Fill_Solid_Cells(dt,time);}

        Postprocess_Substep(dt,time);
        if(!done) Write_Substep("END Substep",substep,0);
        time+=dt;}
}
//#####################################################################
// Simulate_To_Frame
//#####################################################################
template<class TV> void COMPRESSIBLE_DRIVER<TV>::
Simulate_To_Frame(const int target_frame)
{
    example.frame_title=STRING_UTILITIES::string_sprintf("Frame %d",current_frame);
    if(!example.restart) Write_Output_Files(current_frame);

    while(current_frame<target_frame){
        LOG::SCOPE scope("FRAME","Frame %d",++current_frame,1);

        Preprocess_Frame(current_frame);
        Advance_To_Target_Time(example.Time_At_Frame(current_frame));
        Postprocess_Frame(current_frame);

        example.frame_title=STRING_UTILITIES::string_sprintf("Frame %d",current_frame);
        Write_Output_Files(++output_number);
        LOG::cout<<"TIME = "<<time<<std::endl;}
}
//#####################################################################
template class COMPRESSIBLE_DRIVER<VECTOR<float,1> >;
template class COMPRESSIBLE_DRIVER<VECTOR<float,2> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMPRESSIBLE_DRIVER<VECTOR<double,1> >;
template class COMPRESSIBLE_DRIVER<VECTOR<double,2> >;
#endif
