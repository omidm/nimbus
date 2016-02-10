//#####################################################################
// Copyright 2012, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Dynamics/Coupled_Driver/COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Vectors/VECTOR_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT<TV>::
COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT(COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>& example_input)
    : BASE(example_input),example(example_input)
{}
//#####################################################################
// Compute_Dt
//#####################################################################
template<class TV> typename TV::SCALAR COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT<TV>::
Compute_Dt(const T time,const T target_time)
{
    T dt=target_time-time;
    example.Limit_Dt(dt,time);
    return dt;
}
//#####################################################################
// Initialize
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT<TV>::
Initialize()
{
    BASE::Initialize();
    example.Parse_Late_Options();
    example.Log_Parameters();

    example.Initialize_Boundaries();
    example.output_number=output_number;
    example.current_frame=current_frame;

    example.particle_levelset_evolution.Initialize_Domain(example.grid);
    example.particle_levelset_evolution.particle_levelset.Set_Band_Width(6);
    example.collision_bodies_affecting_fluid.Initialize_Grids();
    example.particle_levelset_evolution.Set_Time(time);
    example.particle_levelset_evolution.Set_CFL_Number(example.cfl);
    example.particle_levelset_evolution.particle_levelset.mpi_grid=example.mpi_grid;

    if(example.mpi_grid){example.mpi_grid->Initialize(example.domain_walls);
        example.compressible_boundary=new BOUNDARY_MPI<T_GRID,TV_DIMENSION>(example.mpi_grid,*example.compressible_boundary_scalar);
        example.boundary=new BOUNDARY_MPI<GRID<TV> >(example.mpi_grid,example.boundary_scalar);
        example.phi_boundary=new BOUNDARY_MPI<GRID<TV> >(example.mpi_grid,example.phi_boundary_water);
        example.particle_levelset_evolution.particle_levelset.last_unique_particle_id=example.mpi_grid->rank*30000000;}
    else{example.compressible_boundary=example.compressible_boundary_scalar;
        example.boundary=&example.boundary_scalar;
        example.phi_boundary=&example.phi_boundary_water;}
    example.phi_boundary_water.Set_Velocity_Pointer(example.face_velocities);
    example.projection.Initialize_Grid(example.grid);

    VECTOR<VECTOR<bool,2>,TV::dimension> domain_open_boundaries=VECTOR_UTILITIES::Complement(example.domain_walls);
    example.phi_boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    example.boundary->Set_Constant_Extrapolation(domain_open_boundaries);
    example.particle_levelset_evolution.Levelset_Advection(1).Set_Custom_Advection(example.advection_scalar);

    example.particle_levelset_evolution.Set_Number_Particles_Per_Cell(example.pls_particles_per_cell);
    example.particle_levelset_evolution.Set_Levelset_Callbacks(example);
    example.particle_levelset_evolution.Initialize_FMM_Initialization_Iterative_Solver(true);

    example.particle_levelset_evolution.particle_levelset.levelset.Set_Custom_Boundary(*example.phi_boundary);
    example.particle_levelset_evolution.Bias_Towards_Negative_Particles(false);
    example.particle_levelset_evolution.particle_levelset.Use_Removed_Positive_Particles();
    example.particle_levelset_evolution.particle_levelset.Use_Removed_Negative_Particles();
    example.particle_levelset_evolution.particle_levelset.Store_Unique_Particle_Id();
    example.particle_levelset_evolution.Use_Particle_Levelset(true);
    example.particle_levelset_evolution.particle_levelset.levelset.Set_Collision_Body_List(example.collision_bodies_affecting_fluid);
    example.particle_levelset_evolution.particle_levelset.levelset.Set_Face_Velocities_Valid_Mask(&example.valid_mask);
    example.particle_levelset_evolution.particle_levelset.Set_Collision_Distance_Factors((T).1,1);

    if(!example.restart){example.Initialize_Fluid_State();
        example.collision_bodies_affecting_fluid.Update_Intersection_Acceleration_Structures(false);
        example.collision_bodies_affecting_fluid.Rasterize_Objects();
        example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*example.grid.Minimum_Edge_Length(),5);
        example.particle_levelset_evolution.Make_Signed_Distance((T)example.number_of_cells_to_extrapolate);
        example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time);}
    else{example.Read_Output_Files(example.restart_frame);
        example.collision_bodies_affecting_fluid.Rasterize_Objects();
        example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(false,(T)2*example.grid.Minimum_Edge_Length(),5);} // compute grid visibility (for advection later)

    example.collision_bodies_affecting_fluid.Compute_Grid_Visibility();
    example.particle_levelset_evolution.Set_Seed(2606);
    if(!example.restart) example.particle_levelset_evolution.Seed_Particles(time);
    example.particle_levelset_evolution.Delete_Particles_Outside_Grid();

    ARRAY<T,TV_INT> exchanged_phi_ghost(example.grid.Domain_Indices(8));
    example.particle_levelset_evolution.particle_levelset.levelset.boundary->Fill_Ghost_Cells(example.grid,example.particle_levelset_evolution.phi,exchanged_phi_ghost,0,time,8);
    example.Extrapolate_Velocity_Across_Interface(example.face_velocities,exchanged_phi_ghost,(T)example.number_of_cells_to_extrapolate);
    example.Update_Psi_And_Extrapolate_State_Into_Incompressible_Flow(0,time);
}
//#####################################################################
// Run
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT<TV>::
Run(RANGE<TV_INT>& domain,const T dt,const T time)  // TODO(aanjneya): Modify when using gravity/buoyancy
{
    T_FACE_ARRAYS_SCALAR face_velocities_ghost;face_velocities_ghost.Resize(example.grid,3,false);
    example.boundary->Fill_Ghost_Cells_Face(example.grid,example.face_velocities,face_velocities_ghost,time+dt,3);
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,TV> interpolation;
    PARTICLE_LEVELSET_UNIFORM<GRID<TV> >& pls=example.particle_levelset_evolution.particle_levelset;
    if(pls.use_removed_positive_particles) for(typename GRID<TV>::NODE_ITERATOR iterator(example.grid,domain);iterator.Valid();iterator.Next()) if(pls.removed_positive_particles(iterator.Node_Index())){
        PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>& particles=*pls.removed_positive_particles(iterator.Node_Index());
        for(int p=1;p<=particles.array_collection->Size();p++){
            TV X=particles.X(p),V=interpolation.Clamped_To_Array_Face(example.grid,face_velocities_ghost,X);
            particles.V(p)=V;}}
}
//#####################################################################
// Advance_Levelset
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT<TV>::
Advance_Levelset(const T dt,const T time)
{
    LOG::Time("Compute Occupied Blocks");
    T maximum_fluid_speed=example.face_velocities.Maxabs().Max();
    T max_particle_collision_distance=example.particle_levelset_evolution.particle_levelset.max_collision_distance_factor*example.grid.dX.Max();
    example.collision_bodies_affecting_fluid.Compute_Occupied_Blocks(true,dt*maximum_fluid_speed+2*max_particle_collision_distance+(T).5*example.grid.dX.Max(),10);

    LOG::Time("Advect Phi");
    example.phi_boundary_water.Use_Extrapolation_Mode(false);
    example.particle_levelset_evolution.Advance_Levelset(dt);
    example.phi_boundary_water.Use_Extrapolation_Mode(true);

    LOG::Time("Step Particles");
    T_FACE_ARRAYS_SCALAR face_velocities_ghost;face_velocities_ghost.Resize(example.grid,3,false);
    example.boundary->Fill_Ghost_Cells_Face(example.grid,example.face_velocities,face_velocities_ghost,time+dt,3);
    example.particle_levelset_evolution.particle_levelset.Euler_Step_Particles(face_velocities_ghost,dt,time,true,true,false,false);

    LOG::Time("Advect Removed Particles");
    example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time);
    RANGE<TV_INT> domain(example.grid.Domain_Indices());domain.max_corner+=TV_INT::All_Ones_Vector();
    Run(domain,dt,time);

    LOG::Time("Modify Levelset");
    example.particle_levelset_evolution.particle_levelset.Exchange_Overlap_Particles();
    example.particle_levelset_evolution.Modify_Levelset_And_Particles(&face_velocities_ghost,(T)example.number_of_cells_to_extrapolate);

    LOG::Time("Delete Particles");
    example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time+dt);
    example.particle_levelset_evolution.Delete_Particles_Outside_Grid();
    example.particle_levelset_evolution.particle_levelset.Delete_Particles_In_Local_Maximum_Phi_Cells(1);
    example.particle_levelset_evolution.particle_levelset.Delete_Particles_Far_From_Interface(); // uses visibility
    example.particle_levelset_evolution.particle_levelset.Identify_And_Remove_Escaped_Particles(face_velocities_ghost,1.5,time+dt);
        
    LOG::Time("Reincorporate Particles");
    if(example.particle_levelset_evolution.particle_levelset.use_removed_positive_particles || example.particle_levelset_evolution.particle_levelset.use_removed_negative_particles)
        example.particle_levelset_evolution.particle_levelset.Reincorporate_Removed_Particles(1,1,0,true);
    example.particle_levelset_evolution.Fill_Levelset_Ghost_Cells(time+dt);
}
//#####################################################################
// Advance_One_Time_Step_Explicit_Part
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT<TV>::
Advance_One_Time_Step_Explicit_Part(const T dt,const T time)
{
    LOG::Time("Advect V");
    T_FACE_ARRAYS_SCALAR face_velocities_ghost;face_velocities_ghost.Resize(example.grid,3,false);
    example.boundary->Fill_Ghost_Cells_Face(example.grid,example.face_velocities,face_velocities_ghost,time+dt,3);
    example.eno_advection.Update_Advection_Equation_Face_Lookup(example.grid,example.face_velocities,T_FACE_LOOKUP(face_velocities_ghost),T_FACE_LOOKUP(face_velocities_ghost),*example.boundary,dt,time,0,0,0,0);

    LOG::Time("Advect U");
    example.compressible_boundary->Fill_Ghost_Cells(example.grid,example.U,example.U_ghost,dt,time,3);
    example.conservation_law_solver->Update_Conservation_Law(example.grid,example.U,example.U_ghost,example.psi,dt,example.cfl_eigensystem,example.cfl_eigensystem,example.psi_N,example.compressible_face_velocities);
    example.compressible_boundary->Apply_Boundary_Condition(example.grid,example.U,time+dt);
}
//#####################################################################
// Advance_One_Time_Step_Implicit_Part
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT<TV>::
Advance_One_Time_Step_Implicit_Part(const T dt,const T time)
{
    LOG::Time("Project");
    example.Project(dt,time);

    example.Update_Psi_And_Extrapolate_State_Into_Incompressible_Flow(dt,time);
    ARRAY<T,TV_INT> exchanged_phi_ghost(example.grid.Domain_Indices(8));
    example.particle_levelset_evolution.particle_levelset.levelset.boundary->Fill_Ghost_Cells(example.grid,example.particle_levelset_evolution.phi,exchanged_phi_ghost,0,time,8);
    example.Extrapolate_Velocity_Across_Interface(example.face_velocities,exchanged_phi_ghost,(T)example.number_of_cells_to_extrapolate);
}
//#####################################################################
// Advance_To_Target_Time
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT<TV>::
Advance_To_Target_Time(const T target_time)
{
    bool done=false;
    for(int substep=1;!done;substep++){
        LOG::SCOPE scope("SUBSTEP","substep %d",substep,1);
        T dt=Compute_Dt(time,target_time);dt=min(dt,example.particle_levelset_evolution.CFL(false,false));
        if(example.mpi_grid) example.mpi_grid->Synchronize_Dt(dt);
        EXAMPLE<TV>::Clamp_Time_Step_With_Target_Time(time,target_time,dt,done);

        Advance_Levelset(dt,time);

        RUNGEKUTTA<T_ARRAYS_DIMENSION_SCALAR> rungekutta_U(example.U);RUNGEKUTTA<T_FACE_ARRAYS_SCALAR> rungekutta_velocities(example.face_velocities);
        RUNGEKUTTA<T_ARRAYS_SCALAR> rungekutta_p(example.pressure);

        rungekutta_U.Set_Grid_And_Boundary_Condition(example.grid,*example.compressible_boundary);
        rungekutta_U.Set_Order(example.rk_order);rungekutta_U.Set_Time(time);rungekutta_U.Start(dt);T rk_time=time;
        rungekutta_velocities.Set_Grid_And_Boundary_Condition(example.grid,*example.boundary);
        rungekutta_velocities.Set_Order(example.rk_order);rungekutta_velocities.Set_Time(time);rungekutta_velocities.Start(dt);
        rungekutta_p.Set_Grid_And_Boundary_Condition(example.grid,*example.boundary);
        rungekutta_p.Set_Order(example.rk_order);rungekutta_p.Set_Time(time);rungekutta_p.Start(dt);

        for(int rk_substep=1;rk_substep<=example.rk_order;++rk_substep){
            Advance_One_Time_Step_Explicit_Part(dt,rk_time);
            Advance_One_Time_Step_Implicit_Part(dt,rk_time);
            rk_time=rungekutta_U.Main();rungekutta_velocities.Main();rungekutta_p.Main();}

        if(!done) example.Write_Substep("END Substep",substep,0);
        time+=dt;}
}
//#####################################################################
// Simulate_To_Frame
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT<TV>::
Simulate_To_Frame(const int target_frame)
{
    example.frame_title=STRING_UTILITIES::string_sprintf("Frame %d",example.current_frame);
    if(!example.restart) Write_Output_Files(example.current_frame);

    while(example.current_frame<target_frame){
        LOG::SCOPE scope("FRAME","Frame %d",++example.current_frame,1);

        Advance_To_Target_Time(example.Time_At_Frame(example.current_frame));
        //LOG::Time("Reseed");
        //if((current_frame-example.first_frame)%10==0){
        //    example.particle_levelset_evolution.Reseed_Particles(time);
        //    example.particle_levelset_evolution.Delete_Particles_Outside_Grid();}

        example.frame_title=STRING_UTILITIES::string_sprintf("Frame %d",example.current_frame);
        Write_Output_Files(++example.output_number);
        LOG::cout<<"TIME = "<<time<<std::endl;}
}
//#####################################################################
template class COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT<VECTOR<float,1> >;
template class COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT<VECTOR<float,2> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT<VECTOR<double,1> >;
template class COMPRESSIBLE_INCOMPRESSIBLE_DRIVER_EXPLICIT<VECTOR<double,2> >;
#endif
