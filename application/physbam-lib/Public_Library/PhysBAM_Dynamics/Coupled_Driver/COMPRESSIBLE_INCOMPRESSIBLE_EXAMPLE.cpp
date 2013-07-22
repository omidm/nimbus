//#####################################################################
// Copyright 2012, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Dynamics/Coupled_Driver/COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Equations_Of_State/EOS_GAMMA.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UTILITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_1D_EIGENSYSTEM_F_ADVECTION_ONLY.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_2D_EIGENSYSTEM_F_ADVECTION_ONLY.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_2D_EIGENSYSTEM_G_ADVECTION_ONLY.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_STRUCTURE.h>
#include <PhysBAM_Dynamics/Particles/PARTICLES_FORWARD.h>
#include <PhysBAM_Dynamics/Read_Write/Particles/READ_WRITE_PARTICLES.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_FAST_LEVELSET.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_1D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Tools/Read_Write/Point_Clouds/READ_WRITE_POINT_CLOUD.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE(const STREAM_TYPE& stream_type): 
    BASE(stream_type),mpi_grid(0),rho_incompressible(1e3),viscosity(0),surface_tension(0),scale(0),number_of_cells_to_extrapolate(7),apply_viscosity(false),
    eos(new EOS_GAMMA<T>()),boundary(0),phi_boundary(0),compressible_boundary(0),compressible_boundary_scalar(0),conservation_law_solver(0),particle_levelset_evolution(grid,3),
    projection(grid,pressure,one_over_rho_c_squared,&particle_levelset_evolution.phi),collision_bodies_affecting_fluid(grid)
{
    LOG::Initialize_Logging(false,false,1<<30,true,1);
    Initialize_Particles();
    Initialize_Read_Write_Structures();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
~COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE()
{
    if(mpi_grid) delete compressible_boundary;
    for(int i=1;i<=TV::dimension;i++) delete cfl_eigensystem[i];
    delete conservation_law_solver;
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Write_Output_Files(const int frame) const
{
    std::string f=STRING_UTILITIES::string_sprintf("/%d",frame);
    FILE_UTILITIES::Create_Directory(output_directory+STRING_UTILITIES::string_sprintf("/%d",frame));
    FILE_UTILITIES::Write_To_Text_File(output_directory+STRING_UTILITIES::string_sprintf("/%d/frame_title",frame),frame_title);
    if(write_last_frame) FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/last_frame",frame,"\n");

    // Restart data
    if(frame==first_frame){FILE_UTILITIES::Create_Directory(output_directory);
        FILE_UTILITIES::Create_Directory(output_directory+"/common");
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",grid);
        FILE_UTILITIES::Write_To_Text_File(output_directory+"/common/first_frame",frame,"\n");
        if(mpi_grid) FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/global_grid",mpi_grid->global_grid);}
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+f+"/euler_U",U);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+f+"/euler_psi",psi);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure",pressure);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_N",projection.psi_N);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/psi_D",projection.psi_D);
    const T_PARTICLE_LEVELSET& particle_levelset=particle_levelset_evolution.particle_levelset;
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/levelset",particle_levelset.levelset);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"positive_particles"),particle_levelset.positive_particles);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"negative_particles"),particle_levelset.negative_particles);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_positive_particles"),particle_levelset.removed_positive_particles);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_negative_particles"),particle_levelset.removed_negative_particles);
    FILE_UTILITIES::Write_To_Text_File(output_directory+"/"+f+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);

    // Detailed output data
    T_ARRAYS_SCALAR rho(grid.Domain_Indices()),energy(grid.Domain_Indices()),internal_energy(grid.Domain_Indices());
    T_ARRAYS_VECTOR velocity(grid.Domain_Indices());
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(particle_levelset_evolution.phi(cell)>0){rho(cell)=U(cell)(1);
            internal_energy(cell)=EULER<T_GRID>::e(U(cell));
            energy(cell)=EULER<T_GRID>::Get_Total_Energy(U,cell);
            for(int axis=1;axis<=TV::dimension;axis++) velocity(cell)(axis)=EULER<T_GRID>::Get_Velocity_Component(U(cell),axis);}
        else{rho(cell)=(T)rho_incompressible;
            for(int axis=1;axis<=TV::dimension;axis++){FACE_INDEX<TV::dimension> first_face=iterator.Full_First_Face_Index(axis),second_face=iterator.Full_Second_Face_Index(axis);
                velocity(cell)(axis)=(T).5*(face_velocities(first_face)+face_velocities(second_face));}}}

    FILE_UTILITIES::Write_To_File(stream_type,output_directory+f+"/density",rho);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+f+"/centered_velocities",velocity);
    if(write_debug_data){
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/compressible_implicit_pressure",compressible_pressure);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+f+"/energy",energy);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+f+"/internal_energy",internal_energy);}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Write_Substep(const std::string& title,const int substep,const int level)
{
    if(level<=write_substeps_level){
        frame_title=title;
        std::stringstream ss;
        ss<<"Writing substep ["<<frame_title<<"]: output_number="<<output_number+1<<", time="<<time<<", frame="<<current_frame<<", substep="<<substep<<std::endl;
        LOG::filecout(ss.str());
        Write_Output_Files(++output_number);frame_title="";}
}
//#####################################################################
// Initialize_Grids
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Initialize_Grids()
{
    U.Resize(grid.Domain_Indices());
    U_ghost.Resize(grid.Domain_Indices(3));
    face_velocities.Resize(grid);
    compressible_face_velocities.Resize(grid.Domain_Indices(3));
    entropy.Resize(grid.Domain_Indices());
    pressure.Resize(grid.Domain_Indices(1));
    psi.Resize(grid.Domain_Indices());
    one_over_rho_c_squared.Resize(grid.Domain_Indices());
    compressible_pressure.Resize(grid.Domain_Indices(1));
    psi_N.Resize(grid.Domain_Indices(3),true,false);
    valid_mask.Resize(grid.Domain_Indices(3),true,true,true);
}
//#####################################################################
// Read_Output_Files
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Read_Output_Files(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("/%d",frame);

    // Restart data
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/common/grid",grid);
    if(mpi_grid) FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/common/global_grid",mpi_grid->global_grid);
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+f+"/euler_U",U);
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+f+"/euler_psi",psi);
    
    T_PARTICLE_LEVELSET& particle_levelset=particle_levelset_evolution.particle_levelset;
    FILE_UTILITIES::Read_From_File(stream_type,output_directory+"/"+f+"/levelset",particle_levelset.levelset);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"positive_particles"),particle_levelset.positive_particles);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"negative_particles"),particle_levelset.negative_particles);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_positive_particles"),particle_levelset.removed_positive_particles);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/%s",output_directory.c_str(),frame,"removed_negative_particles"),particle_levelset.removed_negative_particles);
    FILE_UTILITIES::Read_From_Text_File(output_directory+"/"+f+"/last_unique_particle_id",particle_levelset.last_unique_particle_id);
    std::string filename;
    filename=output_directory+"/"+f+"/pressure";
    if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading pressure "<<filename<<std::endl;LOG::filecout(ss.str());FILE_UTILITIES::Read_From_File(stream_type,filename,pressure);}
    filename=output_directory+"/"+f+"/mac_velocities";
    if(FILE_UTILITIES::File_Exists(filename)){std::stringstream ss;ss<<"Reading mac_velocities "<<filename<<std::endl;LOG::filecout(ss.str());FILE_UTILITIES::Read_From_File(stream_type,filename,face_velocities);}
}
//#####################################################################
// Compute_Divergence
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Compute_Divergence(const T_FACE_LOOKUP& face_lookup)
{
    TV one_over_dx=grid.one_over_dX;
    // Compute divergence at cell centers
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        const typename T_FACE_LOOKUP::LOOKUP& lookup=face_lookup.Starting_Point_Cell(cell);
        T divergence=0;for(int axis=1;axis<=T_GRID::dimension;axis++){
            divergence+=(lookup(axis,iterator.Second_Face_Index(axis))-lookup(axis,iterator.First_Face_Index(axis)))*one_over_dx[axis];}
        projection.f(cell)=divergence;}
}
//#####################################################################
// Compute_One_Over_Rho_C_Squared
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Compute_One_Over_Rho_C_Squared()
{
    one_over_rho_c_squared.Fill((T)0.);
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(particle_levelset_evolution.phi(cell)>0){T sound_speed=eos->c(U(cell)(1),EULER<T_GRID>::e(U(cell)));
            one_over_rho_c_squared(cell)=(T)1./(U(cell)(1)*sound_speed*sound_speed);}}
}
//#####################################################################
// Compute_Right_Hand_Side
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Compute_Right_Hand_Side(const T dt)
{
    // compute combined velocities.
    T_FACE_ARRAYS_SCALAR combined_velocities(face_velocities);
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        TV_INT first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
        if(particle_levelset_evolution.phi(first_cell)>0 && particle_levelset_evolution.phi(second_cell)>0)
            combined_velocities(face)=compressible_face_velocities(face);}

    T one_over_dt=(T)1./dt;
    Compute_Divergence(T_FACE_LOOKUP(combined_velocities));
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(particle_levelset_evolution.phi(cell)>0){T p_star=EULER<T_GRID>::p(eos,U(cell));
            projection.f(cell)-=p_star*one_over_rho_c_squared(cell)*one_over_dt;}}
}
//#####################################################################
// Fill_Face_Weights_For_Projection
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Fill_Face_Weights_For_Projection(const T dt,const T time)
{
    compressible_boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,3);
    projection.Set_Variable_beta(true);
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> full_index=iterator.Full_Index();
        TV_INT first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
        if(particle_levelset_evolution.phi(first_cell)<=0 && particle_levelset_evolution.phi(second_cell)<=0) projection.beta_face(full_index)=(T)1./rho_incompressible;
        else if(particle_levelset_evolution.phi(first_cell)>0 && particle_levelset_evolution.phi(second_cell)>0){
            T rho_first_cell=U_ghost(first_cell)(1),rho_second_cell=U_ghost(second_cell)(1);
            T rho_face=(T).5*(rho_first_cell+rho_second_cell);
            projection.beta_face(full_index)=(T)1./rho_face;}
        else{T theta=abs(particle_levelset_evolution.phi(first_cell))/(abs(particle_levelset_evolution.phi(first_cell))+abs(particle_levelset_evolution.phi(second_cell)));
            if(particle_levelset_evolution.phi(first_cell)<=0 && particle_levelset_evolution.phi(second_cell)>0){
                T beta_minus=(T)1./rho_incompressible,beta_plus=(T)1/U_ghost(second_cell)(1);
                projection.beta_face(full_index)=(beta_minus*beta_plus)/((1-theta)*beta_minus+theta*beta_plus);}
            else{T beta_minus=(T)1/U_ghost(first_cell)(1),beta_plus=(T)1./rho_incompressible;
                projection.beta_face(full_index)=(beta_minus*beta_plus)/((1-theta)*beta_minus+theta*beta_plus);}}}
}
//#####################################################################
// Compute_Pressure
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Compute_Pressure(const T dt,const T time)
{
    pressure*=dt;
    Compute_Right_Hand_Side(dt);
    Fill_Face_Weights_For_Projection(dt,time);
    projection.Set_Dt(dt);
    projection.Find_Solution_Regions();     // flood fill
    projection.Solve(time,true);            // solve all regions
    boundary->Apply_Boundary_Condition(grid,pressure,time);
    pressure*=((T)1./dt);
}
//#####################################################################
// Project
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Project(const T dt,const T time)
{
    // compute viscosity
    if(apply_viscosity){T_ARRAYS_SCALAR divergence_velocities(grid.Domain_Indices());
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();T divergence=(T)0.;
            for(int axis=1;axis<=TV::dimension;axis++){FACE_INDEX<TV::dimension> first_face=iterator.Full_First_Face_Index(axis),second_face=iterator.Full_Second_Face_Index(axis);
                divergence+=(face_velocities(second_face)-face_velocities(first_face));}
            divergence_velocities(cell)=divergence;}
        T_ARRAYS_SCALAR divergence_velocities_ghost(grid.Domain_Indices(1));T_FACE_ARRAYS_SCALAR grad_divergence_velocities(grid);
        boundary->Fill_Ghost_Cells(grid,divergence_velocities,divergence_velocities_ghost,dt,time,1);
        ARRAYS_UTILITIES<T_GRID,T>::Compute_Gradient_At_Faces_From_Cell_Data(grid,grad_divergence_velocities,divergence_velocities_ghost);

        // apply viscosity
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
            face_velocities(face)-=dt*viscosity*grad_divergence_velocities(face);}}

    Compute_Density_Weighted_Face_Velocities(dt,time);
    Set_Boundary_Conditions(dt,time);
    Compute_Pressure(dt,time);

    T_ARRAYS_SCALAR p_ghost(grid.Domain_Indices(1));
    boundary->Fill_Ghost_Cells(grid,pressure,p_ghost,dt,time,1);

    T_FACE_ARRAYS_SCALAR p_face(grid);
    Compute_Face_Pressure_From_Cell_Pressures(grid,p_face,p_ghost,dt,time);
    Apply_Pressure(p_ghost,p_face,dt,time);
}
//#####################################################################
// Apply_Pressure
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Apply_Pressure(T_ARRAYS_SCALAR& p_ghost,T_FACE_ARRAYS_SCALAR& p_face,const T dt,const T time)
{
    TV one_over_dx=grid.one_over_dX;
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> full_index=iterator.Full_Index();
        TV_INT first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();int axis=iterator.Axis();
        T phi_face=(T).5*(particle_levelset_evolution.phi(first_cell)+particle_levelset_evolution.phi(second_cell));
        if(phi_face<=0) face_velocities(full_index)-=dt*projection.beta_face(full_index)*(p_ghost(second_cell)-p_ghost(first_cell))*one_over_dx[axis];}

    T_ARRAYS_SCALAR p_hat(p_ghost);T_FACE_ARRAYS_SCALAR p_hat_face(p_face);
    p_hat*=dt;p_hat_face*=dt;
    // update momentum
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(psi(cell)){for(int axis=1;axis<=T_GRID::dimension;axis++){
            TV_INT first_face_index=iterator.First_Face_Index(axis),second_face_index=iterator.Second_Face_Index(axis);
            T p_hat_first_face=p_hat_face.Component(axis)(first_face_index);
            T p_hat_second_face=p_hat_face.Component(axis)(second_face_index);
            T grad_p_hat=(p_hat_second_face-p_hat_first_face)*one_over_dx[axis];
            U(cell)(axis+1)-=grad_p_hat;}}}
    T_FACE_ARRAYS_SCALAR face_velocities_np1(compressible_face_velocities);

    // update face velocities for energy update
    T_FACE_ARRAYS_SCALAR grad_p_hat_face(grid);
    ARRAYS_UTILITIES<T_GRID,T>::Compute_Gradient_At_Faces_From_Cell_Data(grid,grad_p_hat_face,p_hat);
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
        TV_INT first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
        if((!psi.Valid_Index(first_cell) || psi(first_cell)) && (!psi.Valid_Index(second_cell) || psi(second_cell))){
            T rho_face=(U_ghost(first_cell)(1)+U_ghost(second_cell)(1))*(T).5;
            face_velocities_np1.Component(axis)(face_index)-=grad_p_hat_face.Component(axis)(face_index)/rho_face;}}

    // extrapolate incompressible velocities
    ARRAY<T,TV_INT> exchanged_phi_ghost(grid.Domain_Indices(8));T_FACE_ARRAYS_SCALAR incompressible_velocities(face_velocities);
    particle_levelset_evolution.particle_levelset.levelset.boundary->Fill_Ghost_Cells(grid,particle_levelset_evolution.phi,exchanged_phi_ghost,dt,time,8);
    Extrapolate_Velocity_Across_Interface(face_velocities,exchanged_phi_ghost,(T)number_of_cells_to_extrapolate);

    T_ARRAYS_SCALAR& phi=particle_levelset_evolution.phi;
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        TV_INT first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();
        if(phi(first_cell)*phi(second_cell)<0) face_velocities_np1(face)=incompressible_velocities(face);}

    // update energy
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(psi(cell)){T div_p_hat_u=0;
            for(int axis=1;axis<=T_GRID::dimension;axis++){
                TV_INT first_face_index=iterator.First_Face_Index(axis),second_face_index=iterator.Second_Face_Index(axis);
                T p_hat_first_face=p_hat_face.Component(axis)(first_face_index);
                T p_hat_second_face=p_hat_face.Component(axis)(second_face_index);
                div_p_hat_u+=(p_hat_second_face*face_velocities_np1.Component(axis)(second_face_index)-p_hat_first_face*face_velocities_np1.Component(axis)(first_face_index))*one_over_dx[axis];}
            U(cell)(T_GRID::dimension+2)-=div_p_hat_u;}}
    compressible_boundary->Apply_Boundary_Condition(grid,U,time+dt);
}
//#####################################################################
// Adjust_Particle_For_Domain_Boundaries
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Adjust_Particle_For_Domain_Boundaries(PARTICLE_LEVELSET_PARTICLES<TV>& particles,const int index,TV& V,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time)
{
    if(particle_type==PARTICLE_LEVELSET_POSITIVE || particle_type==PARTICLE_LEVELSET_REMOVED_POSITIVE) return;

    TV& X=particles.X(index);TV X_new=X+dt*V;
    T max_collision_distance=particle_levelset_evolution.particle_levelset.Particle_Collision_Distance(particles.quantized_collision_distance(index));
    T min_collision_distance=particle_levelset_evolution.particle_levelset.min_collision_distance_factor*max_collision_distance;
    TV min_corner=grid.domain.Minimum_Corner(),max_corner=grid.domain.Maximum_Corner();
    for(int axis=1;axis<=TV::dimension;axis++){
        if(domain_walls[axis][1] && X_new[axis]<min_corner[axis]+max_collision_distance){
            T collision_distance=X[axis]-min_corner[axis];
            if(collision_distance>max_collision_distance)collision_distance=X_new[axis]-min_corner[axis];
            collision_distance=max(min_collision_distance,collision_distance);
            X_new[axis]+=max((T)0,min_corner[axis]-X_new[axis]+collision_distance);
            V[axis]=max((T)0,V[axis]);X=X_new-dt*V;}
        if(domain_walls[axis][2] && X_new[axis]>max_corner[axis]-max_collision_distance){
            T collision_distance=max_corner[axis]-X[axis];
            if(collision_distance>max_collision_distance) collision_distance=max_corner[axis]-X_new[axis];
            collision_distance=max(min_collision_distance,collision_distance);
            X_new[axis]-=max((T)0,X_new[axis]-max_corner[axis]+collision_distance);
            V[axis]=min((T)0,V[axis]);X=X_new-dt*V;}}
}
//#####################################################################
// Compute_Density_Weighted_Face_Velocities
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Compute_Density_Weighted_Face_Velocities(const T dt,const T time)
{
    compressible_boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,3);
    compressible_face_velocities.Fill(TV());
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT face=iterator.Face_Index();
        TV_INT first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();int axis=iterator.Axis();
        if((!psi.Valid_Index(first_cell) || psi(first_cell)) && (!psi.Valid_Index(second_cell) || psi(second_cell))){
            T rho_first_cell=U_ghost(first_cell)(1),rho_second_cell=U_ghost(second_cell)(1);
            compressible_face_velocities.Component(axis)(face)=(rho_first_cell*EULER<T_GRID>::Get_Velocity_Component(U_ghost,first_cell,axis)+
                    rho_second_cell*EULER<T_GRID>::Get_Velocity_Component(U_ghost,second_cell,axis))/(rho_first_cell+rho_second_cell);}
        else if(!psi.Valid_Index(first_cell) || psi(first_cell))
            compressible_face_velocities.Component(axis)(face)=EULER<T_GRID>::Get_Velocity_Component(U_ghost,first_cell,axis);
        else if(!psi.Valid_Index(second_cell) || psi(second_cell))
            compressible_face_velocities.Component(axis)(face)=EULER<T_GRID>::Get_Velocity_Component(U_ghost,second_cell,axis);}
}
//#####################################################################
// Compute_Face_Pressure_From_Cell_Pressures
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Compute_Face_Pressure_From_Cell_Pressures(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& p_face,const T_ARRAYS_SCALAR& p_cell,const T dt,const T time)
{
    T_ARRAYS_SCALAR& phi=particle_levelset_evolution.phi;
    compressible_boundary->Fill_Ghost_Cells(grid,U,U_ghost,dt,time,3);
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        TV_INT first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();int axis=iterator.Axis();
        if(phi(first_cell)>0 && phi(second_cell)>0){
            T rho_first_cell=U_ghost(first_cell)(1),rho_second_cell=U_ghost(second_cell)(1);
            p_face.Component(axis)(iterator.Face_Index())=(rho_second_cell*p_cell(first_cell)+rho_first_cell*p_cell(second_cell))/(rho_first_cell+rho_second_cell);}
        else if(phi(first_cell)<=0 && phi(second_cell)<=0) p_face(iterator.Full_Index())=(T).5*(p_cell(first_cell)+p_cell(second_cell));
        else{T theta=abs(particle_levelset_evolution.phi(first_cell))/(abs(particle_levelset_evolution.phi(first_cell))+abs(particle_levelset_evolution.phi(second_cell)));
            T rho_first_cell=rho_incompressible,rho_second_cell=rho_incompressible;
            if(phi(first_cell)>0) rho_first_cell=U_ghost(first_cell)(1);
            else rho_second_cell=U_ghost(second_cell)(1);
            T beta_minus=(T)1./rho_first_cell,beta_plus=(T)1./rho_second_cell;
            T p_interface=(theta*beta_plus*p_cell(second_cell)+((T)1.-theta)*beta_minus*p_cell(first_cell))/(theta*beta_plus+((T)1.-theta)*beta_minus);
           
            T tolerance=(T)1e-5;
            if(phi(first_cell)>0){
                if(theta<tolerance){TV_INT left_neighbor=first_cell-TV_INT::Axis_Vector(axis);
                    p_face(face)=p_cell(first_cell)+(T).5*(p_cell(first_cell)-p_cell(left_neighbor));}
                else p_face(face)=p_cell(first_cell)+((T).5/theta)*(p_interface-p_cell(first_cell));}
            else{
                if((T)1.-theta<tolerance){TV_INT right_neighbor=second_cell+TV_INT::Axis_Vector(axis);
                    p_face(face)=p_cell(second_cell)-(T).5*(p_cell(right_neighbor)-p_cell(second_cell));}
                else p_face(face)=p_cell(second_cell)-((T).5/((T)1.-theta))*(p_cell(second_cell)-p_interface);}}}
}
//#####################################################################
// Extrapolate_Velocity_Across_Interface
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Extrapolate_Velocity_Across_Interface(T_FACE_ARRAYS_SCALAR& face_velocities_input,T_ARRAYS_SCALAR& phi_ghost,const T band_width)
{
    // TODO(aanjneya): Fix for MPI.
    T delta=band_width*grid.dX.Max();const int ghost_cells=2*(int)ceil(band_width)+1;
    for(int axis=1;axis<=T_GRID::dimension;axis++){
        T_GRID face_grid=grid.Get_Face_Grid(axis);T_ARRAYS_SCALAR phi_face(face_grid.Domain_Indices(),false);T_ARRAYS_BASE& face_velocity=face_velocities.Component(axis);
        T_ARRAYS_BOOL fixed_face=T_ARRAYS_BOOL(face_grid.Domain_Indices());
        for(FACE_ITERATOR iterator(grid,0,T_GRID::WHOLE_REGION,0,axis);iterator.Valid();iterator.Next()){
            TV_INT index=iterator.Face_Index();phi_face(index)=(T).5*(phi_ghost(iterator.First_Cell_Index())+phi_ghost(iterator.Second_Cell_Index()));
            if(phi_face(index)<=0) fixed_face(index)=true;if(phi_face(index) >= delta && !fixed_face(index)) face_velocity(index)=(T)0;}
        T_EXTRAPOLATION_SCALAR extrapolate(face_grid,phi_face,face_velocity,ghost_cells);extrapolate.Set_Band_Width(band_width);extrapolate.Set_Custom_Seed_Done(&fixed_face);
        extrapolate.Extrapolate();}
}
//#####################################################################
// Extrapolate_State_One_Sided
//#####################################################################
template<class TV> template<class T_TYPE> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Extrapolate_State_One_Sided(T_TYPE& Z,T_ARRAYS_SCALAR& phi)
{
    typedef EXTRAPOLATION_UNIFORM<T_GRID,typename T_TYPE::ELEMENT> T2_EXTRAPOLATION_SCALAR;
    T band_width=(T)number_of_cells_to_extrapolate;const int ghost_cells=2*(int)ceil(band_width)+1;
    T_ARRAYS_SCALAR phi_negated(phi);phi_negated*=(T)-1;
    T2_EXTRAPOLATION_SCALAR extrapolate(grid,phi_negated,Z,ghost_cells);extrapolate.Set_Band_Width(band_width);
    extrapolate.Extrapolate();
}
//#####################################################################
// Extrapolate_State_Into_Incompressible_Flow
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Update_Psi_And_Extrapolate_State_Into_Incompressible_Flow(const T dt,const T time)
{
    T_ARRAYS_SCALAR exchanged_phi_ghost(grid.Domain_Indices(8));T_FACE_ARRAYS_SCALAR incompressible_velocities(face_velocities);
    boundary->Fill_Ghost_Cells(grid,particle_levelset_evolution.phi,exchanged_phi_ghost,0,time,8);
    Extrapolate_Velocity_Across_Interface(incompressible_velocities,exchanged_phi_ghost,(T)number_of_cells_to_extrapolate);

    // compute entropy, pressure and velocity from the compressible region
    T_ARRAYS_VECTOR velocity(grid.Domain_Indices());
    compressible_pressure.Fill((T)0.);entropy.Fill((T)0.);
    particle_levelset_evolution.particle_levelset.levelset.Compute_Normals(0);

    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
       if(psi(cell)){compressible_pressure(cell)=eos->p(U(cell)(1),EULER<T_GRID>::e(U(cell)));
           entropy(cell)=eos->S(U(cell)(1),EULER<T_GRID>::e(U(cell)));
           for(int axis=1;axis<=TV::dimension;axis++) velocity(cell)(axis)=U(cell)(axis+1)/U(cell)(1);}
       bool is_compressible=false;
       if(particle_levelset_evolution.phi(cell)>0) is_compressible=true;
       else{for(int neighbor=1;neighbor<=2*TV::dimension;neighbor++){TV_INT cell_neighbor=iterator.Cell_Neighbor(neighbor);
           if(particle_levelset_evolution.phi(cell_neighbor)>0){is_compressible=true;break;}}}
       psi(cell)=is_compressible;}

    // extrapolate entropy, pressure and velocity from the compressible region
    Extrapolate_State_One_Sided(compressible_pressure,particle_levelset_evolution.phi);
    Extrapolate_State_One_Sided(entropy,particle_levelset_evolution.phi);
    Extrapolate_State_One_Sided(velocity,particle_levelset_evolution.phi);
   
    // make normal component of velocity incompressible
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(particle_levelset_evolution.phi(cell)<=0 && particle_levelset_evolution.phi(cell)>(T)-number_of_cells_to_extrapolate*grid.min_dX){
            TV cell_normal=particle_levelset_evolution.particle_levelset.levelset.Normal(iterator.Location()).Normalized();TV incompressible_velocity;
            for(int axis=1;axis<=TV::dimension;axis++){
                TV_INT first_face=iterator.First_Face_Index(axis),second_face=iterator.Second_Face_Index(axis);
                incompressible_velocity(axis)=(T).5*(incompressible_velocities.Component(axis)(first_face)+incompressible_velocities.Component(axis)(second_face));}
            velocity(cell)=TV::Dot_Product(cell_normal,incompressible_velocity)*cell_normal+velocity(cell)-TV::Dot_Product(cell_normal,velocity(cell))*cell_normal;}}

    // reconstruct the extrapolated state vector
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(particle_levelset_evolution.phi(cell)>(T)-number_of_cells_to_extrapolate*grid.min_dX){
            T rho=eos->rho_From_p_And_S(compressible_pressure(cell),entropy(cell));
            T internal_energy=eos->e_From_p_And_S(compressible_pressure(cell),entropy(cell));
            U(cell)=EULER<T_GRID>::Get_Euler_State_From_rho_velocity_And_internal_energy(rho,velocity(cell),internal_energy);}}
}
//#####################################################################
// Limit_Dt
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Limit_Dt(T& dt,const T time)
{
    TV one_over_dx=grid.one_over_dX,max_lambdas,max_grad_p_over_rho,max_grad_p;
    Update_Psi_And_Extrapolate_State_Into_Incompressible_Flow(dt,time);

    T_FACE_ARRAYS_SCALAR p_approx_face(grid);
    Compute_Face_Pressure_From_Cell_Pressures(grid,p_approx_face,compressible_pressure,dt,time);
    T_ARRAYS_VECTOR grad_p_approx(grid.Domain_Indices());ARRAYS_UTILITIES<T_GRID,T>::Compute_Gradient_At_Cells_From_Face_Data(grid,grad_p_approx,p_approx_face);

    T dt_incompressible=0;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){const TV_INT cell=iterator.Cell_Index();TV X=iterator.Location();
        if(particle_levelset_evolution.phi(cell)>0){for(int axis=1;axis<=T_GRID::dimension;axis++){
            max_lambdas[axis]=max(max_lambdas[axis],cfl_eigensystem[axis]->Maximum_Magnitude_Eigenvalue(U(cell)));
            max_grad_p[axis]=maxabs(max_grad_p[axis],grad_p_approx(cell)[axis]);
            max_grad_p_over_rho[axis]=maxabs(max_grad_p_over_rho[axis],grad_p_approx(cell)[axis]/U(cell)(1));}}
        else{T local_V_norm=0,local_viscosity_norm=(T)0.;
            for(int axis=1;axis<=T_GRID::dimension;axis++){ 
                local_V_norm+=grid.one_over_dX[axis]*maxabs(face_velocities(axis,grid.First_Face_Index_In_Cell(axis,cell)),face_velocities(axis,grid.Second_Face_Index_In_Cell(axis,cell)));
                local_viscosity_norm+=grid.one_over_dX[axis]*grid.one_over_dX[axis];}
            local_viscosity_norm*=((T)2.*viscosity/rho_incompressible);
            T local_surface_tension_norm_squared=-surface_tension*particle_levelset_evolution.particle_levelset.levelset.Compute_Curvature(X)/(rho_incompressible*grid.min_dX*grid.min_dX);
            T local_dt=(T).5*(local_V_norm+local_viscosity_norm+sqrt((local_V_norm+local_viscosity_norm)*(local_V_norm+local_viscosity_norm)+(T)4.*local_surface_tension_norm_squared));
            dt_incompressible=max(dt_incompressible,local_dt);}}
    if(dt_incompressible!=0) dt_incompressible=cfl/dt_incompressible;

    T dt_compressible=0;
    for(int axis=1;axis<=T_GRID::dimension;axis++){T max_u_over_dx=max_lambdas[axis]*one_over_dx[axis];
        dt_compressible+=max_u_over_dx+sqrt(max_u_over_dx*max_u_over_dx+((T)4*max_grad_p_over_rho[axis])*one_over_dx[axis]);}
    dt_compressible*=(T).5;
    if(dt_compressible!=0) dt_compressible=cfl/dt_compressible;
   
    LOG::cout<<"Precomputed Dt: "<<dt<<std::endl;
    if(dt_compressible!=0) LOG::cout<<"Compressible Dt: "<<dt_compressible<<std::endl;
    if(dt_incompressible!=0) LOG::cout<<"Incompressible Dt: "<<dt_incompressible<<std::endl;

    if(dt_incompressible!=0 && dt_compressible!=0) dt=min(dt,min(dt_compressible,dt_incompressible));
    else if(dt_incompressible!=0) dt=min(dt,dt_incompressible);
    else if(dt_compressible!=0) dt=min(dt,dt_compressible);
}
//#####################################################################
// Log_Parameters
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Log_Parameters() const
{
    BASE::Log_Parameters();
    LOG::SCOPE scope("COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE parameters");
    LOG::cout<<"CFL Number = "<<cfl<<std::endl;
    LOG::cout<<"Scale = "<<scale<<std::endl;
    LOG::cout<<"ENO Order = "<<eno_order<<std::endl;
    LOG::cout<<"Use ENO RF = "<<use_rf<<std::endl;
    LOG::cout<<"R-K Order = "<<rk_order<<std::endl;
}
//#####################################################################
// Register_Options
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Register_Options()
{
    BASE::Register_Options();
    if(!parse_args) return;

    // Custom Stuff
    parse_args->Add_Double_Argument("-cfl",(T).5,"Compressible Flow CFL number.");
    parse_args->Add_Integer_Argument("-scale",100,"Resolution of the fluid mesh.");
    parse_args->Add_Integer_Argument("-eno_order",3,"Spatial order of accuracy for flux calculation.");
    parse_args->Add_Integer_Argument("-rk_order",3,"Temporal order of accuracy for time stepping.");
    parse_args->Add_Option_Argument("-write_debug_data","Write debug data.");
    parse_args->Add_Option_Argument("-rf","Use ENO RF.");
}
//#####################################################################
// Parse_Options
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Parse_Options()
{
    BASE::Parse_Options();
    scale=parse_args->Get_Integer_Value("-scale");
    eno_order=parse_args->Get_Integer_Value("-eno_order");
    rk_order=parse_args->Get_Integer_Value("-rk_order");
    cfl=(T)parse_args->Get_Double_Value("-cfl");
    write_debug_data=parse_args->Is_Value_Set("-write_debug_data");
    use_rf=parse_args->Is_Value_Set("-rf");
}
//#####################################################################
// Set_Eigensystems
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Set_Eigensystems()
{
    Set_Eigensystems_Helper(cfl_eigensystem);
}
//#####################################################################
// Function Set_Eigensystems_Helper
//#####################################################################
template<class T> void Set_Eigensystems_Helper(VECTOR<EIGENSYSTEM<T,VECTOR<T,3> >*,1>& eigensystems)
{
    if(eigensystems[1]) delete eigensystems[1];eigensystems[1]=new EULER_1D_EIGENSYSTEM_F_ADVECTION_ONLY<T>();
}
template<class T> void Set_Eigensystems_Helper(VECTOR<EIGENSYSTEM<T,VECTOR<T,4> >*,2>& eigensystems)
{
    if(eigensystems[1]) delete eigensystems[1];eigensystems[1]=new EULER_2D_EIGENSYSTEM_F_ADVECTION_ONLY<T>();
    if(eigensystems[2]) delete eigensystems[2];eigensystems[2]=new EULER_2D_EIGENSYSTEM_G_ADVECTION_ONLY<T>();
}
//#####################################################################
// Parse_Late_Options
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<TV>::
Parse_Late_Options()
{
    BASE::Parse_Late_Options();
    Set_Eigensystems();

    if(use_rf) conservation_law_solver=new CONSERVATION_ENO_RF<GRID<TV>,TV::dimension+2>();
    else conservation_law_solver=new CONSERVATION_ENO_LLF<GRID<TV>,TV::dimension+2>();
    conservation_law_solver->Set_Order(eno_order);
    eno_advection.Set_Order(eno_order);
    Initialize_Grids();
}
//#####################################################################
template class COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<VECTOR<float,1> >;
template class COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<VECTOR<float,2> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<VECTOR<double,1> >;
template class COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE<VECTOR<double,2> >;
#endif
