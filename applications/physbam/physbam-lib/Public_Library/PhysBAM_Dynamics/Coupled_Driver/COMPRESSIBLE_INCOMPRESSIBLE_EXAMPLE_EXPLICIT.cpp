//#####################################################################
// Copyright 2012, Mridul Aanjaneya.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Dynamics/Coupled_Driver/COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_1D_EIGENSYSTEM_F.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_2D_EIGENSYSTEM_F.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_2D_EIGENSYSTEM_G.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_RF.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/IMPLICIT_VISCOSITY_UNIFORM.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_STRUCTURE.h>
#include <PhysBAM_Dynamics/Particles/PARTICLES_FORWARD.h>
#include <PhysBAM_Dynamics/Read_Write/Particles/READ_WRITE_PARTICLES.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UTILITIES.h>
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
template<class TV> COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT(const STREAM_TYPE& stream_type): 
    BASE(stream_type),mpi_grid(0),rho_incompressible(1e3),viscosity(0),surface_tension(0),scale(0),number_of_cells_to_extrapolate(7),apply_viscosity(false),
    boundary(0),phi_boundary(0),compressible_boundary(0),compressible_boundary_scalar(0),conservation_law_solver(0),particle_levelset_evolution(grid,3),
    projection(grid,pressure,false,false,false),collision_bodies_affecting_fluid(grid)
{
    LOG::Initialize_Logging(false,false,1<<30,true,1);
    Initialize_Particles();
    Initialize_Read_Write_Structures();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
~COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT()
{
    if(mpi_grid) delete compressible_boundary;
    for(int i=1;i<=TV::dimension;i++) delete cfl_eigensystem[i];
    delete conservation_law_solver;
}
//#####################################################################
// Write_Output_Files
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
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
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+f+"/entropy",entropy);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+f+"/internal_energy",internal_energy);}
}
//#####################################################################
// Function Write_Substep
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
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
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
Initialize_Grids()
{
    U.Resize(grid.Domain_Indices());
    U_ghost.Resize(grid.Domain_Indices(3));
    face_velocities.Resize(grid);
    compressible_face_velocities.Resize(grid.Domain_Indices(3));
    compressible_pressure.Resize(grid.Domain_Indices());
    entropy.Resize(grid.Domain_Indices());
    pressure.Resize(grid.Domain_Indices(1));
    psi.Resize(grid.Domain_Indices());
    psi_N.Resize(grid.Domain_Indices(3),true,false);
    valid_mask.Resize(grid.Domain_Indices(3),true,true,true);
}
//#####################################################################
// Read_Output_Files
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
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
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
Compute_Divergence(const T_FACE_LOOKUP& face_lookup)
{
    TV one_over_dx=grid.one_over_dX;
    // Compute divergence at cell centers
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        const typename T_FACE_LOOKUP::LOOKUP& lookup=face_lookup.Starting_Point_Cell(cell);T divergence=0;
        if(particle_levelset_evolution.phi(cell)<=0){for(int axis=1;axis<=T_GRID::dimension;axis++){
            divergence+=(lookup(axis,iterator.Second_Face_Index(axis))-lookup(axis,iterator.First_Face_Index(axis)))*one_over_dx[axis];}}
        projection.f(cell)=divergence;}
}
//#####################################################################
// Compute_Right_Hand_Side
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
Compute_Right_Hand_Side(const T dt)
{
    Compute_Divergence(T_FACE_LOOKUP(face_velocities));
}
//#####################################################################
// Fill_Face_Weights_For_Projection
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
Fill_Face_Weights_For_Projection(const T dt,const T time)
{
    projection.Set_Variable_beta(true);
    for(FACE_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> full_index=iterator.Full_Index();
        projection.beta_face(full_index)=(T)1./rho_incompressible;}
}
//#####################################################################
// Compute_Pressure
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
Compute_Pressure(const T dt,const T time)
{
    pressure*=dt;
    Compute_Right_Hand_Side(dt);
    Fill_Face_Weights_For_Projection(dt,time);
    projection.Find_Solution_Regions();     // flood fill
    projection.Solve(time,true);            // solve all regions
    boundary->Apply_Boundary_Condition(grid,pressure,time);
    pressure*=((T)1./dt);
}
//#####################################################################
// Project
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
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

    Set_Boundary_Conditions(time);
    Compute_Pressure(dt,time);

    T_ARRAYS_SCALAR p_ghost(grid.Domain_Indices(1));
    boundary->Fill_Ghost_Cells(grid,pressure,p_ghost,dt,time,1);
    Apply_Pressure(p_ghost,dt,time);
}
//#####################################################################
// Apply_Pressure
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
Apply_Pressure(T_ARRAYS_SCALAR& p_ghost,const T dt,const T time)
{
    TV one_over_dx=grid.one_over_dX;
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> full_index=iterator.Full_Index();
        TV_INT first_cell=iterator.First_Cell_Index(),second_cell=iterator.Second_Cell_Index();int axis=iterator.Axis();
        T phi_face=(T).5*(particle_levelset_evolution.phi(first_cell)+particle_levelset_evolution.phi(second_cell));
        if(phi_face<=0) face_velocities(full_index)-=dt*(p_ghost(second_cell)-p_ghost(first_cell))*one_over_dx[axis]/rho_incompressible;}
}
//#####################################################################
// Extrapolate_Velocity_Across_Interface
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
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
template<class TV> template<class T_TYPE> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
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
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
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
            U(cell)=EULER<T_GRID>::Get_Euler_State_From_rho_velocity_And_internal_energy(rho,velocity(cell),internal_energy);}
        if(particle_levelset_evolution.phi(cell)>0) pressure(cell)=compressible_pressure(cell);}
}
//#####################################################################
// Limit_Dt
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
Limit_Dt(T& dt,const T time)
{
    TV max_lambdas=TV::Constant_Vector((T)0.);T dt_compressible=(T)0.,dt_incompressible=(T)0.;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();TV X=iterator.Location();
        if(psi(cell)){for(int k=1;k<=T_GRID::dimension;++k){
            max_lambdas(k)=max(max_lambdas(k),cfl_eigensystem(k)->Maximum_Magnitude_Eigenvalue(U(cell)));}}
        else{T local_V_norm=(T)0.,local_viscosity_norm=(T)0.;
            for(int axis=1;axis<=T_GRID::dimension;axis++){
                local_V_norm+=grid.one_over_dX[axis]*maxabs(face_velocities(axis,grid.First_Face_Index_In_Cell(axis,cell)),face_velocities(axis,grid.Second_Face_Index_In_Cell(axis,cell)));
                local_viscosity_norm+=grid.one_over_dX[axis]*grid.one_over_dX[axis];}
            local_viscosity_norm*=((T)2.*viscosity/rho_incompressible);
            T local_surface_tension_norm_squared=-surface_tension*particle_levelset_evolution.particle_levelset.levelset.Compute_Curvature(X)/(rho_incompressible*grid.min_dX*grid.min_dX);
            T local_dt=(T).5*(local_V_norm+local_viscosity_norm+sqrt((local_V_norm+local_viscosity_norm)*(local_V_norm+local_viscosity_norm)+(T)4.*local_surface_tension_norm_squared));
            dt_incompressible=max(dt_incompressible,local_dt);}}
    if(dt_incompressible!=0) dt_incompressible=cfl/dt_incompressible;
    TV max_lambda_over_dx=max_lambdas*grid.one_over_dX;
    if(max_lambda_over_dx.Sum()!=0) dt_compressible=cfl/max_lambda_over_dx.Sum();

    LOG::cout<<"Precomputed Dt: "<<dt<<std::endl;
    if(dt_compressible!=0) LOG::cout<<"Compressible Dt: "<<dt_compressible<<std::endl;
    if(dt_incompressible!=0) LOG::cout<<"Incompressible Dt: "<<dt_incompressible<<std::endl;

    if(dt_incompressible!=0 && dt_compressible!=0) dt=min(dt,min(dt_compressible,dt_incompressible));
    else if(dt_compressible==0) dt=min(dt,dt_incompressible);
    else if(dt_incompressible==0) dt=min(dt,dt_compressible);
}
//#####################################################################
// Adjust_Particle_For_Domain_Boundaries
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
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
            if(collision_distance>max_collision_distance) collision_distance=X_new[axis]-min_corner[axis];
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
// Log_Parameters
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
Log_Parameters() const
{
    BASE::Log_Parameters();
    LOG::SCOPE scope("COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT parameters");
    LOG::cout<<"CFL Number = "<<cfl<<std::endl;
    LOG::cout<<"Scale = "<<scale<<std::endl;
    LOG::cout<<"ENO Order = "<<eno_order<<std::endl;
    LOG::cout<<"Use ENO RF = "<<use_rf<<std::endl;
    LOG::cout<<"R-K Order = "<<rk_order<<std::endl;
}
//#####################################################################
// Register_Options
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
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
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
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
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
Set_Eigensystems()
{
    Set_Eigensystems_Helper(cfl_eigensystem);
}
//#####################################################################
// Function Set_Eigensystems_Helper
//#####################################################################
template<class T> void Set_Eigensystems_Helper(VECTOR<EIGENSYSTEM<T,VECTOR<T,3> >*,1>& eigensystems)
{
    if(eigensystems[1]) delete eigensystems[1];eigensystems[1]=new EULER_1D_EIGENSYSTEM_F<T>();
}
template<class T> void Set_Eigensystems_Helper(VECTOR<EIGENSYSTEM<T,VECTOR<T,4> >*,2>& eigensystems)
{
    if(eigensystems[1]) delete eigensystems[1];eigensystems[1]=new EULER_2D_EIGENSYSTEM_F<T>();
    if(eigensystems[2]) delete eigensystems[2];eigensystems[2]=new EULER_2D_EIGENSYSTEM_G<T>();
}
//#####################################################################
// Parse_Late_Options
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<TV>::
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
template class COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<VECTOR<float,1> >;
template class COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<VECTOR<float,2> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<VECTOR<double,1> >;
template class COMPRESSIBLE_INCOMPRESSIBLE_EXAMPLE_EXPLICIT<VECTOR<double,2> >;
#endif
