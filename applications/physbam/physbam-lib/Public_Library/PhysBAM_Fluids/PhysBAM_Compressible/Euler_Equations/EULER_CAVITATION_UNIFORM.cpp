//#####################################################################
// Copyright 2010, Mridul Aanjaneya, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UTILITIES.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/LAPLACE_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_CAVITATION_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS.h>
using namespace PhysBAM;

//#####################################################################
// Constructor
//#####################################################################
template<class TV> EULER_CAVITATION_UNIFORM<TV>::
EULER_CAVITATION_UNIFORM(EULER_UNIFORM<T_GRID>& euler_input, const bool clamp_density_input, const T epsilon_input)
:euler(euler_input),epsilon(epsilon_input),clamp_density(clamp_density_input)
{
    elliptic_solver=new LAPLACE_COLLIDABLE_UNIFORM<T_GRID>(euler.grid,p_cavitation,false,false,true);
    Initialize_Grid();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> EULER_CAVITATION_UNIFORM<TV>::
~EULER_CAVITATION_UNIFORM()
{
    delete elliptic_solver;
}
//#####################################################################
// Initialize_Grid
//#####################################################################
template<class TV> void EULER_CAVITATION_UNIFORM<TV>::
Initialize_Grid()
{
    elliptic_solver->Initialize_Grid(euler.grid);
    p_cavitation.Resize(euler.grid.Domain_Indices(1));
    clamped_momentum_divergence.Resize(euler.grid.Domain_Indices(0));
    clamped_internal_energy_divergence.Resize(euler.grid.Domain_Indices(0));
}
//#####################################################################
// Fill_Ghost_Pressures_Along_Neumann_Boundaries
//#####################################################################
template<class TV> void EULER_CAVITATION_UNIFORM<TV>::
Fill_Ghost_Pressures_Along_Neumann_Boundaries()
{
    // iterate over boundary faces, and if neumann face, copy from second_cell to first_cell
    for(FACE_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){
        int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
        if(elliptic_solver->psi_N.Component(axis)(face_index)){
            TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
            bool first_cell_inside_object=elliptic_solver->psi_D(first_cell_index);
            bool second_cell_inside_object=elliptic_solver->psi_D(second_cell_index);
            bool first_cell_inside_domain=euler.grid.Domain_Indices().Lazy_Inside(first_cell_index);
            bool second_cell_inside_domain=euler.grid.Domain_Indices().Lazy_Inside(second_cell_index);

            if((first_cell_inside_object&&(!second_cell_inside_object))||(second_cell_inside_domain&&(!first_cell_inside_domain)))
                p_cavitation(first_cell_index)=p_cavitation(second_cell_index);
            else if((second_cell_inside_object&&(!first_cell_inside_object))||(first_cell_inside_domain&&(!second_cell_inside_domain)))
                p_cavitation(second_cell_index)=p_cavitation(first_cell_index);
            else PHYSBAM_FATAL_ERROR("Error while copying pressures across Neumann boundaries!");}}
}
//#####################################################################
// Compute_Clamped_Momentum_Divergence
//#####################################################################
template<class TV> void EULER_CAVITATION_UNIFORM<TV>::
Compute_Clamped_Momentum_Divergence(const T dt)
{
    T one_over_dt = (T)1/dt;

    // Compute indices of those cells where density is greater than 2*epsilon
    T_ARRAYS_BOOL sufficient_density_cells(euler.grid.Domain_Indices(0));
    T total_deficient_density=(T)0;
    T total_donor_density=(T)0;

    for(CELL_ITERATOR iterator(euler.grid,0);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        clamped_momentum_divergence(cell_index)=0;
        if(elliptic_solver->psi_D(cell_index)) continue;
        sufficient_density_cells(cell_index)=false;
        if(euler.U(cell_index)(1)>5*epsilon){
            sufficient_density_cells(cell_index)=true;
            total_donor_density+=(euler.U(cell_index)(1)-5*epsilon);}
        else if(euler.U(cell_index)(1)<epsilon){
            clamped_momentum_divergence(cell_index)=min((T)0,(euler.U(cell_index)(1) - epsilon)*one_over_dt);
            std::stringstream ss;ss<<"clamping density at cell_index="<<cell_index<<", density="<<euler.U(cell_index)(1)<<", epsilon="<<epsilon<<std::endl;LOG::filecout(ss.str());
            total_deficient_density+=(epsilon-euler.U(cell_index)(1));}}

    T fractional_contribution_from_each_cell=total_deficient_density/total_donor_density;
    std::stringstream ss;ss<<"Fractional contribution from each cell:"<<fractional_contribution_from_each_cell<<std::endl;LOG::filecout(ss.str());

    if(fractional_contribution_from_each_cell >= 1) PHYSBAM_FATAL_ERROR();

    for(CELL_ITERATOR iterator(euler.grid,0);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(sufficient_density_cells(cell_index))
            clamped_momentum_divergence(cell_index)=fractional_contribution_from_each_cell*(euler.U(cell_index)(1)-5*epsilon)*one_over_dt;}
}
//#####################################################################
// Compute_Clamped_Internal_Energy_Divergence
//#####################################################################
template<class TV> void EULER_CAVITATION_UNIFORM<TV>::
Compute_Clamped_Internal_Energy_Divergence(const T dt)
{
    T one_over_dt = (T)1/dt;

    // Compute indices of those cells where density greater than 2*epsilon
    T_ARRAYS_BOOL sufficient_internal_energy_cells(euler.grid.Domain_Indices(0));
    T total_deficient_energy=(T)0;
    T total_donor_energy=(T)0;

    for(CELL_ITERATOR iterator(euler.grid,0);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        clamped_internal_energy_divergence(cell_index)=0;
        if(elliptic_solver->psi_D(cell_index)) continue;
        sufficient_internal_energy_cells(cell_index)=false;
        if(euler.U(cell_index)(1)*EULER<T_GRID>::e(euler.U,cell_index)>5*epsilon){
            sufficient_internal_energy_cells(cell_index)=true;
            total_donor_energy+=(euler.U(cell_index)(1)*EULER<T_GRID>::e(euler.U,cell_index)-5*epsilon);}
        else if(euler.U(cell_index)(1)*EULER<T_GRID>::e(euler.U,cell_index)<epsilon){
            clamped_internal_energy_divergence(cell_index)=min((T)0,(euler.U(cell_index)(1)*EULER<T_GRID>::e(euler.U,cell_index)-epsilon)*one_over_dt);
            std::stringstream ss;ss<<"clamping energy at cell_index="<<cell_index<<", internal energy="<<euler.U(cell_index)(1)*EULER<T_GRID>::e(euler.U,cell_index)<<", epsilon="<<epsilon<<std::endl;LOG::filecout(ss.str());
            total_deficient_energy+=(epsilon-euler.U(cell_index)(1)*EULER<T_GRID>::e(euler.U,cell_index));}}

    T fractional_contribution_from_each_cell=total_deficient_energy/total_donor_energy;
    std::stringstream ss;ss<<"Fractional contribution from each cell:"<<fractional_contribution_from_each_cell<<std::endl;LOG::filecout(ss.str());

    if(fractional_contribution_from_each_cell >= 1) PHYSBAM_FATAL_ERROR();

    for(CELL_ITERATOR iterator(euler.grid,0);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(sufficient_internal_energy_cells(cell_index))
            clamped_internal_energy_divergence(cell_index)=fractional_contribution_from_each_cell*(euler.U(cell_index)(1)*EULER<T_GRID>::e(euler.U,cell_index)-5*epsilon)*one_over_dt;}
}
//#####################################################################
// Compute_Right_Hand_Side
//#####################################################################
template<class TV> void EULER_CAVITATION_UNIFORM<TV>::
Compute_Right_Hand_Side(const T dt)
{
    if(clamp_density){
        Compute_Clamped_Momentum_Divergence(dt);
        for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){
            TV_INT cell_index=iterator.Cell_Index();
            elliptic_solver->f(cell_index) = -clamped_momentum_divergence(cell_index);}}
    else{
        Compute_Clamped_Internal_Energy_Divergence(dt);
        for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){
            TV_INT cell_index=iterator.Cell_Index();
            elliptic_solver->f(cell_index) = -clamped_internal_energy_divergence(cell_index);}}
}
//#####################################################################
// Compute_Pressure
//#####################################################################
template<class TV> void EULER_CAVITATION_UNIFORM<TV>::
Compute_Pressure(const T dt,const T time)
{
    Compute_Right_Hand_Side(dt);

    elliptic_solver->Find_Solution_Regions(); // flood fill
    elliptic_solver->Solve(time,true); // solve all regions
    
    Fill_Ghost_Pressures_Along_Neumann_Boundaries();
}
//#####################################################################
// Apply_Pressure_To_Density
//#####################################################################
template<class TV> void EULER_CAVITATION_UNIFORM<TV>::
Apply_Pressure_To_Density(const T dt)
{
    T_FACE_ARRAYS_SCALAR grad_p_cavitation_face(euler.grid);
    ARRAYS_UTILITIES<T_GRID,T>::Compute_Gradient_At_Faces_From_Cell_Data(euler.grid,grad_p_cavitation_face,p_cavitation);

    T_ARRAYS_SCALAR laplacian_p_cavitation_cell(euler.grid.Domain_Indices(0));
    ARRAYS_UTILITIES<T_GRID,T>::Compute_Divergence_At_Cells_From_Face_Data(euler.grid,laplacian_p_cavitation_cell,grad_p_cavitation_face);

    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        euler.U(cell_index)(1) += laplacian_p_cavitation_cell(cell_index)*dt;}

    euler.Invalidate_Ghost_Cells();
}
//#####################################################################
// Is_Density_Clamped
//#####################################################################
template<class TV> bool EULER_CAVITATION_UNIFORM<TV>::
Is_Density_Clamped()
{
    return clamp_density;
}
//#####################################################################
// Apply_Pressure_To_Internal_Energy
//#####################################################################
template<class TV> void EULER_CAVITATION_UNIFORM<TV>::
Apply_Pressure_To_Internal_Energy(const T dt)
{
    T_FACE_ARRAYS_SCALAR grad_p_cavitation_face(euler.grid);
    ARRAYS_UTILITIES<T_GRID,T>::Compute_Gradient_At_Faces_From_Cell_Data(euler.grid,grad_p_cavitation_face,p_cavitation);

    T_ARRAYS_SCALAR laplacian_p_cavitation_cell(euler.grid.Domain_Indices(0));
    ARRAYS_UTILITIES<T_GRID,T>::Compute_Divergence_At_Cells_From_Face_Data(euler.grid,laplacian_p_cavitation_cell,grad_p_cavitation_face);

    for(CELL_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        euler.U(cell_index)(T_GRID::dimension+2) += laplacian_p_cavitation_cell(cell_index)*dt;}

    euler.Invalidate_Ghost_Cells();
}
template<class TV> void EULER_CAVITATION_UNIFORM<TV>::
Compute_Face_Pressure_From_Cell_Pressures(const T_GRID& face_grid,T_FACE_ARRAYS_SCALAR& p_face,const T_ARRAYS_SCALAR& p_cell)
{
    TV_INT first_cell_index,second_cell_index;int axis;
    for(FACE_ITERATOR iterator(face_grid);iterator.Valid();iterator.Next()){
        first_cell_index=iterator.First_Cell_Index();second_cell_index=iterator.Second_Cell_Index();axis=iterator.Axis();
        p_face.Component(axis)(iterator.Face_Index())=(p_cell(first_cell_index)+p_cell(second_cell_index))*(T).5;}
}
//#####################################################################
// Apply_Pressure
//#####################################################################
template<class TV> void EULER_CAVITATION_UNIFORM<TV>::
Apply_Pressure(const T dt,const T time, T_FACE_ARRAYS_SCALAR& face_velocities)
{
    if(clamp_density) 
    {
        // Store time star density values
        T_ARRAYS_SCALAR rho_star(euler.grid.Domain_Indices(0));
        for(CELL_ITERATOR iterator(euler.grid,0);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            rho_star(cell_index)=euler.U_ghost(cell_index)(1);}

        Apply_Pressure_To_Density(dt);

        euler.Fill_Ghost_Cells(dt,time,1);

        T_FACE_ARRAYS_SCALAR p_face;
        p_face.Resize(euler.grid);
        Compute_Face_Pressure_From_Cell_Pressures(euler.grid,p_face,p_cavitation);

        T_FACE_ARRAYS_SCALAR grad_p_cavitation_face(euler.grid);
        ARRAYS_UTILITIES<T_GRID,T>::Compute_Gradient_At_Faces_From_Cell_Data(euler.grid,grad_p_cavitation_face,p_cavitation);
        for(FACE_ITERATOR iterator(euler.grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face_index=iterator.Full_Index();
            TV_INT first_cell_index=iterator.First_Cell_Index(), second_cell_index=iterator.Second_Cell_Index();
            T rho_star_face=(T).5*(rho_star(first_cell_index)+rho_star(second_cell_index));
            T rho_np1_face=(T).5*(euler.U_ghost(first_cell_index)(1)+euler.U_ghost(second_cell_index)(1));
            face_velocities(face_index)=(rho_star_face*face_velocities(face_index) - dt*grad_p_cavitation_face(face_index))/rho_np1_face;}


        //T_FACE_ARRAYS_SCALAR face_velocities;
        //face_velocities.Resize(euler.grid);
        //EULER_PROJECTION_UNIFORM<T_GRID>::Compute_Density_Weighted_Face_Velocities(euler.grid,face_velocities,euler.U_ghost,elliptic_solver->psi_N);

        T_ARRAYS_SCALAR density_scaling(euler.grid.Domain_Indices(1));density_scaling.Fill((T)1);
        EULER_PROJECTION_UNIFORM<T_GRID>::Apply_Pressure(p_cavitation,p_face,face_velocities,elliptic_solver->psi_D,elliptic_solver->psi_N,dt,time,density_scaling,0,&euler);
    }
    else
        Apply_Pressure_To_Internal_Energy(dt);
}
//#####################################################################
// Apply_Pressure
//#####################################################################
template<class TV> void EULER_CAVITATION_UNIFORM<TV>::
Apply_Cavitation_Correction(const T dt,const T time, T_FACE_ARRAYS_SCALAR& face_velocities)
{
    Compute_Pressure(dt,time);
    Apply_Pressure(dt,time,face_velocities);
}
//#####################################################################
// Log_Parameters
//#####################################################################
template<class TV> void EULER_CAVITATION_UNIFORM<TV>::
Log_Parameters() const
{
    LOG::SCOPE scope("EULER_CAVITATION_UNIFORM parameters");
    std::stringstream ss;ss<<"epsilon="<<epsilon<<std::endl;LOG::filecout(ss.str());
}

//#####################################################################
template class EULER_CAVITATION_UNIFORM<VECTOR<float,1> >;
template class EULER_CAVITATION_UNIFORM<VECTOR<float,2> >;
template class EULER_CAVITATION_UNIFORM<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EULER_CAVITATION_UNIFORM<VECTOR<double,1> >;
template class EULER_CAVITATION_UNIFORM<VECTOR<double,2> >;
template class EULER_CAVITATION_UNIFORM<VECTOR<double,3> >;
#endif
