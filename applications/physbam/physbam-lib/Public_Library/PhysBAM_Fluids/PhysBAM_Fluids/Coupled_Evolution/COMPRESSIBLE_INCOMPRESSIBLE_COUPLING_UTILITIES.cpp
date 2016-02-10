//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COUPLING_UTILITIES
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Fluids/Coupled_Evolution/COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES.h>
namespace PhysBAM{
//#####################################################################
// Function Constructor
//#####################################################################
template<class TV> COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<TV>::
COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR* incompressible_face_velocities,const T* incompressible_density,
    const T_ARRAYS_SCALAR* incompressible_phi):incompressible_face_velocities_(*incompressible_face_velocities),
    incompressible_density_(*incompressible_density),incompressible_phi_(*incompressible_phi)
{
    p_dirichlet_incompressible.Resize(grid.Domain_Indices(1));
}
//#####################################################################
// Function Extrapolate_Compressible_State_Into_Incompressible_Region
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<TV>::
Extrapolate_Compressible_State_Into_Incompressible_Region(const T dt,const T time,const T bandwidth,const int ghost_cells,const EOS<T>& eos,
    const T_GRID& grid,const T_ARRAYS_SCALAR& phi_ghost,const T_FACE_ARRAYS_SCALAR& incompressible_face_velocities,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,T_ARRAYS_DIMENSION_SCALAR& U)
{
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    T_ARRAYS_SCALAR phi_ghost_negated(phi_ghost),entropy(phi_ghost,false),pressure(phi_ghost,false);T_ARRAYS_VECTOR velocity(phi_ghost.Domain_Indices());

    for(CELL_ITERATOR iterator(grid,ghost_cells);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();T density=U_ghost(cell_index)(1);TV vel=EULER<T_GRID>::Get_Velocity(U_ghost,cell_index);
        T e=EULER<T_GRID>::e(U_ghost(cell_index));
        phi_ghost_negated(cell_index)*=-1;
        pressure(cell_index)=eos.p(density,e);
        entropy(cell_index)=eos.S(density,e);
        velocity(iterator.Cell_Index())=vel;}

    T_EXTRAPOLATION_SCALAR extrapolate_entropy(grid,phi_ghost_negated,entropy,ghost_cells);extrapolate_entropy.Set_Band_Width(bandwidth);extrapolate_entropy.Extrapolate(time,false);
    T_EXTRAPOLATION_SCALAR extrapolate_pressure(grid,phi_ghost_negated,pressure,ghost_cells);extrapolate_pressure.Set_Band_Width(bandwidth);extrapolate_pressure.Extrapolate(time,false);
    T_EXTRAPOLATION_VECTOR extrapolate_velocity(grid,phi_ghost_negated,velocity,ghost_cells);extrapolate_velocity.Set_Band_Width(bandwidth);extrapolate_velocity.Extrapolate(time,false);
    // Find the tangential component and combine it with the normal component from the incompressible region
    TV one_over_two_dx=(T).5*grid.one_over_dX;
    for(CELL_ITERATOR iterator(grid,ghost_cells);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(phi_ghost_negated(cell_index)>0 && phi_ghost_negated(cell_index)<=bandwidth){
            TV v_ext,v_i,grad_phi,N,v_total;
            v_ext=velocity(cell_index);
            for(int axis=1;axis<=T_GRID::dimension;axis++){
                TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                v_i(axis)=(incompressible_face_velocities.Component(axis)(iterator.First_Face_Index(axis))+
                    incompressible_face_velocities.Component(axis)(iterator.Second_Face_Index(axis)))*(T).5;
                grad_phi(axis)=(phi_ghost(cell_index+axis_vector)-phi_ghost(cell_index-axis_vector))*one_over_two_dx(axis);}
            N=grad_phi.Normalized();
            v_total=(TV::Dot_Product(v_i,N)*N)+(v_ext-(TV::Dot_Product(v_ext,N))*N);
            // Get the conserved variables from entropy, velocity, and pressure
            T rho=eos.rho_From_p_And_S(pressure(cell_index),entropy(cell_index));
            EULER<T_GRID>::Set_Euler_State_From_rho_velocity_And_internal_energy(U,cell_index,rho,v_total,eos.e_From_S_And_rho(entropy(cell_index),rho));}}
}
//#####################################################################
// Function Get_Dirichlet_Boundary_Conditions_For_Incompressible_Region
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<TV>::
Get_Dirichlet_Boundary_Conditions_For_Incompressible_Region(const T_GRID& grid,const T_ARRAYS_DIMENSION_SCALAR& U_dirichlet,const EOS<T>& euler_eos,const T incompressible_density,const T dt)
{
    TV_INT cell_index;
    for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){
        cell_index=iterator.Cell_Index();
        const EOS_GAMMA<T> *eos_gamma=dynamic_cast<const EOS_GAMMA<T>*>(&euler_eos);
        // TODO(kwatra): This does not look right. It seems to be calculating p=(gamma-1)e
        p_dirichlet_incompressible(cell_index)=dt*(1/incompressible_density)*(eos_gamma->gamma-1)*
            (U_dirichlet(cell_index)(TV::dimension+2)-((T).5*U_dirichlet(cell_index)(1)*EULER<T_GRID>::Get_Velocity(U_dirichlet,cell_index).Magnitude_Squared()));}
}
//#####################################################################
// Function Compute_Compressible_Incompressible_Face_Velocities
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<TV>::
Compute_Compressible_Incompressible_Face_Velocities(const T_GRID& face_grid,const T_FACE_ARRAYS_SCALAR& incompressible_face_velocities,const T incompressible_density,
    const T_ARRAYS_SCALAR& incompressible_phi,const T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_BOOL& euler_psi,T_FACE_ARRAYS_SCALAR& compressible_face_velocities)
{
    TV_INT first_cell_index,second_cell_index,compressible_cell_index;
    for(FACE_ITERATOR iterator(face_grid);iterator.Valid();iterator.Next()){
        first_cell_index=iterator.First_Cell_Index();second_cell_index=iterator.Second_Cell_Index();
        bool first_cell_euler=face_grid.Inside_Domain(first_cell_index)&&euler_psi(first_cell_index),
             second_cell_euler=face_grid.Inside_Domain(second_cell_index)&&euler_psi(second_cell_index),
             first_cell_incompressible=face_grid.Inside_Domain(first_cell_index)&&incompressible_phi(first_cell_index)<0,
             second_cell_incompressible=face_grid.Inside_Domain(second_cell_index)&&incompressible_phi(second_cell_index)<0;

        TV_INT face_index=iterator.Face_Index();int axis=iterator.Axis();
        if((first_cell_euler&&second_cell_incompressible) || (second_cell_euler&&first_cell_incompressible)){
            if(first_cell_euler) compressible_cell_index=first_cell_index; else compressible_cell_index=second_cell_index;

            T rho_compressible=U(compressible_cell_index)(1),rho_incompressible=incompressible_density;
            compressible_face_velocities.Component(axis)(face_index)=(rho_incompressible*incompressible_face_velocities.Component(axis)(face_index)+
                rho_compressible*EULER<T_GRID>::Get_Velocity_Component(U,compressible_cell_index,axis))/(rho_compressible+rho_incompressible);}
        else if(first_cell_incompressible && second_cell_incompressible) // TODO(jontg): Move this out.
            compressible_face_velocities.Component(axis)(face_index)=incompressible_face_velocities.Component(axis)(face_index);
    }
}
//#####################################################################
// Function Compute_Compressible_Incompressible_Face_Pressures_From_Cell_Pressures
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<TV>::
Compute_Compressible_Incompressible_Face_Pressures_From_Cell_Pressures(const T_GRID& face_grid,
    const T_FACE_ARRAYS_SCALAR& incompressible_face_velocities,const T incompressible_density,const T_ARRAYS_SCALAR& incompressible_phi,
    const T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_BOOL& euler_psi,const T_ARRAYS_SCALAR& p_cell,T_FACE_ARRAYS_SCALAR& p_face)
{
    TV_INT first_cell_index,second_cell_index,compressible_cell_index;
    for(FACE_ITERATOR iterator(face_grid);iterator.Valid();iterator.Next()){
        first_cell_index=iterator.First_Cell_Index();second_cell_index=iterator.Second_Cell_Index();
        bool first_cell_euler=face_grid.Inside_Domain(first_cell_index)&&euler_psi(first_cell_index),
             second_cell_euler=face_grid.Inside_Domain(second_cell_index)&&euler_psi(second_cell_index),
             first_cell_incompressible=face_grid.Inside_Domain(first_cell_index)&&incompressible_phi(first_cell_index)<0,
             second_cell_incompressible=face_grid.Inside_Domain(second_cell_index)&&incompressible_phi(second_cell_index)<0;

        TV_INT face_index=iterator.Face_Index();int axis=iterator.Axis();
        if((first_cell_euler&&second_cell_incompressible) || (second_cell_euler&&first_cell_incompressible)){
            T rho_first_cell,rho_second_cell;
            T u_star_first_cell,u_star_second_cell;
            if(first_cell_euler){
                rho_first_cell=U(first_cell_index)(1);
                u_star_first_cell=EULER<T_GRID>::Get_Velocity_Component(U,first_cell_index,axis);
                rho_second_cell=incompressible_density;
                u_star_second_cell=incompressible_face_velocities.Component(axis)(face_index);}
            else{
                rho_first_cell=incompressible_density;
                u_star_first_cell=incompressible_face_velocities.Component(axis)(face_index);
                rho_second_cell=U(second_cell_index)(1);
                u_star_second_cell=EULER<T_GRID>::Get_Velocity_Component(U,second_cell_index,axis);}

            T correction_term=(T).5*rho_first_cell*rho_second_cell*face_grid.dX[axis]*(u_star_second_cell-u_star_first_cell)/(rho_first_cell+rho_second_cell);
            p_face.Component(axis)(iterator.Face_Index())=(rho_second_cell*p_cell(first_cell_index)+rho_first_cell*p_cell(second_cell_index))/(rho_first_cell+rho_second_cell)-correction_term;}}
}
//#####################################################################
// Function Compute_Compressible_Incompressible_Face_Pressures_From_Cell_Pressures
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<TV>::
Compute_Compressible_Incompressible_Face_Pressures_From_Cell_Pressures(const T_GRID& face_grid,const T_ARRAYS_DIMENSION_SCALAR& U,
    const T_ARRAYS_BOOL& euler_psi,const T_ARRAYS_SCALAR& p_cell,T_FACE_ARRAYS_SCALAR& p_face) const
{
    Compute_Compressible_Incompressible_Face_Pressures_From_Cell_Pressures(face_grid,incompressible_face_velocities_,incompressible_density_,
        incompressible_phi_,U,euler_psi,p_cell,p_face);
}
//#####################################################################
// Function Fill_Incompressible_Beta_Face
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<TV>::
Fill_Incompressible_Beta_Face(const T_GRID& grid,const T incompressible_density,const T_ARRAYS_SCALAR& incompressible_phi,
    T_FACE_ARRAYS_SCALAR& beta_face)
{
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
        TV_INT face_index=iterator.Face_Index(),first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
        if((incompressible_phi(first_cell_index)+incompressible_phi(second_cell_index))*(T).5 < 0){
            beta_face.Component(axis)(face_index)=(T)1/incompressible_density;}}
}
//#####################################################################
// Function Fill_Incompressible_Beta_Face
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<TV>::
Fill_Incompressible_Beta_Face(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& beta_face) const
{
    Fill_Incompressible_Beta_Face(grid,incompressible_density_,incompressible_phi_,beta_face);
}
//#####################################################################
// Function Apply_Pressure_At_Incompressible_Faces
//#####################################################################
template<class TV> void COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<TV>::
Apply_Pressure_At_Incompressible_Faces(const T_GRID& face_grid,const T incompressible_density,const T_ARRAYS_SCALAR& incompressible_phi,
    const T_FACE_ARRAYS_BOOL& psi_N,const T_ARRAYS_SCALAR& p_hat,T_FACE_ARRAYS_SCALAR& incompressible_face_velocities)
{
    TV one_over_dx=face_grid.one_over_dX;
    // TODO(kwatra): check if we use incompressible->projection.p anywhere?
    //for(CELL_ITERATOR iterator(face_grid);iterator.Valid();iterator.Next())
     //   if(!euler->psi(iterator.Cell_Index()))
      //      incompressible->projection.p(iterator.Cell_Index())=p_hat(iterator.Cell_Index())*(1/incompressible_density);
    for(FACE_ITERATOR iterator(face_grid);iterator.Valid();iterator.Next()){
        int axis=iterator.Axis();
        TV_INT face_index=iterator.Face_Index(),first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
        if(!psi_N.Component(axis)(face_index) && (incompressible_phi(first_cell_index)<0 || incompressible_phi(second_cell_index)<0)){
            T orig_face_vel=incompressible_face_velocities.Component(axis)(face_index);
            T first_cell_pressure=p_hat(first_cell_index),second_cell_pressure=p_hat(second_cell_index);
            T new_face_vel=orig_face_vel-((second_cell_pressure-first_cell_pressure)*one_over_dx[axis]*incompressible_density);
            incompressible_face_velocities.Component(axis)(face_index)=new_face_vel;}}
}
//#####################################################################
template class COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<VECTOR<float,1> >;
template class COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<VECTOR<float,2> >;
template class COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<VECTOR<double,1> >;
template class COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<VECTOR<double,2> >;
template class COMPRESSIBLE_INCOMPRESSIBLE_COUPLING_UTILITIES<VECTOR<double,3> >;
#endif
}
