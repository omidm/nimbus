//#####################################################################
// Copyright 2007, Jon Gretarsson, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UTILITIES.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_LAPLACE.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_PROJECTION_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/INCOMPRESSIBLE_COMPRESSIBLE_COUPLING_CALLBACKS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> EULER_PROJECTION_UNIFORM<T_GRID>::
EULER_PROJECTION_UNIFORM(EULER_UNIFORM<T_GRID>* euler_input)
    :is_pressure_scaled(false),save_fluxes(true),use_exact_neumann_face_location(false),use_neumann_condition_for_outflow_boundaries(false),
    pressure_object_boundary(false),hj_eno_order(3),incompressible_coupling_callbacks(0),transition_to_using_implicit_pressure(false)
{
    euler=euler_input;
    fluxes=new T_FACE_ARRAYS_DIMENSION_SCALAR(euler->grid);
    elliptic_solver=new EULER_LAPLACE<POISSON_COLLIDABLE_UNIFORM<T_GRID> >(euler->grid,p,one_over_rho_c_squared);
    elliptic_solver->Set_Variable_beta(true);
    poisson=elliptic_solver;
    pressure_advection_HJ.Set_Order(hj_eno_order);
    Set_Use_Exact_Neumann_Face_Location(use_exact_neumann_face_location);
    pressure_boundary=&pressure_boundary_default;
    Set_Constant_Extrapolation_Pressure_Boundary();
    Initialize_Grid();
}
//#####################################################################
// Destructor 
//#####################################################################
template<class T_GRID> EULER_PROJECTION_UNIFORM<T_GRID>::
~EULER_PROJECTION_UNIFORM()
{
    delete fluxes;delete elliptic_solver;
}
//#####################################################################
// Function Initialize_Grid
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Initialize_Grid()
{
    fluxes->Resize(euler->grid);
    elliptic_solver->Initialize_Grid(euler->grid);
    p.Resize(euler->grid.Domain_Indices(1));p_save_for_projection.Resize(euler->grid.Domain_Indices(1));
    face_velocities.Resize(euler->grid);face_velocities_save.Resize(euler->grid);
    p_advected.Resize(euler->grid.Domain_Indices());
    density_scaling.Resize(euler->grid.Domain_Indices(1));density_scaling.Fill((T)1);
    one_over_rho_c_squared.Resize(euler->grid.Domain_Indices(1)); // TODO(jontg): This only needs a ghost ring in the mpi case.
    p_dirichlet.Resize(euler->grid.Domain_Indices(1));
}
//#####################################################################
// Function p
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Get_Pressure(T_ARRAYS_SCALAR& pressure) const
{
    if(transition_to_using_implicit_pressure){
        std::stringstream ss;ss<<"using implicit pressure, is_pressure_scaled="<<is_pressure_scaled<<std::endl;LOG::filecout(ss.str());
        T scaling=(T)1;if(is_pressure_scaled) scaling=(T)1/dt_scale_pressure;
        for(CELL_ITERATOR iterator(euler->grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            pressure(cell_index)=p(cell_index)*scaling;}}
    else{ 
        for(CELL_ITERATOR iterator(euler->grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            pressure(cell_index)=euler->eos->p(euler->U(cell_index)(1),euler->e(euler->U,cell_index));}}
}
//#####################################################################
// Function Fill_Face_Weights_For_Projection
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Fill_Face_Weights_For_Projection(const T dt,const T time,T_FACE_ARRAYS_SCALAR& beta_face)
{
    euler->Fill_Ghost_Cells(dt,time,1);
    for(FACE_ITERATOR iterator(euler->grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
        T rho_first_cell=euler->U_ghost(iterator.First_Cell_Index())(1),rho_second_cell=euler->U_ghost(iterator.Second_Cell_Index())(1);
        T rho_face=(rho_first_cell+rho_second_cell)*(T).5;
        beta_face.Component(axis)(face_index)=(T)1/rho_face;}

    if(incompressible_coupling_callbacks){
        incompressible_coupling_callbacks->Fill_Incompressible_Beta_Face(euler->grid,beta_face);}
}
//#####################################################################
// Function Get_Ghost_Density
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Get_Ghost_Density(const T dt,const T time,const int number_of_ghost_cells,T_ARRAYS_SCALAR& density_ghost) const
{
    euler->Fill_Ghost_Cells(dt,time,number_of_ghost_cells);
    for(CELL_ITERATOR iterator(euler->grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        density_ghost(cell_index)=euler->U_ghost(cell_index)(1);}
}
//#####################################################################
// Function Get_Ghost_Density
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Get_Ghost_Centered_Velocity(const T dt,const T time,const int number_of_ghost_cells,T_ARRAYS_VECTOR& centered_velocity_ghost) const
{
    euler->Fill_Ghost_Cells(dt,time,number_of_ghost_cells);
    for(CELL_ITERATOR iterator(euler->grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){
        TV_INT cell_index=iterator.Cell_Index();
        if(!euler->psi.Valid_Index(cell_index) || euler->psi(cell_index))
            centered_velocity_ghost(cell_index)=EULER<T_GRID>::Get_Velocity(euler->U_ghost(cell_index));}
}
//#####################################################################
// Function Make_Boundary_Faces_Neumann
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Make_Boundary_Faces_Neumann(T_FACE_ARRAYS_BOOL& psi_N)
{
    for(int side=1;side<=2*TV::dimension;++side)
        if(!euler->mpi_grid || euler->mpi_grid->side_neighbor_ranks(side)>0)
            for(FACE_ITERATOR iterator(euler->grid,0,T_GRID::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){
                TV_INT face_index=iterator.Face_Index();int axis=iterator.Axis();
                psi_N.Component(axis)(face_index)=true;}
}
//#####################################################################
// Function Project 
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Project(const T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    Compute_Pressure(face_velocities,dt,time);

    T_ARRAYS_SCALAR p_ghost;
    Get_Ghost_Pressures(dt,time,elliptic_solver->psi_D,elliptic_solver->psi_N,p,p_ghost);
    T_FACE_ARRAYS_SCALAR p_face;
    Get_Pressure_At_Faces(dt,time,p_ghost,p_face);

    Apply_Pressure(p_ghost,p_face,face_velocities,elliptic_solver->psi_D,elliptic_solver->psi_N,dt,time);
}
//#####################################################################
// Function Compute_Pressure 
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Compute_Pressure(const T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    Scale_Pressure_By_Dt(dt);

    Compute_Right_Hand_Side(face_velocities,dt,time);

    T_FACE_ARRAYS_BOOL psi_N_copy=elliptic_solver->psi_N;
    if(use_neumann_condition_for_outflow_boundaries) Make_Boundary_Faces_Neumann(elliptic_solver->psi_N);
    assert(Consistent_Boundary_Conditions());
    Fill_Face_Weights_For_Projection(dt,time,elliptic_solver->beta_face);
    elliptic_solver->Set_Dt(dt);
    elliptic_solver->Find_Solution_Regions(); // flood fill
    elliptic_solver->Compute_beta_And_Add_Jumps_To_b(dt,time); // only does something for poisson solver
    elliptic_solver->Solve(time,true); // solve all regions
    pressure_boundary->Apply_Boundary_Condition(euler->grid,p,time);
    if(use_neumann_condition_for_outflow_boundaries) elliptic_solver->psi_N=psi_N_copy;
    Unscale_Pressure_By_Dt(dt);
}
//#####################################################################
// Function Get_Dirichlet_Boundary_Conditions
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Get_Dirichlet_Boundary_Conditions(const T_ARRAYS_DIMENSION_SCALAR& U_dirichlet)
{
    TV_INT cell_index;
    for(CELL_ITERATOR iterator(euler->grid,1,T_GRID::GHOST_REGION);iterator.Valid();iterator.Next()){
        cell_index=iterator.Cell_Index();
        p_dirichlet(cell_index)=euler->eos->p(U_dirichlet(cell_index)(1),euler->e(U_dirichlet,cell_index));}
}
//#####################################################################
// Function Set_Dirichlet_Boundary_Conditions
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Set_Dirichlet_Boundary_Conditions(const T time)
{
    assert(!is_pressure_scaled);
    TV_INT cell_index;
    for(CELL_ITERATOR iterator(euler->grid,1,T_GRID::GHOST_REGION);iterator.Valid();iterator.Next()){
        cell_index=iterator.Cell_Index();
        p(cell_index)=p_dirichlet(cell_index);}

    // TODO: should we set unsolved cells to be dirichlet
    //for(CELL_ITERATOR iterator(euler->grid);iterator.Valid();iterator.Next()){cell_index=iterator.Cell_Index();
    //   if(!euler->psi(cell_index)) elliptic_solver->psi_D(cell_index)=true;}

    if(euler->mpi_grid){
        euler->mpi_grid->Exchange_Boundary_Cell_Data(elliptic_solver->psi_D,1,false);
        euler->mpi_grid->Exchange_Boundary_Cell_Data(p,1,false);}
}
//#####################################################################
// Function Compute_Divergence
//#####################################################################
template<class T_GRID> template<class FACE_LOOKUP> void EULER_PROJECTION_UNIFORM<T_GRID>::
Compute_Divergence(const FACE_LOOKUP& face_lookup)
{
    TV one_over_dx=euler->grid.one_over_dX;
    // Compute divergence at cell centers
    for(CELL_ITERATOR iterator(euler->grid);iterator.Valid();iterator.Next()){
        const typename FACE_LOOKUP::LOOKUP& lookup=face_lookup.Starting_Point_Cell(iterator.Cell_Index());
        T divergence=0;for(int axis=1;axis<=T_GRID::dimension;axis++){
            divergence+=(lookup(axis,iterator.Second_Face_Index(axis))-lookup(axis,iterator.First_Face_Index(axis)))*one_over_dx[axis];}
        elliptic_solver->f(iterator.Cell_Index())=divergence;}
}
//#####################################################################
// Function Compute_Advected_Pressure
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Compute_Advected_Pressure(const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const T_FACE_ARRAYS_SCALAR face_velocities_for_solid_faces,const T dt)
{
    TV one_over_dx=euler->grid.one_over_dX;

    //Advect pressure using Hamilton Jacobi
    T_ARRAYS_SCALAR p_ghost(euler->grid.Domain_Indices(3));
    Get_Pressure(p_ghost);
    pressure_boundary->Fill_Ghost_Cells(euler->grid,p_ghost,p_ghost,dt,0,3);

#if 1
    T_ARRAYS_VECTOR v_cell(euler->grid.Domain_Indices(3));
    for(CELL_ITERATOR iterator(euler->grid,3);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        v_cell(cell_index)=euler->Get_Velocity(U_ghost,cell_index);}

    FLOOD_FILL_1D find_connected_components;
    T_ARRAYS_SCALAR rhs(euler->grid.Domain_Indices());ARRAY<T,VECTOR<int,1> > p_1d,u,u_px;
    for(int axis=1;axis<=T_GRID::dimension;axis++){
        GRID<VECTOR<T,1> > grid_1d=euler->grid.Get_1D_Grid(axis);
        p_1d.Resize(grid_1d.Domain_Indices(3));u.Resize(grid_1d.Domain_Indices(3));u_px.Resize(grid_1d.Domain_Indices(3));
        T_GRID_LOWER_DIM lower_dimension_grid=euler->grid.Remove_Dimension(axis);
        ARRAY<int,VECTOR<int,1> > colors(1,grid_1d.counts.x);
        ARRAY<bool,VECTOR<int,1> > psi_N_axis(1,grid_1d.counts.x+1);
        VECTOR<int,2> region_boundary;
        for(CELL_ITERATOR_LOWER_DIM iterator(lower_dimension_grid);iterator.Valid();iterator.Next()){TV_INT_LOWER_DIM cell_index=iterator.Cell_Index();
            for(int i=-2;i<=grid_1d.counts.x+3;i++) u(i)=v_cell(cell_index.Insert(i,axis))[axis];
            for(int i=1;i<=grid_1d.counts.x;i++) colors(i)=euler->psi(cell_index.Insert(i,axis))?0:-1;
            for(int i=1;i<=grid_1d.counts.x+1;i++) psi_N_axis(i)=elliptic_solver->psi_N(axis,cell_index.Insert(i,axis));
            int number_of_regions=find_connected_components.Flood_Fill(colors,psi_N_axis);
            for(int color=1;color<=number_of_regions;color++){
                for(int i=-2;i<=grid_1d.counts.x+3;i++) p_1d(i)=p_ghost(cell_index.Insert(i,axis));
                region_boundary=find_connected_components.region_boundaries(color);
                VECTOR<bool,2> psi_N_boundary=VECTOR<bool,2>(psi_N_axis(region_boundary.x),psi_N_axis(region_boundary.y+1));
                pressure_object_boundary.Fill_Ghost_Cells_Neumann(grid_1d,p_1d,face_velocities_for_solid_faces,cell_index,axis,hj_eno_order,
                    use_exact_neumann_face_location,VECTOR<int,2>(1,grid_1d.counts.x),region_boundary,psi_N_boundary,0);
                pressure_advection_HJ.Advection_Solver(region_boundary.x,region_boundary.y,grid_1d.dX.x,p_1d,u,u_px);
                for(int k=region_boundary.x;k<=region_boundary.y;k++) rhs(cell_index.Insert(k,axis))+=u_px(k);}}}
    for(CELL_ITERATOR iterator(euler->grid);iterator.Valid();iterator.Next())
        p_advected(iterator.Cell_Index())=p_ghost(iterator.Cell_Index())-dt*rhs(iterator.Cell_Index());
#else
    FACE_LOOKUP_UNIFORM<T_GRID> face_lookup(face_velocities);
    pressure_advection_HJ.Update_Advection_Equation_Cell_Lookup(euler->grid,p_advected,p_ghost,face_lookup,*pressure_boundary,dt,(const T)0);
#endif
}
//#####################################################################
// Function Compute_Right_Hand_Side
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Compute_Right_Hand_Side(const T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time)
{
    Compute_Divergence(T_FACE_LOOKUP(face_velocities));
    for(CELL_ITERATOR iterator(euler->grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(euler->psi(cell_index)) elliptic_solver->f(cell_index)-=p_advected(cell_index)*one_over_rho_c_squared(cell_index)*(1/dt);}
}
//#####################################################################
// Function Compute_One_Over_rho_C_Squared
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Compute_One_Over_rho_c_Squared()
{
    for(CELL_ITERATOR iterator(euler->grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        if(euler->psi(cell_index)){
            T rho=euler->U_ghost(cell_index)(1),e=euler->e(euler->U_ghost,cell_index);
            T c_squared=euler->eos->c_squared(rho,e);
            if(c_squared <= 0) LOG::cerr<<"Cell "<<cell_index<<" with state "<<euler->U_ghost(cell_index)<<" has negative sound speed!"<<std::endl;
            one_over_rho_c_squared(cell_index)=(T)1/(rho*c_squared);}
        else one_over_rho_c_squared(cell_index)=0;}
}
//#####################################################################
// Function Compute_Density_Weighted_Face_Velocities
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Compute_Density_Weighted_Face_Velocities(const T dt,const T time,const T_FACE_ARRAYS_BOOL& psi_N)
{
    euler->Fill_Ghost_Cells(dt,time,1);
    Compute_Density_Weighted_Face_Velocities(euler->grid,face_velocities,euler->U_ghost,euler->psi,psi_N);
}
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Compute_Density_Weighted_Face_Velocities(const T_GRID& face_grid,T_FACE_ARRAYS_SCALAR& face_velocities,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const T_ARRAYS_BOOL& psi,const T_FACE_ARRAYS_BOOL& psi_N)
{
    TV_INT first_cell_index,second_cell_index;int axis;
    for(FACE_ITERATOR iterator(face_grid);iterator.Valid();iterator.Next()){
        first_cell_index=iterator.First_Cell_Index();second_cell_index=iterator.Second_Cell_Index();axis=iterator.Axis();
        if(!psi_N.Component(axis)(iterator.Face_Index()) && ((!psi.Valid_Index(first_cell_index) || psi(first_cell_index)) && (!psi.Valid_Index(second_cell_index) || psi(second_cell_index)))){
            T rho_first_cell=U_ghost(first_cell_index)(1),rho_second_cell=U_ghost(second_cell_index)(1);
            face_velocities.Component(axis)(iterator.Face_Index())=(rho_first_cell*EULER<T_GRID>::Get_Velocity_Component(U_ghost,first_cell_index,axis)+
                    rho_second_cell*EULER<T_GRID>::Get_Velocity_Component(U_ghost,second_cell_index,axis))/(rho_first_cell+rho_second_cell);}
        else face_velocities(iterator.Full_Index())=(T)0;}
}
//#####################################################################
// Function Compute_Face_Pressure_From_Cell_Pressures
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Compute_Face_Pressure_From_Cell_Pressures(const T_GRID& face_grid,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const T_ARRAYS_BOOL& psi,T_FACE_ARRAYS_SCALAR& p_face,const T_ARRAYS_SCALAR& p_cell)
{
    TV_INT first_cell_index,second_cell_index;int axis;
    for(FACE_ITERATOR iterator(face_grid);iterator.Valid();iterator.Next()){
        first_cell_index=iterator.First_Cell_Index();second_cell_index=iterator.Second_Cell_Index();axis=iterator.Axis();
        if((!psi.Valid_Index(first_cell_index) || psi(first_cell_index)) && (!psi.Valid_Index(second_cell_index) || psi(second_cell_index))){
            T rho_first_cell=U_ghost(first_cell_index)(1),rho_second_cell=U_ghost(second_cell_index)(1);
            p_face.Component(axis)(iterator.Face_Index())=(rho_second_cell*p_cell(first_cell_index)+rho_first_cell*p_cell(second_cell_index))/(rho_first_cell+rho_second_cell);}
        else if(!psi.Valid_Index(first_cell_index) || psi(first_cell_index))
            p_face(iterator.Full_Index()) = p_cell(first_cell_index);
        else if(!psi.Valid_Index(second_cell_index) || psi(second_cell_index))
            p_face(iterator.Full_Index()) = p_cell(second_cell_index);}
}
//#####################################################################
// Function Fill_Ghost_Pressures
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Get_Ghost_Pressures(const T dt,const T time,const T_ARRAYS_BOOL& psi_D,const T_FACE_ARRAYS_BOOL& psi_N,
    const T_ARRAYS_SCALAR& pressure,T_ARRAYS_SCALAR& p_ghost)
{
    assert(!is_pressure_scaled);
        
    p_ghost.Resize(euler->grid.Domain_Indices(1));
    pressure_boundary->Fill_Ghost_Cells(euler->grid,pressure,p_ghost,dt,time,1);
    // Restore dirichlet pressures at non-periodic dirichlet cells which are extrapolated by pressure_boundary
    if(!use_neumann_condition_for_outflow_boundaries){
        for(int axis=1;axis<=T_GRID::dimension;axis++) if(!elliptic_solver->periodic_boundary[axis]){
            for(int axis_side=1;axis_side<=2;axis_side++){int side=2*(axis-1)+axis_side;
                for(CELL_ITERATOR iterator(euler->grid,1,T_GRID::GHOST_REGION,side);iterator.Valid();iterator.Next()){
                    TV_INT cell_index=iterator.Cell_Index();
                    TV_INT boundary_face_index=side&1?iterator.Second_Face_Index(axis):iterator.First_Face_Index(axis);
                    if(psi_D(cell_index) && !psi_N.Component(axis)(boundary_face_index)){
                        p_ghost(cell_index)=p_dirichlet(cell_index);}}}}}
    if(euler->mpi_grid) euler->mpi_grid->Exchange_Boundary_Cell_Data(p_ghost,1,false);
}
//#####################################################################
// Function Get_Pressure_At_Faces
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Get_Pressure_At_Faces(const T dt,const T time,const T_ARRAYS_SCALAR& p_ghost,T_FACE_ARRAYS_SCALAR& p_face)
{
    euler->Fill_Ghost_Cells(dt,time,1);
    p_face.Resize(euler->grid);
    Compute_Face_Pressure_From_Cell_Pressures(euler->grid,euler->U_ghost,euler->psi,p_face,p_ghost);

    if(incompressible_coupling_callbacks){
        assert(!is_pressure_scaled);
        incompressible_coupling_callbacks->Compute_Compressible_Incompressible_Face_Pressures_From_Cell_Pressures(euler->grid,euler->U_ghost,euler->psi,p_ghost,p_face);}
}
//#####################################################################
// Function Apply_Pressure
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Apply_Pressure(const T_ARRAYS_SCALAR& p_ghost,const T_FACE_ARRAYS_SCALAR& p_face,const T_FACE_ARRAYS_SCALAR& face_velocities_star,
    const T_ARRAYS_BOOL& psi_D,const T_FACE_ARRAYS_BOOL& psi_N,const T dt,const T time)
{
    assert(!is_pressure_scaled);
    Apply_Pressure(p_ghost,p_face,face_velocities,elliptic_solver->psi_D,elliptic_solver->psi_N,dt,time,density_scaling,
            save_fluxes?fluxes:0,euler);
}
//#####################################################################
// Function Apply_Pressure
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Apply_Pressure(const T_ARRAYS_SCALAR& p_ghost,const T_FACE_ARRAYS_SCALAR& p_face,const T_FACE_ARRAYS_SCALAR& face_velocities_star,
    const T_ARRAYS_BOOL& psi_D,const T_FACE_ARRAYS_BOOL& psi_N,const T dt,const T time, const T_ARRAYS_SCALAR& density_scaling,
    T_FACE_ARRAYS_DIMENSION_SCALAR *fluxes,EULER_UNIFORM<T_GRID>* euler)
{
    T_ARRAYS_SCALAR p_hat(p_ghost);
    T_FACE_ARRAYS_SCALAR p_hat_face(p_face);
    p_hat*=dt;
    p_hat_face*=dt;

    euler->Fill_Ghost_Cells(dt,time,1);
    
    TV one_over_dx=euler->grid.one_over_dX;T one_over_dt=1/dt;
    // update momentum
    for(CELL_ITERATOR iterator(euler->grid);iterator.Valid();iterator.Next()) if(euler->psi(iterator.Cell_Index())){
        TV_INT cell_index=iterator.Cell_Index();
        for(int axis=1;axis<=T_GRID::dimension;axis++){
            TV_INT first_face_index=iterator.First_Face_Index(axis),second_face_index=iterator.Second_Face_Index(axis);
            T p_hat_first_face=p_hat_face.Component(axis)(first_face_index);
            T p_hat_second_face=p_hat_face.Component(axis)(second_face_index);
            T grad_p_hat=(p_hat_second_face-p_hat_first_face)*one_over_dx[axis];
            euler->U(cell_index)(axis+1)-=grad_p_hat/density_scaling(cell_index);
            if(fluxes){
                (*fluxes).Component(axis)(first_face_index)(axis+1)=p_hat_first_face*one_over_dt;
                (*fluxes).Component(axis)(second_face_index)(axis+1)=p_hat_second_face*one_over_dt;}}}

    T_FACE_ARRAYS_SCALAR face_velocities_np1(face_velocities_star);

    if(!euler->apply_cavitation_correction){
        // update face velocities for energy update
        T_FACE_ARRAYS_SCALAR grad_p_hat_face(euler->grid);
        ARRAYS_UTILITIES<T_GRID,T>::Compute_Gradient_At_Faces_From_Cell_Data(euler->grid,grad_p_hat_face,p_hat);
        for(FACE_ITERATOR iterator(euler->grid);iterator.Valid();iterator.Next()){int axis=iterator.Axis();TV_INT face_index=iterator.Face_Index();
            TV_INT first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
            if(!psi_N.Component(axis)(face_index) && (!euler->psi.Valid_Index(first_cell_index) || euler->psi(first_cell_index)) && (!euler->psi.Valid_Index(second_cell_index) || euler->psi(second_cell_index))){
                T rho_face=(euler->U_ghost(first_cell_index)(1)/density_scaling(first_cell_index)+euler->U_ghost(second_cell_index)(1)/density_scaling(second_cell_index))*(T).5;
                face_velocities_np1.Component(axis)(face_index)-=grad_p_hat_face.Component(axis)(face_index)/rho_face;}}}

   // update energy
   for(CELL_ITERATOR iterator(euler->grid);iterator.Valid();iterator.Next()) if(euler->psi(iterator.Cell_Index())){TV_INT cell_index=iterator.Cell_Index();
       T div_p_hat_u=0;
       for(int axis=1;axis<=T_GRID::dimension;axis++){
           TV_INT first_face_index=iterator.First_Face_Index(axis),second_face_index=iterator.Second_Face_Index(axis);
           T p_hat_first_face=p_hat_face.Component(axis)(first_face_index);
           T p_hat_second_face=p_hat_face.Component(axis)(second_face_index);
           div_p_hat_u+=(p_hat_second_face*face_velocities_np1.Component(axis)(second_face_index)-p_hat_first_face*face_velocities_np1.Component(axis)(first_face_index))*one_over_dx[axis];
           if(fluxes){
               (*fluxes).Component(axis)(first_face_index)(T_GRID::dimension+2)=p_hat_first_face*one_over_dt*face_velocities_np1.Component(axis)(first_face_index);
               (*fluxes).Component(axis)(second_face_index)(T_GRID::dimension+2)=p_hat_second_face*one_over_dt*face_velocities_np1.Component(axis)(second_face_index);}}
       euler->U(cell_index)(T_GRID::dimension+2)-=div_p_hat_u/density_scaling(cell_index);}

    //if(incompressible_coupling_callbacks) incompressible_coupling_callbacks->Apply_Pressure_At_Incompressible_Faces();

    euler->boundary->Apply_Boundary_Condition(euler->grid,euler->U,time+dt);
    euler->Invalidate_Ghost_Cells();
}
//#####################################################################
// Function Consistent_Boundary_Conditions
//#####################################################################
template<class T_GRID> bool EULER_PROJECTION_UNIFORM<T_GRID>::
Consistent_Boundary_Conditions() const
{
    if(euler->mpi_grid){
        std::stringstream ss;ss<<"checking for consistent mpi boundaries"<<std::endl;LOG::filecout(ss.str());
        T_ARRAYS_BOOL psi_D_ghost(elliptic_solver->psi_D);T_FACE_ARRAYS_SCALAR psi_N_ghost(euler->grid);
        euler->mpi_grid->Exchange_Boundary_Cell_Data(psi_D_ghost,1);
        for(int axis=1;axis<=T_GRID::dimension;axis++)for(int axis_side=1;axis_side<=2;axis_side++){int side=2*(axis-1)+axis_side;
            if(euler->mpi_grid->Neighbor(axis,axis_side)){TV_INT exterior_cell_offset=axis_side==1?-TV_INT::Axis_Vector(axis):TV_INT();
                for(FACE_ITERATOR iterator(euler->grid,0,T_GRID::BOUNDARY_REGION,side);iterator.Valid();iterator.Next()){
                    TV_INT face=iterator.Face_Index(),cell=face+exterior_cell_offset;int axis=iterator.Axis();
                    psi_N_ghost(axis,face)=(T)elliptic_solver->psi_N(axis,face);
                    assert(elliptic_solver->psi_D(cell)==psi_D_ghost(cell));}}}
        euler->mpi_grid->Assert_Common_Face_Data(psi_N_ghost);}
    return true;
}
//#####################################################################
// Function Log_Parameters
//#####################################################################
template<class T_GRID> void EULER_PROJECTION_UNIFORM<T_GRID>::
Log_Parameters() const
{
    LOG::SCOPE scope("EULER_PROJECTION_UNIFORM parameters");
    std::stringstream ss;ss<<"save_fluxes="<<save_fluxes<<std::endl;
    ss<<"use_exact_neumann_face_location="<<use_exact_neumann_face_location<<std::endl;
    ss<<"use_neumann_condition_for_outflow_boundaries="<<use_neumann_condition_for_outflow_boundaries<<std::endl;
    ss<<"hj_eno_order="<<hj_eno_order<<std::endl;LOG::filecout(ss.str());
}
//#####################################################################
template class EULER_PROJECTION_UNIFORM<GRID<VECTOR<float,1> > >;
template class EULER_PROJECTION_UNIFORM<GRID<VECTOR<float,2> > >;
template class EULER_PROJECTION_UNIFORM<GRID<VECTOR<float,3> > >;
template void EULER_PROJECTION_UNIFORM<GRID<VECTOR<float,1> > >::Compute_Divergence(const FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > &face_lookup);
template void EULER_PROJECTION_UNIFORM<GRID<VECTOR<float,2> > >::Compute_Divergence(const FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > &face_lookup);
template void EULER_PROJECTION_UNIFORM<GRID<VECTOR<float,3> > >::Compute_Divergence(const FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > &face_lookup);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EULER_PROJECTION_UNIFORM<GRID<VECTOR<double,1> > >;
template class EULER_PROJECTION_UNIFORM<GRID<VECTOR<double,2> > >;
template class EULER_PROJECTION_UNIFORM<GRID<VECTOR<double,3> > >;
template void EULER_PROJECTION_UNIFORM<GRID<VECTOR<double,1> > >::Compute_Divergence(const FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > &face_lookup);
template void EULER_PROJECTION_UNIFORM<GRID<VECTOR<double,2> > >::Compute_Divergence(const FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > &face_lookup);
template void EULER_PROJECTION_UNIFORM<GRID<VECTOR<double,3> > >::Compute_Divergence(const FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > &face_lookup);
#endif
