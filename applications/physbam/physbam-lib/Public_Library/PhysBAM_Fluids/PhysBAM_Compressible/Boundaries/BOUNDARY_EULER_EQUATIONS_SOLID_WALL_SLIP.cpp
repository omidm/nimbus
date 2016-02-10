//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Boundaries/BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>::
BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP(EULER_UNIFORM<T_GRID>* euler_input,const T_FACE_VECTOR rho_far_field,
    const T_FACE_VECTOR p_far_field,const TV_FACE_VECTOR velocity_far_field,const T inflow_attenuation_input,
    const TV_SIDES& constant_extrapolation,const bool always_attenuate_input,const T_FACE_VECTOR linear_attenuations_input,
    const T_FACE_VECTOR_BOOL linear_attenuation_faces_input)
{
    euler=euler_input;
    Set_Constant_Extrapolation(constant_extrapolation);
    inflow_attenuation=inflow_attenuation_input;
    always_attenuate=always_attenuate_input;
    linear_attenuations=linear_attenuations_input;
    linear_attenuation_faces=linear_attenuation_faces_input;

    T gamma=dynamic_cast<EOS_GAMMA<T>*>(euler->eos)->gamma;
    for(int side=1;side<=2*T_GRID::dimension;side++){
        T e_far_field=euler->eos->e_From_p_And_rho(p_far_field(side),rho_far_field(side));
        T c_far_field=euler->eos->c(rho_far_field(side),e_far_field);
        S_far_field(side)=euler->eos->S(rho_far_field(side),e_far_field);
        iL_far_field(side)=-velocity_far_field(side)[1]+2*c_far_field/(gamma-1);
        iR_far_field(side)=velocity_far_field(side)[1]+2*c_far_field/(gamma-1);    
        U_far_field(side)=EULER<T_GRID>::Get_Euler_State_From_rho_velocity_And_internal_energy(rho_far_field(side),velocity_far_field(side),e_far_field);}
}
//#####################################################################
// Function Attenuate_To_Far_Field_Values_Using_Riemann_Invariants
//#####################################################################
template<class T_GRID> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>::
Attenuate_To_Far_Field_Values_Using_Riemann_Invariants(const T_ARRAYS_DIMENSION_BASE& u_ghost,const TV_INT& node_index,const int side,TV_DIMENSION &U,const T dt) const
{
    // TODO: Implement for multi-dimension.
    T gamma=dynamic_cast<EOS_GAMMA<T>*>(euler->eos)->gamma;
    T net_inflow_attenuation=exp((inflow_attenuation-1)*dt);

    T rho=u_ghost(node_index)(1);
    T u_velocity=EULER<T_GRID>::Get_Velocity_Component(u_ghost,node_index,1);
    T e=EULER<T_GRID>::e(u_ghost,node_index);

    T c=euler->eos->c(rho,e);
    T S=euler->eos->S(rho,e);
    T iL=-u_velocity+2*c/(gamma-1);
    T iR=u_velocity+2*c/(gamma-1);

    int flip=side&1?1:-1;
    if(flip*u_velocity > 0) S=S_far_field(side)+net_inflow_attenuation*(S-S_far_field(side));
    if(flip*(u_velocity-c) > 0) iL=iL_far_field(side)+net_inflow_attenuation*(iL-iL_far_field(side));
    if(flip*(u_velocity+c) > 0) iR=iR_far_field(side)+net_inflow_attenuation*(iR-iR_far_field(side));

    u_velocity=(iR-iL)/2;
    c=(iR-u_velocity)*(gamma-1)*(T).5;
    rho=pow(sqr(c)/(gamma*S),1/(gamma-1));
    e=euler->eos->e_From_S_And_rho(S,rho);

    U[1]=rho;
    U[2]=rho*u_velocity;
    U[3]=rho*(e+sqr(u_velocity)/2);
}
//#####################################################################
// Function Attenuate_To_Far_Field_Values_Using_Characteristics
//#####################################################################
template<class T_GRID> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>::
Attenuate_To_Far_Field_Values_Using_Characteristics(const T_ARRAYS_DIMENSION_BASE& u_ghost,const TV_INT& node_index,const int side,TV_DIMENSION &U,const T dt) const
{
    T net_inflow_attenuation=exp((inflow_attenuation-1)*dt);
    if(linear_attenuation_faces[side]){net_inflow_attenuation=(1-linear_attenuations(side));}

    if(always_attenuate){
        U=u_ghost(node_index);
        U=U_far_field(side)+net_inflow_attenuation*(U-U_far_field(side));
        return;}
   
    int axis=(side+1)/2;
    MATRIX<T,d,d> L,R;
    TV_DIMENSION LU,LU_far_field;
    ARRAY<T,VECTOR<int,1> > lambda(1,d),lambda_left(1,d),lambda_right(1,d);

    // extract a 1d array
    int min_index=u_ghost.Domain_Indices().min_corner[axis],max_index=u_ghost.Domain_Indices().max_corner[axis];
    TV_INT current_index=node_index;
    ARRAY<TV_DIMENSION,VECTOR<int,1> > U_1d(min_index,max_index);
    for(int i=min_index;i<=max_index;i++){current_index[axis]=i;U_1d(i)=u_ghost(current_index);}

    // apply L
    euler->eigensystems_default[axis]->Eigenvectors(U_1d,node_index[axis],L,R);
    euler->eigensystems_default[axis]->Eigenvalues(U_1d,node_index[axis],lambda,lambda_left,lambda_right);
    for(int k=1;k<=d;k++){LU(k)=0;LU_far_field(k)=0;
        for(int kk=1;kk<=d;kk++){LU(k)+=L(k,kk)*u_ghost(node_index)(kk);LU_far_field(k)+=L(k,kk)*U_far_field(side)(kk);}}

    // attenuate
    int flip=side&1?1:-1;
    for(int k=1;k<=d;k++) if(flip*lambda_left(k)>0) LU(k)=LU_far_field(k)+net_inflow_attenuation*(LU(k)-LU_far_field(k));

    // apply R
    for(int k=1;k<=d;k++){U(k)=0;for(int kk=1;kk<=d;kk++) U(k)+=LU(kk)*R(kk,k);}
}
//#####################################################################
// Function Attenuate_To_Far_Field_Values
//#####################################################################
template<class T_GRID> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>::
Attenuate_To_Far_Field_Values(const T_ARRAYS_DIMENSION_BASE& u_ghost,const TV_INT& node_index,const int side,TV_DIMENSION &U,const T dt) const
{
    U=u_ghost(node_index);
    T net_inflow_attenuation=exp((inflow_attenuation-1)*dt);
    U=U_far_field(side)+net_inflow_attenuation*(U-U_far_field(side));
}
//#####################################################################
// Function Fill_Single_Ghost_Region
//#####################################################################
template<class T_GRID> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>::
Fill_Single_Ghost_Region(const T_GRID& grid,T_ARRAYS_DIMENSION_BASE& u_ghost,const RANGE<TV_INT>& region,const int side,const T dt,const T time,const int number_of_ghost_cells) const
{
    int total_number_of_ghost_cells=u_ghost.Domain_Indices().max_corner.x-grid.Domain_Indices().max_corner.x;
    int ghost_cells_to_fill=number_of_ghost_cells;
    bool attenuate_dirichlet_cells=true;
    int axis=(side+1)/2;
    if(Constant_Extrapolation(side)){
        RANGE<TV_INT> region_new(region);
        if(!attenuate_dirichlet_cells){
            int already_filled_cells=total_number_of_ghost_cells-ghost_cells_to_fill;
            side&1?region_new.max_corner(axis)-=already_filled_cells:region_new.min_corner(axis)+=already_filled_cells;}
        int boundary=Boundary(side,region_new);
        if(attenuate_dirichlet_cells){
            Fill_Single_Ghost_Region(grid,u_ghost,side,region_new); // extrapolate values first, so that eigenvectors calculated at boundary flux values are calculated correctly.
            for(CELL_ITERATOR iterator(grid,region_new);iterator.Valid();iterator.Next()){
                TV_INT cell=iterator.Cell_Index(),boundary_node=cell;boundary_node[axis]=boundary;TV_DIMENSION U_boundary;
                Attenuate_To_Far_Field_Values_Using_Characteristics(u_ghost,boundary_node,side,U_boundary,dt);
                u_ghost(cell)=U_boundary;}}
        else{
            for(CELL_ITERATOR iterator(grid,region_new);iterator.Valid();iterator.Next()){
                TV_INT cell=iterator.Cell_Index(),boundary_node=cell;boundary_node[axis]=boundary;
                u_ghost(cell)=u_ghost(boundary_node);}}}
    else{
        int boundary=Boundary(side,region);
        int reflection_times_two;
        if(grid.Is_MAC_Grid())
            reflection_times_two=2*boundary+(side&1?-1:1);
        else reflection_times_two=2*boundary;
        for(CELL_ITERATOR iterator(grid,region);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            TV_INT reflected_node=cell;reflected_node[axis]=reflection_times_two-cell[axis];
            T rho=u_ghost(reflected_node)(1);
            TV velocity=EULER<T_GRID>::Get_Velocity(u_ghost,reflected_node);velocity(axis)*=-1;
            T e=EULER<T_GRID>::e(u_ghost,reflected_node);
            EULER<T_GRID>::Set_Euler_State_From_rho_velocity_And_internal_energy(u_ghost,cell,rho,velocity,e);}}
}
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T_GRID> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>::
Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_DIMENSION_BASE& u,T_ARRAYS_DIMENSION_BASE& u_ghost,const T dt,const T time,const int number_of_ghost_cells)
{
    T_ARRAYS_DIMENSION_BASE::Put(u,u_ghost); // interior
    ARRAY<RANGE<TV_INT> > regions;Find_Ghost_Regions(grid,regions,number_of_ghost_cells);
    for(int side=1;side<=2*T_GRID::dimension;side++){
        Fill_Single_Ghost_Region(grid,u_ghost,regions(side),side,dt,time,number_of_ghost_cells);}
}
//#####################################################################
// Function Apply_Boundary_Condition_Single_Side
//#####################################################################
template<class T_GRID> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>::
Apply_Boundary_Condition_Single_Side(const T_GRID& grid,T_ARRAYS_DIMENSION_BASE& u,const int side,const T time) const 
{
    if(grid.Is_MAC_Grid()) return;
    int axis=(side+1)/2;
    int axis_side=2-side%2;
    if(!euler->mpi_grid||!euler->mpi_grid->Neighbor(axis,axis_side)){
        if(!Constant_Extrapolation(side)) for(CELL_ITERATOR iterator(grid,0,T_GRID::BOUNDARY_INTERIOR_REGION,side);iterator.Valid();iterator.Next()){
            TV_INT boundary_node=iterator.Cell_Index();
            T rho=u(boundary_node)(1);
            TV velocity=EULER<T_GRID>::Get_Velocity(u,boundary_node);
            T e=EULER<T_GRID>::e(u,boundary_node);
            velocity[axis]=0; // Boundary condition
            EULER<T_GRID>::Set_Euler_State_From_rho_velocity_And_internal_energy(u,boundary_node,rho,velocity,e);}}
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class T_GRID> void BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<T_GRID>::
Apply_Boundary_Condition(const T_GRID& grid,T_ARRAYS_DIMENSION_BASE& u,const T time)
{
    if(grid.Is_MAC_Grid()) return;
    for(int side=1;side<=2*T_GRID::dimension;side++) Apply_Boundary_Condition_Single_Side(grid,u,side,time);
}
//#####################################################################
template class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<GRID<VECTOR<float,1> > >;
template class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<GRID<VECTOR<float,2> > >;
template class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<GRID<VECTOR<double,1> > >;
template class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<GRID<VECTOR<double,2> > >;
template class BOUNDARY_EULER_EQUATIONS_SOLID_WALL_SLIP<GRID<VECTOR<double,3> > >;
#endif
