//#####################################################################
// Copyright 2002-2009, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Avi Robinson-Mosher, Tamar Shinar, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEPARABLE_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/AVERAGING_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Dynamics/Level_Sets/FAST_LEVELSET_ADVECTION.h>
using namespace PhysBAM;
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T_GRID> void FAST_LEVELSET_ADVECTION<T_GRID>::
Euler_Step(const T_ARRAYS_VECTOR& V,const T dt,const T time,const int number_of_ghost_cells)
{
    T_GRID& grid=((T_FAST_LEVELSET*)levelset)->grid;
    T_ARRAYS_SCALAR& phi=((T_FAST_LEVELSET*)levelset)->phi;
    T_ARRAYS_SCALAR phi_ghost(grid.Domain_Indices(number_of_ghost_cells));((T_FAST_LEVELSET*)levelset)->boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,number_of_ghost_cells);

    if(local_semi_lagrangian_advection){
        LINEAR_INTERPOLATION_UNIFORM<T_GRID,T> interpolation;
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();if(abs(phi_ghost(cell)) <= ((T_FAST_LEVELSET*)levelset)->half_band_width)
            phi(cell)=interpolation.Clamped_To_Array(grid,phi_ghost,iterator.Location()-dt*V(cell));}}
    else if(local_advection_spatial_order){
        T_ARRAYS_SCALAR rhs(grid.Domain_Indices());
        Euler_Step_High_Order_Helper(grid,V,phi,phi_ghost,rhs,local_advection_spatial_order,((T_FAST_LEVELSET*)levelset)->half_band_width);
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();if(abs(phi(cell)) <= ((T_FAST_LEVELSET*)levelset)->half_band_width) phi(cell)-=dt*rhs(cell);}}
    else // use the advection routine in the level set base class
        advection->Update_Advection_Equation_Node(grid,phi,phi_ghost,V,*((T_FAST_LEVELSET*)levelset)->boundary,dt,time);

    ((T_FAST_LEVELSET*)levelset)->boundary->Apply_Boundary_Condition(grid,phi,time+dt);
}
//#####################################################################
// Function Euler_Step
//#####################################################################
template<class T_GRID> void FAST_LEVELSET_ADVECTION<T_GRID>::
Euler_Step(const T_FACE_ARRAYS_SCALAR& V,const T dt,const T time,const int number_of_ghost_cells)
{
    T_GRID& grid=((T_FAST_LEVELSET*)levelset)->grid;
    T_ARRAYS_SCALAR& phi=((T_FAST_LEVELSET*)levelset)->phi;
    
    assert(grid.Is_MAC_Grid() && advection); // for now use advection in base class
    T_ARRAYS_SCALAR phi_ghost(grid.Domain_Indices(number_of_ghost_cells));((T_FAST_LEVELSET*)levelset)->boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,number_of_ghost_cells);
    advection->Update_Advection_Equation_Cell(grid,phi,phi_ghost,V,*((T_FAST_LEVELSET*)levelset)->boundary,dt,time);
    ((T_FAST_LEVELSET*)levelset)->boundary->Apply_Boundary_Condition(grid,phi,time+dt);
}
//#####################################################################
// Function Euler_Step_High_Order_Helper
//#####################################################################
template<class T> static void
Euler_Step_High_Order_Helper(const GRID<VECTOR<T,1> >& grid,const ARRAY<VECTOR<T,1> ,VECTOR<int,1> >& V,const ARRAY<T,VECTOR<int,1> >& phi,const ARRAY<T,VECTOR<int,1> >& phi_ghost,ARRAY<T,VECTOR<int,1> >& rhs,const int spatial_order,
    const T half_band_width)
{
    int m=grid.counts.x;T dx=grid.dX.x;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,1> > phi_1d_x(1-ghost_cells,m+ghost_cells),u_1d(1,m),distance_1d_x(1,m);
    for(int i=1-ghost_cells;i<=m+ghost_cells;i++) phi_1d_x(i)=phi_ghost(i);
    for(int i=1;i<=m;i++){u_1d(i)=V(i).x;distance_1d_x(i)=phi(i);}
    if(spatial_order == 5) Local_WENO_Advect(m,dx,phi_1d_x,u_1d,distance_1d_x,rhs,half_band_width);
    else Local_ENO_Advect(spatial_order,m,dx,phi_1d_x,u_1d,distance_1d_x,rhs,half_band_width);
}
template<class T> static void
Euler_Step_High_Order_Helper(const GRID<VECTOR<T,2> >& grid,const ARRAY<VECTOR<T,2> ,VECTOR<int,2> >& V,const ARRAY<T,VECTOR<int,2> >& phi,const ARRAY<T,VECTOR<int,2> >& phi_ghost,ARRAY<T,VECTOR<int,2> >& rhs,const int spatial_order,
    const T half_band_width)
{
    int m=grid.counts.x,n=grid.counts.y;T dx=grid.dX.x,dy=grid.dX.y;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,1> > phi_1d_x(1-ghost_cells,m+ghost_cells),u_1d(1,m),distance_1d_x(1,m),u_phix_1d(1,m);
    for(int j=1;j<=n;j++){
        for(int i=1-ghost_cells;i<=m+ghost_cells;i++) phi_1d_x(i)=phi_ghost(i,j);
        for(int i=1;i<=m;i++){u_1d(i)=V(i,j).x;distance_1d_x(i)=phi(i,j);}
        if(spatial_order == 5) Local_WENO_Advect(m,dx,phi_1d_x,u_1d,distance_1d_x,u_phix_1d,half_band_width);
        else Local_ENO_Advect(spatial_order,m,dx,phi_1d_x,u_1d,distance_1d_x,u_phix_1d,half_band_width);
        for(int i=1;i<=m;i++) rhs(i,j)=u_phix_1d(i);}
    ARRAY<T,VECTOR<int,1> > phi_1d_y(1-ghost_cells,n+ghost_cells),v_1d(1,n),distance_1d_y(1,n),v_phiy_1d(1,n);
    for(int i=1;i<=m;i++){
        for(int j=1-ghost_cells;j<=n+ghost_cells;j++) phi_1d_y(j)=phi_ghost(i,j);
        for(int j=1;j<=n;j++){v_1d(j)=V(i,j).y;distance_1d_y(j)=phi(i,j);}
        if(spatial_order == 5) Local_WENO_Advect(n,dy,phi_1d_y,v_1d,distance_1d_y,v_phiy_1d,half_band_width);
        else Local_ENO_Advect(spatial_order,n,dy,phi_1d_y,v_1d,distance_1d_y,v_phiy_1d,half_band_width);
        for(int j=1;j<=n;j++) rhs(i,j)+=v_phiy_1d(j);}
}
template<class T> static void
Euler_Step_High_Order_Helper(const GRID<VECTOR<T,3> >& grid,const ARRAY<VECTOR<T,3> ,VECTOR<int,3> >& V,const ARRAY<T,VECTOR<int,3> >& phi,const ARRAY<T,VECTOR<int,3> >& phi_ghost,ARRAY<T,VECTOR<int,3> >& rhs,const int spatial_order,
    const T half_band_width)
{
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;T dx=grid.dX.x,dy=grid.dX.y,dz=grid.dX.z;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,1> > phi_1d_x(1-ghost_cells,m+ghost_cells),u_1d(1,m),distance_1d_x(1,m),u_phix_1d(1,m);
    for(int j=1;j<=n;j++) for(int ij=1;ij<=mn;ij++){
        for(int i=1-ghost_cells;i<=m+ghost_cells;i++) phi_1d_x(i)=phi_ghost(i,j,ij);
        for(int i=1;i<=m;i++){u_1d(i)=V(i,j,ij).x;distance_1d_x(i)=phi(i,j,ij);}
        if(spatial_order == 5) Local_WENO_Advect(m,dx,phi_1d_x,u_1d,distance_1d_x,u_phix_1d,half_band_width);
        else Local_ENO_Advect(spatial_order,m,dx,phi_1d_x,u_1d,distance_1d_x,u_phix_1d,half_band_width);
        for(int i=1;i<=m;i++) rhs(i,j,ij)=u_phix_1d(i);}
    ARRAY<T,VECTOR<int,1> > phi_1d_y(1-ghost_cells,n+ghost_cells),v_1d(1,n),distance_1d_y(1,n),v_phiy_1d(1,n);
    for(int i=1;i<=m;i++) for(int ij=1;ij<=mn;ij++){
        for(int j=1-ghost_cells;j<=n+ghost_cells;j++) phi_1d_y(j)=phi_ghost(i,j,ij);
        for(int j=1;j<=n;j++){v_1d(j)=V(i,j,ij).y;distance_1d_y(j)=phi(i,j,ij);}
        if(spatial_order == 5) Local_WENO_Advect(n,dy,phi_1d_y,v_1d,distance_1d_y,v_phiy_1d,half_band_width);
        else Local_ENO_Advect(spatial_order,n,dy,phi_1d_y,v_1d,distance_1d_y,v_phiy_1d,half_band_width);
        for(int j=1;j<=n;j++) rhs(i,j,ij)+=v_phiy_1d(j);}
    ARRAY<T,VECTOR<int,1> > phi_1d_z(1-ghost_cells,mn+ghost_cells),w_1d(1,mn),distance_1d_z(1,mn),w_phiz_1d(1,mn);
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){
        for(int ij=1-ghost_cells;ij<=mn+ghost_cells;ij++) phi_1d_z(ij)=phi_ghost(i,j,ij);
        for(int ij=1;ij<=mn;ij++){w_1d(ij)=V(i,j,ij).z;distance_1d_z(ij)=phi(i,j,ij);}
        if(spatial_order == 5) Local_WENO_Advect(mn,dz,phi_1d_z,w_1d,distance_1d_z,w_phiz_1d,half_band_width);
        else Local_ENO_Advect(spatial_order,mn,dz,phi_1d_z,w_1d,distance_1d_z,w_phiz_1d,half_band_width);
        for(int ij=1;ij<=mn;ij++) rhs(i,j,ij)+=w_phiz_1d(ij);}
}
//#####################################################################
// Functions Reinitialize
//#####################################################################
template<class T_GRID> void FAST_LEVELSET_ADVECTION<T_GRID>::
Reinitialize(const int time_steps,const T time)
{
    T_GRID& grid=((T_FAST_LEVELSET*)levelset)->grid;
    T_ARRAYS_SCALAR& phi=((T_FAST_LEVELSET*)levelset)->phi;
    
    T large_band=((T_FAST_LEVELSET*)levelset)->half_band_width+grid.dX.Max()*(1+min(3,local_advection_spatial_order));
    T_ARRAYS_SCALAR signed_distance(grid.Domain_Indices());
    ((T_FAST_LEVELSET*)levelset)->Get_Signed_Distance_Using_FMM(signed_distance,time,large_band);

    T_ARRAYS_SCALAR sign_phi(grid.Domain_Indices()); // smeared out sign function
    T epsilon=sqr(grid.dX.Max());
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        sign_phi(cell)=phi(cell)/sqrt(sqr(phi(cell))+epsilon);}

    T dt=reinitialization_cfl*grid.dX.Min();
    RUNGEKUTTA<T_ARRAYS_SCALAR> rungekutta(phi);
    rungekutta.Set_Grid_And_Boundary_Condition(grid,*((T_FAST_LEVELSET*)levelset)->boundary);
    rungekutta.Set_Order(reinitialization_runge_kutta_order);
    rungekutta.Set_Time(time);
    rungekutta.Pseudo_Time();
    for(int k=1;k<=time_steps;k++){
        rungekutta.Start(dt);
        for(int kk=1;kk<=rungekutta.order;kk++){Euler_Step_Of_Reinitialization(signed_distance,sign_phi,dt,time);rungekutta.Main();}
        }

    T min_DX=grid.dX.Min();
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();if(abs(signed_distance(cell)) <= large_band){
        if(abs(signed_distance(cell)) > ((T_FAST_LEVELSET*)levelset)->half_band_width) phi(cell)=signed_distance(cell); // outer band - use the FMM solution
        else if(abs(signed_distance(cell)-phi(cell)) > min_DX) phi(cell)=signed_distance(cell);}} // inner band - use FMM if errors look big
}
//#####################################################################
// Functions Euler_Step_Of_Reinitialization
//#####################################################################
template<class T_GRID> void FAST_LEVELSET_ADVECTION<T_GRID>::
Euler_Step_Of_Reinitialization(const T_ARRAYS_SCALAR& signed_distance,const T_ARRAYS_SCALAR& sign_phi,const T dt,const T time)
{
    T_GRID& grid=((T_FAST_LEVELSET*)levelset)->grid;
    T_ARRAYS_SCALAR& phi=((T_FAST_LEVELSET*)levelset)->phi;
    
    int ghost_cells=3;
    T_ARRAYS_SCALAR phi_ghost(grid.Domain_Indices(ghost_cells));((T_FAST_LEVELSET*)levelset)->boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,dt,time,ghost_cells);
    T_ARRAYS_SCALAR rhs(grid.Domain_Indices());

    Euler_Step_Of_Reinitialization_High_Order_Helper(grid,signed_distance,phi,phi_ghost,rhs,reinitialization_spatial_order,((T_FAST_LEVELSET*)levelset)->half_band_width);

    T min_DX=grid.dX.Min();
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();if(abs(signed_distance(cell)) <= ((T_FAST_LEVELSET*)levelset)->half_band_width){
        phi(cell)-=dt*sign_phi(cell)*(sqrt(rhs(cell))-1);
        if(LEVELSET_UTILITIES<T>::Interface(phi_ghost(cell),phi(cell))) phi(cell)=LEVELSET_UTILITIES<T>::Sign(phi_ghost(cell))*((T_FAST_LEVELSET*)levelset)->small_number*min_DX;}}

    ((T_FAST_LEVELSET*)levelset)->boundary->Apply_Boundary_Condition(grid,phi,time); // time not incremented - pseudo-time
}
//#####################################################################
// Functions Euler_Step_Of_Reinitialization_High_Order_Helper
//#####################################################################
template<class T> static void
Euler_Step_Of_Reinitialization_High_Order_Helper(const GRID<VECTOR<T,1> >& grid,const ARRAY<T,VECTOR<int,1> >& signed_distance,const ARRAY<T,VECTOR<int,1> >& phi,const ARRAY<T,VECTOR<int,1> >& phi_ghost,ARRAY<T,VECTOR<int,1> >& rhs,
    const int spatial_order,const T half_band_width)
{
    int m=grid.counts.x;T dx=grid.dX.x;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,1> > phi_1d_x(1-ghost_cells,m+3),distance_1d_x(1,m),phix_minus(1,m),phix_plus(1,m);
    for(int i=1-ghost_cells;i<=m+3;i++) phi_1d_x(i)=phi_ghost(i);
    for(int i=1;i<=m;i++) distance_1d_x(i)=signed_distance(i);
    if(spatial_order == 5) Local_WENO_Reinitialize(m,dx,phi_1d_x,distance_1d_x,phix_minus,phix_plus,half_band_width);
    else Local_ENO_Reinitialize(spatial_order,m,dx,phi_1d_x,distance_1d_x,phix_minus,phix_plus,half_band_width);
    for(int i=1;i<=m;i++) if(abs(signed_distance(i)) <= half_band_width){
        if(LEVELSET_UTILITIES<T>::Sign(phi(i)) < 0) rhs(i)=sqr(max(-phix_minus(i),phix_plus(i),(T)0));
        else rhs(i)=sqr(max(phix_minus(i),-phix_plus(i),(T)0));}
}
template<class T> static void
Euler_Step_Of_Reinitialization_High_Order_Helper(const GRID<VECTOR<T,2> >& grid,const ARRAY<T,VECTOR<int,2> >& signed_distance,const ARRAY<T,VECTOR<int,2> >& phi,const ARRAY<T,VECTOR<int,2> >& phi_ghost,ARRAY<T,VECTOR<int,2> >& rhs,
    const int spatial_order,const T half_band_width)
{
    int m=grid.counts.x,n=grid.counts.y;T dx=grid.dX.x,dy=grid.dX.y;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,1> > phi_1d_x(1-ghost_cells,m+3),distance_1d_x(1,m),phix_minus(1,m),phix_plus(1,m);
    for(int j=1;j<=n;j++){
        for(int i=1-ghost_cells;i<=m+3;i++) phi_1d_x(i)=phi_ghost(i,j);
        for(int i=1;i<=m;i++) distance_1d_x(i)=signed_distance(i,j);
        if(spatial_order == 5) Local_WENO_Reinitialize(m,dx,phi_1d_x,distance_1d_x,phix_minus,phix_plus,half_band_width);
        else Local_ENO_Reinitialize(spatial_order,m,dx,phi_1d_x,distance_1d_x,phix_minus,phix_plus,half_band_width);
        for(int i=1;i<=m;i++) if(abs(signed_distance(i,j)) <= half_band_width){
            if(LEVELSET_UTILITIES<T>::Sign(phi(i,j)) < 0) rhs(i,j)=sqr(max(-phix_minus(i),phix_plus(i),(T)0));
            else rhs(i,j)=sqr(max(phix_minus(i),-phix_plus(i),(T)0));}}

    ARRAY<T,VECTOR<int,1> > phi_1d_y(1-ghost_cells,n+3),distance_1d_y(1,n),phiy_minus(1,n),phiy_plus(1,n);
    for(int i=1;i<=m;i++){
        for(int j=1-ghost_cells;j<=n+3;j++) phi_1d_y(j)=phi_ghost(i,j);
        for(int j=1;j<=n;j++) distance_1d_y(j)=signed_distance(i,j);
        if(spatial_order == 5) Local_WENO_Reinitialize(n,dy,phi_1d_y,distance_1d_y,phiy_minus,phiy_plus,half_band_width);
        else Local_ENO_Reinitialize(spatial_order,n,dy,phi_1d_y,distance_1d_y,phiy_minus,phiy_plus,half_band_width);
        for(int j=1;j<=n;j++) if(abs(signed_distance(i,j)) <= half_band_width){
            if(LEVELSET_UTILITIES<T>::Sign(phi(i,j)) < 0) rhs(i,j)+=sqr(max(-phiy_minus(j),phiy_plus(j),(T)0));
            else rhs(i,j)+=sqr(max(phiy_minus(j),-phiy_plus(j),(T)0));}}
}
template<class T> static void
Euler_Step_Of_Reinitialization_High_Order_Helper(const GRID<VECTOR<T,3> >& grid,const ARRAY<T,VECTOR<int,3> >& signed_distance,const ARRAY<T,VECTOR<int,3> >& phi,const ARRAY<T,VECTOR<int,3> >& phi_ghost,ARRAY<T,VECTOR<int,3> >& rhs,
    const int spatial_order,const T half_band_width)
{
    int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;T dx=grid.dX.x,dy=grid.dX.y,dz=grid.dX.z;
    int ghost_cells=3;
    ARRAY<T,VECTOR<int,1> > phi_1d_x(1-ghost_cells,m+3),distance_1d_x(1,m),phix_minus(1,m),phix_plus(1,m);
    for(int j=1;j<=n;j++) for(int ij=1;ij<=mn;ij++){
        for(int i=1-ghost_cells;i<=m+3;i++) phi_1d_x(i)=phi_ghost(i,j,ij);
        for(int i=1;i<=m;i++) distance_1d_x(i)=signed_distance(i,j,ij);
        if(spatial_order == 5) Local_WENO_Reinitialize(m,dx,phi_1d_x,distance_1d_x,phix_minus,phix_plus,half_band_width);
        else Local_ENO_Reinitialize(spatial_order,m,dx,phi_1d_x,distance_1d_x,phix_minus,phix_plus,half_band_width);
        for(int i=1;i<=m;i++) if(abs(signed_distance(i,j,ij)) <= half_band_width){
            if(LEVELSET_UTILITIES<T>::Sign(phi(i,j,ij)) < 0) rhs(i,j,ij)=sqr(max(-phix_minus(i),phix_plus(i),(T)0));
            else rhs(i,j,ij)=sqr(max(phix_minus(i),-phix_plus(i),(T)0));}}

    ARRAY<T,VECTOR<int,1> > phi_1d_y(1-ghost_cells,n+3),distance_1d_y(1,n),phiy_minus(1,n),phiy_plus(1,n);
    for(int i=1;i<=m;i++) for(int ij=1;ij<=mn;ij++){
        for(int j=1-ghost_cells;j<=n+3;j++) phi_1d_y(j)=phi_ghost(i,j,ij);
        for(int j=1;j<=n;j++) distance_1d_y(j)=signed_distance(i,j,ij);
        if(spatial_order == 5) Local_WENO_Reinitialize(n,dy,phi_1d_y,distance_1d_y,phiy_minus,phiy_plus,half_band_width);
        else Local_ENO_Reinitialize(spatial_order,n,dy,phi_1d_y,distance_1d_y,phiy_minus,phiy_plus,half_band_width);
        for(int j=1;j<=n;j++) if(abs(signed_distance(i,j,ij)) <= half_band_width){
            if(LEVELSET_UTILITIES<T>::Sign(phi(i,j,ij)) < 0) rhs(i,j,ij)+=sqr(max(-phiy_minus(j),phiy_plus(j),(T)0));
            else rhs(i,j,ij)+=sqr(max(phiy_minus(j),-phiy_plus(j),(T)0));}}

    ARRAY<T,VECTOR<int,1> > phi_1d_z(1-ghost_cells,mn+3),distance_1d_z(1,mn),phiz_minus(1,mn),phiz_plus(1,mn);
    for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){
        for(int ij=1-ghost_cells;ij<=mn+3;ij++) phi_1d_z(ij)=phi_ghost(i,j,ij);
        for(int ij=1;ij<=mn;ij++) distance_1d_z(ij)=signed_distance(i,j,ij);
        if(spatial_order == 5) Local_WENO_Reinitialize(mn,dz,phi_1d_z,distance_1d_z,phiz_minus,phiz_plus,half_band_width);
        else Local_ENO_Reinitialize(spatial_order,mn,dz,phi_1d_z,distance_1d_z,phiz_minus,phiz_plus,half_band_width);
        for(int ij=1;ij<=mn;ij++) if(abs(signed_distance(i,j,ij)) <= half_band_width){
            if(LEVELSET_UTILITIES<T>::Sign(phi(i,j,ij)) < 0) rhs(i,j,ij)+=sqr(max(-phiz_minus(ij),phiz_plus(ij),(T)0));
            else rhs(i,j,ij)+=sqr(max(phiz_minus(ij),-phiz_plus(ij),(T)0));}}
}
//#####################################################################
// Function Local_WENO_Advect
//#####################################################################
// phi is (-2,m+3), u, distance and u_phix are (1,m)
template<class T> static void
Local_WENO_Advect(const int m,const T dx,const ARRAY<T,VECTOR<int,1> >& phi,const ARRAY<T,VECTOR<int,1> >& u,const ARRAY<T,VECTOR<int,1> >& distance,ARRAY<T,VECTOR<int,1> >& u_phix,const T half_band_width)
{
    T epsilon=(T)1e-6*sqr(dx); // 1e-6 works since phi is a distance function - sqr(dx) since undivided differences are used
    T one_over_dx=1/dx;
    for(int i=1;i<=m;i++) if(abs(distance(i)) <= half_band_width){ // one_over_dx since undivided differences are used
        if(u(i) > 0) u_phix(i)=u(i)*ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<T,1> >,T>::WENO(phi(i-2)-phi(i-3),phi(i-1)-phi(i-2),phi(i)-phi(i-1),phi(i+1)-phi(i),phi(i+2)-phi(i+1),epsilon)*one_over_dx;
        else u_phix(i)=u(i)*ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<T,1> >,T>::WENO(phi(i+3)-phi(i+2),phi(i+2)-phi(i+1),phi(i+1)-phi(i),phi(i)-phi(i-1),phi(i-1)-phi(i-2),epsilon)*one_over_dx;}
}
//#####################################################################
// Function Local_ENO_Advect
//#####################################################################
// order = 1, 2 or 3, phi is (-2,m_3), u, distance and u_phix are (1,m)
template<class T> static void
Local_ENO_Advect(const int order,const int m,const T dx,const ARRAY<T,VECTOR<int,1> >& phi,const ARRAY<T,VECTOR<int,1> >& u,const ARRAY<T,VECTOR<int,1> >& distance,ARRAY<T,VECTOR<int,1> >& u_phix,const T half_band_width)
{
    T one_over_dx=1/dx;
    if(order == 1){for(int i=1;i<=m;i++) if(abs(distance(i)) <= half_band_width){
        if(u(i) > 0) u_phix(i)=u(i)*(phi(i)-phi(i-1))*one_over_dx;
        else u_phix(i)=u(i)*(phi(i+1)-phi(i))*one_over_dx;}}
    else if(order == 2){for(int i=1;i<=m;i++) if(abs(distance(i)) <= half_band_width){
        if(u(i) > 0) u_phix(i)=u(i)*(phi(i)-phi(i-1)+(T).5*minmag(phi(i)-2*phi(i-1)+phi(i-2),phi(i+1)-2*phi(i)+phi(i-1)))*one_over_dx;
        else u_phix(i)=u(i)*(phi(i+1)-phi(i)-(T).5*minmag(phi(i+2)-2*phi(i+1)+phi(i),phi(i+1)-2*phi(i)+phi(i-1)))*one_over_dx;}}
    else if(order == 3){for(int i=1;i<=m;i++) if(abs(distance(i)) <= half_band_width){
        if(u(i) > 0){
            T D2_left=phi(i)-2*phi(i-1)+phi(i-2),D2_right=phi(i+1)-2*phi(i)+phi(i-1); // both scaled (multiplied by) 2*dx*dx
            if(abs(D2_left) <= abs(D2_right))
                u_phix(i)=u(i)*(phi(i)-phi(i-1)+(T).5*D2_left+(T)one_third*minmag(phi(i)-3*(phi(i-1)-phi(i-2))-phi(i-3),phi(i+1)-3*(phi(i)-phi(i-1))-phi(i-2)))*one_over_dx;
            else u_phix(i)=u(i)*(phi(i)-phi(i-1)+(T).5*D2_right-(T)one_sixth*minmag(phi(i+2)-3*(phi(i+1)-phi(i))-phi(i-1),phi(i+1)-3*(phi(i)-phi(i-1))-phi(i-2)))*one_over_dx;}
        else{
            T D2_left=phi(i+2)-2*phi(i+1)+phi(i),D2_right=phi(i+1)-2*phi(i)+phi(i-1); // both scaled (multiplied by) -2*dx*dx
            if(abs(D2_left) <= abs(D2_right))
                u_phix(i)=u(i)*(phi(i+1)-phi(i)-(T).5*D2_left+(T)one_third*minmag(phi(i+3)-3*(phi(i+2)-phi(i+1))+phi(i),phi(i+2)-3*(phi(i+1)-phi(i))+phi(i-1)))*one_over_dx;
            else u_phix(i)=u(i)*(phi(i+1)-phi(i)-(T).5*D2_right-(T)one_sixth*minmag(phi(i+1)-3*(phi(i)-phi(i-1))+phi(i-2),phi(i+2)-3*(phi(i+1)-phi(i))+phi(i-1)))*one_over_dx;}}}
}
//#####################################################################
// Function Local_WENO_Reinitialize
//#####################################################################
// phi is (-2,m+3), distance and phix_minus and phix_plus are (1,m)
template<class T> static void
Local_WENO_Reinitialize(const int m,const T dx,const ARRAY<T,VECTOR<int,1> >& phi,const ARRAY<T,VECTOR<int,1> >& distance,ARRAY<T,VECTOR<int,1> >& phix_minus,ARRAY<T,VECTOR<int,1> >& phix_plus,const T half_band_width)
{
    T epsilon=(T)1e-6*sqr(dx); // 1e-6 works since phi is a distance function - sqr(dx) since undivided differences are used
    T one_over_dx=1/dx;
    for(int i=1;i<=m;i++) if(abs(distance(i)) <= half_band_width){ // one_over_dx since undivided differences are used
        phix_minus(i)=ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<T,1> >,T>::WENO(phi(i-2)-phi(i-3),phi(i-1)-phi(i-2),phi(i)-phi(i-1),phi(i+1)-phi(i),phi(i+2)-phi(i+1),epsilon)*one_over_dx;
        phix_plus(i)=ADVECTION_SEPARABLE_UNIFORM<GRID<VECTOR<T,1> >,T>::WENO(phi(i+3)-phi(i+2),phi(i+2)-phi(i+1),phi(i+1)-phi(i),phi(i)-phi(i-1),phi(i-1)-phi(i-2),epsilon)*one_over_dx;}
}
//#####################################################################
// Function Local_ENO_Reinitialize
//#####################################################################
// order = 1, 2 or 3, phi is (-2,m_3), phix_minus and phix_plus are (1,m)
template<class T> static void
Local_ENO_Reinitialize(const int order,const int m,const T dx,const ARRAY<T,VECTOR<int,1> >& phi,const ARRAY<T,VECTOR<int,1> >& distance,ARRAY<T,VECTOR<int,1> >& phix_minus,ARRAY<T,VECTOR<int,1> >& phix_plus,const T half_band_width)
{
    T one_over_dx=1/dx;
    if(order == 1){for(int i=1;i<=m;i++) if(abs(distance(i)) <= half_band_width){
        phix_minus(i)=(phi(i)-phi(i-1))*one_over_dx;
        phix_plus(i)=(phi(i+1)-phi(i))*one_over_dx;}}
    else if(order == 2){for(int i=1;i<=m;i++) if(abs(distance(i)) <= half_band_width){
        phix_minus(i)=(phi(i)-phi(i-1)+(T).5*minmag(phi(i)-2*phi(i-1)+phi(i-2),phi(i+1)-2*phi(i)+phi(i-1)))*one_over_dx;
        phix_plus(i)=(phi(i+1)-phi(i)-(T).5*minmag(phi(i+2)-2*phi(i+1)+phi(i),phi(i+1)-2*phi(i)+phi(i-1)))*one_over_dx;}}
    else if(order == 3){for(int i=1;i<=m;i++) if(abs(distance(i)) <= half_band_width){
        T D2_left=phi(i)-2*phi(i-1)+phi(i-2),D2_right=phi(i+1)-2*phi(i)+phi(i-1); // both scaled (multiplied by) 2*dx*dx
        if(abs(D2_left) <= abs(D2_right))
            phix_minus(i)=(phi(i)-phi(i-1)+(T).5*D2_left+(T)one_third*minmag(phi(i)-3*(phi(i-1)-phi(i-2))-phi(i-3),phi(i+1)-3*(phi(i)-phi(i-1))-phi(i-2)))*one_over_dx;
        else phix_minus(i)=(phi(i)-phi(i-1)+(T).5*D2_right-(T)one_sixth*minmag(phi(i+2)-3*(phi(i+1)-phi(i))-phi(i-1),phi(i+1)-3*(phi(i)-phi(i-1))-phi(i-2)))*one_over_dx;
        D2_left=phi(i+2)-2*phi(i+1)+phi(i);D2_right=phi(i+1)-2*phi(i)+phi(i-1); // both scaled (multiplied by) -2*dx*dx
        if(abs(D2_left) <= abs(D2_right))
            phix_plus(i)=(phi(i+1)-phi(i)-(T).5*D2_left+(T)one_third*minmag(phi(i+3)-3*(phi(i+2)-phi(i+1))+phi(i),phi(i+2)-3*(phi(i+1)-phi(i))+phi(i-1)))*one_over_dx;
        else phix_plus(i)=(phi(i+1)-phi(i)-(T).5*D2_right-(T)one_sixth*minmag(phi(i+1)-3*(phi(i)-phi(i-1))+phi(i-2),phi(i+2)-3*(phi(i+1)-phi(i))+phi(i-1)))*one_over_dx;}}
}
//#####################################################################
template class FAST_LEVELSET_ADVECTION<GRID<VECTOR<float,1> > >;
template class FAST_LEVELSET_ADVECTION<GRID<VECTOR<float,2> > >;
template class FAST_LEVELSET_ADVECTION<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FAST_LEVELSET_ADVECTION<GRID<VECTOR<double,1> > >;
template class FAST_LEVELSET_ADVECTION<GRID<VECTOR<double,2> > >;
template class FAST_LEVELSET_ADVECTION<GRID<VECTOR<double,3> > >;
#endif
