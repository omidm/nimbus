#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_HORIZONTAL.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_Y.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_TRANSFER_ITERATOR.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_RLE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> PROJECTION_RLE<T_GRID>::
PROJECTION_RLE(const T_GRID& grid_input,ARRAY<T>& V_input)
    :PROJECTION<T>(false),grid(grid_input),V(V_input),face_velocities(V),laplace(grid,p)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> PROJECTION_RLE<T_GRID>::
~PROJECTION_RLE()
{}
//#####################################################################
// Function Make_Divergence_Free
//#####################################################################
template<class T_GRID> void PROJECTION_RLE<T_GRID>::
Make_Divergence_Free(const T time,const ARRAY<T>* phi_ghost)
{
    // find f - total flux out of each cell
    laplace.f.Fill(0);
    T_GRID::template Horizontal_Face_Loop<Compute_Horizontal_Divergence>(*this);
    Compute_Vertical_Divergence();
    if(use_non_zero_divergence) PHYSBAM_NOT_IMPLEMENTED();

    // find the pressure
    laplace.Find_Solution_Regions(); // flood fill
    if(laplace.solve_neumann_regions) Enforce_Velocity_Compatibility(); // make all the right hand sides compatible
    else T_GRID::template Face_Loop<Zero_Out_Neumann_Regions>(*this); // zero out neumann region velocities to prevent gravity accumulation
    laplace.Solve(time,true,phi_ghost); // solve all regions

    // find divergence free u and v
    T_GRID::template Horizontal_Face_Loop<Apply_Horizontal_Pressure_Gradients>(*this,phi_ghost);
    Apply_Vertical_Pressure_Gradients(phi_ghost);
}
//#####################################################################
// Function Compute_Horizontal_Divergence
//#####################################################################
template<class T_GRID> template<class T_FACE> void PROJECTION_RLE<T_GRID>::
Compute_Horizontal_Divergence::Apply(PROJECTION_RLE<T_GRID>& projection)
{
    T area=projection.grid.uniform_grid.Face_Size(T_FACE::Axis());
    const ARRAY<T>& V=projection.V;LAPLACE_COLLIDABLE_RLE<T_GRID>& laplace=projection.laplace;
    if(projection.grid.long_run_faces_horizontal==1){
        for(T_FACE face(projection.grid,0);face;face++){int f=face.Face(),c1=face.cell1.Cell(),c2=face.cell2.Cell();
            T flux=face.Length()*area*V(f),j_mid=(T).5*(face.j()+face.jmax()-1);
            if(face.cell1.Short()) laplace.f(c1)+=flux;else{T w=(j_mid-face.cell1.j)/(face.cell1.length-1);laplace.f(c1)+=(1-w)*flux;laplace.f(c1+1)+=w*flux;}
            if(face.cell2.Short()) laplace.f(c2)-=flux;else{T w=(j_mid-face.cell2.j)/(face.cell2.length-1);laplace.f(c2)-=(1-w)*flux;laplace.f(c2+1)-=w*flux;}}}
    else{
        for(T_FACE face(projection.grid,0);face;face++){int f=face.Face(),c1=face.cell1.Cell(),c2=face.cell2.Cell();
            if(face.Short()){
                T flux=area*V(f);
                if(face.cell1.Short()) laplace.f(c1)+=flux;else{T w=(T)(face.j()-face.cell1.j)/(face.cell1.length-1);laplace.f(c1)+=(1-w)*flux;laplace.f(c1+1)+=w*flux;}
                if(face.cell2.Short()) laplace.f(c2)-=flux;else{T w=(T)(face.j()-face.cell2.j)/(face.cell2.length-1);laplace.f(c2)-=(1-w)*flux;laplace.f(c2+1)-=w*flux;}}
            else{
                int len_f=face.Length()-1,len_c1=face.cell1.length-1,len_c2=face.cell2.length-1,dj1=face.cell1.j-face.j(),dj2=face.cell2.j-face.j();
                T scale_1=area*(len_f+1)/(6*len_c1),scale_2=area*(len_f+1)/(6*len_c2),Vh_lo=V(f),Vh_hi=V(f+1);
                T half_length=(T).5*area*(len_f+1);
                T a1_hi_hi=scale_1*(2*len_f-3*dj1+1),a1_hi_lo=scale_1*(len_f-3*dj1-1),a1_lo_hi=half_length-a1_hi_hi,a1_lo_lo=half_length-a1_hi_lo;
                T a2_hi_hi=scale_2*(2*len_f-3*dj2+1),a2_hi_lo=scale_2*(len_f-3*dj2-1),a2_lo_hi=half_length-a2_hi_hi,a2_lo_lo=half_length-a2_hi_lo;
                laplace.f(c1)+=a1_lo_lo*Vh_lo+a1_lo_hi*Vh_hi;laplace.f(c1+1)+=a1_hi_lo*Vh_lo+a1_hi_hi*Vh_hi;
                laplace.f(c2)-=a2_lo_lo*Vh_lo+a2_lo_hi*Vh_hi;laplace.f(c2+1)-=a2_hi_lo*Vh_lo+a2_hi_hi*Vh_hi;}}}
}
//#####################################################################
// Function Compute_Vertical_Divergence
//#####################################################################
template<class T_GRID> void PROJECTION_RLE<T_GRID>::
Compute_Vertical_Divergence()
{
    T area=grid.uniform_grid.Face_Size(2);
    for(FACE_Y_ITERATOR face(grid,0,false);face;face++){int f=face.Face(),c2=face.cell2.Cell(),c1=c2-1;
        T flux=area*V(f);
        laplace.f(c1)+=flux;laplace.f(c2)-=flux;
        if(face.cell2.Long()){T flux2=area*V(f+1);laplace.f(c2)+=flux2;laplace.f(c2+1)-=flux2;}}
}
//#####################################################################
// Function Apply_Horizontal_Pressure_Gradients
//#####################################################################
template<class T_GRID> template<class T_FACE> void PROJECTION_RLE<T_GRID>::
Apply_Horizontal_Pressure_Gradients::Apply(const PROJECTION_RLE<T_GRID>& projection,const ARRAY<T>* phi_ghost)
{
    T one_over_length=1/projection.grid.uniform_grid.dX[T_FACE::Axis()];
    ARRAY<T>& V=projection.V;
    const LAPLACE_COLLIDABLE_RLE<T_GRID>& laplace=projection.laplace;
    const ARRAY<bool> &psi_N=laplace.psi_N,&psi_D=laplace.psi_D;const ARRAY<T>& p=projection.p;
    for(T_FACE face(projection.grid,0);face;face++){int f=face.Face(),c1=face.cell1.Cell(),c2=face.cell2.Cell();if(!(psi_N(f)||(psi_D(c1)&&psi_D(c2)))){
        if(laplace.second_order_cut_cell_method && LEVELSET_UTILITIES<T>::Interface((*phi_ghost)(c1),(*phi_ghost)(c2))){
            if(!face.Both_Cells_Short()) PHYSBAM_FATAL_ERROR();
            if(psi_D(c1)&&!psi_D(c2)) V(f)-=(p(c2)-laplace.u_interface(f))*one_over_length/max(LEVELSET_UTILITIES<T>::Theta((*phi_ghost)(c2),(*phi_ghost)(c1)),laplace.second_order_cut_cell_threshold);
            else if(!psi_D(c1)&&psi_D(c2)) V(f)-=(laplace.u_interface(f)-p(c1))*one_over_length/max(LEVELSET_UTILITIES<T>::Theta((*phi_ghost)(c1),(*phi_ghost)(c2)),laplace.second_order_cut_cell_threshold);
            else V(f)-=(p(c2)-p(c1))*one_over_length;}
        else if(projection.grid.long_run_faces_horizontal==2 && face.Long()){
            T p_slope_1=(p(c1+1)-p(c1))/(face.cell1.length-1),p_slope_2=(p(c2+1)-p(c2))/(face.cell2.length-1);int len_f=face.Length()-1;
            T dp_lo=(p(c2)+(face.j()-face.cell2.j)*p_slope_2)-(p(c1)+(face.j()-face.cell1.j)*p_slope_1);
            V(f)-=dp_lo*one_over_length;V(f+1)-=(dp_lo+len_f*(p_slope_2-p_slope_1))*one_over_length;}
        else{
            T dp=p(c2)-p(c1),j_mid=(T).5*(face.j()+face.jmax()-1);
            if(face.cell1.Long()) dp-=(p(c1+1)-p(c1))*(j_mid-face.cell1.j)/(face.cell1.length-1);
            if(face.cell2.Long()) dp+=(p(c2+1)-p(c2))*(j_mid-face.cell2.j)/(face.cell2.length-1);
            V(f)-=dp*one_over_length;}}}
}
//#####################################################################
// Function Apply_Vertical_Pressure_Gradients
//#####################################################################
template<class T_GRID> void PROJECTION_RLE<T_GRID>::
Apply_Vertical_Pressure_Gradients(const ARRAY<T>* phi_ghost)
{
    T one_over_length=1/grid.uniform_grid.dX.y;
    ARRAY<bool> &psi_N=laplace.psi_N,&psi_D=laplace.psi_D;
    for(FACE_Y_ITERATOR face(grid,0,false);face;face++){int f=face.Face(),c2=face.cell2.Cell(),c1=c2-1;if(!(psi_N(f)||(psi_D(c1)&&psi_D(c2)))){
        if(laplace.second_order_cut_cell_method && LEVELSET_UTILITIES<T>::Interface((*phi_ghost)(c1),(*phi_ghost)(c2))){
            if(!face.Both_Cells_Short()) PHYSBAM_FATAL_ERROR();
            if(psi_D(c1)&&!psi_D(c2)) V(f)-=(p(c2)-laplace.u_interface(f))*one_over_length/max(LEVELSET_UTILITIES<T>::Theta((*phi_ghost)(c2),(*phi_ghost)(c1)),laplace.second_order_cut_cell_threshold);
            else if(!psi_D(c1)&&psi_D(c2)) V(f)-=(laplace.u_interface(f)-p(c1))*one_over_length/max(LEVELSET_UTILITIES<T>::Theta((*phi_ghost)(c1),(*phi_ghost)(c2)),laplace.second_order_cut_cell_threshold);
            else V(f)-=(p(c2)-p(c1))*one_over_length;}
        else V(f)-=(p(c2)-p(c1))*one_over_length;}}
    for(CELL_ITERATOR cell(grid,0);cell;cell++)if(cell.Long()){int c=cell.Cell(),f=cell.Face_Y()+1;if(!(psi_N(f)||psi_D(c))){
        V(f)-=(p(c+1)-p(c))*one_over_length/(cell.length-1);}}
}
//#####################################################################
// Function Enforce_Velocity_Compatibility
//#####################################################################
// For each Neumann region, modifies divergence in cells touching its boundary so that the sum of f in the region is zero. Also modifies each bounding face velocity to the desired velocity
// (or the average of the desired velocity if the face joins two different Neumann regions).
template<class T_GRID> void PROJECTION_RLE<T_GRID>::
Enforce_Velocity_Compatibility()
{
    int first_neumann_region=laplace.filled_region_touches_dirichlet.Find(false);
    if(!first_neumann_region) return; // don't need to do any of the following
    for(CELL_ITERATOR cell(grid,0);cell;cell++){int c=cell.Cell(),color=laplace.filled_region_colors(c);
        if(color>0 && !laplace.filled_region_touches_dirichlet(color)){std::stringstream ss;ss<<"found neumann region containing cell "<<c<<", color "<<color<<", I "<<cell.I()<<std::endl;break;LOG::filecout(ss.str());}}
    if(grid.long_run_faces_horizontal!=1) PHYSBAM_NOT_IMPLEMENTED();
    if(laplace.mpi_grid && first_neumann_region<=laplace.laplace_mpi->filled_region_ranks.m) PHYSBAM_NOT_IMPLEMENTED();

    ARRAY<T> compatibility_fraction(laplace.number_of_regions),boundary_area(laplace.number_of_regions);
    // compute compatibility errors
    for(CELL_ITERATOR cell(grid,0);cell;cell++){
        int c=cell.Cell(),color=laplace.filled_region_colors(c);
        if(color>0 && !laplace.filled_region_touches_dirichlet(color)) compatibility_fraction(color)+=laplace.f(c);}
    T_GRID::template Face_Loop<Compute_Boundary_Areas>(*this,boundary_area);
    for(int i=1;i<=laplace.number_of_regions;i++)if(boundary_area(i)) compatibility_fraction(i)/=boundary_area(i);
    // adjust the compatibility error to zero
    T_GRID::template Face_Loop<Zero_Out_Compatibility_Errors>(*this,compatibility_fraction);
}
//#####################################################################
// Function Compute_Boundary_Areas
//#####################################################################
template<class T_GRID> template<class T_FACE> void PROJECTION_RLE<T_GRID>::
Compute_Boundary_Areas::Apply(const PROJECTION_RLE<T_GRID>& projection,ARRAY<T>& boundary_area)
{
    T short_area=projection.grid.uniform_grid.Face_Size(T_FACE::Axis());
    const ARRAY<int>& filled_region_colors=projection.laplace.filled_region_colors;
    const ARRAY<bool>& filled_region_touches_dirichlet=projection.laplace.filled_region_touches_dirichlet;
    for(T_FACE face(projection.grid,0,T_FACE::Axis()!=2);face;face++){
        int c1=face.Cell(0),c2=face.Cell(1),color1=filled_region_colors(c1),color2=filled_region_colors(c2);
        if(color1!=color2){T area=short_area*face.Length();
            if(color1>0 && !filled_region_touches_dirichlet(color1)) boundary_area(color1)+=area;
            if(color2>0 && !filled_region_touches_dirichlet(color2)) boundary_area(color2)+=area;}}
}
//#####################################################################
// Function Zero_Out_Compatibility_Errors
//#####################################################################
template<class T_GRID> template<class T_FACE> void PROJECTION_RLE<T_GRID>::
Zero_Out_Compatibility_Errors::Apply(PROJECTION_RLE<T_GRID>& projection,const ARRAY<T>& compatibility_fraction)
{
    T short_area=projection.grid.uniform_grid.Face_Size(T_FACE::Axis());
    ARRAY<T>& V=projection.V;
    ARRAY<int>& filled_region_colors=projection.laplace.filled_region_colors;
    ARRAY<bool>& filled_region_touches_dirichlet=projection.laplace.filled_region_touches_dirichlet;
    for(T_FACE face(projection.grid,0,T_FACE::Axis()!=2);face;face++){
        int f=face.Face(),c1=face.Cell(0),c2=face.Cell(1);
        int color1=filled_region_colors(c1),color2=filled_region_colors(c2);
        if(color1<=0 || filled_region_touches_dirichlet(color1)) color1=0;
        if(color2<=0 || filled_region_touches_dirichlet(color2)) color2=0;
        if(color1!=color2){
            T V_delta1=0,V_delta2=0,area=short_area*face.Length();
            if(color1){V_delta1=-compatibility_fraction(color1);projection.laplace.f(c1)+=area*V_delta1;}
            if(color2){V_delta2=compatibility_fraction(color2);projection.laplace.f(c2)-=area*V_delta2;}
            T V_delta=V_delta1+V_delta2;if(color1&&color2) V_delta*=(T).5;V(f)+=V_delta;}}
}
//#####################################################################
// Function Zero_OOut_Neumann_Regions
//#####################################################################
template<class T_GRID> template<class T_FACE> void PROJECTION_RLE<T_GRID>::
Zero_Out_Neumann_Regions::Apply(PROJECTION_RLE<T_GRID>& projection)
{
    ARRAY<T>& V=projection.V;
    ARRAY<int>& filled_region_colors=projection.laplace.filled_region_colors;
    for(T_FACE face(projection.grid,0,T_FACE::Axis()!=2);face;face++){
        int f=face.Face(),c1=face.Cell(0),c2=face.Cell(1);
        if(filled_region_colors(c1)==-2 && filled_region_colors(c2)==-2){V(f)=0;if(face.Long()) V(f+1)=0;}} // unsolved neumann regions are colored -2
}
//#####################################################################
// Function Transfer_Pressure
//#####################################################################
template<class T_GRID> void PROJECTION_RLE<T_GRID>::
Transfer_Pressure(const T_GRID& new_grid)
{
    ARRAY<T> new_p(new_grid.number_of_cells);
    for(RLE_GRID_TRANSFER_ITERATOR<T_GRID,CELL_ITERATOR> cells(grid,new_grid,3);cells;cells++){
        int destination=cells.destination.Cell(),jmin_destination=cells.destination.j,jmax_destination=cells.destination.jmax();
        if(cells.destination.Short()){
            new_p(destination)=cells.source.Cell_Value(p,jmin_destination);}
        else{
            if(cells.source.j<=jmin_destination) new_p(destination)=cells.source.Cell_Value(p,jmin_destination);
            if(cells.source.jmax()>=jmax_destination) new_p(destination+1)=cells.source.Cell_Value(p,jmax_destination-1);}}
    ARRAY<T>::Exchange_Arrays(p,new_p);
}
//#####################################################################
template class PROJECTION_RLE<RLE_GRID_2D<float> >;
template class PROJECTION_RLE<RLE_GRID_3D<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PROJECTION_RLE<RLE_GRID_2D<double> >;
template class PROJECTION_RLE<RLE_GRID_3D<double> >;
#endif
#endif
