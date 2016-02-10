#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_OCTREE_HELPER.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_OCTREE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_QUADTREE.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/PROJECTION_DYADIC.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID> PROJECTION_DYADIC<T_GRID>::
PROJECTION_DYADIC(T_GRID& grid_input,ARRAY<T>* flame_speed_input,const bool flame_input)
    :PROJECTION<T>(flame_input),grid(grid_input),flame_speed(flame_speed_input)
{
    if(!flame){laplace=new LAPLACE_COLLIDABLE_DYADIC<T_GRID>(grid,p);elliptic_solver=laplace;poisson=0;}
    else{poisson=new POISSON_COLLIDABLE_DYADIC<T_GRID>(grid,p);elliptic_solver=poisson;laplace=0;}
    Initialize_Grid();
    if(flame) poisson->Use_Internal_Level_Set();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID> PROJECTION_DYADIC<T_GRID>::
~PROJECTION_DYADIC()
{
    delete elliptic_solver;
}
//#####################################################################
// Function Make_Divergence_Free
//#####################################################################
template<class T_GRID> void PROJECTION_DYADIC<T_GRID>::
Make_Divergence_Free(const T dt,const T time)
{
    // TODO: this should be a modified enforce compatibility call to ensure the t-junctions don't blow up - refine near objects works for now
    // average the t-junction face velocities together to satisfy the second order pressure solver
    if(!flame){
        MAP_MESH::template Map_Faces_For_Reduction<T>(grid.uniform_grid,grid.cells,1,NULL,face_velocities,INCOMPRESSIBLE::Compute_Aggregate_Face_Velocities);
        for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
            if(iterator.Deepest_Cell()->Depth_Of_This_Cell()>iterator.Other_Cell()->Depth_Of_This_Cell())
                face_velocities(iterator.Face_Index())=face_velocities(iterator.Other_Face_Index());}}

    // find f - divergence of the velocity
    elliptic_solver->f.Fill(0);
    if(!flame){
        for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
            if(elliptic_solver->psi_D(iterator.Deepest_Cell_Index())&&elliptic_solver->psi_D(iterator.Other_Cell_Index())) continue;
            T face_velocity=face_velocities(iterator.Face_Index());T small_size=iterator.Face_DX(),large_size=iterator.Other_Face_DX();T ratio=iterator.Face_Size()/iterator.Other_Face_Size();
            if(iterator.Deepest_Cell_Is_First_Cell()){
                elliptic_solver->f(iterator.Deepest_Cell_Index())+=face_velocity/small_size;elliptic_solver->f(iterator.Other_Cell_Index())-=ratio*face_velocity/large_size;}
            else{elliptic_solver->f(iterator.Deepest_Cell_Index())-=face_velocity/small_size;elliptic_solver->f(iterator.Other_Cell_Index())+=ratio*face_velocity/large_size;}}}
    else{
        for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
            if(elliptic_solver->psi_D(iterator.Deepest_Cell_Index())&&elliptic_solver->psi_D(iterator.Other_Cell_Index())) continue;
            T small_size=iterator.Face_DX(),large_size=iterator.Other_Face_DX();T ratio=iterator.Face_Size()/iterator.Other_Face_Size();
            CELL* deepest_cell=iterator.Deepest_Cell();T small_phi=0;for(int i=0;i<T_GRID::number_of_nodes_per_cell;i++)small_phi+=poisson->levelset->phi(deepest_cell->Node(i));
            CELL* other_cell=iterator.Other_Cell();T large_phi=0;for(int i=0;i<T_GRID::number_of_nodes_per_cell;i++)large_phi+=poisson->levelset->phi(other_cell->Node(i));
            T small_velocity=Face_Velocity_With_Ghost_Value(face_velocities,iterator.Face_Index(),small_phi);
            T large_velocity=Face_Velocity_With_Ghost_Value(face_velocities,iterator.Other_Face_Index(),large_phi);
            if(iterator.Deepest_Cell_Is_First_Cell()){
                elliptic_solver->f(iterator.Deepest_Cell_Index())+=small_velocity/small_size;
                elliptic_solver->f(iterator.Other_Cell_Index())-=ratio*large_velocity/large_size;}
            else{
                elliptic_solver->f(iterator.Deepest_Cell_Index())-=small_velocity/small_size;
                elliptic_solver->f(iterator.Other_Cell_Index())+=ratio*large_velocity/large_size;}}}
    if(use_non_zero_divergence) for(int i=1;i<=grid.number_of_cells;i++) elliptic_solver->f(i)-=divergence(i);

    // find the pressure
    elliptic_solver->Find_Solution_Regions(); // flood fill
    elliptic_solver->Compute_beta_And_Add_Jumps_To_b(dt,time); // only does something for poisson solver
    if(elliptic_solver->solve_neumann_regions) Enforce_Velocity_Compatibility(); // make all the right hand sides compatible
    elliptic_solver->Solve(time,true); // solve all regions

    // find divergence free face velocities
    if(!flame){
        for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
            int face_index=iterator.Face_Index();CELL* cell1=iterator.First_Cell();CELL* cell2=iterator.Second_Cell();
            if(elliptic_solver->psi_N(face_index))continue; // if the face separating the two cells is Neumann, don't touch it
            T distance=cell2->Center()[iterator.Axis()+1]-cell1->Center()[iterator.Axis()+1];
            if(poisson) face_velocities(face_index)-=poisson->beta_face(face_index)*(p(cell2->Cell())-p(cell1->Cell()))/distance;
            else if(!laplace->second_order_cut_cell_method) face_velocities(face_index)-=(p(cell2->Cell())-p(cell1->Cell()))/distance;
            else{ // second order cut cell method
                ARRAY<bool>& psi_D=laplace->psi_D;ARRAY<T>& cell_centered_phi_ghost=laplace->levelset->phi;
                if(psi_D(cell1->Cell()) && !psi_D(cell2->Cell())){
                    distance*=max(LEVELSET_UTILITIES<T>::Theta(cell_centered_phi_ghost(cell2->Cell()),cell_centered_phi_ghost(cell1->Cell())),laplace->second_order_cut_cell_threshold);
                    face_velocities(face_index)-=(p(cell2->Cell())-laplace->u_interface(face_index))/distance;}
                else if(!psi_D(cell1->Cell()) && psi_D(cell2->Cell())){
                    distance*=max(LEVELSET_UTILITIES<T>::Theta(cell_centered_phi_ghost(cell1->Cell()),cell_centered_phi_ghost(cell2->Cell())),laplace->second_order_cut_cell_threshold);
                    face_velocities(face_index)-=(laplace->u_interface(face_index)-p(cell1->Cell()))/distance;}
                else face_velocities(face_index)-=(p(cell2->Cell())-p(cell1->Cell()))/distance;}}}
    else{
        for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
            int face_index=iterator.Face_Index();CELL* cell1=iterator.First_Cell();CELL* cell2=iterator.Second_Cell();
            if(elliptic_solver->psi_N(face_index))continue; // if the face separating the two cells is Neumann, don't touch it
            T update=poisson->beta_face(face_index)*(p(cell2->Cell())-p(cell1->Cell()))/(cell2->Center()[iterator.Axis()+1]-cell1->Center()[iterator.Axis()+1]);
            face_velocities(face_index)-=update;}}

    // average the t-junction face velocities together to satisfy the second order pressure solver
    // note that we still use the first order scheme, but fix it in this postprocess (area weighting face velocities = area weighting pressure derivatives)
    if(!flame){
        MAP_MESH::template Map_Faces_For_Reduction<T>(grid.uniform_grid,grid.cells,1,NULL,face_velocities,INCOMPRESSIBLE::Compute_Aggregate_Face_Velocities);
        for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
            if(iterator.Deepest_Cell()->Depth_Of_This_Cell()>iterator.Other_Cell()->Depth_Of_This_Cell())
                face_velocities(iterator.Face_Index())=face_velocities(iterator.Other_Face_Index());}}

    // fix the jump in pressure - interior only
    if(poisson && poisson->u_jumps){
        ARRAY<T> phi_ghost(grid.number_of_cells);poisson->levelset->boundary->Fill_Ghost_Cells_Cell(grid,poisson->levelset->phi,phi_ghost,time);
        for(FACE_ITERATOR iterator(grid,grid.Map_Interior_Faces());iterator.Valid();iterator.Next()){
            int face_index=iterator.Face_Index();CELL* cell1=iterator.First_Cell();CELL* cell2=iterator.Second_Cell();
            if(!LEVELSET_UTILITIES<T>::Interface(phi_ghost(cell1->Cell()),phi_ghost(cell2->Cell()))||(poisson->psi_D(cell1->Cell())&&poisson->psi_D(cell2->Cell()))) continue;
            if(poisson->psi_N(face_index)) continue;
            T dx=(T).5*(cell1->DX()[iterator.Axis()+1]+cell2->DX()[iterator.Axis()+1]);
            T update=poisson->beta_face(face_index)/dx*LEVELSET_UTILITIES<T>::Sign(phi_ghost(cell2->Cell()))*
                LEVELSET_UTILITIES<T>::Average(phi_ghost(cell2->Cell()),poisson->u_jump(cell2->Cell()),phi_ghost(cell1->Cell()),poisson->u_jump(cell1->Cell()));
            face_velocities(face_index)+=update;}}
}
//#####################################################################
// Function Enforce_Velocity_Compatibility
//#####################################################################
namespace PhysBAM{ template<class T_GRID> struct ENFORCE_VELOCITY_COMPATIBILITY_HELPER{PROJECTION_DYADIC<T_GRID>* projection;ARRAY<typename T_GRID::SCALAR> *boundary_size,*compatibility_fraction;};}
template<class T_GRID> static void Calculate_Compatibility_Fraction_Helper(void* data,const typename T_GRID::CELL* cell)
{
    ENFORCE_VELOCITY_COMPATIBILITY_HELPER<T_GRID>* helper=(ENFORCE_VELOCITY_COMPATIBILITY_HELPER<T_GRID>*)data;
    int color=helper->projection->elliptic_solver->filled_region_colors(cell->Cell());
    if(color>0&&!helper->projection->elliptic_solver->filled_region_touches_dirichlet(color)){
        (*helper->compatibility_fraction)(color)+=helper->projection->elliptic_solver->f(cell->Cell())*cell->Cell_Size();}
}
template<class T_GRID> void PROJECTION_DYADIC<T_GRID>::
Enforce_Velocity_Compatibility()
{
    ARRAY<T> compatibility_fraction(elliptic_solver->number_of_regions),boundary_size(elliptic_solver->number_of_regions);
    ENFORCE_VELOCITY_COMPATIBILITY_HELPER<T_GRID> helper;helper.projection=this;helper.boundary_size=&boundary_size;helper.compatibility_fraction=&compatibility_fraction;
    MAP_MESH::Map_Cells(grid.uniform_grid,grid.cells,0,&helper,Calculate_Compatibility_Fraction_Helper<T_GRID>);

    for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
        CELL* cell1=iterator.First_Cell();CELL* cell2=iterator.Second_Cell();
        int color1=elliptic_solver->filled_region_colors(cell1->Cell()),color2=elliptic_solver->filled_region_colors(cell2->Cell());
        if(color1==color2)continue;T face_size=iterator.Face_Size();
        if(color1>0&&!elliptic_solver->filled_region_touches_dirichlet(color1)) boundary_size(color1)+=face_size;
        if(color2>0&&!elliptic_solver->filled_region_touches_dirichlet(color2)) boundary_size(color2)+=face_size;}

    for(int i=1;i<=elliptic_solver->number_of_regions;i++) if(boundary_size(i)) compatibility_fraction(i)/=boundary_size(i);

    for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
        CELL* cell1=iterator.First_Cell();CELL* cell2=iterator.Second_Cell();
        int color1=elliptic_solver->filled_region_colors(cell1->Cell()),color2=elliptic_solver->filled_region_colors(cell2->Cell());
        if(color1<=0||elliptic_solver->filled_region_touches_dirichlet(color1))color1=0;if(color2<=0||elliptic_solver->filled_region_touches_dirichlet(color2))color2=0;
        if(color1!=color2){T u_delta1=0,u_delta2=0;
            if(color1){u_delta1=-compatibility_fraction(color1);elliptic_solver->f(cell1->Cell())+=u_delta1/cell1->DX()[iterator.Axis()+1];}
            if(color2){u_delta2=compatibility_fraction(color2);elliptic_solver->f(cell2->Cell())-=u_delta2/cell2->DX()[iterator.Axis()+1];}
            T u_delta=u_delta1+u_delta2;if(color1&&color2)u_delta*=(T).5;
            face_velocities(iterator.Face_Index())+=u_delta;}}
}
//#####################################################################
// Function Update_Phi_And_Move_Velocity_Discontinuity
//#####################################################################
template<class T_GRID> void PROJECTION_DYADIC<T_GRID>::
Update_Phi_And_Move_Velocity_Discontinuity(const T_LEVELSET& levelset,const T time,const bool update_phi_only)
{
    assert(flame);assert(flame_speed);
    ARRAY<T> phi_face_new(grid.number_of_faces);T_LINEAR_INTERPOLATION_DYADIC_HELPER::Interpolate_From_Cells_To_Faces(grid,levelset.phi,phi_face_new);
    if(!update_phi_only){
        for(FACE_ITERATOR iterator(grid,grid.Map_Regular_Faces());iterator.Valid();iterator.Next()){
            if(LEVELSET_UTILITIES<T>::Interface(phi_face_for_flame(iterator.Face_Index()),phi_face_new(iterator.Face_Index()))){
                if(phi_face_new(iterator.Face_Index())<=0)face_velocities(iterator.Face_Index())-=Face_Jump(iterator.Face_Index());
                else face_velocities(iterator.Face_Index())+=Face_Jump(iterator.Face_Index());}}}
    phi_face_for_flame=phi_face_new;poisson->levelset->phi=levelset.phi;poisson->levelset->Compute_Normals();
}
//#####################################################################
template class PROJECTION_DYADIC<OCTREE_GRID<float> >;
template class PROJECTION_DYADIC<QUADTREE_GRID<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PROJECTION_DYADIC<OCTREE_GRID<double> >;
template class PROJECTION_DYADIC<QUADTREE_GRID<double> >;
#endif
#endif
