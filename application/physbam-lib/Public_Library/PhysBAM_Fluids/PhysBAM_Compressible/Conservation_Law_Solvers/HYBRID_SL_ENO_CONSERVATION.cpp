//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HYBRID_SL_ENO_CONSERVATION  
//##################################################################### 
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UTILITIES.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/BOUNDARY_OBJECT.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/HYBRID_SL_ENO_CONSERVATION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/BOUNDARY_OBJECT_REFLECTION.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER.h>
using namespace PhysBAM;
//#####################################################################
// Function Update_Conservation_Law
//#####################################################################
template<class T_GRID,int d> void HYBRID_SL_ENO_CONSERVATION<T_GRID,d>::
Update_Conservation_Law(T_GRID& grid,T_ARRAYS_DIMENSION_SCALAR& U,const T_ARRAYS_DIMENSION_SCALAR& U_ghost,const T_ARRAYS_BOOL& psi,const T dt,
    VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>& eigensystems,VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>& eigensystems_explicit,const T_FACE_ARRAYS_BOOL& psi_N,
    const T_FACE_ARRAYS_SCALAR& face_velocities,const bool thinshell,const TV_BOOL& outflow_boundaries,VECTOR<EIGENSYSTEM<T,TV_DIMENSION>*,T_GRID::dimension>* eigensystems_auxiliary,
    T_FACE_ARRAYS_DIMENSION_SCALAR* fluxes_auxiliary)
{
    const T cell_volume=grid.Cell_Size();
    const T one_over_cell_volume=(T)1/cell_volume;

    T_ARRAYS_BOOL regular_cell(psi);
    T_ARRAYS_BOOL cell_near_interface(psi),cell_near_interface_tmp(psi);
    T_FACE_ARRAYS_DIMENSION_SCALAR& face_fluxes(conservation->fluxes);
    T_ARRAYS_DIMENSION_SCALAR rhs(U.Domain_Indices(),true);

    {
        LOG::SCOPE scope("Regular Update for Hybrid scheme.");
        for(CELL_ITERATOR iterator(grid,3);iterator.Valid();iterator.Next()){
            if(regular_cell(iterator.Cell_Index())){
                bool compute_self_weight=true;for(int dim=1;dim<=TV::dimension;++dim) compute_self_weight &= flux_face(iterator.Full_First_Face_Index(dim)) & flux_face(iterator.Full_Second_Face_Index(dim));
                regular_cell(iterator.Cell_Index())=compute_self_weight;
                cell_near_interface_tmp(iterator.Cell_Index())=!compute_self_weight;}}
        for(CELL_ITERATOR iterator(grid,2);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            cell_near_interface(cell_index)=cell_near_interface_tmp(cell_index);
            if(cell_near_interface_tmp(cell_index)){
                for(int dim=1;dim<=TV::dimension;++dim) {
                    if(psi(cell_index+TV_INT::Axis_Vector(dim))) cell_near_interface(cell_index+TV_INT::Axis_Vector(dim))=true;
                    if(psi(cell_index-TV_INT::Axis_Vector(dim))) cell_near_interface(cell_index-TV_INT::Axis_Vector(dim))=true;}
#if 0
                if(TV::dimension==2){
                    if(psi(cell_index+VECTOR<int,2>(-1,-1))) cell_near_interface(cell_index+VECTOR<int,2>(-1,-1))=true; if(psi(cell_index+VECTOR<int,2>( 1,-1))) cell_near_interface(cell_index+VECTOR<int,2>( 1,-1))=true;
                    if(psi(cell_index+VECTOR<int,2>(-1, 1))) cell_near_interface(cell_index+VECTOR<int,2>(-1, 1))=true; if(psi(cell_index+VECTOR<int,2>( 1, 1))) cell_near_interface(cell_index+VECTOR<int,2>( 1, 1))=true;}
                if(TV::dimension==3){
                    if(psi(cell_index+VECTOR<int,3>(-1,-1,-1))) cell_near_interface(cell_index+VECTOR<int,3>(-1,-1,-1))=true; if(psi(cell_index+VECTOR<int,3>( 1,-1,-1))) cell_near_interface(cell_index+VECTOR<int,3>( 1,-1,-1))=true;
                    if(psi(cell_index+VECTOR<int,3>(-1, 1,-1))) cell_near_interface(cell_index+VECTOR<int,3>(-1, 1,-1))=true; if(psi(cell_index+VECTOR<int,3>( 1, 1,-1))) cell_near_interface(cell_index+VECTOR<int,3>( 1, 1,-1))=true;
                    if(psi(cell_index+VECTOR<int,3>(-1,-1, 1))) cell_near_interface(cell_index+VECTOR<int,3>(-1,-1, 1))=true; if(psi(cell_index+VECTOR<int,3>( 1,-1, 1))) cell_near_interface(cell_index+VECTOR<int,3>( 1,-1, 1))=true;
                    if(psi(cell_index+VECTOR<int,3>(-1, 1, 1))) cell_near_interface(cell_index+VECTOR<int,3>(-1, 1, 1))=true; if(psi(cell_index+VECTOR<int,3>( 1, 1, 1))) cell_near_interface(cell_index+VECTOR<int,3>( 1, 1, 1))=true;}
#endif
            }
        }

        conservation->fluxes.Resize(grid.Domain_Indices(),true,false);
        conservation->Compute_Flux(grid,U,U_ghost,psi,dt,eigensystems,eigensystems_explicit,psi_N,face_velocities,outflow_boundaries,rhs,thinshell,eigensystems_auxiliary,fluxes_auxiliary);
    }

    {
        LOG::SCOPE scope("Irregular update for hybrid scheme (no collision bodies present).");
        U.Fill(TV_DIMENSION());
        for(int dim=1;dim<=TV_DIMENSION::dimension;++dim){
            ARRAY<PAIR<T,INDEX> > weights;
            ARRAY<ARRAY<int>,INDEX> donors;donors.Resize(grid.Domain_Indices(3));
            ARRAY<ARRAY<int>,INDEX> receivers;receivers.Resize(grid.Domain_Indices(3));
            T_ARRAYS_SCALAR sigma;sigma.Resize(grid.Domain_Indices(3));

            for(FACE_ITERATOR iterator(grid,2);iterator.Valid();iterator.Next()){
                FACE_INDEX<TV::dimension> face_index=iterator.Full_Index();
                if(flux_face(face_index) && (cell_near_interface(iterator.First_Cell_Index()) || cell_near_interface(iterator.Second_Cell_Index()))){
                    T weight=dt*face_fluxes(face_index)(dim);
                    INDEX donor_cell    = (weight >= 0 ? iterator.First_Cell_Index() : iterator.Second_Cell_Index());
                    INDEX receiver_cell = (weight >= 0 ? iterator.Second_Cell_Index() : iterator.First_Cell_Index());
                    weight = abs(weight);
                    int index=weights.Append(PAIR<T, INDEX>(weight,receiver_cell));
                    sigma(donor_cell)+=weight; donors(donor_cell).Append(index);receivers(receiver_cell).Append(index);}}

            // "Finish" cells that are entirely updated via the high order method.
            for(CELL_ITERATOR iterator(grid,2);iterator.Valid();iterator.Next()){
                INDEX cell_index=iterator.Cell_Index();
                if(regular_cell(cell_index) && cell_near_interface(cell_index)){
                    T weight=U_ghost(cell_index)(dim)*cell_volume - sigma(cell_index);
                    int index=weights.Append(PAIR<T,INDEX>(weight,cell_index));
                    donors(cell_index).Append(index); receivers(cell_index).Append(index);
                    sigma(cell_index) = U_ghost(cell_index)(dim)*cell_volume;}}

            // Compute backward-cast weights
            LINEAR_INTERPOLATION_UNIFORM<T_GRID,T> linear;
            for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){
                INDEX cell_index=iterator.Cell_Index();
                if(!regular_cell(cell_index)){
                    RANGE<TV> cell_preimage=iterator.Bounding_Box() - dt*linear.Clamped_To_Array_Face(grid,face_velocities,iterator.Location());
                    RANGE<TV_INT> affected_cells(grid.Index(cell_preimage.min_corner),grid.Index(cell_preimage.max_corner)+TV_INT::All_Ones_Vector());
                    for(CELL_ITERATOR intersecting_iter(grid,affected_cells);intersecting_iter.Valid();intersecting_iter.Next()){
                        INDEX donor_cell=intersecting_iter.Cell_Index();
                        if(regular_cell(donor_cell)) continue;
                        T weight=U_ghost(donor_cell)(dim) * cell_preimage.Intersection_Area(intersecting_iter.Bounding_Box());
                        if(weight > (T)1e-14){
                            int index=weights.Append(PAIR<T,INDEX>(weight,cell_index));
                            sigma(donor_cell) += weight; donors(donor_cell).Append(index); receivers(cell_index).Append(index);}}}}

            // Clamp and forward-cast weights
            for(CELL_ITERATOR iterator(grid,1);iterator.Valid();iterator.Next()){
                INDEX cell_index=iterator.Cell_Index();
                if(!regular_cell(cell_index)){
                    T cell_stuff = U_ghost(cell_index)(dim) * cell_volume;
                    if(abs(sigma(cell_index)) > abs(cell_stuff)) {
                        T one_over_sigma = cell_stuff / sigma(cell_index);
                        for(int i=1;i<= donors(cell_index).Size();++i){
                            weights(donors(cell_index)(i)).x *= one_over_sigma;}}
                    else {
                        T remainder = (cell_stuff - sigma(cell_index))/cell_volume;
                        RANGE<TV> cell_postimage=iterator.Bounding_Box() + dt*linear.Clamped_To_Array_Face(grid,face_velocities,iterator.Location());
                        RANGE<TV_INT> affected_cells(grid.Index(cell_postimage.min_corner),grid.Index(cell_postimage.max_corner)+TV_INT::All_Ones_Vector());
                        for(CELL_ITERATOR intersecting_iter(grid,affected_cells);intersecting_iter.Valid();intersecting_iter.Next()){
                            INDEX receiver_cell=intersecting_iter.Cell_Index();
                            T weight= remainder * cell_postimage.Intersection_Area(intersecting_iter.Bounding_Box());
                            int index=weights.Append(PAIR<T,INDEX>(weight,receiver_cell));
                            donors(cell_index).Append(index); receivers(receiver_cell).Append(index);}}}}

            for(int i=1;i<=weights.Size();++i){
                PAIR<T,INDEX>& weight_pair(weights(i));
                if(U.Valid_Index(weight_pair.y)) U(weight_pair.y)(dim) += weight_pair.x;}}
    }

    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){INDEX cell_index=iterator.Cell_Index();
        if(cell_near_interface(cell_index)) U(cell_index) *= one_over_cell_volume;
        else if(psi(cell_index)) U(cell_index) = U_ghost(cell_index) - dt*rhs(cell_index);}
}
template class HYBRID_SL_ENO_CONSERVATION<GRID<VECTOR<float,1> >,3>;
template class HYBRID_SL_ENO_CONSERVATION<GRID<VECTOR<float,2> >,4>;
template class HYBRID_SL_ENO_CONSERVATION<GRID<VECTOR<float,3> >,5>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class HYBRID_SL_ENO_CONSERVATION<GRID<VECTOR<double,1> >,3>;
template class HYBRID_SL_ENO_CONSERVATION<GRID<VECTOR<double,2> >,4>;
template class HYBRID_SL_ENO_CONSERVATION<GRID<VECTOR<double,3> >,5>;
#endif
