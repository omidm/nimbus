//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Dynamics/Incompressible_Flows/DIFFUSION.h>
#define DEBUG 1

using namespace PhysBAM;

template<class T_GRID,class T2> DIFFUSION_UNIFORM<T_GRID,T2>::
DIFFUSION_UNIFORM(T_MPI_GRID* mpi_grid_input)
    :diffuse_weights(false),diffuse_errors(false),mpi_grid(mpi_grid_input),num_diffusion_iterations(10),evenodd(0),evenodd_cell(0),max_value(0),min_value(0),hard_max(0),hard_min(-1)
{
    for(int i=1;i<=TV::dimension;i++){solid_walls(i)(1)=true;solid_walls(i)(2)=true;}        
    if(mpi_grid) mpi_grid->Initialize(solid_walls);
    for(int i=1;i<=TV::dimension;i++){mpi_boundary(i)(1)=!solid_walls(i)(1);mpi_boundary(i)(2)=!solid_walls(i)(2);}
}

template<class T_GRID,class T2> DIFFUSION_UNIFORM<T_GRID,T2>::
~DIFFUSION_UNIFORM()
{
}

template<class T_GRID,class T2> void DIFFUSION_UNIFORM<T_GRID,T2>::
Cell_Diffusion_Value_Helper(FACE_ITERATOR& iterator,T_ARRAYS_T2& Z,ARRAY<bool,TV_INT>* inside)
{
    if(inside && !((*inside)(iterator.First_Cell_Index()) && (*inside)(iterator.Second_Cell_Index()))) return;
    T2 Z_diff=(Z(iterator.First_Cell_Index())-Z(iterator.Second_Cell_Index()));Z_diff*=0.5;
    Z(iterator.Second_Cell_Index())+=Z_diff;Z(iterator.First_Cell_Index())-=Z_diff;
}

template<class T_GRID,class T2> void DIFFUSION_UNIFORM<T_GRID,T2>::
Cell_Diffusion_Sum_Helper(FACE_ITERATOR& iterator,ARRAY<T,TV_INT>& sum_jc_cell,T_ARRAYS_T2& Z,ARRAY<bool,TV_INT>* inside)
{
    if(inside && !((*inside)(iterator.First_Cell_Index()) && (*inside)(iterator.Second_Cell_Index()))) return;
    T wjc_diff=(sum_jc_cell(iterator.First_Cell_Index())-sum_jc_cell(iterator.Second_Cell_Index()))/(T)2.;
    T Z_diff=wjc_diff>0?wjc_diff/sum_jc_cell(iterator.First_Cell_Index()):wjc_diff/sum_jc_cell(iterator.Second_Cell_Index());
    sum_jc_cell(iterator.Second_Cell_Index())+=wjc_diff;sum_jc_cell(iterator.First_Cell_Index())-=wjc_diff;
    T2 local_Z=wjc_diff>0?Z(iterator.First_Cell_Index()):Z(iterator.Second_Cell_Index());
    Z(iterator.Second_Cell_Index())+=local_Z*(Z_diff);Z(iterator.First_Cell_Index())-=local_Z*(Z_diff);
}

template<class T_GRID,class T2> void DIFFUSION_UNIFORM<T_GRID,T2>::
Cell_Diffusion_Error_Helper(FACE_ITERATOR& iterator,ARRAY<T,TV_INT>& error_cell,T_ARRAYS_T2& Z,ARRAY<bool,TV_INT>* inside)
{
    if(inside && !((*inside)(iterator.First_Cell_Index()) && (*inside)(iterator.Second_Cell_Index()))) return;
    T wjc_diff=(error_cell(iterator.First_Cell_Index())-error_cell(iterator.Second_Cell_Index()))/(T)2.;
    if(hard_min>=0 && Z(iterator.Second_Cell_Index())-wjc_diff<hard_min) wjc_diff=Z(iterator.Second_Cell_Index())-hard_min;
    else if(hard_min>=0 && Z(iterator.First_Cell_Index())+wjc_diff<hard_min) wjc_diff=hard_min-Z(iterator.First_Cell_Index());
    if(hard_max && Z(iterator.Second_Cell_Index())-wjc_diff>hard_max && Z(iterator.First_Cell_Index())+wjc_diff>hard_max){
        if(Z(iterator.Second_Cell_Index())<hard_max || Z(iterator.First_Cell_Index())<hard_max){
            if(Z(iterator.Second_Cell_Index())<hard_max) wjc_diff=Z(iterator.Second_Cell_Index())-hard_max;
            else if(Z(iterator.First_Cell_Index())<hard_max) wjc_diff=hard_max-Z(iterator.First_Cell_Index());
            error_cell(iterator.Second_Cell_Index())+=wjc_diff;error_cell(iterator.First_Cell_Index())-=wjc_diff;
            Z(iterator.Second_Cell_Index())-=wjc_diff;Z(iterator.First_Cell_Index())+=wjc_diff;
            wjc_diff=(error_cell(iterator.First_Cell_Index())-error_cell(iterator.Second_Cell_Index()))/(T)2.;}}
    else if(hard_max && Z(iterator.Second_Cell_Index())-wjc_diff>hard_max) wjc_diff=Z(iterator.Second_Cell_Index())-hard_max;
    else if(hard_max && Z(iterator.First_Cell_Index())+wjc_diff>hard_max) wjc_diff=hard_max-Z(iterator.First_Cell_Index());
    error_cell(iterator.Second_Cell_Index())+=wjc_diff;error_cell(iterator.First_Cell_Index())-=wjc_diff;
    Z(iterator.Second_Cell_Index())-=wjc_diff;Z(iterator.First_Cell_Index())+=wjc_diff;
    if(max_value){
        if(Z(iterator.First_Cell_Index())>max_value) error_cell(iterator.First_Cell_Index())=max_value-Z(iterator.First_Cell_Index());
        else if(Z(iterator.First_Cell_Index())<min_value) error_cell(iterator.First_Cell_Index())=min_value-Z(iterator.First_Cell_Index());
        else error_cell(iterator.First_Cell_Index())=0;
        if(Z(iterator.Second_Cell_Index())>max_value) error_cell(iterator.Second_Cell_Index())=max_value-Z(iterator.Second_Cell_Index());
        else if(Z(iterator.Second_Cell_Index())<min_value) error_cell(iterator.Second_Cell_Index())=min_value-Z(iterator.Second_Cell_Index());
        else error_cell(iterator.Second_Cell_Index())=0;}
}

template<class T_GRID,class T2> void DIFFUSION_UNIFORM<T_GRID,T2>::
Cell_Diffusion_Helper(FACE_ITERATOR& iterator,ARRAY<T,TV_INT>* sum_jc_cell,T_ARRAYS_T2& Z,ARRAY<bool,TV_INT>* inside)
{
    if(!sum_jc_cell) Cell_Diffusion_Value_Helper(iterator,Z,inside);
    else if(diffuse_weights) Cell_Diffusion_Sum_Helper(iterator,*sum_jc_cell,Z,inside);
    else if(diffuse_errors) Cell_Diffusion_Error_Helper(iterator,*sum_jc_cell,Z,inside);
    else PHYSBAM_FATAL_ERROR();
}

template<class T_GRID,class T2> void DIFFUSION_UNIFORM<T_GRID,T2>::
Face_Diffusion_Sum_Helper(const GRID<TV>& grid,FACE_INDEX<TV::dimension>& first_face_index,FACE_INDEX<TV::dimension>& second_face_index,ARRAY<T,FACE_INDEX<TV::dimension> >& sum_jc,T_FACE_ARRAYS_SCALAR& Z,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside)
{
    if(inside && !((*inside)(first_face_index) && (*inside)(second_face_index))) return;
    if(sum_jc(first_face_index)==0 && sum_jc(second_face_index)==0) return;
    T wjc_diff=(sum_jc(first_face_index)-sum_jc(second_face_index))/(T)2.;
    T Z_diff=wjc_diff>0?wjc_diff/sum_jc(first_face_index):wjc_diff/sum_jc(second_face_index);
    sum_jc(second_face_index)+=wjc_diff;sum_jc(first_face_index)-=wjc_diff;
    T local_Z=wjc_diff>0?Z(first_face_index):Z(second_face_index);
    Z(second_face_index)+=local_Z*(Z_diff);Z(first_face_index)-=local_Z*(Z_diff);
}

template<class T_GRID,class T2> void DIFFUSION_UNIFORM<T_GRID,T2>::
Face_Diffusion_Helper(FACE_ITERATOR& iterator,int axis,ARRAY<T,FACE_INDEX<TV::dimension> >* sum_jc,T_FACE_ARRAYS_SCALAR& Z,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside)
{
    assert(axis!=iterator.Axis());
    FACE_INDEX<TV::dimension> first_face_index=FACE_INDEX<TV::dimension>(axis,iterator.First_Cell_Index()),second_face_index=FACE_INDEX<TV::dimension>(axis,iterator.Second_Cell_Index());
    assert(diffuse_weights);Face_Diffusion_Sum_Helper(iterator.grid,first_face_index,second_face_index,*sum_jc,Z,inside);
}

template<class T_GRID,class T2> void DIFFUSION_UNIFORM<T_GRID,T2>::
Face_Diffusion_Helper(CELL_ITERATOR& iterator,int axis,ARRAY<T,FACE_INDEX<TV::dimension> >* sum_jc,T_FACE_ARRAYS_SCALAR& Z,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside)
{
    FACE_INDEX<TV::dimension> first_face_index=FACE_INDEX<TV::dimension>(axis,iterator.First_Face_Index(axis)),second_face_index=FACE_INDEX<TV::dimension>(axis,iterator.Second_Face_Index(axis));
    assert(diffuse_weights);Face_Diffusion_Sum_Helper(iterator.grid,first_face_index,second_face_index,*sum_jc,Z,inside);
}

template<class T_GRID,class T2> void DIFFUSION_UNIFORM<T_GRID,T2>::
Cell_Diffusion(const T_GRID& grid,T_ARRAYS_T2& Z,T_BOUNDARY_T2& boundary,ARRAY<T,TV_INT>* sum_jc_cell,BOUNDARY_UNIFORM<T_GRID,T>* boundary_sum,ARRAY<bool,TV_INT>* inside)
{
    T_ARRAYS_T2 Z_ghost_local(grid.Domain_Indices(1));
    ARRAY<T,TV_INT>* sum_jc_cell_ghost=0;if(mpi_grid && sum_jc_cell) sum_jc_cell_ghost=new ARRAY<T,TV_INT>(grid.Domain_Indices(1));
    for(int iter=1;iter<=num_diffusion_iterations;iter++){
        for(int axis=1;axis<=TV::dimension;axis++){
            if(evenodd_cell==0){
                RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner-=TV_INT::All_Ones_Vector();domain.min_corner+=TV_INT::All_Ones_Vector();domain.max_corner(axis)+=1;
                for(FACE_ITERATOR iterator(grid,domain,axis);iterator.Valid();iterator.Next()) Cell_Diffusion_Helper(iterator,sum_jc_cell,Z,inside);}
            else{
                ARRAY<FACE_ITERATOR*> faces;
                RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner-=TV_INT::All_Ones_Vector();domain.min_corner+=TV_INT::All_Ones_Vector();domain.max_corner(axis)+=1;
                for(FACE_ITERATOR iterator(grid,domain,axis);iterator.Valid();iterator.Next()) faces.Append(new FACE_ITERATOR(iterator));
                for(int i=faces.m;i>=1;i--){FACE_ITERATOR& iterator=*faces(i);Cell_Diffusion_Helper(iterator,sum_jc_cell,Z,inside);}
                for(int i=faces.m;i>=1;i--) delete faces(i);}
            if(mpi_grid){
                boundary.Fill_Ghost_Cells(grid,Z,Z_ghost_local,0,0,1); //Sync for mpi boundaries
                if(sum_jc_cell) boundary_sum->Fill_Ghost_Cells(grid,*sum_jc_cell,*sum_jc_cell_ghost,0,0,1); //Sync for mpi boundaries
                for(int side=1;side<=2;side++) if(mpi_boundary(axis)(side)){
                    if(evenodd_cell==0){
                        for(FACE_ITERATOR iterator(grid,0,GRID<TV>::BOUNDARY_REGION,2*(axis-1)+side);iterator.Valid();iterator.Next()) Cell_Diffusion_Helper(iterator,sum_jc_cell_ghost,Z_ghost_local,inside);}
                    else{
                        ARRAY<FACE_ITERATOR*> faces;
                        for(FACE_ITERATOR iterator(grid,0,GRID<TV>::BOUNDARY_REGION,2*(axis-1)+side);iterator.Valid();iterator.Next()) faces.Append(new FACE_ITERATOR(iterator));
                        for(int i=faces.m;i>=1;i--){FACE_ITERATOR& iterator=*faces(i);Cell_Diffusion_Helper(iterator,sum_jc_cell_ghost,Z_ghost_local,inside);}
                        for(int i=faces.m;i>=1;i--) delete faces(i);}}
                for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) Z(iterator.Cell_Index())=Z_ghost_local(iterator.Cell_Index());
                for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) (*sum_jc_cell)(iterator.Cell_Index())=(*sum_jc_cell_ghost)(iterator.Cell_Index());}}
        evenodd_cell++;evenodd_cell=evenodd_cell%2;}
    delete sum_jc_cell_ghost;
}

template<class T_GRID,class T2> void DIFFUSION_UNIFORM<T_GRID,T2>::
Face_Diffusion(const T_GRID& grid,ARRAY<T,FACE_INDEX<TV::dimension> >* sum_jc,T_FACE_ARRAYS_SCALAR& Z,T_BOUNDARY& boundary,BOUNDARY_UNIFORM<T_GRID,T>* boundary_sum,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside)
{
    T_FACE_ARRAYS_SCALAR Z_ghost_local(grid,1);
    ARRAY<T,FACE_INDEX<TV::dimension> >* sum_jc_ghost=0;if(sum_jc) sum_jc_ghost=new ARRAY<T,FACE_INDEX<TV::dimension> >(grid,1);
    for(int iter=1;iter<=num_diffusion_iterations;iter++){
        for(int axis=1;axis<=TV::dimension;axis++){
            if(evenodd==0){
                for(int axis2=1;axis2<=TV::dimension;axis2++){if(axis==axis2) continue;//handled above
                    RANGE<TV_INT> domain=grid.Domain_Indices();
                    domain.max_corner(axis)+=1;domain.min_corner(axis2)+=1;
                    for(FACE_ITERATOR iterator(grid,domain,axis2);iterator.Valid();iterator.Next()) Face_Diffusion_Helper(iterator,axis,sum_jc,Z,inside);}
                RANGE<TV_INT> domain=grid.Domain_Indices();
                if((mpi_grid && mpi_boundary(axis)(1))) domain.min_corner(axis)+=1; 
                if((mpi_grid && mpi_boundary(axis)(2))) domain.max_corner(axis)-=1; 
                for(CELL_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()) Face_Diffusion_Helper(iterator,axis,sum_jc,Z,inside);}
            else{
                ARRAY<FACE_ITERATOR*> faces;ARRAY<CELL_ITERATOR*> cells;
                for(int axis2=1;axis2<=TV::dimension;axis2++){if(axis==axis2) continue;//handled above
                    RANGE<TV_INT> domain=grid.Domain_Indices();
                    domain.max_corner(axis)+=1;domain.min_corner(axis2)+=1;
                    for(FACE_ITERATOR iterator(grid,domain,axis2);iterator.Valid();iterator.Next()) faces.Append(new FACE_ITERATOR(iterator));}
                RANGE<TV_INT> domain=grid.Domain_Indices();
                if((mpi_grid && mpi_boundary(axis)(1))) domain.min_corner(axis)+=1; 
                if((mpi_grid && mpi_boundary(axis)(2))) domain.max_corner(axis)-=1; 
                for(CELL_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()) cells.Append(new CELL_ITERATOR(iterator));
                for(int i=faces.m;i>=1;i--){FACE_ITERATOR& iterator=*faces(i);Face_Diffusion_Helper(iterator,axis,sum_jc,Z,inside);}
                for(int i=cells.m;i>=1;i--){CELL_ITERATOR& iterator=*cells(i);Face_Diffusion_Helper(iterator,axis,sum_jc,Z,inside);}
                for(int i=faces.m;i>=1;i--) delete faces(i);
                for(int i=cells.m;i>=1;i--) delete cells(i);}
            if(mpi_grid){
                boundary.Apply_Boundary_Condition_Face(grid,Z,0);
                boundary.Fill_Ghost_Cells_Face(grid,Z,Z_ghost_local,0,1); //Sync for mpi boundaries
                if(sum_jc){
                    boundary_sum->Apply_Boundary_Condition_Face(grid,*sum_jc,0);
                    boundary_sum->Fill_Ghost_Cells_Face(grid,*sum_jc,*sum_jc_ghost,0,1);} //Sync for mpi boundaries
                for(int side=1;side<=2;side++){
                    if(evenodd==0){
                        for(int axis2=1;axis2<=TV::dimension;axis2++){if(axis==axis2 || !mpi_boundary(axis2)(side)) continue;//handled above
                            RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner(axis)++;
                            if(side==1) domain.max_corner(axis2)=domain.min_corner(axis2);
                            else{domain.max_corner(axis2)++;domain.min_corner(axis2)=domain.max_corner(axis2);}
                            for(FACE_ITERATOR iterator(grid,domain,axis2);iterator.Valid();iterator.Next()) Face_Diffusion_Helper(iterator,axis,sum_jc_ghost,Z_ghost_local,inside);}
                        RANGE<TV_INT> domain=grid.Domain_Indices(1);
                        if(side==1) domain.max_corner(axis)=domain.min_corner(axis)+1;
                        else domain.min_corner(axis)=domain.max_corner(axis)-1;
                        if(mpi_boundary(axis)(side)) for(CELL_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()) Face_Diffusion_Helper(iterator,axis,sum_jc_ghost,Z_ghost_local,inside);}
                    else{
                        ARRAY<FACE_ITERATOR*> faces;ARRAY<CELL_ITERATOR*> cells;
                        for(int axis2=1;axis2<=TV::dimension;axis2++){if(axis==axis2 || !mpi_boundary(axis2)(side)) continue;//handled above
                            RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner(axis)++;
                            if(side==1) domain.max_corner(axis2)=domain.min_corner(axis2);
                            else{domain.max_corner(axis2)++;domain.min_corner(axis2)=domain.max_corner(axis2);}
                            for(FACE_ITERATOR iterator(grid,domain,axis2);iterator.Valid();iterator.Next()) faces.Append(new FACE_ITERATOR(iterator));}
                        RANGE<TV_INT> domain=grid.Domain_Indices(1);
                        if(side==1) domain.max_corner(axis)=domain.min_corner(axis)+1;
                        else domain.min_corner(axis)=domain.max_corner(axis)-1;
                        if(mpi_boundary(axis)(side)) for(CELL_ITERATOR iterator(grid,domain);iterator.Valid();iterator.Next()) cells.Append(new CELL_ITERATOR(iterator));
                        for(int i=faces.m;i>=1;i--){FACE_ITERATOR& iterator=*faces(i);Face_Diffusion_Helper(iterator,axis,sum_jc_ghost,Z_ghost_local,inside);}
                        for(int i=cells.m;i>=1;i--){CELL_ITERATOR& iterator=*cells(i);Face_Diffusion_Helper(iterator,axis,sum_jc_ghost,Z_ghost_local,inside);}
                        for(int i=faces.m;i>=1;i--) delete faces(i);
                        for(int i=cells.m;i>=1;i--) delete cells(i);}}
                for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) Z(iterator.Full_Index())=Z_ghost_local(iterator.Full_Index());
                for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) (*sum_jc)(iterator.Full_Index())=(*sum_jc_ghost)(iterator.Full_Index());
                boundary.Apply_Boundary_Condition_Face(grid,Z,0);
                if(sum_jc) boundary_sum->Apply_Boundary_Condition_Face(grid,*sum_jc,0);}}
        evenodd++;evenodd=evenodd%2;}
    delete sum_jc_ghost;
}

template<class T_GRID,class T2> bool DIFFUSION_UNIFORM<T_GRID,T2>::
Is_MPI_Boundary(const RANGE<TV_INT>& inside_domain,const FACE_INDEX<TV::dimension>& face){
    PHYSBAM_NOT_IMPLEMENTED();
}

template<class T_GRID,class T2> bool DIFFUSION_UNIFORM<T_GRID,T2>::
Is_MPI_Boundary(const RANGE<TV_INT>& inside_domain,const TV_INT& index){
    PHYSBAM_NOT_IMPLEMENTED();
}

template class DIFFUSION_UNIFORM<GRID<VECTOR<float,1> >,float>;
template class DIFFUSION_UNIFORM<GRID<VECTOR<float,2> >,float>;
template class DIFFUSION_UNIFORM<GRID<VECTOR<float,3> >,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DIFFUSION_UNIFORM<GRID<VECTOR<double,1> >,double>;
template class DIFFUSION_UNIFORM<GRID<VECTOR<double,2> >,double>;
template class DIFFUSION_UNIFORM<GRID<VECTOR<double,3> >,double>;
#endif


