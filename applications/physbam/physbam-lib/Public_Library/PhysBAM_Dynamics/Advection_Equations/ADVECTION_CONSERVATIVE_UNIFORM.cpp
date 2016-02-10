//#####################################################################
// Copyright 2010, Mridul Aanjaneya, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Parallel_Computation/BOUNDARY_MPI.h>
#include <PhysBAM_Tools/Parallel_Computation/MPI_UNIFORM_GRID.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_UNIFORM.h>
#include <PhysBAM_Dynamics/Advection_Equations/ADVECTION_CONSERVATIVE_UNIFORM.h>
//#define DEBUG 1

using namespace PhysBAM;

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
ADVECTION_CONSERVATIVE_UNIFORM(const T_GRID& grid,T_MPI_GRID* mpi_grid_input)
    :use_second_order(false),clamp_weights(true),number_of_ghost_cells(3),num_iterations(0),num_diffusion_iterations(0),evenodd(0),evenodd_cell(0),num_steps(1),collision_objects(0),elasticity(0),mpi_grid(mpi_grid_input),
    pls(0),density(1000),num_cells(1),find_closest_point(false),use_pls(false),diffuse_quantity(false)
{
    BOUNDARY_UNIFORM<T_GRID,T>* boundary_scalar=new BOUNDARY_UNIFORM<T_GRID,T>();boundary_scalar->Set_Fixed_Boundary(true,1);
    boundary_sum=boundary_scalar;
    if(mpi_grid) boundary_sum=new BOUNDARY_MPI<T_GRID>(mpi_grid,*boundary_scalar);
    sum_jc_cell.Resize(grid.Domain_Indices());sum_jc.Resize(grid);
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) sum_jc_cell(iterator.Cell_Index())=1;
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) sum_jc(iterator.Full_Index())=1;
    weights_to.Resize(grid,number_of_ghost_cells),weights_from.Resize(grid,number_of_ghost_cells),sum.Resize(grid),sum_jc.Resize(grid);
    weights_to_cell.Resize(grid.Domain_Indices(number_of_ghost_cells)),weights_from_cell.Resize(grid.Domain_Indices(number_of_ghost_cells)),sum_cell.Resize(grid.Domain_Indices());
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
~ADVECTION_CONSERVATIVE_UNIFORM()
{
    delete boundary_sum;
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Node(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
    const T_ARRAYS_VECTOR& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
    const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max)
{
    PHYSBAM_FATAL_ERROR();
}
template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Cell(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
    const T_ARRAYS_VECTOR& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
    const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max)
{
    PHYSBAM_FATAL_ERROR();
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Clamp_Weights(const GRID<TV>& grid,const RANGE<TV_INT>& inside_domain,ARRAY<PAIR<TV_INT,T> >& weights)
{
    T delta=(T)1e-2;
    for(int i=1;i<=weights.m;i++) if(!inside_domain.Lazy_Inside(weights(i).x)) for(int j=1;j<=TV::dimension;j++){
        TV_INT index=TV_INT::All_Ones_Vector()*2;index(j)=weights(i).x(j);
        if(!inside_domain.Lazy_Inside(index)){
            int side=1;if(index(j)>inside_domain.max_corner(j)) side=2;else assert(index(j)<inside_domain.min_corner(j));
            if(solid_walls(j)(side)) weights(i).y=0;}}
    if(collision_objects){COLLISION_GEOMETRY_ID id;
        for(int i=1;i<=weights.m;i++) if(collision_objects->Implicit_Geometry_Lazy_Inside_Any_Body(grid.Center(weights(i).x),id)) weights(i).y=0;}
    T sum=0;for(int i=1;i<=weights.m;i++) sum+=weights(i).y;
    if(sum>delta && sum!=1) for(int i=1;i<=weights.m;i++) weights(i).y/=sum;
    assert(sum==0 || sum>delta);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Clamp_Weights(const GRID<TV>& grid,const RANGE<TV_INT>& inside_domain,ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& weights)
{
    T delta=(T)1e-2;
    for(int i=1;i<=weights.m;i++) if(!inside_domain.Lazy_Inside(weights(i).x.index)) for(int j=1;j<=TV::dimension;j++){
        TV_INT index=TV_INT::All_Ones_Vector()*2;index(j)=weights(i).x.index(j);
        if(!inside_domain.Lazy_Inside(index)){
            int side=1;if(index(j)>inside_domain.max_corner(j)) side=2;else assert(index(j)<inside_domain.min_corner(j));
            if(solid_walls(j)(side)) weights(i).y=0;}}
    if(collision_objects){COLLISION_GEOMETRY_ID id;
        for(int i=1;i<=weights.m;i++) if(collision_objects->Implicit_Geometry_Lazy_Inside_Any_Body(grid.Axis_X_Face(weights(i).x),id)) weights(i).y=0;}
    T sum=0;for(int i=1;i<=weights.m;i++) sum+=weights(i).y;
    if(sum>delta && sum!=1) for(int i=1;i<=weights.m;i++) weights(i).y/=sum;
    else if(sum<delta){sum=0;for(int i=1;i<=weights.m;i++) weights(i).y=0;}
    assert(sum==0 || sum>delta);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Clamp_X(const GRID<TV>& grid,TV& X)
{
    X=grid.Clamp(X,number_of_ghost_cells);
    bool outside=false;RANGE<TV> range=grid.Domain();
    if(!range.Lazy_Inside(X)) for(int i=1;i<=TV::dimension;i++){TV point=range.min_corner;point(i)=X(i);
        if(!range.Lazy_Inside(point)){
            int side=1;if(point(i)>range.max_corner(i)) side=2;else assert(point(i)<range.min_corner(i));
            if(solid_walls(i)(side)) outside=true;}}
    if(outside) X=grid.Clamp(X);
    if(collision_objects){COLLISION_GEOMETRY_ID id;
        if(collision_objects->Implicit_Geometry_Lazy_Inside_Any_Body(X,id)){TV pre_X=X;
            X=collision_objects->bodies(id)->Implicit_Geometry_Closest_Point_On_Boundary(X,(number_of_ghost_cells+1)*grid.dX.Max());
            if(elasticity){T phi_value;
                T magnitude=(X-pre_X).Magnitude()*elasticity;
                TV normal=collision_objects->bodies(id)->Implicit_Geometry_Extended_Normal(X,phi_value).Normalized();
                X=X+magnitude*normal;}}}
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> bool ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Is_Outside(const GRID<TV>& grid,const TV_INT& cell)
{
    COLLISION_GEOMETRY_ID id;
    return (collision_objects && collision_objects->Inside_Any_Body(grid.Center(cell),0,id));
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> bool ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Is_Outside(const GRID<TV>& grid,const RANGE<TV_INT>& inside_domain,const FACE_INDEX<TV::dimension>& face)
{
    if(!inside_domain.Lazy_Inside(face.index)) for(int i=1;i<=TV::dimension;i++){
        TV_INT index=TV_INT::All_Ones_Vector()*2;index(i)=face.index(i);
        if(!inside_domain.Lazy_Inside(index)){
            int side=1;if(index(i)>inside_domain.max_corner(i)) side=2;else assert(index(i)<inside_domain.min_corner(i));
            if(solid_walls(i)(side)) return true;}}
    COLLISION_GEOMETRY_ID id;
    return (collision_objects && collision_objects->Inside_Any_Body(grid.Axis_X_Face(face),0,id));
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> bool ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Is_MPI_Boundary(const RANGE<TV_INT>& inside_domain,const FACE_INDEX<TV::dimension>& face)
{
    return Is_MPI_Boundary(inside_domain,face.index);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> bool ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Is_MPI_Boundary(const RANGE<TV_INT>& inside_domain,const TV_INT& index)
{
    if(!inside_domain.Lazy_Inside(index)) for(int i=1;i<=TV::dimension;i++){
        TV_INT tmp_index=TV_INT::All_Ones_Vector()*2;tmp_index(i)=index(i);
        if(!inside_domain.Lazy_Inside(tmp_index)){
            int side=1;if(tmp_index(i)>inside_domain.max_corner(i)) side=2;else assert(tmp_index(i)<inside_domain.min_corner(i));
            if(mpi_boundary(i)(side)) return true;}}
    return false;
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Clean_Weights(ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& weights)
{
    bool rescale=false;T delta=(T)1e-4;
    for(int i=1;i<=weights.m;i++){
        assert(weights(i).y>-delta);
        if(weights(i).y<delta){weights(i).y=0;rescale=true;}}
    if(rescale){
        T sum=0;for(int i=1;i<=weights.m;i++) sum+=weights(i).y;
        if(sum>0) for(int i=1;i<=weights.m;i++) weights(i).y/=sum;}
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Cell_Diffusion_Helper(FACE_ITERATOR& iterator,ARRAY<T,TV_INT>& sum_jc_cell,T_ARRAYS_T2& Z,ARRAY<bool,TV_INT>* inside)
{
    if(inside && !((*inside)(iterator.First_Cell_Index()) && (*inside)(iterator.Second_Cell_Index()))) return;
    if(Is_Outside(iterator.grid,iterator.First_Cell_Index()) || Is_Outside(iterator.grid,iterator.Second_Cell_Index())) return;
    if(diffuse_quantity){
        T2 Z_diff=(Z(iterator.First_Cell_Index())-Z(iterator.Second_Cell_Index()));(*reinterpret_cast<T*>(&Z_diff))*=0.5; //TODO: UNSAFE REMOVE THIS
        Z(iterator.Second_Cell_Index())+=Z_diff;Z(iterator.First_Cell_Index())-=Z_diff;
        return;}
    T wjc_diff=(sum_jc_cell(iterator.First_Cell_Index())-sum_jc_cell(iterator.Second_Cell_Index()))/(T)2.;
    T Z_diff=wjc_diff>0?wjc_diff/sum_jc_cell(iterator.First_Cell_Index()):wjc_diff/sum_jc_cell(iterator.Second_Cell_Index());
    sum_jc_cell(iterator.Second_Cell_Index())+=wjc_diff;sum_jc_cell(iterator.First_Cell_Index())-=wjc_diff;
    T2 local_Z=wjc_diff>0?Z(iterator.First_Cell_Index()):Z(iterator.Second_Cell_Index());
    Z(iterator.Second_Cell_Index())+=local_Z*(Z_diff);Z(iterator.First_Cell_Index())-=local_Z*(Z_diff);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Face_Diffusion_Helper(const GRID<TV>& grid,FACE_INDEX<TV::dimension>& first_face_index,FACE_INDEX<TV::dimension>& second_face_index,ARRAY<T,FACE_INDEX<TV::dimension> >& sum_jc,T_FACE_ARRAYS_SCALAR& Z,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside)
{
    if(inside && !((*inside)(first_face_index) && (*inside)(second_face_index))) return;
    if(Is_Outside(grid,grid.Domain_Indices(),first_face_index) || Is_Outside(grid,grid.Domain_Indices(),second_face_index)) return;
    if(sum_jc(first_face_index)==0 && sum_jc(second_face_index)==0) return;
    T wjc_diff=(sum_jc(first_face_index)-sum_jc(second_face_index))/(T)2.;
    T Z_diff=wjc_diff>0?wjc_diff/sum_jc(first_face_index):wjc_diff/sum_jc(second_face_index);
    sum_jc(second_face_index)+=wjc_diff;sum_jc(first_face_index)-=wjc_diff;
    T local_Z=wjc_diff>0?Z(first_face_index):Z(second_face_index);
    Z(second_face_index)+=local_Z*(Z_diff);Z(first_face_index)-=local_Z*(Z_diff);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Face_Diffusion_Helper(FACE_ITERATOR& iterator,int axis,ARRAY<T,FACE_INDEX<TV::dimension> >& sum_jc,T_FACE_ARRAYS_SCALAR& Z,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside)
{
    assert(axis!=iterator.Axis());
    FACE_INDEX<TV::dimension> first_face_index=FACE_INDEX<TV::dimension>(axis,iterator.First_Cell_Index()),second_face_index=FACE_INDEX<TV::dimension>(axis,iterator.Second_Cell_Index());
    Face_Diffusion_Helper(iterator.grid,first_face_index,second_face_index,sum_jc,Z,inside);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Face_Diffusion_Helper(CELL_ITERATOR& iterator,int axis,ARRAY<T,FACE_INDEX<TV::dimension> >& sum_jc,T_FACE_ARRAYS_SCALAR& Z,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside)
{
    FACE_INDEX<TV::dimension> first_face_index=FACE_INDEX<TV::dimension>(axis,iterator.First_Face_Index(axis)),second_face_index=FACE_INDEX<TV::dimension>(axis,iterator.Second_Face_Index(axis));
    Face_Diffusion_Helper(iterator.grid,first_face_index,second_face_index,sum_jc,Z,inside);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Cell_Diffusion(const T_GRID& grid,ARRAY<T,TV_INT>& sum_jc_cell,T_ARRAYS_T2& Z,T_BOUNDARY_T2& boundary,BOUNDARY_UNIFORM<T_GRID,T>* boundary_sum,ARRAY<bool,TV_INT>* inside)
{
    T_ARRAYS_T2 Z_ghost_local(grid.Domain_Indices(1));
    ARRAY<T,TV_INT> sum_jc_cell_ghost(grid.Domain_Indices(1));
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
                boundary_sum->Fill_Ghost_Cells(grid,sum_jc_cell,sum_jc_cell_ghost,0,0,1); //Sync for mpi boundaries
                for(int side=1;side<=2;side++) if(mpi_boundary(axis)(side)){
                    if(evenodd_cell==0){
                        for(FACE_ITERATOR iterator(grid,0,GRID<TV>::BOUNDARY_REGION,2*(axis-1)+side);iterator.Valid();iterator.Next()) Cell_Diffusion_Helper(iterator,sum_jc_cell_ghost,Z_ghost_local,inside);}                        
                    else{
                        ARRAY<FACE_ITERATOR*> faces;
                        for(FACE_ITERATOR iterator(grid,0,GRID<TV>::BOUNDARY_REGION,2*(axis-1)+side);iterator.Valid();iterator.Next()) faces.Append(new FACE_ITERATOR(iterator));
                        for(int i=faces.m;i>=1;i--){FACE_ITERATOR& iterator=*faces(i);Cell_Diffusion_Helper(iterator,sum_jc_cell_ghost,Z_ghost_local,inside);}
                        for(int i=faces.m;i>=1;i--) delete faces(i);}}
                for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) Z(iterator.Cell_Index())=Z_ghost_local(iterator.Cell_Index());
                for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) sum_jc_cell(iterator.Cell_Index())=sum_jc_cell_ghost(iterator.Cell_Index());}}
        evenodd_cell++;evenodd_cell=evenodd_cell%2;}
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Face_Diffusion(const T_GRID& grid,ARRAY<T,FACE_INDEX<TV::dimension> >& sum_jc,T_FACE_ARRAYS_SCALAR& Z,T_BOUNDARY& boundary,BOUNDARY_UNIFORM<T_GRID,T>* boundary_sum,ARRAY<bool,FACE_INDEX<TV::dimension> >* inside)
{
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Before diffusion",0,0);
    T_FACE_ARRAYS_SCALAR Z_ghost_local(grid,1);
    ARRAY<T,FACE_INDEX<TV::dimension> > sum_jc_ghost(grid,1);
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
                boundary_sum->Apply_Boundary_Condition_Face(grid,sum_jc,0);
                boundary_sum->Fill_Ghost_Cells_Face(grid,sum_jc,sum_jc_ghost,0,1); //Sync for mpi boundaries
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
                for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) sum_jc(iterator.Full_Index())=sum_jc_ghost(iterator.Full_Index());
                boundary.Apply_Boundary_Condition_Face(grid,Z,0);
                boundary_sum->Apply_Boundary_Condition_Face(grid,sum_jc,0);}}
        evenodd++;evenodd=evenodd%2;}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After diffusion",0,0);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
    const T_FACE_LOOKUP& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
    const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max)
{
    PHYSBAM_ASSERT(!Z_min || !Z_max); //we don't support extrema yet
    weights_to_cell.Resize(grid.Domain_Indices(number_of_ghost_cells)),weights_from_cell.Resize(grid.Domain_Indices(number_of_ghost_cells));
    for(CELL_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        weights_to_cell(cell).Resize(0);weights_from_cell(cell).Resize(0);}
     
    LOG::Time("Backwards");
    COLLISION_GEOMETRY_ID id;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(Is_Outside(grid,cell)) continue;
        int local_num_steps=(collision_objects && collision_objects->Inside_Any_Body(grid.Center(cell),(T)number_of_ghost_cells*grid.dX.Max()/(T)2.,id))?num_steps:1;
        TV X=iterator.Location()-(dt/local_num_steps)*averaging.Face_To_Cell_Vector(grid,cell,face_velocities);
        for(int i=2;i<=local_num_steps;i++) X-=(dt/local_num_steps)*interpolation.Clamped_To_Array_Face(grid,face_velocities,X);
        if(use_second_order) X+=(iterator.Location()-(X+dt*interpolation.Clamped_To_Array_Face(grid,face_velocities,X)))/2.;
        Clamp_X(grid,X);
        ARRAY<PAIR<TV_INT,T> > backwards_weights=interpolation.Clamped_To_Array_Weights(grid,Z_ghost,X);
        Clamp_Weights(grid,grid.Domain_Indices(),backwards_weights);
        for(int i=1;i<=backwards_weights.m;i++){assert(backwards_weights(i).y>-1e-2);if(backwards_weights(i).y<0) backwards_weights(i).y=0;}
        for(int i=1;i<=backwards_weights.m;i++){if(backwards_weights(i).y==0) continue;
            weights_to_cell(backwards_weights(i).x).Append(PAIR<TV_INT,T>(cell,backwards_weights(i).y));
            weights_from_cell(cell).Append(PAIR<TV_INT,int>(backwards_weights(i).x,weights_to_cell(backwards_weights(i).x).m));}}
    LOG::Time("Backwards MPI");
    if(mpi_grid) mpi_grid->Sync_Common_Cell_Weights_To(weights_to_cell,weights_from_cell,number_of_ghost_cells);

    LOG::Time("Sum");
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        ARRAY<PAIR<TV_INT,T> >& local_weights=weights_to_cell(cell);
        sum_cell(cell)=0;for(int i=1;i<=local_weights.m;i++) sum_cell(cell)+=local_weights(i).y;}

    LOG::Time("Forwards");
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        if(Is_Outside(grid,cell)) continue;
        if(clamp_weights) if(sum_cell(cell)>=1) continue;
        T remaining=(1-sum_cell(cell));
        int local_num_steps=(collision_objects && collision_objects->Inside_Any_Body(grid.Center(cell),(T)number_of_ghost_cells*grid.dX.Max()/(T)2.,id))?num_steps:1;
        TV X=iterator.Location()+(dt/local_num_steps)*averaging.Face_To_Cell_Vector(grid,cell,face_velocities);
        for(int i=2;i<=local_num_steps;i++) X+=(dt/local_num_steps)*interpolation.Clamped_To_Array_Face(grid,face_velocities,X);
        if(use_second_order) X+=(iterator.Location()-(X-dt*interpolation.Clamped_To_Array_Face(grid,face_velocities,X)))/2.;
        Clamp_X(grid,X);
        ARRAY<PAIR<TV_INT,T> > forward_weights=interpolation.Clamped_To_Array_Weights(grid,Z_ghost,X);
        Clamp_Weights(grid,grid.Domain_Indices(),forward_weights);
        for(int i=1;i<=forward_weights.m;i++){assert(forward_weights(i).y>-1e-2);if(forward_weights(i).y<0) forward_weights(i).y=0;}
        for(int i=1;i<=forward_weights.m;i++){if(forward_weights(i).y==0) continue;
            int index=0;for(int j=1;j<=weights_to_cell(cell).m;j++) if(weights_to_cell(cell)(j).x==forward_weights(i).x) index=j;
            if(index) weights_to_cell(cell)(index).y+=forward_weights(i).y*remaining;
            else{
                weights_to_cell(cell).Append(PAIR<TV_INT,T>(forward_weights(i).x,forward_weights(i).y*remaining));
                weights_from_cell(forward_weights(i).x).Append(PAIR<TV_INT,int>(cell,weights_to_cell(cell).m));}}}
    if(mpi_grid) mpi_grid->Sync_Common_Cell_Weights_From(weights_to_cell,weights_from_cell,number_of_ghost_cells);

    LOG::Time("Clamp");
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        ARRAY<PAIR<TV_INT,T> >& local_weights=weights_to_cell(cell);
        T sum=0;for(int i=1;i<=local_weights.m;i++) sum+=local_weights(i).y;
        if(sum) for(int i=1;i<=local_weights.m;i++){assert(sum>1e-6);
            local_weights(i).y=local_weights(i).y/sum;}}
    if(mpi_grid) mpi_grid->Sync_Common_Cell_Weights_From(weights_to_cell,weights_from_cell,number_of_ghost_cells);

#ifdef DEBUG
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        ARRAY<PAIR<TV_INT,T> >& local_weights=weights_to_cell(cell);
        T sum=0;for(int i=1;i<=local_weights.m;i++) sum+=local_weights(i).y;}
#endif

#ifdef DEBUG
    T_ARRAYS_T2 face_back(grid.Domain_Indices());
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT face=iterator.Cell_Index();
        face_back(face)=Z(face);}
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT face=iterator.Cell_Index();
        Z(face)=T2();Z(face)+=sum_jc_cell(face);}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Sum jc cell 1",0,0);
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT face=iterator.Cell_Index();
        Z(face)=face_back(face);}
#endif

    LOG::Time("Iterations");
    ARRAY<T,TV_INT> sum_ic(grid.Domain_Indices(number_of_ghost_cells));
    if(num_iterations || num_diffusion_iterations) boundary_sum->Fill_Ghost_Cells(grid,sum_jc_cell,sum_ic,dt,time,number_of_ghost_cells);
    for(int i=1;i<=num_iterations;i++){
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            ARRAY<PAIR<TV_INT,int> >& local_weights=weights_from_cell(cell);
            T sum=0;for(int i=1;i<=local_weights.m;i++) sum+=sum_ic(local_weights(i).x)*weights_to_cell(local_weights(i).x)(local_weights(i).y).y;
            for(int i=1;i<=local_weights.m;i++){assert(sum_jc_cell(cell)>1e-2);
                weights_to_cell(local_weights(i).x)(local_weights(i).y).y=weights_to_cell(local_weights(i).x)(local_weights(i).y).y/sum;}}
        if(mpi_grid) mpi_grid->Sync_Common_Cell_Weights_To(weights_to_cell,weights_from_cell,number_of_ghost_cells);
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
            ARRAY<PAIR<TV_INT,T> >& local_weights=weights_to_cell(cell);
            T sum=0;for(int i=1;i<=local_weights.m;i++) sum+=local_weights(i).y;
            for(int i=1;i<=local_weights.m;i++){assert(sum>1e-2);
                local_weights(i).y=local_weights(i).y/sum;}}
        if(mpi_grid) mpi_grid->Sync_Common_Cell_Weights_From(weights_to_cell,weights_from_cell,number_of_ghost_cells);}

    LOG::Time("Mapping");
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){Z(iterator.Cell_Index())=T2();}
    for(CELL_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        ARRAY<PAIR<TV_INT,T> >& local_weights=weights_to_cell(cell);
        for(int i=1;i<=local_weights.m;i++)
            if(grid.Domain_Indices().Lazy_Inside(local_weights(i).x)) Z(local_weights(i).x)+=local_weights(i).y*Z_ghost(cell);}

    LOG::Time("Set up diffusion");
    if(num_diffusion_iterations) for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell=iterator.Cell_Index();
        ARRAY<PAIR<TV_INT,int> >& local_weights=weights_from_cell(cell);
        sum_jc_cell(cell)=0;for(int i=1;i<=local_weights.m;i++){
            if(!use_pls || phi1(local_weights(i).x)<=0) sum_jc_cell(cell)+=sum_ic(local_weights(i).x)*weights_to_cell(local_weights(i).x)(local_weights(i).y).y;}}

#ifdef DEBUG
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT face=iterator.Cell_Index();
        face_back(face)=Z(face);}
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT face=iterator.Cell_Index();
        Z(face)=T2();Z(face)+=sum_jc_cell(face);}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Sum jc cell",0,0);
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT face=iterator.Cell_Index();
        Z(face)=face_back(face);}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Before Diffusion",0,0);
#endif

    LOG::Time("Diffusion");
    if(!use_pls && num_diffusion_iterations) Cell_Diffusion(grid,sum_jc_cell,Z,boundary,boundary_sum);

#ifdef DEBUG
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT face=iterator.Cell_Index();
        face_back(face)=Z(face);}
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT face=iterator.Cell_Index();
        Z(face)=T2();Z(face)+=sum_jc_cell(face);}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Sum jc cell final",0,0);
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT face=iterator.Cell_Index();
        Z(face)=face_back(face);}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After Diffusion",0,0);
#endif
}
 
template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Face_Lookup_PLS(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
    const T_FACE_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
    const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max)
{
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Before Advection",0,0);
    const T delta=(T)1e-2;
    T_LEVELSET lsv1(const_cast<T_GRID&>(grid),phi1),lsv2(const_cast<T_GRID&>(grid),phi2);
    T_FACE_ARRAYS_BOOL inside1(grid,number_of_ghost_cells),inside2(grid,number_of_ghost_cells);
    T_FACE_ARRAYS_SCALAR phiface1(grid,number_of_ghost_cells),phiface2(grid,number_of_ghost_cells);
    for(FACE_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){
        phiface1(iterator.Full_Index())=((T)0.5)*((phi1.Valid_Index(iterator.First_Cell_Index())?phi1(iterator.First_Cell_Index()):phi1(iterator.Second_Cell_Index()))+
                                            (phi1.Valid_Index(iterator.Second_Cell_Index())?phi1(iterator.Second_Cell_Index()):phi1(iterator.First_Cell_Index())));
        phiface2(iterator.Full_Index())=((T)0.5)*((phi2.Valid_Index(iterator.First_Cell_Index())?phi2(iterator.First_Cell_Index()):phi2(iterator.Second_Cell_Index()))+
                                            (phi2.Valid_Index(iterator.Second_Cell_Index())?phi2(iterator.Second_Cell_Index()):phi2(iterator.First_Cell_Index())));}
    for(FACE_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){inside1(iterator.Full_Index())=phiface1(iterator.Full_Index())<0;inside2(iterator.Full_Index())=phiface2(iterator.Full_Index())<0;}
//    for(CELL_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){
//        for(int axis=1;axis<=TV::dimension;axis++){
//            if(!inside1(FACE_INDEX<TV::dimension>(axis,iterator.First_Face_Index(axis)))) inside1(FACE_INDEX<TV::dimension>(axis,iterator.First_Face_Index(axis)))=lsv1.Phi(iterator.Location())<=0;
//            if(!inside1(FACE_INDEX<TV::dimension>(axis,iterator.Second_Face_Index(axis)))) inside1(FACE_INDEX<TV::dimension>(axis,iterator.Second_Face_Index(axis)))=lsv1.Phi(iterator.Location())<=0;
//            if(!inside2(FACE_INDEX<TV::dimension>(axis,iterator.First_Face_Index(axis)))) inside2(FACE_INDEX<TV::dimension>(axis,iterator.First_Face_Index(axis)))=lsv2.Phi(iterator.Location())<=0;
//            if(!inside2(FACE_INDEX<TV::dimension>(axis,iterator.Second_Face_Index(axis)))) inside2(FACE_INDEX<TV::dimension>(axis,iterator.Second_Face_Index(axis)))=lsv2.Phi(iterator.Location())<=0;}}
    PHYSBAM_ASSERT(!Z_min || !Z_max); //we don't support extrema yet
    momentum_lost.Resize(grid.Domain_Indices(number_of_ghost_cells+3));
    for(CELL_ITERATOR iterator(grid,number_of_ghost_cells+3);iterator.Valid();iterator.Next()) momentum_lost(iterator.Cell_Index())=0;
    weights_to.Resize(grid,number_of_ghost_cells),weights_from.Resize(grid,number_of_ghost_cells);
    for(FACE_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        weights_to(face).Resize(0);weights_from(face).Resize(0);}
    LOG::Time("Face forwards");
    COLLISION_GEOMETRY_ID id;
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner+=TV_INT::Axis_Vector(face.axis);
        if(!inside2(iterator.Full_Index())) continue; //Don't get velocity if not inside lsv
        if(Is_Outside(grid,domain,face)) continue;
        int local_num_steps=(collision_objects && collision_objects->Inside_Any_Body(grid.Axis_X_Face(face),(T)number_of_ghost_cells*grid.dX.Max()/(T)2.,id))?num_steps:1;
        TV X=iterator.Location()-(dt/local_num_steps)*averaging.Face_To_Face_Vector(grid,face.axis,face.index,face_velocities);
        for(int i=2;i<=local_num_steps;i++) X-=(dt/local_num_steps)*interpolation.Clamped_To_Array_Face(grid,face_velocities,X);
        if(use_second_order) X+=(iterator.Location()-(X+dt*interpolation.Clamped_To_Array_Face(grid,face_velocities,X)))/2.;
        Clamp_X(grid,X);
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> > backwards_weights=interpolation.Clamped_To_Array_Face_Component_Weights(face.axis,grid,Z_ghost.Starting_Point_Face(face.axis,face.index),X);
        for(int i=1;i<=backwards_weights.m;i++) if(!inside1(backwards_weights(i).x)) backwards_weights(i).y=0;
        Clamp_Weights(grid,domain,backwards_weights);
        Clean_Weights(backwards_weights);
        T sum=0;for(int i=1;i<=backwards_weights.m;i++) sum+=backwards_weights(i).y;
        assert(sum<=1+delta);
        if(sum<1 && sum>=delta) for(int i=1;i<=backwards_weights.m;i++) backwards_weights(i).y=backwards_weights(i).y/sum;
        else if(sum<delta){
            TV axis_vector=TV::Axis_Vector(iterator.Axis());
            TV_INT cell1=iterator.First_Cell_Index(),cell2=iterator.Second_Cell_Index();
            TV X_center1=grid.Center(cell1)-dt*averaging.Face_To_Cell_Vector(grid,cell1,face_velocities),X_center2=grid.Center(cell2)-dt*averaging.Face_To_Cell_Vector(grid,cell2,face_velocities);
            TV X_cell1=X_center1+grid.dX/2.*axis_vector,X_cell2=X_center2-grid.dX/2.*axis_vector;
            assert(lsv2.Phi(grid.Center(cell1))<=0||lsv2.Phi(grid.Center(cell2))<=0);
            //assert(lsv1.Phi(X_center1)<=0||lsv1.Phi(X_center2)<=0);
            if(lsv1.Phi(X_cell1)>0){
                if(lsv1.Phi(X_center1)<0) X_cell1=X_center1+grid.dX/2.*axis_vector*lsv1.Phi(X_center1)/(lsv1.Phi(X_center1)-lsv1.Phi(X_cell1));
                else if(lsv1.Phi(X_center1)<lsv1.Phi(X_cell1)) X_cell1=X_center1;}
            if(lsv1.Phi(X_cell2)>0){
                if(lsv1.Phi(X_center2)<0) X_cell2=X_center2-grid.dX/2.*axis_vector*lsv1.Phi(X_center2)/(lsv1.Phi(X_center2)-lsv1.Phi(X_cell2));
                else if(lsv1.Phi(X_center2)<lsv1.Phi(X_cell2)) X_cell2=X_center2;}
            assert(!use_second_order);
            Clamp_X(grid,X_cell1);
            Clamp_X(grid,X_cell2);
            backwards_weights=interpolation.Clamped_To_Array_Face_Component_Weights(face.axis,grid,Z_ghost.Starting_Point_Face(face.axis,face.index),X_cell1);
            ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> > backwards_weights2=interpolation.Clamped_To_Array_Face_Component_Weights(face.axis,grid,Z_ghost.Starting_Point_Face(face.axis,face.index),X_cell2);
            for(int i=1;i<=backwards_weights.m;i++) backwards_weights(i).y*=0.5;
            for(int i=1;i<=backwards_weights2.m;i++) backwards_weights2(i).y*=0.5;
            for(int i=1;i<=backwards_weights2.m;i++){
                int index=0;for(int j=1;j<=backwards_weights.m;j++) if(backwards_weights(j).x==backwards_weights2(i).x) index=j;
                if(index) backwards_weights(index).y+=backwards_weights2(i).y;
                else backwards_weights.Append(backwards_weights2(i));}
            for(int i=1;i<=backwards_weights.m;i++) if(!inside1(backwards_weights(i).x)) backwards_weights(i).y=0;
            Clamp_Weights(grid,domain,backwards_weights); 
            sum=0;for(int i=1;i<=backwards_weights.m;i++) sum+=backwards_weights(i).y;
            assert(sum<=1+delta);
            if(sum<1 && sum>=delta) for(int i=1;i<=backwards_weights.m;i++) backwards_weights(i).y=backwards_weights(i).y/sum;
            else if(sum<delta){
                for(int i=1;i<=backwards_weights.m;i++) backwards_weights(i).y=0;
                PHYSBAM_WARNING("Cannot find backward weights for this cell");}}
        //sum=0;for(int i=1;i<=backwards_weights.m;i++) sum+=backwards_weights(i).y;
        //assert((sum>delta));
        for(int i=1;i<=backwards_weights.m;i++){assert(backwards_weights(i).y>-delta);if(backwards_weights(i).y<0) backwards_weights(i).y=0;}
        for(int i=1;i<=backwards_weights.m;i++){if(backwards_weights(i).y==0) continue;
            weights_to(backwards_weights(i).x).Append(PAIR<FACE_INDEX<TV::dimension>,T>(face,backwards_weights(i).y));
            weights_from(face).Append(PAIR<FACE_INDEX<TV::dimension>,int>(backwards_weights(i).x,weights_to(backwards_weights(i).x).m));}}
    if(mpi_grid) mpi_grid->Sync_Common_Face_Weights_To(weights_to,weights_from,number_of_ghost_cells);

    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& local_weights=weights_to(face);
        sum(face)=0;for(int i=1;i<=local_weights.m;i++) sum(face)+=local_weights(i).y;}

    LOG::Time("Face backwards");
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        if(clamp_weights) if(sum(face)>=1) continue;
        if(!inside1(iterator.Full_Index())) continue; //Don't get velocity if not inside lsv
        RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner+=TV_INT::Axis_Vector(face.axis);
        T remaining=(1-sum(face));
        if(Is_Outside(grid,domain,face)) continue;
        int local_num_steps=(collision_objects && collision_objects->Inside_Any_Body(grid.Axis_X_Face(face),(T)number_of_ghost_cells*grid.dX.Max()/(T)2.,id))?num_steps:1;
        TV X=iterator.Location()+(dt/local_num_steps)*averaging.Face_To_Face_Vector(grid,face.axis,face.index,face_velocities);
        for(int i=2;i<=local_num_steps;i++) X+=(dt/local_num_steps)*interpolation.Clamped_To_Array_Face(grid,face_velocities,X);
        if(use_second_order) X+=(iterator.Location()-(X-dt*interpolation.Clamped_To_Array_Face(grid,face_velocities,X)))/2.;
        Clamp_X(grid,X);
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> > forward_weights=interpolation.Clamped_To_Array_Face_Component_Weights(face.axis,grid,Z_ghost.Starting_Point_Face(face.axis,face.index),X);
        for(int i=1;i<=forward_weights.m;i++) if(!inside2(forward_weights(i).x)) forward_weights(i).y=0;
        Clamp_Weights(grid,domain,forward_weights);
        Clean_Weights(forward_weights);
        T sum=0;for(int i=1;i<=forward_weights.m;i++) sum+=forward_weights(i).y;
        assert(sum<=1+delta);
        if(sum<1 && sum>=delta) for(int i=1;i<=forward_weights.m;i++) forward_weights(i).y=forward_weights(i).y/sum;
        else if(sum<delta){int iterations=100,iter=0;
            for(int i=1;i<=forward_weights.m;i++) forward_weights(i).y=0;
            TV X_start=X;
            if(!(lsv2.Phi(X)>num_cells*grid.dX(face.axis) && !find_closest_point && pls && pls->Fix_Momentum_With_Escaped_Particles(X,face_velocities.Starting_Point_Face(face.axis,face.index).V_face,density*grid.dX.Product()*remaining,1.5,1,time,false))){
                while(lsv2.Phi(X)>0 && iter<iterations){iter++;X=X-(T)(lsv2.Phi(X))*lsv2.Extended_Normal(X).Normalized();Clamp_X(grid,X);}
                forward_weights=interpolation.Clamped_To_Array_Face_Component_Weights(face.axis,grid,Z_ghost.Starting_Point_Face(face.axis,face.index),X);
                for(int i=1;i<=forward_weights.m;i++) if(!inside2(forward_weights(i).x)) forward_weights(i).y=0;
                sum=0;for(int i=1;i<=forward_weights.m;i++) sum+=forward_weights(i).y;
                assert(sum<=1+delta);
                if(sum<1 && sum>=delta) for(int i=1;i<=forward_weights.m;i++) forward_weights(i).y=forward_weights(i).y/sum;
                else if(sum<delta){X=X_start;iter=0;
                    while(lsv2.Phi(X)>0 && iter<iterations){iter++;X=X+(T)(lsv2.Phi(X))*lsv2.Extended_Normal(X).Normalized();Clamp_X(grid,X);}
                    forward_weights=interpolation.Clamped_To_Array_Face_Component_Weights(face.axis,grid,Z_ghost.Starting_Point_Face(face.axis,face.index),X);
                    for(int i=1;i<=forward_weights.m;i++) if(!inside2(forward_weights(i).x)) forward_weights(i).y=0;
                    sum=0;for(int i=1;i<=forward_weights.m;i++) sum+=forward_weights(i).y;
                    assert(sum<=1+delta);
                    if(sum<1 && sum>=delta) for(int i=1;i<=forward_weights.m;i++) forward_weights(i).y=forward_weights(i).y/sum;
                    if(sum<delta){
                        weights_to(face).Append(PAIR<FACE_INDEX<TV::dimension>,T>(FACE_INDEX<TV::dimension>(),remaining));
                        PHYSBAM_WARNING("Cannot find LSV or Particle for Momentum");
                        for(int i=1;i<=forward_weights.m;i++) forward_weights(i).y=0;
                        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> > forward_weights_local=interpolation.Clamped_To_Array_Face_Component_Weights(face.axis,grid,Z_ghost.Starting_Point_Face(face.axis,face.index),X_start);
                        for(int i=1;i<=forward_weights_local.m;i++){
                            momentum_lost(forward_weights_local(i).x.index)+=density*grid.dX.Product()*forward_weights_local(i).y*remaining/(T)2.;
                            momentum_lost(forward_weights_local(i).x.index+TV_INT::Axis_Vector(forward_weights_local(i).x.axis))+=density*grid.dX.Product()*forward_weights_local(i).y*remaining/(T)2.;}}}}
            else weights_to(face).Append(PAIR<FACE_INDEX<TV::dimension>,T>(FACE_INDEX<TV::dimension>(),remaining));}
        Clamp_Weights(grid,domain,forward_weights);
        for(int i=1;i<=forward_weights.m;i++){if(forward_weights(i).y<0) forward_weights(i).y=0;}
        for(int i=1;i<=forward_weights.m;i++){if(forward_weights(i).y==0) continue;
            int index=0;for(int j=1;j<=weights_to(face).m;j++) if(weights_to(face)(j).x==forward_weights(i).x) index=j;
            if(index) weights_to(face)(index).y+=forward_weights(i).y*remaining;
            else{
                weights_to(face).Append(PAIR<FACE_INDEX<TV::dimension>,T>(forward_weights(i).x,forward_weights(i).y*remaining));
                weights_from(forward_weights(i).x).Append(PAIR<FACE_INDEX<TV::dimension>,int>(face,weights_to(face).m));}}}
    if(mpi_grid) mpi_grid->Sync_Common_Face_Weights_From(weights_to,weights_from,number_of_ghost_cells);

    LOG::Time("Face Clamp");
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        if(!inside1(iterator.Full_Index())) continue; //Don't clamp weights for outside cells
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& local_weights=weights_to(face);
        T sum=0;for(int i=1;i<=local_weights.m;i++) sum+=local_weights(i).y;
        if(sum) for(int i=1;i<=local_weights.m;i++){assert(sum>delta);local_weights(i).y=local_weights(i).y/sum;}}
    if(mpi_grid){
        mpi_grid->ignore_boundary_faces=true;
        mpi_grid->Sync_Common_Face_Weights_From(weights_to,weights_from,number_of_ghost_cells);
        mpi_grid->ignore_boundary_faces=false;}

#ifdef DEBUG
    T_FACE_ARRAYS_SCALAR face_back(grid);
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        face_back(face)=Z(face);}
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        Z(face)=phiface1(face);}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Phi1",0,0);
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        Z(face)=phiface2(face);}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Phi2",0,0);
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        Z(face)=face_back(face);}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After Weights",0,0);
#endif
    
#ifdef DEBUG
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        face_back(face)=Z(face);}
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& local_weights=weights_to(face);
        T sum=0;for(int i=1;i<=local_weights.m;i++) sum+=local_weights(i).y;
        Z(face)=sum;}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("Weights w",0,0);
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        Z(face)=face_back(face);}
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After Weights",0,0);
#endif

    LOG::Time("Face Iter");
    ARRAY<T,FACE_INDEX<TV::dimension> > sum_ic(grid,number_of_ghost_cells,false);
    if(num_iterations || num_diffusion_iterations) boundary_sum->Fill_Ghost_Cells_Face(grid,sum_jc,sum_ic,time,number_of_ghost_cells);
    for(int iter=1;iter<=num_iterations;iter++){
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
            if(!inside2(face)) continue; //Don't clamp weights for outside cells
            ARRAY<PAIR<FACE_INDEX<TV::dimension>,int> >& local_weights=weights_from(face);
            T sum=0;for(int i=1;i<=local_weights.m;i++) sum+=sum_ic(local_weights(i).x)*weights_to(local_weights(i).x)(local_weights(i).y).y;
            for(int i=1;i<=local_weights.m;i++){assert(sum>delta);
                weights_to(local_weights(i).x)(local_weights(i).y).y=weights_to(local_weights(i).x)(local_weights(i).y).y/sum;}}
        if(mpi_grid){
            mpi_grid->ignore_boundary_faces=true;
            mpi_grid->Sync_Common_Face_Weights_To(weights_to,weights_from,number_of_ghost_cells);
            mpi_grid->ignore_boundary_faces=false;}
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
            if(!inside1(face)) continue; //Don't clamp weights for outside cells
            ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& local_weights=weights_to(face);
            T sum=0;for(int i=1;i<=local_weights.m;i++) sum+=local_weights(i).y;
            for(int i=1;i<=local_weights.m;i++){assert(sum>delta);
                local_weights(i).y=local_weights(i).y/sum;}}
        if(mpi_grid){
            mpi_grid->ignore_boundary_faces=true;
            mpi_grid->Sync_Common_Face_Weights_From(weights_to,weights_from,number_of_ghost_cells);
            mpi_grid->ignore_boundary_faces=false;}}

    LOG::Time("Face Diffusion setup");
    if(num_diffusion_iterations) for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        if(!inside2(face)){sum_jc(face)=1;continue;} //Don't clamp weights for outside cells
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,int> >& local_weights=weights_from(face);
        sum_jc(face)=0;for(int i=1;i<=local_weights.m;i++) sum_jc(face)+=sum_ic(local_weights(i).x)*weights_to(local_weights(i).x)(local_weights(i).y).y;}

    LOG::Time("Face Map");
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){if(!inside2(iterator.Full_Index())) continue;Z(iterator.Full_Index())=T();}
    for(FACE_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        if(!inside1(iterator.Full_Index())) continue; //Don't get velocity if not inside lsv
        RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner+=TV_INT::Axis_Vector(face.axis);
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& local_weights=weights_to(face);
        for(int i=1;i<=local_weights.m;i++) if(domain.Lazy_Inside(local_weights(i).x.index)) Z(local_weights(i).x)+=local_weights(i).y*Z_ghost(face);}

    LOG::Time("Face Diffusion");
    if(num_diffusion_iterations) Face_Diffusion(grid,sum_jc,Z,boundary,boundary_sum,&inside2);
    PHYSBAM_DEBUG_WRITE_SUBSTEP("After Advection",0,0);
}

template<class T_GRID,class T2,class T_AVERAGING,class T_INTERPOLATION> void ADVECTION_CONSERVATIVE_UNIFORM<T_GRID,T2,T_AVERAGING,T_INTERPOLATION>::
Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
    const T_FACE_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
    const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max)
{
    if(use_pls){Update_Advection_Equation_Face_Lookup_PLS(grid,Z,Z_ghost,face_velocities,boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max);return;}
    PHYSBAM_ASSERT(!Z_min || !Z_max); //we don't support extrema yet
    weights_to.Resize(grid,number_of_ghost_cells),weights_from.Resize(grid,number_of_ghost_cells);
    for(FACE_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        weights_to(face).Resize(0);weights_from(face).Resize(0);}

    COLLISION_GEOMETRY_ID id;
    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner+=TV_INT::Axis_Vector(face.axis);
        int local_num_steps=(collision_objects && collision_objects->Inside_Any_Body(grid.Axis_X_Face(face),(T)number_of_ghost_cells*grid.dX.Max()/(T)2.,id))?num_steps:1;
        TV X=iterator.Location()-(dt/local_num_steps)*averaging.Face_To_Face_Vector(grid,face.axis,face.index,face_velocities);
        for(int i=2;i<=local_num_steps;i++) X-=(dt/local_num_steps)*interpolation.Clamped_To_Array_Face(grid,face_velocities,X);
        if(use_second_order) X+=(iterator.Location()-(X+dt*interpolation.Clamped_To_Array_Face(grid,face_velocities,X)))/2.;
        Clamp_X(grid,X);
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> > backwards_weights=interpolation.Clamped_To_Array_Face_Component_Weights(face.axis,grid,Z_ghost.Starting_Point_Face(face.axis,face.index),X);
        Clamp_Weights(grid,domain,backwards_weights);
        Clean_Weights(backwards_weights);
        for(int i=1;i<=backwards_weights.m;i++){if(backwards_weights(i).y==0) continue;
            weights_to(backwards_weights(i).x).Append(PAIR<FACE_INDEX<TV::dimension>,T>(face,backwards_weights(i).y));
            weights_from(face).Append(PAIR<FACE_INDEX<TV::dimension>,int>(backwards_weights(i).x,weights_to(backwards_weights(i).x).m));}}
    if(mpi_grid) mpi_grid->Sync_Common_Face_Weights_To(weights_to,weights_from,number_of_ghost_cells);
    LOG::Time("After SL");

    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& local_weights=weights_to(face);
        sum(face)=0;for(int i=1;i<=local_weights.m;i++) sum(face)+=local_weights(i).y;}
    LOG::Time("After Sum");

    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        if(clamp_weights) if(sum(face)>=1) continue;
        RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner+=TV_INT::Axis_Vector(face.axis);
        T remaining=(1-sum(face));
        int local_num_steps=(collision_objects && collision_objects->Inside_Any_Body(grid.Axis_X_Face(face),(T)number_of_ghost_cells*grid.dX.Max()/(T)2.,id))?num_steps:1;
        TV X=iterator.Location()+(dt/local_num_steps)*averaging.Face_To_Face_Vector(grid,face.axis,face.index,face_velocities);
        for(int i=2;i<=local_num_steps;i++) X+=(dt/local_num_steps)*interpolation.Clamped_To_Array_Face(grid,face_velocities,X);
        if(use_second_order) X+=(iterator.Location()-(X-dt*interpolation.Clamped_To_Array_Face(grid,face_velocities,X)))/2.;
        Clamp_X(grid,X);
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> > forward_weights=interpolation.Clamped_To_Array_Face_Component_Weights(face.axis,grid,Z_ghost.Starting_Point_Face(face.axis,face.index),X);
        Clamp_Weights(grid,domain,forward_weights);
        Clean_Weights(forward_weights);
        for(int i=1;i<=forward_weights.m;i++){
            int index=0;for(int j=1;j<=weights_to(face).m;j++) if(weights_to(face)(j).x==forward_weights(i).x) index=j;
            if(index) weights_to(face)(index).y+=forward_weights(i).y*remaining;
            else{if(forward_weights(i).y==0) continue;
                weights_to(face).Append(PAIR<FACE_INDEX<TV::dimension>,T>(forward_weights(i).x,forward_weights(i).y*remaining));
                weights_from(forward_weights(i).x).Append(PAIR<FACE_INDEX<TV::dimension>,int>(face,weights_to(face).m));}}}
    if(mpi_grid) mpi_grid->Sync_Common_Face_Weights_From(weights_to,weights_from,number_of_ghost_cells);
    LOG::Time("After Back");

    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& local_weights=weights_to(face);
        T sum=0;for(int i=1;i<=local_weights.m;i++) sum+=local_weights(i).y;
        for(int i=1;i<=local_weights.m;i++){assert(sum>1e-2);
            local_weights(i).y=local_weights(i).y/sum;}}
    if(mpi_grid){
        mpi_grid->ignore_boundary_faces=true;
        mpi_grid->Sync_Common_Face_Weights_From(weights_to,weights_from,number_of_ghost_cells);
        mpi_grid->ignore_boundary_faces=false;}
    LOG::Time("After Clamp");

    ARRAY<T,FACE_INDEX<TV::dimension> > sum_ic(grid,number_of_ghost_cells,false);
    if(num_iterations || num_diffusion_iterations) boundary_sum->Fill_Ghost_Cells_Face(grid,sum_jc,sum_ic,time,number_of_ghost_cells);
    for(int iter=1;iter<=num_iterations;iter++){
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
            ARRAY<PAIR<FACE_INDEX<TV::dimension>,int> >& local_weights=weights_from(face);
            T sum=0;for(int i=1;i<=local_weights.m;i++) sum+=sum_ic(local_weights(i).x)*weights_to(local_weights(i).x)(local_weights(i).y).y;
            for(int i=1;i<=local_weights.m;i++){assert(sum>1e-2);
                weights_to(local_weights(i).x)(local_weights(i).y).y=weights_to(local_weights(i).x)(local_weights(i).y).y/sum;}}
        if(mpi_grid){
            mpi_grid->ignore_boundary_faces=true;
            mpi_grid->Sync_Common_Face_Weights_To(weights_to,weights_from,number_of_ghost_cells);
            mpi_grid->ignore_boundary_faces=false;}
        for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
            ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& local_weights=weights_to(face);
            T sum=0;for(int i=1;i<=local_weights.m;i++) sum+=local_weights(i).y;
            for(int i=1;i<=local_weights.m;i++){assert(sum>1e-2);
                local_weights(i).y=local_weights(i).y/sum;}}
        if(mpi_grid){
            mpi_grid->ignore_boundary_faces=true;
            mpi_grid->Sync_Common_Face_Weights_From(weights_to,weights_from,number_of_ghost_cells);
            mpi_grid->ignore_boundary_faces=false;}}
    LOG::Time("After Clamp Iter");

    if(num_diffusion_iterations) for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,int> >& local_weights=weights_from(face);
        sum_jc(face)=0;for(int i=1;i<=local_weights.m;i++) sum_jc(face)+=sum_ic(local_weights(i).x)*weights_to(local_weights(i).x)(local_weights(i).y).y;}
    LOG::Time("After Sum");

    for(FACE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){Z(iterator.Full_Index())=T();}
    for(FACE_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){FACE_INDEX<TV::dimension> face=iterator.Full_Index();
        RANGE<TV_INT> domain=grid.Domain_Indices();domain.max_corner+=TV_INT::Axis_Vector(face.axis);
        ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> >& local_weights=weights_to(face);
        for(int i=1;i<=local_weights.m;i++)
            if(domain.Lazy_Inside(local_weights(i).x.index)) Z(local_weights(i).x)+=local_weights(i).y*Z_ghost(face);}
    LOG::Time("After Put");

    if(num_diffusion_iterations) Face_Diffusion(grid,sum_jc,Z,boundary,boundary_sum);
    LOG::Time("After Diffusion");
}

template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,1> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,2> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,3> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,1> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >,CUBIC_MN_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,1> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >,QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,2> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >,CUBIC_MN_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,2> >,float,AVERAGING_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >,QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,1> >,VECTOR<float,3>,AVERAGING_UNIFORM<GRID<VECTOR<float,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,1> >,VECTOR<float,3>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,2> >,VECTOR<float,4>,AVERAGING_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,2> >,VECTOR<float,4>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<float,3> >,VECTOR<float,5>,AVERAGING_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<float,3> >,VECTOR<float,5>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,1> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,2> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,3> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,1> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > >,CUBIC_MN_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,1> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > >,QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,2> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >,CUBIC_MN_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,2> >,double,AVERAGING_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >,QUADRATIC_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,1> >,VECTOR<double,3>,AVERAGING_UNIFORM<GRID<VECTOR<double,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,1> >,VECTOR<double,3>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,2> >,VECTOR<double,4>,AVERAGING_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,2> >,VECTOR<double,4>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > > >;
template class ADVECTION_CONSERVATIVE_UNIFORM<GRID<VECTOR<double,3> >,VECTOR<double,5>,AVERAGING_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >,LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<double,3> >,VECTOR<double,5>,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > > >;
#endif


