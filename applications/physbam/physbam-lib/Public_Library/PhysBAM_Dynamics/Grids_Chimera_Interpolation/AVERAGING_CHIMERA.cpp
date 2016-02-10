//#####################################################################
// Copyright 2005-2010, Eran Guendelman, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Dynamics/Grids_Chimera_Interpolation/AVERAGING_CHIMERA.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Dynamics/Grids_Chimera_PDE_Linear/Parallel_Computation/LAPLACE_CHIMERA_MPI.h>
#include <PhysBAM_Dynamics/Grids_Chimera_PDE_Linear/Parallel_Computation/LAPLACE_CHIMERA_GRID_MPI.h>

using namespace PhysBAM;

//#####################################################################
// Function Face_To_Cell_Vector_Weights
//#####################################################################
template<class T_GRID> ARRAY<TRIPLE<int,FACE_INDEX<T_GRID::VECTOR_T::dimension>,typename T_GRID::VECTOR_T> > AVERAGING_CHIMERA<T_GRID>::Face_To_Cell_Vector_Weights(const GRID_CELL_INDEX& grid_cell_index) const
{
    ARRAY<TRIPLE<int,D_FACE_INDEX,TV> > weights;
    //ARRAY<T> areas;
    
    for(int axis=1;axis<=TV::dimension;axis++)
        for(int side=1;side<=2;side++){
            D_FACE_INDEX face_index(axis,grid_cell_index.y+(side-1)*TV_INT::Axis_Vector(axis));
            if(laplace_grid.Chimera_Face(grid_cell_index.x,face_index)){
                //areas.Append(laplace_grid.Grid(grid_cell_index.x).Face_Size(axis));
                weights.Append(Tuple(0,face_index,laplace_grid.Frame(grid_cell_index.x).r.Rotate(TV::Axis_Vector(axis))));}}
    
    if(laplace_grid.cell_indices_to_incident_voronoi_face_indices.Contains(grid_cell_index)){
        const ARRAY<int>& incident=laplace_grid.cell_indices_to_incident_voronoi_face_indices.Get(grid_cell_index);
        for(int i=1;i<=incident.Size();i++){
            const VECTOR<GRID_CELL_INDEX,2>& grid_cell_indices=laplace_grid.voronoi_faces(incident(i)).x;
            TV location_1=laplace_grid.Frame(grid_cell_indices(1).x)*laplace_grid.Grid(grid_cell_indices(1).x).X(grid_cell_indices(1).y);
            TV location_2=laplace_grid.Frame(grid_cell_indices(2).x)*laplace_grid.Grid(grid_cell_indices(2).x).X(grid_cell_indices(2).y);
            TV normal=(location_2-location_1).Normalized();
            //areas.Append(laplace_grid.voronoi_faces(incident(i)).y);
            weights.Append(Tuple(incident(i),D_FACE_INDEX(),normal));}}
    
    MATRIX_MXN<T> normal_matrix(TV::dimension,weights.Size());
    //MATRIX_MXN<T> area_matrix(weights.Size(),weights.Size());
    for(int i=1;i<=weights.Size();i++){
        //weight_matrix(i,i)=areas(i);
        normal_matrix.Set_Column(i,weights(i).z);}

    MATRIX<T,TV::dimension,TV::dimension> nnti=MATRIX<T,TV::dimension,TV::dimension>(normal_matrix/**area_matrix*/*normal_matrix.Transposed()).Inverse();
    for(int i=1;i<=weights.Size();i++)
        weights(i).z=nnti/**areas(i)*/*weights(i).z;

    //EETODO: work out partial least squares formulation for when we don't have all the face data if possible - for MPI
    //EETODO: use area weights?
    //EETODO: use momentum averaging?
    
    return weights;
}
//#####################################################################
// Function Face_To_Cell_Vector_Weights
//#####################################################################
template<class T_GRID> typename T_GRID::VECTOR_T AVERAGING_CHIMERA<T_GRID>::Face_To_Cell_Vector(LAPLACE_CHIMERA_MPI<T_GRID>& laplace,LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time,const VECTOR_ND<T>& u,const GRID_CELL_INDEX& grid_cell_index) const
{
    ARRAY<TRIPLE<int,D_FACE_INDEX,TV> > weights=Face_To_Cell_Vector_Weights(grid_cell_index);
    TV vector;
    for(int i=1;i<=weights.Size();i++){
        int matrix_face_index;
        if(weights(i).x) matrix_face_index=laplace.Matrix_Face_Index(weights(i).x);
        else matrix_face_index=laplace.Matrix_Face_Index(grid_cell_index.x,weights(i).y);
        if(matrix_face_index)
            vector+=u(matrix_face_index)*weights(i).z;
        else{
            T face_value;
            TV location,normal;
            if(weights(i).x){
                const VECTOR<GRID_CELL_INDEX,2>& grid_cell_indices=laplace_grid.voronoi_faces(weights(i).x).x;
                TV location_1=laplace_grid.Frame(grid_cell_indices(1).x)*laplace_grid.Grid(grid_cell_indices(1).x).X(grid_cell_indices(1).y);
                TV location_2=laplace_grid.Frame(grid_cell_indices(2).x)*laplace_grid.Grid(grid_cell_indices(2).x).X(grid_cell_indices(2).y);
                location=(location_1+location_2)*(T)0.5;
                normal=(location_2-location-1).Normalized();}
            else{
                location=laplace_grid.Frame(grid_cell_index.x)*laplace_grid.Grid(grid_cell_index.x).Face(weights(i).y.axis,weights(i).y.index);
                normal=laplace_grid.Frame(grid_cell_index.x).r.Rotate(TV::Axis_Vector(weights(i).y.axis));}
            callbacks.Get_Neumann_Boundary_Condition(location,normal,face_value,time);
            vector+=face_value*weights(i).z;}}
    return vector;
}
//#####################################################################
// Function Build_Face_To_Cell_Matrix
//#####################################################################
template<class T_GRID> void AVERAGING_CHIMERA<T_GRID>::Build_Face_To_Cell_Matrix(LAPLACE_CHIMERA_MPI<T_GRID>& laplace,LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time)
{
    ARRAY<int> face_to_cell_matrix_row_lengths(laplace.n_matrix_cells*TV::dimension);
    for(int grid_index=1;grid_index<=laplace_grid.n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index))
            for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                int matrix_cell_index=laplace.Matrix_Cell_Index(grid_index,cell_index);
                if(matrix_cell_index){
                    int incident_faces=0;
                    for(int axis=1;axis<=TV::dimension;axis++)
                        for(int side=1;side<=2;side++)
                            if(laplace.Matrix_Face_Index(grid_index,D_FACE_INDEX(axis,cell_index+(side-1)*TV_INT::Axis_Vector(axis))))
                                incident_faces++;
                    GRID_CELL_INDEX grid_cell_index=GRID_CELL_INDEX(grid_index,cell_index);
                    if(laplace_grid.cell_indices_to_incident_voronoi_face_indices.Contains(grid_cell_index)){
                        const ARRAY<int>& incident=laplace_grid.cell_indices_to_incident_voronoi_face_indices.Get(grid_cell_index);
                        for(int i=1;i<=incident.Size();i++)
                            if(laplace.Matrix_Face_Index(incident(i)))
                                incident_faces++;}
                    for(int axis=1;axis<=TV::dimension;axis++)
                        face_to_cell_matrix_row_lengths((axis-1)*laplace.n_matrix_cells+matrix_cell_index)=incident_faces;}}

    face_to_cell_matrix.Set_Row_Lengths(face_to_cell_matrix_row_lengths);
    face_to_cell_matrix.n=laplace.n_matrix_faces;
    for(int grid_index=1;grid_index<=laplace_grid.n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index))
            for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                int matrix_cell_index=laplace.Matrix_Cell_Index(grid_index,cell_index);
                if(matrix_cell_index){
                    ARRAY<TRIPLE<int,D_FACE_INDEX,TV> > weights=Face_To_Cell_Vector_Weights(GRID_CELL_INDEX(grid_index,cell_index));
                    for(int i=1;i<=weights.Size();i++){
                        int matrix_face_index;
                        if(weights(i).x) matrix_face_index=laplace.Matrix_Face_Index(weights(i).x);
                        else matrix_face_index=laplace.Matrix_Face_Index(grid_index,weights(i).y);
                        if(matrix_face_index)
                            for(int axis=1;axis<=TV::dimension;axis++)
                                face_to_cell_matrix.Set_Element((axis-1)*laplace.n_matrix_cells+matrix_cell_index,matrix_face_index,weights(i).z(axis));}}}
}
//#####################################################################
// Function Build_Face_To_Cell_Matrix
//#####################################################################
template<class T_GRID> void AVERAGING_CHIMERA<T_GRID>::Build_Face_To_Cell_Vector(LAPLACE_CHIMERA_MPI<T_GRID>& laplace,LAPLACE_CHIMERA_MPI_CALLBACKS<T_GRID>& callbacks,const T time)
{
    face_to_cell_vector.Resize(laplace.n_matrix_cells*TV::dimension);
    face_to_cell_vector.Fill(0);
    for(int grid_index=1;grid_index<=laplace_grid.n_global_grids;grid_index++)
        if(laplace_grid.Local_Grid(grid_index))
            for(CELL_ITERATOR iterator(laplace_grid.Grid(grid_index));iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                int matrix_cell_index=laplace.Matrix_Cell_Index(grid_index,cell_index);
                if(matrix_cell_index){
                    ARRAY<TRIPLE<int,D_FACE_INDEX,TV> > weights=Face_To_Cell_Vector_Weights(GRID_CELL_INDEX(grid_index,cell_index));
                    for(int i=1;i<=weights.Size();i++){
                        int matrix_face_index;
                        if(weights(i).x) matrix_face_index=laplace.Matrix_Face_Index(weights(i).x);
                        else matrix_face_index=laplace.Matrix_Face_Index(grid_index,weights(i).y);
                        if(!matrix_face_index){
                            TV location,normal;
                            if(weights(i).x){
                                const VECTOR<GRID_CELL_INDEX,2>& grid_cell_indices=laplace_grid.voronoi_faces(weights(i).x).x;
                                TV location_1=laplace_grid.Frame(grid_cell_indices(1).x)*laplace_grid.Grid(grid_cell_indices(1).x).X(grid_cell_indices(1).y);
                                TV location_2=laplace_grid.Frame(grid_cell_indices(2).x)*laplace_grid.Grid(grid_cell_indices(2).x).X(grid_cell_indices(2).y);
                                location=(location_1+location_2)*(T)0.5;
                                normal=(location_2-location_1).Normalized();}
                            else{
                                location=laplace_grid.Frame(grid_index)*laplace_grid.Grid(grid_index).Face(weights(i).y.axis,weights(i).y.index);
                                normal=laplace_grid.Frame(grid_index).r.Rotate(TV::Axis_Vector(weights(i).y.axis));}
                            T value;
                            if(callbacks.Get_Neumann_Boundary_Condition(location,normal,value,time))
                                for(int axis=1;axis<=TV::dimension;axis++)
                                    face_to_cell_vector((axis-1)*laplace.n_matrix_cells+matrix_cell_index)+=value*weights(i).z(axis);}}}}
}
//#####################################################################

template class AVERAGING_CHIMERA<GRID<VECTOR<float,1> > >;
template class AVERAGING_CHIMERA<GRID<VECTOR<float,2> > >;
template class AVERAGING_CHIMERA<GRID<VECTOR<float,3> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class AVERAGING_CHIMERA<GRID<VECTOR<double,1> > >;
template class AVERAGING_CHIMERA<GRID<VECTOR<double,2> > >;
template class AVERAGING_CHIMERA<GRID<VECTOR<double,3> > >;
#endif
