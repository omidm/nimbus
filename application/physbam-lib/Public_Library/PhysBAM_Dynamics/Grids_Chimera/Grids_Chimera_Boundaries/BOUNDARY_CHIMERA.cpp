//#####################################################################
// Copyright 2011, Linhai Qiu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_CHIMERA  
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Dynamics/Grids_Chimera/Grids_Chimera_Boundaries/BOUNDARY_CHIMERA.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Dynamics/Grids_Chimera/Parallel_Computation/CHIMERA_GRID.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T2> BOUNDARY_CHIMERA<T_GRID,T2>::
BOUNDARY_CHIMERA(CHIMERA_GRID<T_GRID,T2>* chimera_grid_input,T_BOUNDARY_T2& boundary_input)
    :chimera_grid(chimera_grid_input),boundary(boundary_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID,class T2> BOUNDARY_CHIMERA<T_GRID,T2>::
~BOUNDARY_CHIMERA()
{}
//#####################################################################
// Function Set_Constant_Extrapolation
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_CHIMERA<T_GRID,T2>::
Set_Constant_Extrapolation(const TV_SIDES& constant_extrapolation_input)
{
    boundary.Set_Constant_Extrapolation(constant_extrapolation_input);
}
//#####################################################################
// Function Constant_Extrapolation
//#####################################################################
template<class T_GRID,class T2> bool BOUNDARY_CHIMERA<T_GRID,T2>::
Constant_Extrapolation(const int side) const
{
    return boundary.Constant_Extrapolation(side);
}
//#####################################################################
// Function Set_Fixed_Boundary
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_CHIMERA<T_GRID,T2>::
Set_Fixed_Boundary(const bool use_fixed_boundary_input,const T2 fixed_boundary_value_input)
{
    boundary.Set_Fixed_Boundary(use_fixed_boundary_input,fixed_boundary_value_input);
}
//#####################################################################
// Function Apply_Boundary_Condition
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_CHIMERA<T_GRID,T2>::
Apply_Boundary_Condition(const T_GRID& grid,T_ARRAYS_T2& u,const T time)
{
}
//#####################################################################
// Function Apply_Boundary_Condition_Face
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_CHIMERA<T_GRID,T2>::
Apply_Boundary_Condition_Face(const T_GRID& grid,T_FACE_ARRAYS_T2& u,const T time)
{
}
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_CHIMERA<T_GRID,T2>::
Fill_Ghost_Cells(const T_GRID& grid,const T_ARRAYS_T2& u,T_ARRAYS_T2& u_ghost,const T dt,const T time,const int number_of_ghost_cells_input)
{
    int grid_local_index=chimera_grid->Find_Current_Grid_Local_Index(grid);
    T_ARRAYS_T2::Put(chimera_grid->scalar_field_buffer(grid_local_index),u_ghost);
}
//#####################################################################
// Function Fill_Ghost_Cells_Face
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_CHIMERA<T_GRID,T2>::
Fill_Ghost_Cells_Face(const T_GRID& grid,const T_FACE_ARRAYS_T2& u,T_FACE_ARRAYS_T2& u_ghost,const T time,const int number_of_ghost_cells_input)
{
    int grid_local_index=chimera_grid->Find_Current_Grid_Local_Index(grid);
    T_FACE_ARRAYS_T2::Put(chimera_grid->face_scalar_field_buffer(grid_local_index),u_ghost);
}
//################################################################################
// Function Interpolation_From_Delaunay_Triangle_Vertices
//################################################################################
template<class T_GRID,class T2> typename T_GRID::SCALAR BOUNDARY_CHIMERA<T_GRID,T2>::
Interpolation_From_Delaunay_Triangle_Vertices(const VECTOR<T,1> location,const ARRAY<VECTOR<T,1> > points,const ARRAY<T> values,bool& selected,const T tolerance,const ARRAY<T>* const phi_values,const T interface_value)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
//################################################################################
// Function Interpolation_From_Delaunay_Triangle_Vertices
//################################################################################
template<class T_GRID,class T2> typename T_GRID::SCALAR BOUNDARY_CHIMERA<T_GRID,T2>::
Interpolation_From_Delaunay_Triangle_Vertices(const VECTOR<T,2> location,const ARRAY<VECTOR<T,2> > points,const ARRAY<T> values,bool& selected,const T tolerance,const ARRAY<T>* const phi_values,const T interface_value)
{
    typedef VECTOR<T,2> TV;typedef typename BASIC_SIMPLEX_POLICY<TV,2>::SIMPLEX T_SIMPLEX;
    typedef VECTOR<T,3> TV_COOR;
    VECTOR<TV,3> selected_points;ARRAY<T> selected_values(3);VECTOR<T,3> selected_phi_values;T interpolated_value=(T)0;
    selected=false;
    for(int p1=1;p1<=points.Size();p1++){
        for(int p2=p1+1;p2<=points.Size();p2++){
            for(int p3=p2+1;p3<=points.Size();p3++){
                selected_points(1)=points(p1);selected_points(2)=points(p2);selected_points(3)=points(p3);
                T_SIMPLEX triangle(selected_points);triangle.Fix_Orientation();
                T edge_tolerance=tolerance*triangle.Maximum_Edge_Length();
                if(triangle.Minimum_Altitude()>edge_tolerance)
                {
                    TV circumcenter=T_SIMPLEX::Circumcenter(selected_points(1),selected_points(2),selected_points(3));
                    T radius=(circumcenter-points(p1)).Magnitude();
                    bool delaunay=true;
                    for(int p4=1;p4<=points.Size();p4++){
                        if(p4==p1 || p4==p2 || p4==p3) continue;
                        if((circumcenter-points(p4)).Magnitude()<(radius-edge_tolerance)) delaunay=false;}
                    if(delaunay && triangle.Inside(location,edge_tolerance)){
                        selected=true;
                        selected_values(1)=values(p1);selected_values(2)=values(p2);selected_values(3)=values(p3);
                        if(phi_values){
                            selected_phi_values(1)=(*phi_values)(p1);selected_phi_values(2)=(*phi_values)(p2);selected_phi_values(3)=(*phi_values)(p3);}
                        break;}}}
            if(selected) break;}
        if(selected) break;}

    int positive_count=0;
    if(phi_values) for(int p=1;p<=3;p++) if(selected_phi_values(p)>0) positive_count++;
    if(positive_count!=0 && positive_count!=3){
        T interpolated_phi(0);
        TV_COOR bary_coordinates_test_phi=T_SIMPLEX::Barycentric_Coordinates(location,selected_points);
        for(int p=1;p<=3;p++) interpolated_phi+=bary_coordinates_test_phi(p)*selected_phi_values(p);
        if(interpolated_phi>0){
            return interface_value;}

        ARRAY<TV> X(3);X=selected_points;
        ARRAY<VECTOR<int,3> > left_tris,right_tris; 
        T_SIMPLEX triangle(selected_points);
        triangle.Cut_Simplex(X,VECTOR<int,3>(1,2,3),selected_points,selected_phi_values,left_tris,right_tris);

        for(int i=1;i<=left_tris.Size();i++){
            VECTOR<int,3> tri=left_tris(i);
            VECTOR<TV,3> tri_points(X(tri(1)),X(tri(2)),X(tri(3)));
            T_SIMPLEX cut_triangle(tri_points);cut_triangle.Fix_Orientation();
            if(cut_triangle.Inside(location,tolerance*cut_triangle.Maximum_Edge_Length())){
                VECTOR<bool,3> tri_point_selected;tri_point_selected.Fill(false);
                for(int select_index=1;select_index<=3;select_index++)if(!tri.Contains(select_index))
                    for(int j=1;j<=3;j++)if(tri(j)>3 && !tri_point_selected(j)){
                        selected_points(select_index)=tri_points(j);
                        selected_values(select_index)=interface_value;
                        tri_point_selected(j)=true;break;}
                break;
            }
            if(i==left_tris.Size()) return interface_value;
        }
    }
    
    TV_COOR bary_coordinates=T_SIMPLEX::Barycentric_Coordinates(location,selected_points);
    for(int p=1;p<=3;p++) interpolated_value+=bary_coordinates(p)*selected_values(p);
    //assert(selected);
    return interpolated_value;
}
//################################################################################
// Function Interpolation_From_Delaunay_Triangle_Vertices
//################################################################################
template<class T_GRID,class T2> typename T_GRID::SCALAR BOUNDARY_CHIMERA<T_GRID,T2>::
Interpolation_From_Delaunay_Triangle_Vertices(const VECTOR<T,3> location,const ARRAY<VECTOR<T,3> > points,const ARRAY<T> values,bool& selected,const T tolerance,const ARRAY<T>* const phi_values,const T interface_value)
{
    typedef VECTOR<T,3> TV;typedef typename BASIC_SIMPLEX_POLICY<TV,3>::SIMPLEX T_SIMPLEX;
    typedef VECTOR<T,4> TV_COOR;
    VECTOR<TV,4> selected_points;ARRAY<T> selected_values(4);VECTOR<T,4> selected_phi_values;T interpolated_value=(T)0;
    selected=false;
    int p1=1;
    for(int p2=p1+1;p2<=points.Size();p2++){
        for(int p3=p2+1;p3<=points.Size();p3++){
            for(int p4=p3+1;p4<=points.Size();p4++){
                selected_points(1)=points(p1);selected_points(2)=points(p2);selected_points(3)=points(p3);selected_points(4)=points(p4);
                T_SIMPLEX tetrahedron(selected_points);
                T edge_tolerance=tolerance*tetrahedron.Maximum_Edge_Length();
                if(tetrahedron.Volume()>edge_tolerance*edge_tolerance*edge_tolerance/6){
                    TV circumcenter=T_SIMPLEX::Circumcenter(selected_points(1),selected_points(2),selected_points(3),selected_points(4));
                    T radius=(circumcenter-points(p1)).Magnitude();
                    bool delaunay=true;
                    for(int p5=1;p5<=points.Size();p5++){
                        if(p5==p1 || p5==p2 || p5==p3 || p5==p4) continue;
                        if((circumcenter-points(p5)).Magnitude()<radius-edge_tolerance) delaunay=false;}
                    if(delaunay && T_SIMPLEX::Barycentric_Inside(location,selected_points(1),selected_points(2),selected_points(3),selected_points(4),edge_tolerance)){
                        selected=true;
                        selected_values(1)=values(p1);selected_values(2)=values(p2);selected_values(3)=values(p3);selected_values(4)=values(p4);
                        if(phi_values){
                            selected_phi_values(1)=(*phi_values)(p1);selected_phi_values(2)=(*phi_values)(p2);
                            selected_phi_values(3)=(*phi_values)(p3);selected_phi_values(4)=(*phi_values)(p4);}
                        break;}}}
            if(selected) break;}
        if(selected) break;}

    int positive_count=0;
    if(phi_values) for(int p=1;p<=4;p++) if(selected_phi_values(p)>0) positive_count++;
    if(positive_count!=0 && positive_count!=4){
        T interpolated_phi(0);
        TV_COOR bary_coordinates_test_phi=T_SIMPLEX::Barycentric_Coordinates(location,selected_points);
        for(int p=1;p<=4;p++) interpolated_phi+=bary_coordinates_test_phi(p)*selected_phi_values(p);
        if(interpolated_phi>0) return interface_value;

        ARRAY<TV> X(4);X=selected_points;
        ARRAY<VECTOR<int,4> > left_tets,right_tets; 
        T_SIMPLEX tetrahedron(selected_points);
        tetrahedron.Cut_Simplex(X,VECTOR<int,4>(1,2,3,4),selected_points,selected_phi_values,left_tets,right_tets);
        for(int i=1;i<=left_tets.Size();i++){
            VECTOR<int,4> tet=left_tets(i);
            VECTOR<TV,4> tet_points(X(tet(1)),X(tet(2)),X(tet(3)),X(tet(4)));
            T_SIMPLEX cut_tetrahedron(tet_points);
            if(cut_tetrahedron.Inside(location,-tolerance*cut_tetrahedron.Maximum_Edge_Length())){
                VECTOR<bool,4> tet_point_selected;tet_point_selected.Fill(false);
                for(int select_index=1;select_index<=4;select_index++)if(!tet.Contains(select_index))
                    for(int j=1;j<=4;j++)if(tet(j)>4 && !tet_point_selected(j)){
                        selected_points(select_index)=tet_points(j);
                        selected_values(select_index)=interface_value;
                        tet_point_selected(j)=true;break;}
                break;
            }
            if(i==left_tets.Size()){
                std::cout<<"no tet enclosing!"<<std::endl;
                PHYSBAM_FATAL_ERROR();
            }
        }
    }
    
    TV_COOR bary_coordinates=T_SIMPLEX::Barycentric_Coordinates(location,selected_points);
    for(int p=1;p<=4;p++){
        //std::cout<<"debug selected points:"<<selected_points(p)<<",weights"<<bary_coordinates(p)<<",selected_value:"<<selected_values(p)<<std::endl;
        interpolated_value+=bary_coordinates(p)*selected_values(p);}
    //assert(selected);
    return interpolated_value;
}
//#####################################################################
template class BOUNDARY_CHIMERA<GRID<VECTOR<float,1> >,float>;
template class BOUNDARY_CHIMERA<GRID<VECTOR<float,2> >,float>;
template class BOUNDARY_CHIMERA<GRID<VECTOR<float,3> >,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BOUNDARY_CHIMERA<GRID<VECTOR<double,1> >,double>;
template class BOUNDARY_CHIMERA<GRID<VECTOR<double,2> >,double>;
template class BOUNDARY_CHIMERA<GRID<VECTOR<double,3> >,double>;
#endif
