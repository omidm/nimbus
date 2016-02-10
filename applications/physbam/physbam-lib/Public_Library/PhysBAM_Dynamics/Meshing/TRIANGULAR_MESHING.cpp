//#####################################################################
// Copyright 2003-2005, Christopher Allocco, Ronald Fedkiw, Geoffrey Irving
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGULAR_MESHING
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Dynamics/Meshing/TRIANGULAR_MESHING.h>
using namespace PhysBAM;
//#####################################################################
// Function Snap_Nodes_To_Level_Set_Boundary
//#####################################################################
template<class T>
void TRIANGULAR_MESHING<T>::
Snap_Nodes_To_Level_Set_Boundary(const int output_frame,const int iterations)
{
    triangle_mesh.Initialize_Boundary_Nodes();
    for(int t=1;t<=triangle_mesh.boundary_nodes->m;t++) for(int k=1;k<=iterations;k++){
        int node=(*triangle_mesh.boundary_nodes)(t);VECTOR<T,2> X=particles.X(node);
        particles.X(node)-=implicit_curve.Extended_Phi(X)*implicit_curve.Extended_Normal(X);}
    Write_Tri_File_Format(output_frame,output_directory);
}
//#####################################################################
// Function Initialize_Optimization
//#####################################################################
template<class T> void TRIANGULAR_MESHING<T>::
Initialize_Optimization(const bool verbose)
{
    triangle_mesh.Initialize_Incident_Elements();triangle_mesh.Initialize_Boundary_Mesh();
    triangle_mesh.boundary_mesh->Initialize_Neighbor_Nodes();triangle_mesh.boundary_mesh->Initialize_Incident_Elements();
    triangle_mesh.Initialize_Boundary_Nodes();

    int i,j,k;VECTOR<T,2> xi,xj,xk;
    triangle_mesh.elements(1).Get(i,j,k);xi=particles.X(i);xj=particles.X(j);xk=particles.X(k);
    initial_min_altitude=TRIANGLE_2D<T>::Minimum_Altitude(xi,xj,xk);
    initial_area=TRIANGLE_2D<T>::Area(xi,xj,xk);
    map_from_nodes_to_boundary_list.Resize(1,triangle_mesh.number_nodes);
    for(i=1;i<=triangle_mesh.boundary_nodes->m;i++) map_from_nodes_to_boundary_list((*triangle_mesh.boundary_nodes)(i))=i;
    for(i=1;i<=layers.m;i++) delete layers(i);layers.Resize(1);layers(1)=triangle_mesh.boundary_nodes;
    triangle_mesh.boundary_nodes=0; // we don't need it hanging off the mesh object any more
    if(verbose) {std::stringstream ss;ss<<"boundary layer has "<<layers(1)->m<<" nodes"<<std::endl;LOG::filecout(ss.str());}
    ARRAY<bool,VECTOR<int,1> > marked(1,triangle_mesh.number_nodes);for(i=1;i<=layers(1)->m;i++) marked((*layers(1))(i))=true;
    for(int l=2;;l++){
        layers.Append(new ARRAY<int>);
        for(int i=1;i<=layers(l-1)->m;i++){
            j=(*layers(l-1))(i);
            for(k=1;k<=(*triangle_mesh.incident_elements)(j).m;k++) for(int a=1;a<=3;a++){
                int b=triangle_mesh.elements((*triangle_mesh.incident_elements)(j)(k))(a);
                if(!marked(b)){layers(l)->Append(b);marked(b)=true;}}}
        if(layers(l)->m == 0){delete layers(l);layers.Remove_End();break;}
        if(verbose) {std::stringstream ss;ss<<"layer "<<l<<" has "<<layers(l)->m<<" nodes"<<std::endl;LOG::filecout(ss.str());}}
    boundary_mesh_normals.Resize(1,layers(1)->m);
    Compute_Boundary_Mesh_Normals();
}
//#####################################################################
// Function Create_Final_Mesh_With_Optimization
//#####################################################################
template<class T> void TRIANGULAR_MESHING<T>::
Create_Final_Mesh_With_Optimization(const int number_of_initial_steps,const int number_of_final_steps,int* frame_number,const bool verbose)
{
    int local_frame=0;if(!frame_number) frame_number=&local_frame;
    int i;for(i=1;i<=number_of_initial_steps;i++){
        worst_boundary_quality=FLT_MAX;worst_interior_quality=FLT_MAX;
        if(verbose) {std::stringstream ss;ss<<"Working on initial iteration "<<i<<" of "<<number_of_initial_steps<<"+"<<number_of_final_steps<<std::endl;LOG::filecout(ss.str());}
        Optimize_Boundary_Layer((T).05);if(verbose) LOG::filecout(".");
        int j;for(j=2;j<=layers.m;j++){Optimize_Interior_Layer(j);if(verbose) LOG::filecout(".");}
        for(j=layers.m;j>=2;j--){Optimize_Interior_Layer(j,true);if(verbose) LOG::filecout(".");}
        Optimize_Boundary_Layer((T).05,true);if(verbose) LOG::filecout(".\n");
        Write_Tri_File_Format(*frame_number,output_directory);++*frame_number;}
    for(i=1;i<=number_of_final_steps;i++){
        worst_boundary_quality=FLT_MAX;worst_interior_quality=FLT_MAX;
        if(verbose) {std::stringstream ss;ss<<"Working on iteration "<<i<<" of "<<number_of_final_steps<<" (full step towards boundary)"<<std::endl;LOG::filecout(ss.str());}
        Optimize_Boundary_Layer(1);if(verbose) LOG::filecout(".");
        int j;for(j=2;j<=layers.m;j++){Optimize_Interior_Layer(j);if(verbose) LOG::filecout(".");}
        for(j=layers.m;j>=2;j--){Optimize_Interior_Layer(j,true);if(verbose) LOG::filecout(".");}
        Optimize_Boundary_Layer(1,true);if(verbose) LOG::filecout(".\n");
        Write_Tri_File_Format(*frame_number,output_directory);++*frame_number;}
}
//#####################################################################
// Function Optimize_Boundary_Layer
//#####################################################################
template<class T> void TRIANGULAR_MESHING<T>::
Optimize_Boundary_Layer(const T compression_fraction,const bool reverse)
{
    ARRAY<VECTOR<T,2> > directions(2);ARRAY<int>& nodes=*layers(1);
    int i;for(i=1;i<=nodes.m;i++) 
        particles.X(nodes(i))-=compression_fraction*implicit_curve(particles.X(nodes(i)))*boundary_mesh_normals(map_from_nodes_to_boundary_list(nodes(i)));
    Compute_Boundary_Mesh_Normals();
    for(i=1;i<=nodes.m;i++){
        int j;if(reverse) j=nodes(nodes.m+1-i);else j=nodes(i);
        VECTOR<T,2> normal=boundary_mesh_normals(map_from_nodes_to_boundary_list(j));
        directions(1)=VECTOR<T,2>(normal.y,-normal.x);
        directions(1).Normalize();
        directions(2)=-directions(1);
        Search_For_Best_Position(j,directions,true);} 
}
//#####################################################################
// Function Optimize_Interior_Layer
//#####################################################################
template<class T> void TRIANGULAR_MESHING<T>::
Optimize_Interior_Layer(const int layer,const bool reverse)
{
    ARRAY<VECTOR<T,2> > directions(5);
    directions(1)=VECTOR<T,2>(1,0);
    directions(2)=VECTOR<T,2>((T).30901699437494742410229341718282,(T).95105651629515357211643933337938);
    directions(3)=VECTOR<T,2>((T)-.80901699437494742410229341718282,(T).58778525229247312916870595463907);
    directions(4)=VECTOR<T,2>((T)-.80901699437494742410229341718282,(T)-.58778525229247312916870595463907);
    directions(5)=VECTOR<T,2>((T).30901699437494742410229341718282,(T)-.95105651629515357211643933337938);
    for(int i=1;i<=layers(layer)->m;i++){int j;if(reverse) j=(*layers(layer))(layers(layer)->m+1-i);else j=(*layers(layer))(i);Search_For_Best_Position(j,directions);}
}
//#####################################################################
// Function Search_For_Best_Position -- Need to modify
//#####################################################################
template<class T> void TRIANGULAR_MESHING<T>::
Search_For_Best_Position(const int node,const ARRAY<VECTOR<T,2> >& directions,bool include_boundary_terms)
{
    T best_quality=Quality_Of_Worst_Incident_Triangle(node);
    if(best_quality<=worst_interior_quality) worst_interior_quality=best_quality;
    VECTOR<T,2> best_x(particles.X(node)),xi,xj,xk;T alpha=FLT_MAX;
    for(int s=1;s<=(*triangle_mesh.incident_elements)(node).m;s++){
        int t=(*triangle_mesh.incident_elements)(node)(s);
        int i,j,k;triangle_mesh.elements(t).Get(i,j,k);xi=particles.X(i);xj=particles.X(j);xk=particles.X(k);
        T area=TRIANGLE_2D<T>::Area(xi,xj,xk);
        if (i==node)alpha=min(alpha,2*area/(xj-xk).Magnitude());
        if (j==node)alpha=min(alpha,2*area/(xi-xk).Magnitude());
        if (k==node)alpha=min(alpha,2*area/(xi-xj).Magnitude());
    }
    alpha*=(T).3;
    int strikes=0,last_direction=1;
    while(strikes<=3){
        T localbest_quality=best_quality;VECTOR<T,2> localbest_x=best_x;
        for(int d=1;d<=directions.m;d++){
            int this_direction;if(d%2) this_direction=last_direction+d/2;else this_direction=last_direction-d/2;
            this_direction=(this_direction+directions.m-1)%directions.m+1;
            particles.X(node)=best_x+alpha*directions(this_direction);
            T q=Quality_Of_Worst_Incident_Triangle(node);
            if(q>localbest_quality){localbest_quality=q;localbest_x=particles.X(node);last_direction=this_direction;break;}}
        if(localbest_quality>best_quality){best_quality=localbest_quality;best_x=localbest_x;}
        else{strikes++;alpha*=(T).45;}}
}
//#####################################################################
// Function Quality_Of_Worst_Incident_Triangle
//#####################################################################
template<class T> T TRIANGULAR_MESHING<T>::
Quality_Of_Worst_Incident_Triangle(const int node)
{
    TRIANGLE_2D<T> triangle;T worst_quality=1;
    for(int s=1;s<=(*triangle_mesh.incident_elements)(node).m;s++){
        int t=(*triangle_mesh.incident_elements)(node)(s);
        int i,j,k;triangle_mesh.elements(t).Get(i,j,k);
        VECTOR<T,2> xi(particles.X(i).x,particles.X(i).y),xj(particles.X(j).x,particles.X(j).y),
                     xk(particles.X(k).x,particles.X(k).y);
        triangle.Specify_Three_Points(xi,xj,xk);
        T min_altitude = triangle.Minimum_Altitude(xi,xj,xk);
        worst_quality=min<T>(worst_quality,1/triangle.Aspect_Ratio()+(T).0/((min_altitude<initial_min_altitude)?((initial_min_altitude/min_altitude)/**(initial_min_altitude/min_altitude)*/):((min_altitude/initial_min_altitude)/**(min_altitude/initial_min_altitude)*/)));}  
    return worst_quality;
}
//#####################################################################
// Function Compute_Boundary_Mesh_Normals
//#####################################################################
template<class T> void TRIANGULAR_MESHING<T>::
Compute_Boundary_Mesh_Normals()
{
    boundary_mesh_normals.Fill(VECTOR<T,2>());
    for(int t=1;t<=triangle_mesh.boundary_mesh->elements.m;t++){
        int i,j;triangle_mesh.boundary_mesh->elements(t).Get(i,j);
        VECTOR<T,2> normal(particles.X(j).y-particles.X(i).y,-particles.X(j).x+particles.X(i).x);
        normal.Normalize();
        boundary_mesh_normals(map_from_nodes_to_boundary_list(i))+=normal;
        boundary_mesh_normals(map_from_nodes_to_boundary_list(j))+=normal;
    }
    for(int i=1;i<=boundary_mesh_normals.counts.x;i++) boundary_mesh_normals(i).Normalize();
}
//#####################################################################
// Function Create_Initial_Mesh
//#####################################################################
template<class T> void TRIANGULAR_MESHING<T>::
Create_Initial_Mesh(const T lattice_cell_size,const bool use_adaptive_refinement,const int max_subdivision_levels,const bool discard_to_get_nice_topology,const bool verbose)
{  
    // initial mesh
    implicit_curve.Update_Box();
    VECTOR<T,2> size=implicit_curve.box.Edge_Lengths();
    T x_cell_size=lattice_cell_size;
    T y_cell_size=lattice_cell_size*sqrt((T).75);
    int m=(int)ceil(size.x/x_cell_size),n=(int)ceil(size.y/y_cell_size);
    GRID<TV> grid(m,n,implicit_curve.box.min_corner.x,implicit_curve.box.min_corner.x+x_cell_size*m,implicit_curve.box.min_corner.y,implicit_curve.box.min_corner.y+y_cell_size*n);
    triangulated_area.Initialize_Equilateral_Mesh_And_Particles(grid);
    Write_Tri_File_Format(999,output_directory);
    triangulated_area.Discard_Triangles_Outside_Implicit_Curve(implicit_curve); 
    triangulated_area.Discard_Valence_Zero_Particles_And_Renumber();
    Write_Tri_File_Format(1000,output_directory);
}
//#####################################################################
// Function Write_Diagnostic_Files
//#####################################################################
template<class T> void TRIANGULAR_MESHING<T>::
Write_Diagnostic_Files(const int stepnumber)
{
//  std::ostream* output=FILE_UTILITIES::Safe_Open_Output(STRING_UTILITIES::string_sprintf("%s/diagnostics/diagnostic_info%.3d",output_directory.c_str(),stepnumber),false);
//  int index=0;
//  *output <<  "max_phi (index) " << triangulated_area.Maximum_Magnitude_Phi_On_Boundary(implicit_curve,&index); input << "  (" << index << ")" << std::endl;
//  *output << "max_velocity (index) " << particles.Maximum_Speed(&index);input << "  (" << index << ")" <<  std::endl;
//  *output << "max_aspect_ratio (index) " << triangulated_area.Maximum_Aspect_Ratio(&index);input << "  (" << index << ")" <<  std::endl;
//  *output << "max_boundary_aspect_ratio (index) " << triangulated_area.Maximum_Boundary_Aspect_Ratio(&index);input << "  (" << index << ")" <<  std::endl;
//  *output << "max_interior_aspect_ratio (index) " << triangulated_area.Maximum_Interior_Aspect_Ratio(&index);input << "  (" << index << ")" <<  std::endl;
//  *output << "avg_boundary_aspect_ratio " << triangulated_area.Average_Boundary_Aspect_Ratio() << std::endl;
//  *output << "avg_interior_aspect_ratio" << triangulated_area.Average_Interior_Aspect_Ratio() << std::endl;
//  *output << "min_area " << triangulated_area.Minimum_Area(&index) << std::endl; 
//  *output << "min_angle " << 180/pi*triangulated_area.Minimum_Angle() << std::endl; 
//  *output << "max_angle " << 180/pi*triangulated_area.Maximum_Angle() << std::endl;
//  *output << "min_dihedral_angle " << 180/pi*triangulated_area.Minimum_Dihedral_Angle() << std::endl;
//  *output << "max_dihedral_angle " << 180/pi*triangulated_area.Maximum_Dihedral_Angle() << std::endl;
//  if(use_masses_and_springs){
//      *output << "min_edge_length " << triangulated_area.Minimum_Edge_Length() << std::endl;
//      *output << "max_edge_length " << triangulated_area.Maximum_Edge_Length() << std::endl;}
//  *output << "min_altitude " << triangulated_area.Minimum_Altitude() << std::endl; 
//  if(edge_springs){
//      *output << "max_edge_compression " << edge_springs->Maximum_Compression_Or_Expansion_Fraction(&index);input << "  (" << index << ")" << std::endl;}
//  delete output;
}
//#####################################################################
// Function Write_Tri_File_Format
//#####################################################################
template<class T> void TRIANGULAR_MESHING<T>::
Write_Tri_File_Format(const int stepnumber,const std::string& output_directory)
{
    FILE_UTILITIES::Write_To_File<T>(STRING_UTILITIES::string_sprintf("%s/triangulated_area_%d.tri2d",output_directory.c_str(),stepnumber),triangulated_area);
}
//#####################################################################
template class TRIANGULAR_MESHING<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TRIANGULAR_MESHING<double>;
#endif
