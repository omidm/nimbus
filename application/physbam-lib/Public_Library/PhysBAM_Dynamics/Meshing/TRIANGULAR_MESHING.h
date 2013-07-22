//#####################################################################
// Copyright 2003-2005, Christopher Allocco, Ronald Fedkiw, Geoffrey Irving
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGULAR_MESHING
//##################################################################### 
#ifndef __TRIANGULAR_MESHING__
#define __TRIANGULAR_MESHING__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
namespace PhysBAM{

template<class TV> class IMPLICIT_OBJECT;

template<class T>
class TRIANGULAR_MESHING:public NONCOPYABLE
{
    typedef VECTOR<T,2> TV;
public:
    IMPLICIT_OBJECT<TV>& implicit_curve;
    TRIANGLE_MESH triangle_mesh;
    TRIANGULATED_AREA<T> triangulated_area;
    COLLISION_GEOMETRY_COLLECTION<TV> collision_body_list;
    PARTICLES<VECTOR<T,2> > particles;
    std::string output_directory;
    bool preserve_altitude;
private:
    T curvature_subdivision_threshold,interpolation_error_subdivision_threshold;
    T initial_min_altitude;
    T initial_area;
    ARRAY<int,VECTOR<int,1> > map_from_nodes_to_boundary_list;
    ARRAY<ARRAY<int>*> layers; // each sublist is a list of nodes in one layer, with the first sublist being the boundary
    ARRAY<VECTOR<T,2> ,VECTOR<int,1> > boundary_mesh_normals;
    T worst_boundary_quality,worst_interior_quality;

public:
    TRIANGULAR_MESHING(IMPLICIT_OBJECT<TV>& implicit_curve_input)
        :implicit_curve(implicit_curve_input),triangulated_area(triangle_mesh,particles)
    {
        particles.Store_Velocity();particles.Store_Mass();
        Set_Output_Directory("meshing_data");
    }

    ~TRIANGULAR_MESHING()
    {for(int i=1;i<=layers.m;i++) delete layers(i);layers.Resize(0);}

    void Set_Curvature_Subdivision_Threshold(const T curvature_subdivision_threshold_input=.7)
    {curvature_subdivision_threshold=curvature_subdivision_threshold_input;}

    void Set_Interpolation_Error_Subdivision_Threshold(const T interpolation_error_subdivision_threshold_input=.08)
    {interpolation_error_subdivision_threshold=interpolation_error_subdivision_threshold_input;}

    void Set_Output_Directory(const std::string& output_directory_input)
    {output_directory=output_directory_input;}

//#####################################################################
    void Snap_Nodes_To_Level_Set_Boundary(const int output_frame=1,const int iterations=100);
    void Initialize_Optimization(const bool verbose=true);
    void Create_Final_Mesh_With_Optimization(const int number_of_initial_steps=4,const int number_of_final_steps=6,int* frame_number=0,const bool verbose=true);
private:
    void Optimize_Boundary_Layer(const T compression_fraction,const bool reverse=false);
    void Optimize_Interior_Layer(const int layer,const bool reverse=false);
    void Search_For_Best_Position(const int node,const ARRAY<VECTOR<T,2> >& direction,bool include_boundary_terms=false);
    T Quality_Of_Worst_Incident_Triangle(const int node);
    void Compute_Boundary_Mesh_Normals();
public:
    void Create_Initial_Mesh(const T lattice_cell_size=0,const bool use_adaptive_refinement=false,const int max_subdivision_levels=7,
        const bool discard_to_get_nice_topology=false,const bool verbose=true);
private:
    void Write_Diagnostic_Files(const int stepnumber=0);
public:
    void Write_Tri_File_Format(const int stepnumber,const std::string& output_directory); 
//#####################################################################
};
}
#endif
