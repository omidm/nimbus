//#####################################################################
// Copyright 2006-2007, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VOF_ADVECTION
//#####################################################################
#ifndef __VOF_ADVECTION__
#define __VOF_ADVECTION__

#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/LAPLACE_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION_POLICY.h>
#include <PhysBAM_Dynamics/Meshing/RED_GREEN_MESHING_POLICY.h>
#include <PhysBAM_Dynamics/Meshing/RED_GREEN_TETRAHEDRA.h>
#include <PhysBAM_Dynamics/Meshing/RED_GREEN_TRIANGLES.h>
#include <PhysBAM_Dynamics/Particles/PARTICLES_FORWARD.h>
namespace PhysBAM{

template<class T_GRID> class FLUIDS_PARAMETERS_UNIFORM;
template<class T_GRID> class PARTICLE_LEVELSET_IMPLICIT_OBJECT;
template<class TV> class GRID;

template<class TV>
class VOF_ADVECTION:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;typedef typename GRID<TV>::VECTOR_INT TV_INT;
    typedef typename LEVELSET_ADVECTION_POLICY<GRID<TV> >::FAST_LEVELSET_ADVECTION_T T_FAST_LEVELSET_ADVECTION;
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET_IMPLICIT_OBJECT T_LEVELSET_IMPLICIT_OBJECT;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension>::SIMPLEX T_SIMPLEX;typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension>::SIMPLEX_FACE T_SIMPLEX_FACE;
    typedef typename BASIC_GEOMETRY_POLICY<TV>::HYPERPLANE T_HYPERPLANE;typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID<TV>::NODE_ITERATOR NODE_ITERATOR;typedef typename GRID<TV>::FACE_ITERATOR FACE_ITERATOR;typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename INTERPOLATION_POLICY<GRID<TV> >::FACE_LOOKUP T_FACE_LOOKUP;typedef VECTOR<int,TV::dimension+1> T_ELEMENT;
    typedef VECTOR<int,TV::dimension> T_FACE_ELEMENT;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;typedef typename T_ARRAYS_SCALAR::template REBIND<ARRAY<int> >::TYPE T_ARRAYS_ARRAY_INT;
    typedef typename T_ARRAYS_SCALAR::template REBIND<ARRAY<T_ELEMENT> >::TYPE T_ARRAYS_ARRAY_OF_SIMPLICES;
    typedef typename T_ARRAYS_SCALAR::template REBIND<ARRAY<int>*>::TYPE T_ARRAYS_ARRAY_OF_INDICES;
    typedef typename T_ARRAYS_SCALAR::template REBIND<ARRAY<T> >::TYPE T_ARRAYS_ARRAY_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef TRIPLE<TV_INT,int,T> PARTICLE_VOLUME_REFERENCE;
    typedef typename T_ARRAYS_SCALAR::template REBIND<ARRAY<PARTICLE_VOLUME_REFERENCE> >::TYPE T_ARRAYS_ARRAY_PARTICLE_VOLUME_REFERENCE;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<ARRAY<bool> >::TYPE T_ARRAYS_ARRAY_BOOL;

    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<GRID<TV> >::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FLOOD_FILL T_FLOOD_FILL;
    typedef typename INTERPOLATION_POLICY<GRID<TV> >::LINEAR_INTERPOLATION_MAC_HELPER T_LINEAR_INTERPOLATION_MAC_HELPER;
    typedef typename INTERPOLATION_POLICY<GRID<TV> >::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
    typedef typename T_LINEAR_INTERPOLATION_SCALAR::template REBIND<TV>::TYPE T_LINEAR_INTERPOLATION_VECTOR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<PARTICLE_LEVELSET_PARTICLES<TV>*>::TYPE T_ARRAYS_PARTICLE_LEVELSET_PARTICLES;
    typedef typename LEVELSET_POLICY<GRID<TV> >::PARTICLE_LEVELSET T_PARTICLE_LEVELSET;
    typedef typename LEVELSET_POLICY<GRID<TV> >::FAST_LEVELSET_T T_FAST_LEVELSET;

    typedef typename RED_GREEN_MESHING_POLICY<TV>::RED_GREEN_SIMPLICES T_RED_GREEN_SIMPLICES;

    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::dimension>::OBJECT T_OBJECT;typedef typename T_OBJECT::MESH T_MESH;

    PARTICLE_LEVELSET_IMPLICIT_OBJECT<TV>* implicit_object;
    PARTICLE_LEVELSET_UNIFORM<GRID<TV> >& particle_levelset;
    T_FAST_LEVELSET_ADVECTION& levelset_advection;
    GRID<TV>& grid;

    T_ARRAYS_SCALAR grid_nodal_phis;

public:
    T_ARRAYS_BOOL dirichlet_neighbor;
    // object and helpers
    T_OBJECT object;
    T_MESH object_mesh;
    PARTICLES<TV> object_particles;
    T_OBJECT old_object;
    T_MESH old_object_mesh;
    PARTICLES<TV> old_object_particles;
    ARRAY<bool> fixed_segment_list;
    ARRAY<bool> fixed_particle_list;
    ARRAY<ARRAY<TV_INT> > level_simplex_cells;
    ARRAY<int> simplex_list;
    ARRAY<TV> old_particle_X;
    T_RED_GREEN_SIMPLICES preimage;
    ARRAY<T> original_simplex_volumes;

    // helper variables for computing geometric volume in cell
    ARRAY<T_ELEMENT> cell_refinement_simplices;
    ARRAY<T> cell_phis;
    ARRAY<TV> cell_particle_X;

    // helper variables for creating particles from postimage simplices
    T_ARRAYS_ARRAY_OF_INDICES cell_to_postimage_particle_map;
    int last_simplex_particle;
    ARRAY<TV> simplex_particles;
    ARRAY<T_ELEMENT> real_simplices,junk_simplices;

    // helper variables for preimage advection
    ARRAY<TV> preimage_particles; // invisible particles
    ARRAY<T> preimage_particle_material_volumes;
    ARRAY<int> cell_particles;
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,TV> velocity_interpolation;
    ARRAY<T> simplex_preimage_material_volume;
    T_ARRAYS_SCALAR cell_preimage_material_volume;
    ARRAY<T> simplex_preimage_volume;
    ARRAY<T> lower_dimensional_preimage_volume;
    ARRAY<bool> lower_dimensional_preimage;
    T_ARRAYS_INT nodes_to_particles_map;
    ARRAY<bool> material_particles;
    ARRAY<ARRAY<int> > adjacency_list;
    
    T_ARRAYS_BOOL fixed_cells;
    T_ARRAYS_BOOL borders_non_fixed_material_cell;
    T_ARRAYS_BOOL borders_fixed_material_cell;
    T_FACE_ARRAYS_SCALAR V;
    T_FACE_ARRAYS_SCALAR V_ghost;
    int runge_kutta_order_particles;
    bool use_frozen_velocity;
    T_ARRAYS_ARRAY_PARTICLE_VOLUME_REFERENCE removed_negative_particle_volumes_in_cell;
    T_ARRAYS_BOOL cells_valid;

    // helper variables for postimage rasterization
    T_ARRAYS_ARRAY_SCALAR cut_material_simplex_postimage_volumes;
    T_ARRAYS_ARRAY_OF_SIMPLICES simplex_lists;

    // acceleration structure
    T_ARRAYS_ARRAY_INT postimage_simplices_in_cells;
    
public:
    T full_cell_size;
    T length_scale,epsilon;
    int maximum_refinement_depth;
    int minimum_refinement_depth;
    int maximum_material_refinement_depth;
    int minimum_material_refinement_depth;
    FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >* fluids_parameters;

    ARRAY<T> phis;
    ARRAY<T> old_phis;
    bool debugging;
    ARRAY<T_ELEMENT> debug_all_simplices;

    // advection variables
    T_ARRAYS_INT marked_cells;
    T_ARRAYS_SCALAR p;
    T_FACE_ARRAYS_SCALAR volume_fluxes;
    LAPLACE_COLLIDABLE_UNIFORM<GRID<TV> > laplace;
    T_GRID_BASED_COLLISION_GEOMETRY* body_list;
    MPI_UNIFORM_GRID<GRID<TV> >* mpi_grid;

    // volume of material variables
    T_ARRAYS_SCALAR volume_of_material;
    T_ARRAYS_SCALAR geometric_volume;
    bool volume_of_material_initialized;
    bool geometry_initialized;

    T_ARRAYS_ARRAY_INT simplices_in_cells;

    VOF_ADVECTION(PARTICLE_LEVELSET_UNIFORM<GRID<TV> >& particle_levelset_input,T_FAST_LEVELSET_ADVECTION& levelset_advection_input,FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >* fluids_parameters_input=0);
    ~VOF_ADVECTION();

    static void Refined_Object_Initialization_Helper(ARRAY<VECTOR<int,4> >& tets)
    {// corner tets
    tets.Append(VECTOR<int,4>(2,6,4,1));
    tets.Append(VECTOR<int,4>(3,7,1,4));
    tets.Append(VECTOR<int,4>(5,6,1,7));
    tets.Append(VECTOR<int,4>(8,6,7,4));
    // inside tet
    tets.Append(VECTOR<int,4>(1,7,6,4));}

    static void Refined_Object_Initialization_Helper(ARRAY<VECTOR<int,3> >& tris)
    {tris.Append(VECTOR<int,3>(1,2,3));tris.Append(VECTOR<int,3>(2,4,3));}

    static void Initialize_Face_Neighbors(TRIANGULATED_AREA<T>& object)
    {object.mesh.Initialize_Edge_Triangles();}

    static void Initialize_Face_Neighbors(TETRAHEDRALIZED_VOLUME<T>& object)
    {object.mesh.Initialize_Triangle_Tetrahedrons();}

    static void Initialize_Segment_Neighbors(TRIANGULATED_AREA<T>& object)
    {object.mesh.Initialize_Edge_Triangles();}

    static void Initialize_Segment_Neighbors(TETRAHEDRALIZED_VOLUME<T>& object)
    {object.mesh.Initialize_Edge_Tetrahedrons();}

    static int Faces(const TRIANGULATED_AREA<T>& object)
    {return object.mesh.edge_triangles->m;}

    static int Faces(const TETRAHEDRALIZED_VOLUME<T>& object)
    {return object.mesh.triangle_tetrahedrons->m;}

    static VECTOR<int,2> Face_Neighbors(const TRIANGULATED_AREA<T>& object,const int index)
    {VECTOR<int,2> face;ARRAY<int>& neighbors=(*object.mesh.edge_triangles)(index);
    face[1]=neighbors(1);face[2]=(neighbors.m>1?neighbors(2):0);return face;}

    static VECTOR<int,2> Face_Neighbors(const TETRAHEDRALIZED_VOLUME<T>& object,const int index)
    {return (*object.mesh.triangle_tetrahedrons)(index);}

    static ARRAY<int>& Segment_Neighbors(const TRIANGULATED_AREA<T>& object,const int index)
    {return (*object.mesh.edge_triangles)(index);}

    static ARRAY<int>& Segment_Neighbors(const TETRAHEDRALIZED_VOLUME<T>& object,const int index)
    {return (*object.mesh.edge_tetrahedrons)(index);}

    static VECTOR<int,2> Face(const TRIANGULATED_AREA<T>& object,const int index)
    {return (*object.mesh.segment_mesh).elements(index);}

    static VECTOR<int,3> Face(const TETRAHEDRALIZED_VOLUME<T>& object,const int index)
    {return (*object.mesh.triangle_mesh).elements(index);}

    static bool Refinement_Condition(const INDIRECT_ARRAY<ARRAY<T>,T_ELEMENT&>& phis,const INDIRECT_ARRAY<ARRAY_VIEW<TV>,T_ELEMENT&>& X)
    {T maximum_edge_length_magnitude_squared=0; 
    for(int i=1;i<=GRID<TV>::dimension;i++)for(int j=i+1;j<=GRID<TV>::dimension+1;j++)
        maximum_edge_length_magnitude_squared=max(maximum_edge_length_magnitude_squared,(X(i)-X(j)).Magnitude_Squared());
    T result=sqr(phis(1));for(int i=2;i<=GRID<TV>::dimension+1;i++) result=min(result,sqr(phis(i)));
    return result<=maximum_edge_length_magnitude_squared;}

    static bool Refinement_Condition(const ARRAY_VIEW<TV>& X,const ARRAY<T>& phis,const T_ELEMENT& indices)
    {T maximum_edge_length_magnitude_squared=0; 
    for(int i=1;i<=GRID<TV>::dimension;i++)for(int j=i+1;j<=GRID<TV>::dimension+1;j++)
        maximum_edge_length_magnitude_squared=max(maximum_edge_length_magnitude_squared,(X(indices[i])-X(indices[j])).Magnitude_Squared());
    T result=sqr(phis(indices[1]));for(int i=2;i<=GRID<TV>::dimension+1;i++) result=min(result,sqr(phis(indices[i])));
    return result<=maximum_edge_length_magnitude_squared;}

    static T Signed_Size(const VECTOR<int,4>& tet,const ARRAY<VECTOR<T,3> >& node_locations)
    {return TETRAHEDRON<T>::Signed_Volume(node_locations(tet[1]),node_locations(tet[2]),node_locations(tet[3]),node_locations(tet[4]));}

    static T Signed_Size(const VECTOR<int,3>& tri,const ARRAY<VECTOR<T,2> >& node_locations)
    {return TRIANGLE_2D<T>::Signed_Area(node_locations(tri[1]),node_locations(tri[2]),node_locations(tri[3]));}

    static TRIANGLE_3D<T> Simplex_From_Nodes(const VECTOR<int,3>& triangle,const ARRAY<VECTOR<T,3> >& node_locations)
    {return TRIANGLE_3D<T>(node_locations(triangle[1]),node_locations(triangle[2]),node_locations(triangle[3]));}

    static SEGMENT_2D<T> Simplex_From_Nodes(const VECTOR<int,2>& segment,const ARRAY<VECTOR<T,2> >& node_locations)
    {return SEGMENT_2D<T>(node_locations(segment[1]),node_locations(segment[2]));}

    static T Half_Boundary_Measure(const VECTOR<int,4>& tet,const ARRAY<VECTOR<T,3> >& node_locations)
    {return (T).5*(TRIANGLE_3D<T>::Area(node_locations(tet[1]),node_locations(tet[2]),node_locations(tet[3]))+TRIANGLE_3D<T>::Area(node_locations(tet[1]),node_locations(tet[4]),node_locations(tet[2]))+
    TRIANGLE_3D<T>::Area(node_locations(tet[1]),node_locations(tet[3]),node_locations(tet[4]))+TRIANGLE_3D<T>::Area(node_locations(tet[2]),node_locations(tet[4]),node_locations(tet[3])));}

    static T Half_Boundary_Measure(const VECTOR<int,3>& tris,const ARRAY<VECTOR<T,2> >& node_locations)
    {return (T).5*((node_locations(tris[1])-node_locations(tris[2])).Magnitude()+(node_locations(tris[2])-node_locations(tris[3])).Magnitude()+(node_locations(tris[3])-node_locations(tris[1])).Magnitude());}

    void Set_Full_Cell_Size(const T full_cell_size_input)
    {full_cell_size=full_cell_size_input;}

    void Set_Maximum_Refinement_Depth(const int maximum_refinement_depth_input)
    {maximum_refinement_depth=maximum_refinement_depth_input;}

    void Set_Minimum_Refinement_Depth(const int minimum_refinement_depth_input)
    {minimum_refinement_depth=minimum_refinement_depth_input;}

    T Phi(const TV& location) const
    {return (*implicit_object)(location);}

    void Precompute_Cell_Particle_Influence()
    {PARTICLE_LEVELSET_IMPLICIT_OBJECT<TV>* pls_implicit_object=dynamic_cast<PARTICLE_LEVELSET_IMPLICIT_OBJECT<TV>*>(implicit_object);
    if(pls_implicit_object) pls_implicit_object->Precompute_Cell_Particle_Influence();}

    void Use_Frozen_Velocity(const bool use_frozen_velocity_input=true)
    {use_frozen_velocity=use_frozen_velocity_input;}

    void Set_Runge_Kutta_Order_Particles(const int order=2)
    {runge_kutta_order_particles=order;assert(order >= 1 && order <= 3);}

    TV Velocity(const T_FACE_LOOKUP& face_velocities,const TV& location,const T time,const bool use_analytic_velocities)
    {//if(use_analytic_velocities) return fluids_parameters->callbacks->Get_Analytic_Velocity(location,time);
    return velocity_interpolation.Clamped_To_Array_Face(grid,face_velocities,location);}

    void Initialize_Refinement_Target_Segments(const ARRAY<int>& segments,RED_GREEN_TETRAHEDRA<T>& geometry)
    {for(int i=1;i<=segments.m;i++) geometry.refinement_target_segments.Insert(segments(i),true);}

    void Initialize_Refinement_Target_Segments(const ARRAY<int>& segments,RED_GREEN_TRIANGLES<VECTOR<T,2> >& geometry)
    {}

    void Clean_Refinement_Target_Segments(RED_GREEN_TETRAHEDRA<T>& geometry)
    {geometry.refinement_target_segments.Remove_All();}

    void Clean_Refinement_Target_Segments(RED_GREEN_TRIANGLES<VECTOR<T,2> >& geometry)
    {}

    static void Write_VOF_Object(const STREAM_TYPE stream_type,const std::string& vof_object_file,TRIANGULATED_AREA<T>& triangulated_area)
    {FILE_UTILITIES::Write_To_File(stream_type,vof_object_file,triangulated_area);}

    static void Write_VOF_Object(const STREAM_TYPE stream_type,const std::string& vof_object_file,TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume)
    {delete tetrahedralized_volume.triangulated_surface;tetrahedralized_volume.triangulated_surface=0;
    FILE_UTILITIES::Write_To_File(stream_type,vof_object_file,tetrahedralized_volume.Get_Boundary_Object());}

//#####################################################################
    T Negative_Material(const TV_INT& cell_index,const bool force_full_refinement=false);
    template<class T_ARRAYS_PARTICLES> void Precompute_Particle_Volumes_In_Cells(T_ARRAYS_ARRAY_PARTICLE_VOLUME_REFERENCE& particle_volumes_in_cell,T_ARRAYS_PARTICLES& particles,int sign);
    template<class T_ARRAYS_PARTICLES> T Get_Outside_Particle_Volume_In_Cell(T_ARRAYS_ARRAY_PARTICLE_VOLUME_REFERENCE& particle_volumes_in_cell,T_ARRAYS_PARTICLES& particles,const TV_INT& cell_index);
    template<class T_ARRAYS_PARTICLES> void Modify_Particle_Material_Volume_In_Cell(T_ARRAYS_ARRAY_PARTICLE_VOLUME_REFERENCE& particle_volumes_in_cell,T_ARRAYS_PARTICLES& particles,const TV_INT& cell_index,const T volume_of_material_density);
    T Create_Preimage_Particles_From_Old_Postimage_Simplices(const TV_INT& cell_index,const T cell_volume_of_material);
    void Create_Geometry();
    void Cut_Simplicial_Object_With_Zero_Isocontour();
    TV Iterative_Find_Interface(TV left,TV right,const int interations=7) const;
    void Advect_Material_Preimages(const T_FACE_LOOKUP& face_velocities,const T dt,const T start_time,const bool use_analytic_velocities);
    void Rasterize_Material_Postimages();
    void Adjust_Node_For_Domain_Boundaries(TV& node);
    void Make_Approximately_Incompressible(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    void Refine_Or_Coarsen_Geometry();
    void Adjust_Levelset_With_Material_Volumes();
    template<class T_ARRAYS_PARTICLES> void Second_Order_Runge_Kutta_Step_Particles(const T_FACE_LOOKUP& V_lookup,T_ARRAYS_PARTICLES& particles,T_PARTICLE_LEVELSET& particle_levelset,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time);
    template<class T_ARRAYS_PARTICLES> void Second_Order_Runge_Kutta_Step_Particles(const T_FACE_ARRAYS_SCALAR& V,T_ARRAYS_PARTICLES& particles,const T_ARRAYS_VECTOR& center_velocities,const PARTICLE_LEVELSET_PARTICLE_TYPE particle_type,const T dt,const T time);
    void Perform_Conservative_Advection(const T_FACE_LOOKUP& face_velocities,const T dt,const T time);
    void Create_Initial_Meshing_For_Cell(GRID<TV>& grid,TRIANGULATED_AREA<T>& object,const VECTOR<int,2>& cell_index,const ARRAY<int,VECTOR<int,2> >& nodes_to_particles_map,ARRAY<VECTOR<int,2> >& simplex_cells);
    void Create_Initial_Meshing_For_Cell(GRID<TV>& grid,TETRAHEDRALIZED_VOLUME<T>& object,const VECTOR<int,3>& cell_index,const ARRAY<int,VECTOR<int,3> >& nodes_to_particles_map,ARRAY<VECTOR<int,3> >& simplex_cells);
    static void Refine_Simplex(ARRAY<VECTOR<int,3> >& tris,ARRAY<VECTOR<T,2> >& particle_X,const VECTOR<int,3>& indices);
    static void Refine_Simplex(ARRAY<VECTOR<int,4> >& tets,ARRAY<VECTOR<T,3> >& particle_X,const VECTOR<int,4>& indices);
    bool Find_Fixed_Cell_On_Point(const TV& point,TV_INT& opposing_cell,TV& face_normal,const T length_scale) const;
    void Set_Up_For_Refinement();
    void Read(const STREAM_TYPE stream_type,const std::string& input_directory,const std::string& f);
    void Write(const STREAM_TYPE stream_type,const std::string& output_directory,const std::string& f);
//#####################################################################
};

// 1D specialization
template<class T>
class VOF_ADVECTION<VECTOR<T,1> >
{
    typedef VECTOR<T,1> TV;typedef VECTOR<int,1> TV_INT;
    typedef typename INTERPOLATION_POLICY<GRID<TV> >::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename LEVELSET_ADVECTION_POLICY<GRID<TV> >::FAST_LEVELSET_ADVECTION_T T_FAST_LEVELSET_ADVECTION;
public:
    ARRAY<T,FACE_INDEX<1> > volume_fluxes;
    ARRAY<T,VECTOR<int,1> > volume_of_material;
    bool volume_of_material_initialized;
    LAPLACE_COLLIDABLE_UNIFORM<GRID<TV> > laplace;

    virtual ~VOF_ADVECTION(){}
    VOF_ADVECTION(PARTICLE_LEVELSET_UNIFORM<GRID<TV> >& particle_levelset_input,T_FAST_LEVELSET_ADVECTION& levelset_advection, FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >* fluids_parameters_input=0)
        :laplace(particle_levelset_input.levelset.grid,volume_of_material,false,false,true){}

//#####################################################################
    void Use_Frozen_Velocity(const bool use_frozen_velocity_input=true){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    void Set_Runge_Kutta_Order_Particles(const int order=2){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    void Make_Approximately_Incompressible(ARRAY<T,FACE_INDEX<1> >& face_velocities,const T dt,const T time){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    void Create_Geometry(){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    void Refine_Or_Coarsen_Geometry(){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    void Adjust_Levelset_With_Material_Volumes(){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    void Set_Full_Cell_Size(const T full_cell_size_input){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    void Set_Maximum_Refinement_Depth(const int maximum_refinement_depth_input){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    void Set_Minimum_Refinement_Depth(const int minimum_refinement_depth_input){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    T Negative_Material(const VECTOR<int,1>& cell_index,const bool force_full_refinement=false){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    void Set_Up_For_Refinement(){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    void Rasterize_Material_Postimages(){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    void Perform_Conservative_Advection(const T_FACE_LOOKUP& face_velocities,T dt,const T time){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    void Read(const STREAM_TYPE stream_type,const std::string& input_directory,const std::string& f){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    void Write(const STREAM_TYPE stream_type,const std::string& output_directory,const std::string& f) const{PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};
}
#endif
