//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Craig Schroeder, Andrew Selle, Tamar Shinar, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_REPULSIONS
//#####################################################################
#ifndef __TRIANGLE_REPULSIONS__
#define __TRIANGLE_REPULSIONS__

#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Parallel_Computation/PARTITION_ID.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLES_COLLISIONS_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
namespace PhysBAM{

template<class TV> class MPI_SOLIDS;
template<class TV> class TRIANGLE_COLLISION_PARAMETERS;

template<class TV>
struct POINT_FACE_REPULSION_PAIR
{
    typedef typename TV::SCALAR T;
    typedef T SCALAR;
    enum WORKAROUND {d=TV::m,count=TV::m+1};

    VECTOR<int,d+1> nodes; // point,node1,node2,(node3) (order implies side)
    T distance;
    TV weights;
    TV normal;

    static T Total_Repulsion_Thickness(ARRAY_VIEW<const T> repulsion_thickness,const VECTOR<int,d+1>& nodes)
    {return max(repulsion_thickness(nodes[1]),ARRAYS_COMPUTATIONS::Min(repulsion_thickness.Subset(nodes.Remove_Index(1))));}

    T Total_Repulsion_Thickness(ARRAY_VIEW<const T> repulsion_thickness) const
    {return Total_Repulsion_Thickness(repulsion_thickness,nodes);}
};

template<class T>
struct EDGE_EDGE_REPULSION_PAIR<VECTOR<T,1> >
{
    typedef T SCALAR;
    enum WORKAROUND {d=1,count=1};
    typedef VECTOR<T,1> TV;

    VECTOR<int,1> nodes; // vertex nodes
    TV weights;
    T distance;

    static T Total_Repulsion_Thickness(ARRAY_VIEW<const T> repulsion_thickness,const VECTOR<int,1>& nodes)
    {return repulsion_thickness(nodes[1]);}

    T Total_Repulsion_Thickness(ARRAY_VIEW<const T> repulsion_thickness) const
    {return 0;}
};

template<class T>
struct EDGE_EDGE_REPULSION_PAIR<VECTOR<T,2> >
{
    typedef T SCALAR;
    enum WORKAROUND {d=2,count=2};
    typedef VECTOR<T,2> TV;

    VECTOR<int,2> nodes; // vertex nodes
    TV collision_free_normal;
    T distance;
    static VECTOR<T,2> weights;
    TV normal;

    static T Total_Repulsion_Thickness(ARRAY_VIEW<const T> repulsion_thickness,const VECTOR<int,2>& nodes)
    {return max(repulsion_thickness(nodes[1]),repulsion_thickness(nodes[2]));}

    T Total_Repulsion_Thickness(ARRAY_VIEW<const T> repulsion_thickness) const
    {return Total_Repulsion_Thickness(repulsion_thickness,nodes);}
};

template<class T>
struct EDGE_EDGE_REPULSION_PAIR<VECTOR<T,3> >
{
    typedef T SCALAR;
    enum WORKAROUND {d=3,count=4};
    typedef VECTOR<T,3> TV;

    VECTOR<int,4> nodes; // segment nodes
    TV collision_free_normal;
    T distance;
    VECTOR<T,2> weights;
    TV normal;

    static T Total_Repulsion_Thickness(ARRAY_VIEW<const T> repulsion_thickness,const VECTOR<int,4>& nodes)
    {return max(min(repulsion_thickness(nodes[1]),repulsion_thickness(nodes[2])),min(repulsion_thickness(nodes[3]),repulsion_thickness(nodes[4])));}

    T Total_Repulsion_Thickness(ARRAY_VIEW<const T> repulsion_thickness) const
    {return Total_Repulsion_Thickness(repulsion_thickness,nodes);}
};

template<class TV>
struct PRECOMPUTE_PROJECT_POINT_FACE
{
    typedef typename TV::SCALAR T;

    VECTOR<TV,4> v_scaled_normals;
    TV weights,normal;
    VECTOR<int,4> nodes;

    void Project(INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,4>&> v) const
    {T s=TV::Dot_Product(v(1)-weights(1)*v(2)-weights(2)*v(3)-weights(3)*v(4),normal);
    for(int i=1;i<=4;i++) v(i)+=s*v_scaled_normals(i);}

    void Precompute(const INDIRECT_ARRAY<ARRAY_VIEW<T>,VECTOR<int,4>&> mass,const TV& weights_input,const TV& normal_input);
};

template<class TV>
struct PRECOMPUTE_PROJECT_EDGE_EDGE
{
    typedef typename TV::SCALAR T;

    VECTOR<TV,4> v_scaled_normals;
    TV normal_13,normal_12,normal_34,normal;
    VECTOR<int,4> nodes;
    VECTOR<T,2-(TV::m==1)> weights;

    void Project(INDIRECT_ARRAY<ARRAY_VIEW<TV>,VECTOR<int,4>&> v) const
    {T s=TV::Dot_Product(v(1)-v(3),normal_13)+TV::Dot_Product(v(2)-v(1),normal_12)+TV::Dot_Product(v(4)-v(3),normal_34);
    for(int i=1;i<=4;i++) v(i)+=s*v_scaled_normals(i);}

    void Precompute(const INDIRECT_ARRAY<ARRAY_VIEW<T>,VECTOR<int,4>&> mass,const VECTOR<T,2>& weights_input,const TV& normal_input);
};

template<class TV>
class TRIANGLE_REPULSIONS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    typedef typename BASIC_SIMPLEX_POLICY<TV,d-1>::SIMPLEX T_FACE;
    typedef typename BASIC_SIMPLEX_POLICY<TV,d-1>::SIMPLEX_FACE T_EDGE;
public:
    TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry;
    ARRAY<T> repulsion_thickness; // sized to full_particles
    T youngs_modulus;
    T spring_limiter_fraction;
    T friction_coefficient;
    bool compute_point_face_friction,compute_edge_edge_friction;
    bool compute_point_face_inelastic_collision_repulsion,compute_edge_edge_inelastic_collision_repulsion;
    bool compute_point_face_repulsion,compute_edge_edge_repulsion;
    int point_face_inelastic_collision_repulsion_attempts,edge_edge_inelastic_collision_repulsion_attempts;
    bool output_repulsion_results;
    bool perform_attractions;
    T attractions_threshold;
    T hierarchy_repulsion_thickness_multiplier;
    T repulsion_thickness_detection_multiplier;
    HASHTABLE<VECTOR<int,d+1> > omit_point_face_repulsion_pairs;
    HASHTABLE<VECTOR<int,2*d-2> > omit_edge_edge_repulsion_pairs;
    // repulsion interaction cache
    ARRAY<POINT_FACE_REPULSION_PAIR<TV> > point_face_interaction_pairs;
    ARRAY<EDGE_EDGE_REPULSION_PAIR<TV> > edge_edge_interaction_pairs;
    ARRAY<PRECOMPUTE_PROJECT_POINT_FACE<TV> > internal_point_face_precomputed;
    ARRAY<PRECOMPUTE_PROJECT_EDGE_EDGE<TV> > internal_edge_edge_precomputed;
    // MPI data
    MPI_SOLIDS<TV>* mpi_solids;
    ARRAY<ARRAY<int>,PARTITION_ID> point_face_send_particles,point_face_receive_particles,edge_edge_send_particles,edge_edge_receive_particles;
    ARRAY<POINT_FACE_REPULSION_PAIR<TV> > point_face_boundary_pairs,point_face_internal_pairs;
    ARRAY<EDGE_EDGE_REPULSION_PAIR<TV> > edge_edge_boundary_pairs,edge_edge_internal_pairs;
    bool use_gauss_jacobi;
    ARRAY<TV> impulse_velocities;
    ARRAY<TV> pf_target_impulses;
    ARRAY<T> pf_old_speeds;
    ARRAY<TV> pf_normals;
    ARRAY<TV> ee_target_impulses;
    ARRAY<VECTOR<T,2*d-2> > ee_target_weights;
    ARRAY<T> ee_old_speeds;
    ARRAY<TV> ee_normals;

    TRIANGLE_REPULSIONS(TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry);
    ~TRIANGLE_REPULSIONS();

    // Note that if repulsion thickness is not constant, repulsion will be discontinuous since we pick the min repulsion among repelling elements.
    void Set_Repulsion_Thickness(const T thickness=(T)1e-3)
    {repulsion_thickness.Resize(geometry.deformable_body_collection.particles.array_collection->Size(),false,false);ARRAYS_COMPUTATIONS::Fill(repulsion_thickness,thickness);}

    void Set_Repulsion_Thickness(ARRAY_VIEW<const T> thickness)
    {repulsion_thickness.Resize(geometry.deformable_body_collection.particles.array_collection->Size(),false,false);ARRAY<T>::Copy(thickness,repulsion_thickness);}

    void Clamp_Repulsion_Thickness(ARRAY_VIEW<const T> max_value)
    {repulsion_thickness.Resize(geometry.deformable_body_collection.particles.array_collection->Size());
    for(int k=1;k<=repulsion_thickness.m;k++) repulsion_thickness(k)=min(repulsion_thickness(k),max_value(k));}

    void Clamp_Repulsion_Thickness_With_Meshes(const T scale=(T).4)
    {Clamp_Repulsion_Thickness_With_Meshes(geometry.deformable_body_collection.particles.X,scale);}

    void Set_Friction_Coefficient(const T friction_coefficient_input=(T).1)
    {friction_coefficient=friction_coefficient_input;}

    void Compute_Point_Face_Friction(const bool compute=true)
    {compute_point_face_friction=compute;}

    void Compute_Edge_Edge_Friction(const bool compute=true)
    {compute_edge_edge_friction=compute;}

    void Compute_Point_Face_Inelastic_Collision_Repulsion(const bool compute=true,const int attempts=3)
    {compute_point_face_inelastic_collision_repulsion=compute;Set_Attempts_For_Point_Face_Inelastic_Collision_Repulsion(attempts);}

    void Set_Attempts_For_Point_Face_Inelastic_Collision_Repulsion(const int attempts=3)
    {point_face_inelastic_collision_repulsion_attempts=attempts;}

    void Compute_Edge_Edge_Inelastic_Collision_Repulsion(const bool compute=true,const int attempts=3)
    {compute_edge_edge_inelastic_collision_repulsion=compute;Set_Attempts_For_Edge_Edge_Inelastic_Collision_Repulsion(attempts);}

    void Set_Attempts_For_Edge_Edge_Inelastic_Collision_Repulsion(const int attempts=3)
    {edge_edge_inelastic_collision_repulsion_attempts=attempts;}

    void Compute_Point_Face_Repulsion(const bool compute=true)
    {compute_point_face_repulsion=compute;}

    void Compute_Edge_Edge_Repulsion(const bool compute=true)
    {compute_edge_edge_repulsion=compute;}

    void Output_Repulsion_Results(const bool output=true)
    {output_repulsion_results=output;}

    void Set_Gauss_Jacobi(const bool use_gauss_jacobi_input=false)
    {use_gauss_jacobi=use_gauss_jacobi_input;}

    static void Project_All_Moving_Constraints(const ARRAY<PRECOMPUTE_PROJECT_POINT_FACE<VECTOR<T,1> > >& point_face_precomputed,
        const ARRAY<PRECOMPUTE_PROJECT_EDGE_EDGE<VECTOR<T,1> > >& edge_edge_precomputed,ARRAY_VIEW<VECTOR<T,1> > field){}

    static void Project_All_Moving_Constraints(const ARRAY<PRECOMPUTE_PROJECT_POINT_FACE<VECTOR<T,2> > >& point_face_precomputed,
        const ARRAY<PRECOMPUTE_PROJECT_EDGE_EDGE<VECTOR<T,2> > >& edge_edge_precomputed,ARRAY_VIEW<VECTOR<T,2> > field){}

    void Set_Collision_Pairs(ARRAY<PRECOMPUTE_PROJECT_POINT_FACE<VECTOR<T,1> > >& point_face_precomputed,
        ARRAY<PRECOMPUTE_PROJECT_EDGE_EDGE<VECTOR<T,1> > >& edge_edge_precomputed,ARRAY<POINT_FACE_REPULSION_PAIR<VECTOR<T,1> > >& point_face_pairs,
        ARRAY<EDGE_EDGE_REPULSION_PAIR<VECTOR<T,1> > >& edge_edge_pairs,const T repulsion_thickness_multiplier){}

    void Set_Collision_Pairs(ARRAY<PRECOMPUTE_PROJECT_POINT_FACE<VECTOR<T,2> > >& point_face_precomputed,
        ARRAY<PRECOMPUTE_PROJECT_EDGE_EDGE<VECTOR<T,2> > >& edge_edge_precomputed,ARRAY<POINT_FACE_REPULSION_PAIR<VECTOR<T,2> > >& point_face_pairs,
        ARRAY<EDGE_EDGE_REPULSION_PAIR<VECTOR<T,2> > >& edge_edge_pairs,const T repulsion_thickness_multiplier){}

//#####################################################################
    void Initialize(TRIANGLE_COLLISION_PARAMETERS<TV>& triangle_collision_parameters);
    void Clean_Memory();
    void Clamp_Repulsion_Thickness_With_Meshes(ARRAY_VIEW<const TV> X,const T scale=(T).4);
    void Turn_Off_Repulsions_Based_On_Current_Proximity(const T extra_factor_on_distance=(T)1.5);
    void Update_Faces_And_Hierarchies_With_Collision_Free_Positions(const ARRAY_VIEW<TV>* X_other);
    void Compute_Interaction_Pairs(ARRAY_VIEW<const TV> X_other);
    int Adjust_Velocity_For_Self_Repulsion(const T dt,bool use_saved_pairs);
    template<class T_ARRAY1,class T_ARRAY2> void Update_Repulsion_Pairs_Using_History(T_ARRAY1& point_face_interaction_pairs,T_ARRAY2& edge_edge_interaction_pairs,bool prune_separating);
    int Adjust_Velocity_For_Self_Repulsion_Using_History(const T dt,const bool use_repulsions,bool use_saved_pairs);
    template<class T_ARRAY1,class T_ARRAY2> int Apply_Repulsions_To_Velocities(const T dt,T_ARRAY1& point_face_interaction_pairs,T_ARRAY2& edge_edge_interaction_pairs,
        const bool use_repulsions,bool use_saved_pairs);
    template<class T_ARRAY1,class T_ARRAY2> int Apply_Repulsions_To_Velocities(const T dt,T_ARRAY1& point_face_boundary_pairs,T_ARRAY2& edge_edge_boundary_pairs,
        T_ARRAY1& point_face_internal_pairs,T_ARRAY2& edge_edge_internal_pairs,const bool use_repulsions);
    void Output_Interaction_Pairs(const STREAM_TYPE stream_type,const std::string& filename) const;
    static void Project_All_Moving_Constraints(const ARRAY<PRECOMPUTE_PROJECT_POINT_FACE<VECTOR<T,3> > >& point_face_precomputed,
        const ARRAY<PRECOMPUTE_PROJECT_EDGE_EDGE<VECTOR<T,3> > >& edge_edge_precomputed,ARRAY_VIEW<VECTOR<T,3> >& field);
    template<class TV2>
    void Set_Collision_Pairs(ARRAY<PRECOMPUTE_PROJECT_POINT_FACE<VECTOR<T,3> > >& point_face_precomputed,
        ARRAY<PRECOMPUTE_PROJECT_EDGE_EDGE<VECTOR<T,3> > >& edge_edge_precomputed,ARRAY<POINT_FACE_REPULSION_PAIR<VECTOR<T,3> > >& point_face_pairs,
        ARRAY<EDGE_EDGE_REPULSION_PAIR<TV2> >& edge_edge_pairs,const T repulsion_thickness_multiplier);
    template<class T_ARRAY> void Scale_And_Apply_Point_Face_Impulses(const T_ARRAY& pairs);
    template<class T_ARRAY> void Scale_And_Apply_Edge_Edge_Impulses(const T_ARRAY& pairs);
private:
    int Get_Faces_Near_Points(STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1,STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2,ARRAY_VIEW<const TV> X_other,const bool use_processor_cull);
    int Get_Edges_Near_Edges(STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1,STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2,ARRAY_VIEW<const TV> X_other,const bool use_processor_cull);
    template<class T_PAIR> T Repulsion_Impulse(TV& direction,const T dt,const T_PAIR& pair,const TV& relative_velocity,const bool elastic_repulsion,const bool friction)
        PHYSBAM_ALWAYS_INLINE PHYSBAM_FLATTEN;
    template<class T_ARRAY> void Adjust_Velocity_For_Point_Face_Repulsion(const T dt,const T_ARRAY& pairs,const bool elastic_repulsion,const bool friction,bool use_repulsions);
    template<class T_ARRAY,class S> void Adjust_Velocity_For_Edge_Edge_Repulsion_Helper(const T dt,const T_ARRAY& pairs,const bool elastic_repulsion,const bool friction,const VECTOR<S,2>&,bool use_repulsions);
    template<class T_ARRAY,class S> void Adjust_Velocity_For_Edge_Edge_Repulsion_Helper(const T dt,const T_ARRAY& pairs,const bool elastic_repulsion,const bool friction,const VECTOR<S,3>&,bool use_repulsions);
    template<class T_ARRAY> void Adjust_Velocity_For_Edge_Edge_Repulsion(const T dt,const T_ARRAY& pairs,const bool elastic_repulsion,const bool friction,bool use_repulsions);
//#####################################################################
};
}
#include <PhysBAM_Solids/PhysBAM_Deformables/Read_Write/Collisions_And_Interactions/READ_WRITE_EDGE_EDGE_REPULSION_PAIR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Read_Write/Collisions_And_Interactions/READ_WRITE_POINT_FACE_REPULSION_PAIR.h>
#endif
