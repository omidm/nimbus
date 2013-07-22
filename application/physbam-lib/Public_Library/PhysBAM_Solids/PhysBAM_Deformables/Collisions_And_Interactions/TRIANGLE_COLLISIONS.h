//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_COLLISIONS
//#####################################################################
#ifndef __TRIANGLE_COLLISIONS__
#define __TRIANGLE_COLLISIONS__

#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
namespace PhysBAM{

template<class TV> class MPI_SOLIDS;
template<class TV> class TRIANGLE_COLLISION_PARAMETERS;

template<class TV>
class TRIANGLE_COLLISIONS:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    typedef typename BASIC_SIMPLEX_POLICY<TV,d-1>::SIMPLEX T_FACE;
    typedef typename BASIC_SIMPLEX_POLICY<TV,d-1>::SIMPLEX_FACE T_EDGE;
public:
    TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry;
    const ARRAY<T>& repulsion_thickness;
    T final_repulsion_youngs_modulus;
    T final_repulsion_limiter_fraction;
    T collision_thickness;
    T restitution_coefficient;
    bool compute_point_face_collisions,compute_edge_edge_collisions;
    int nonrigid_collision_attempts,rigid_collision_attempts;
    bool limit_rigid_collision_attempts;
    bool output_collision_results;
    ARRAY<VECTOR<int,d+1> > point_face_pairs_internal, point_face_pairs_external;
    ARRAY<VECTOR<int,2*d-2> > edge_edge_pairs_internal, edge_edge_pairs_external;
    // MPI data
    MPI_SOLIDS<TV>* mpi_solids;
    bool use_gauss_jacobi;

    struct GAUSS_JACOBI_PF_DATA
    {
        GAUSS_JACOBI_PF_DATA(TV& target_impulse_input,VECTOR<T,d+1>& target_weight_input,TV& target_normal_input,T& old_speed_input)
            :target_impulse(target_impulse_input),target_weight(target_weight_input),target_normal(target_normal_input),old_speed(old_speed_input)
        {}
        TV& target_impulse;
        VECTOR<T,d+1>& target_weight;
        TV& target_normal;
        T& old_speed;
    };

    struct GAUSS_JACOBI_EE_DATA
    {
        GAUSS_JACOBI_EE_DATA(TV& target_impulse_input,VECTOR<T,2*d-2>& target_weight_input,TV& target_normal_input,T& old_speed_input)
            :target_impulse(target_impulse_input),target_weight(target_weight_input),target_normal(target_normal_input),old_speed(old_speed_input)
        {}
        TV& target_impulse;
        VECTOR<T,2*d-2>& target_weight;
        TV& target_normal;
        T& old_speed;
    };
private:
    int final_point_face_repulsions,final_point_face_collisions,final_edge_edge_repulsions,final_edge_edge_collisions;
    ARRAY<bool> recently_modified_full;

    ARRAY<TV> impulse_velocities;
    ARRAY<TV> pf_target_impulses;
    ARRAY<VECTOR<T,d+1> > pf_target_weights;
    ARRAY<T> pf_old_speeds;
    ARRAY<TV> pf_normals;
    ARRAY<TV> ee_target_impulses;
    ARRAY<VECTOR<T,2*d-2> > ee_target_weights;
    ARRAY<T> ee_old_speeds;
    ARRAY<TV> ee_normals;
public:

    TRIANGLE_COLLISIONS(TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& geometry,const ARRAY<T>& repulsion_thickness);
    ~TRIANGLE_COLLISIONS();

    void Set_Collision_Thickness(const T thickness=1e-6)
    {collision_thickness=thickness;}

    void Set_Restitution_Coefficient(const T restitution_coefficient_input=0) // default is inelastic
    {restitution_coefficient=restitution_coefficient_input;}

    void Compute_Point_Face_Collisions(const bool compute=true,const int nonrigid_collision_attempts=3)
    {compute_point_face_collisions=compute;Set_Attempts_For_Nonrigid_Collisions(nonrigid_collision_attempts);}

    void Compute_Edge_Edge_Collisions(const bool compute=true,const int nonrigid_collision_attempts=3)
    {compute_edge_edge_collisions=compute;Set_Attempts_For_Nonrigid_Collisions(nonrigid_collision_attempts);}

    void Set_Attempts_For_Nonrigid_Collisions(const int attempts=3) // after this switch to rigid collisions
    {nonrigid_collision_attempts=attempts;}

    void Set_Attempts_For_Rigid_Collisions(const bool limit_attempts=true,const int attempts=10)
    {limit_rigid_collision_attempts=limit_attempts;rigid_collision_attempts=attempts;}

    void Output_Collision_Results(const bool output=true)
    {output_collision_results=output;}

    void Set_Gauss_Jacobi(const bool use_gauss_jacobi_input=false)
    {use_gauss_jacobi=use_gauss_jacobi_input;}

//#####################################################################
    int Adjust_Velocity_For_Self_Collisions(const T dt,const T time,const bool exit_early);
    void Initialize(TRIANGLE_COLLISION_PARAMETERS<TV>& triangle_collision_parameters);
    void Update_Swept_Hierachies_And_Compute_Pairs(ARRAY_VIEW<TV> X,ARRAY_VIEW<TV> X_self_collision_free,ARRAY_VIEW<bool> recently_modified,const T detection_thickness);
private:
    // TODO: rename from "get" to "append"...
    void Scale_And_Apply_Impulses();
    void Get_Moving_Faces_Near_Moving_Points(STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1,STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2,ARRAY<VECTOR<int,d+1> >& pairs_internal,ARRAY<VECTOR<int,d+1> >& pairs_external,const T detection_thickness);
    void Get_Moving_Edges_Near_Moving_Edges(STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_1,STRUCTURE_INTERACTION_GEOMETRY<TV>& structure_2,ARRAY<VECTOR<int,2*d-2> >& pairs_internal,ARRAY<VECTOR<int,2*d-2> >& pairs_external,const T detection_thickness);
    int Adjust_Velocity_For_Point_Face_Collision(const T dt,const bool rigid,ARRAY<ARRAY<int> >& rigid_lists,ARRAY<int>& list_index,const ARRAY<VECTOR<int,d+1> >& pairs,
        const T attempt_ratio,const bool final_repulsion_only,const bool exit_early);
    bool Point_Face_Collision(GAUSS_JACOBI_PF_DATA& pf_data,const VECTOR<int,d+1>& nodes,const T dt,const T repulsion_thickness,T& collision_time,const T attempt_ratio,const bool exit_early);
    bool Point_Face_Final_Repulsion(GAUSS_JACOBI_PF_DATA& pf_data,const VECTOR<int,d+1>& nodes,const T dt,const T repulsion_thickness,T& collision_time,const T attempt_ratio,const bool exit_early);
    int Adjust_Velocity_For_Edge_Edge_Collision(const T dt,const bool rigid,ARRAY<ARRAY<int> >& rigid_lists,ARRAY<int>& list_index,const ARRAY<VECTOR<int,2*d-2> >& pairs,
        const T attempt_ratio,const bool final_repulsion_only,const bool exit_early);
    int Prune_Non_Intersecting_Pairs(const T dt,ARRAY<VECTOR<int,d+1> >& point_face_pairs,ARRAY<VECTOR<int,2*d-2> >& edge_edge_pairs,const T attempt_ratio);
    bool Edge_Edge_Collision(GAUSS_JACOBI_EE_DATA& ee_data,const VECTOR<int,2*d-2>& nodes,const T dt,T& collision_time,const T attempt_ratio,const bool exit_early);
    bool Edge_Edge_Final_Repulsion(GAUSS_JACOBI_EE_DATA& ee_data,const VECTOR<int,2*d-2>& nodes,const T dt,T& collision_time,const T attempt_ratio,const bool exit_early);
    template<int d2> void Add_To_Rigid_Lists(ARRAY<ARRAY<int> >& rigid_lists,ARRAY<int>& list_index,const VECTOR<int,d2>& nodes);
    void Apply_Rigid_Body_Motions(const T dt,ARRAY<ARRAY<int> >& rigid_lists);
public: // subdivision stuff
    void Stop_Nodes_Before_Self_Collision(const T dt);
//#####################################################################
};
}
#endif
