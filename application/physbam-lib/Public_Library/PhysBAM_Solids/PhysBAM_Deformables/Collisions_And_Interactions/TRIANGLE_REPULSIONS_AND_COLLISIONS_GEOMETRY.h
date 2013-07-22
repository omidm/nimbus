//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY
//#####################################################################
#ifndef __TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY__
#define __TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_OBJECT_FORWARD.h>
namespace PhysBAM{

template<class TV> class BINDING;
template<class TV> class STRUCTURE_INTERACTION_GEOMETRY;
template<class TV> class RIGID_BODY_STATE;
template<class TV> class TRIANGLE_COLLISION_PARAMETERS;
template<class TV> class STRUCTURE;
template<class TV> class DEFORMABLE_BODY_COLLECTION;
template<class TV> class MPI_SOLIDS;

template<class TV>
class TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    struct UNUSABLE{};
public:
    MPI_SOLIDS<TV>* mpi_solids;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection;
    ARRAY<STRUCTURE<TV>*> structures;
    ARRAY<STRUCTURE_INTERACTION_GEOMETRY<TV>*> structure_geometries;
    ARRAY<VECTOR<int,2> > interacting_structure_pairs;
    ARRAY<TV> X_self_collision_free,V_self_collision_free;
    //ARRAY<RIGID_BODY_STATE<TV> > rigid_body_particle_state_collision_free;
    ARRAY<bool> modified_full;
    HASHTABLE<VECTOR<int,d+1> > intersecting_point_face_pairs;
    HASHTABLE<VECTOR<int,2*d-2> > intersecting_edge_edge_pairs;

    bool allow_intersections;
    T allow_intersections_tolerance;
    T small_number;
    bool use_gauss_jacobi;

    bool output_number_checked;

    struct MASS_MODIFIER
    {
        virtual ~MASS_MODIFIER(){}
        virtual void Point_Face_Mass(const T attempt_ratio,const VECTOR<int,d+1>& nodes,const VECTOR<T,d>& weights,VECTOR<T,d+1>& one_over_mass)=0;
        virtual void Point_Face_Mass(const T attempt_ratio,const VECTOR<int,d+1>& nodes,const VECTOR<T,d>& weights,ARRAY_VIEW<T>& one_over_mass)=0;
        virtual void Point_Face_Mass_Revert(const VECTOR<int,d+1>& nodes,ARRAY_VIEW<T>& one_over_mass)=0;
        virtual void Edge_Edge_Mass(const T attempt_ratio,const VECTOR<int,2*d-2>& nodes,const VECTOR<T,2>& weights,VECTOR<T,2*d-2>& one_over_mass)=0;
        virtual void Edge_Edge_Mass(const T attempt_ratio,const VECTOR<int,2*d-2>& nodes,const VECTOR<T,2>& weights,ARRAY_VIEW<T>& one_over_mass)=0;
        virtual void Edge_Edge_Mass_Revert(const VECTOR<int,2*d-2>& nodes,ARRAY_VIEW<T>& one_over_mass)=0;
        virtual void Reorder_Pairs(ARRAY<VECTOR<int,2*d-2> >& edge_edge_pairs,ARRAY<VECTOR<int,d+1> >& point_face_pairs)=0;
    };
    MASS_MODIFIER* mass_modifier;

    TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY(DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection);
    ~TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY();

    void Set_Allow_Intersections_Tolerance(const T allow_intersections_tolerance_input=1e-8)
    {allow_intersections_tolerance=allow_intersections_tolerance_input;}

    void Set_Small_Number(const T small_number_input=1e-8)
    {small_number=small_number_input;}

    void Output_Number_Checked(const bool output=true)
    {output_number_checked=output;}

    void Set_Gauss_Jacobi(const bool use_gauss_jacobi_input=false)
    {use_gauss_jacobi=use_gauss_jacobi_input;}

//#####################################################################
    void Initialize(TRIANGLE_COLLISION_PARAMETERS<TV>& triangle_collision_parameters);
    void Build_Collision_Geometry();
    void Build_Topological_Structure_Of_Hierarchies();
    void Allow_Intersections(const bool allow_intersections_input=true);
    bool Check_For_Intersection(const bool grow_thickness_to_find_first_self_intersection=false,const T thickness=0,VECTOR<int,2>* interaction_pair=0) const;
    bool Check_For_Intersection(const bool grow_thickness_to_find_first_self_intersection,const T thickness,ARRAY<VECTOR<int,2> >& interaction_pairs) const;
    void Save_Self_Collision_Free_State(); // assumes mass does not change
    void Restore_Self_Collision_Free_State();
    void Compute_Intersecting_Segment_Face_Pairs();
private:
    template<class S> void Compute_Intersecting_Pairs_Helper(typename IF<NOT<INTS_EQUAL<d,1>::value>::value,const TV*,UNUSABLE>::TYPE input);
    template<class S> void Compute_Intersecting_Pairs_Helper(typename IF<INTS_EQUAL<d,1>::value,const TV*,UNUSABLE>::TYPE input);
//#####################################################################
};
}
#endif
