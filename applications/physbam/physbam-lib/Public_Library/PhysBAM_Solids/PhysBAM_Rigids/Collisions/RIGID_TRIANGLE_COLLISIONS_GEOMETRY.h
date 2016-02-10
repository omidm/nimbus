//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_TRIANGLE_COLLISIONS_GEOMETRY
//#####################################################################
#ifndef __RIGID_TRIANGLE_COLLISIONS_GEOMETRY__
#define __RIGID_TRIANGLE_COLLISIONS_GEOMETRY__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
namespace PhysBAM{

template<class TV> class BINDING;
template<class TV> class RIGID_STRUCTURE_INTERACTION_GEOMETRY;
template<class TV> class RIGID_BODY_STATE;
template<class TV> class TRIANGLE_COLLISION_PARAMETERS;
template<class TV> class STRUCTURE;
template<class TV> class MPI_RIGIDS;

template<class TV>
class RIGID_TRIANGLE_COLLISIONS_GEOMETRY:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    struct UNUSABLE{};
public:
    MPI_RIGIDS<TV>* mpi_solids;
    ARRAY<STRUCTURE<TV>*> structures;
    ARRAY<RIGID_STRUCTURE_INTERACTION_GEOMETRY<TV>*> structure_geometries;
    ARRAY<VECTOR<int,2> > interacting_structure_pairs;

    bool allow_intersections;
    T allow_intersections_tolerance;
    T small_number;

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

    RIGID_TRIANGLE_COLLISIONS_GEOMETRY();
    ~RIGID_TRIANGLE_COLLISIONS_GEOMETRY();

    void Set_Allow_Intersections_Tolerance(const T allow_intersections_tolerance_input=1e-8)
    {allow_intersections_tolerance=allow_intersections_tolerance_input;}

    void Set_Small_Number(const T small_number_input=1e-8)
    {small_number=small_number_input;}

    void Output_Number_Checked(const bool output=true)
    {output_number_checked=output;}

//#####################################################################
    void Build_Collision_Geometry();
    void Build_Topological_Structure_Of_Hierarchies();
    void Allow_Intersections(const bool allow_intersections_input=true);
    void Save_Current_State(const ARRAY_VIEW<TV>& X,const ARRAY_VIEW<ROTATION<TV> >& roation,const ARRAY_VIEW<TV>& V);
    void Save_Self_Collision_Free_State(const ARRAY_VIEW<TV>& X,const ARRAY_VIEW<ROTATION<TV> >& roation,const ARRAY_VIEW<TV>& V); // assumes mass does not change
    void Save_Self_Collision_Free_State(); // assumes mass does not change
    void Restore_Self_Collision_Free_State();
    void Compute_Intersecting_Segment_Face_Pairs();
private:
//#####################################################################
};
}
#endif
