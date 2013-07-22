//#####################################################################
// Copyright 2005-2008, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COLLISION_BODY
//#####################################################################
#ifndef __COLLISION_BODY__
#define __COLLISION_BODY__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_POLICY.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_POLICY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/ONE.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Collisions/COLLISIONS_GEOMETRY_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TETRAHEDRON_COLLISION_BODY.h>

namespace PhysBAM{

template<class T_GRID,class ID> class OBJECTS_IN_CELL;
template<class TV> class BINDING_LIST;
template<class TV> class SOFT_BINDINGS;
template<class TV> class PARTICLES;
template<class TV> class FRAME;
template<class TV> class RANGE;
template<class TV> class COLLISION_PARTICLE_STATE;
template<class TV> class COLLISION_GEOMETRY;
template<class TV> class GRID;

template<class TV>
class COLLISION_BODY_HELPER
{
private:
    typedef typename TV::SCALAR T;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    typedef typename DYADIC_GRID_POLICY<TV>::DYADIC_GRID T_DYADIC_GRID;
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    typedef typename RLE_GRID_POLICY<TV>::RLE_GRID T_RLE_GRID;
#endif
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename REBIND<T_ARRAYS_SCALAR,bool>::TYPE T_ARRAYS_BOOL;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename REBIND<T_FACE_ARRAYS_SCALAR,int>::TYPE T_FACE_ARRAYS_INT;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX T_SIMPLEX;
    typedef typename IF<TV::dimension==2,T,typename IF<TV::dimension==1,ONE,TV>::TYPE>::TYPE T_WEIGHTS;
public:
    static int Adjust_Nodes_For_Collisions(COLLISION_GEOMETRY<TV>& body,ARRAY_VIEW<const TV> X_old,PARTICLES<TV>& collision_particles,SOFT_BINDINGS<TV>& soft_bindings,
        const ARRAY<int>& nodes_to_check,const ARRAY<bool>& particle_on_surface,const T collision_tolerance,ARRAY<COLLISION_PARTICLE_STATE<TV> >& collision_particle_state,
        ARRAY<COLLISION_GEOMETRY_ID>& particle_to_collision_geometry_id,const T max_relative_velocity,const T dt,const HASHTABLE<int,T> *friction_table,
        const HASHTABLE<int,T> *thickness_table);
    static void Adjust_Point_For_Collision(COLLISION_GEOMETRY<TV>& body,const TV& X_old,TV& X,TV& V,const T point_mass,const T penetration_depth,const T dt,const T one_over_dt,
        const T max_relative_velocity,COLLISION_PARTICLE_STATE<TV>& collision,T local_coefficient_of_friction);
    static void Adjust_Point_For_Collision(COLLISION_GEOMETRY<TV>& body,TV& X,const T penetration_depth);
};

template<class TV>
class COLLISION_BODY
{
public:
    typedef typename TV::SCALAR T;
    static int Adjust_Nodes_For_Collisions(COLLISION_GEOMETRY<TV>& body, ARRAY_VIEW<const TV> X_old,PARTICLES<TV>& collision_particles,SOFT_BINDINGS<TV>& soft_bindings,
        const ARRAY<int>& nodes_to_check,const ARRAY<bool>& particle_on_surface,const T collision_tolerance,ARRAY<COLLISION_PARTICLE_STATE<TV> >& collision_particle_state,
        ARRAY<COLLISION_GEOMETRY_ID>& particle_to_collision_geometry_id,const T max_relative_velocity,const T dt,const HASHTABLE<int,T> *friction_table,
        const HASHTABLE<int,T> *thickness_table)
    {
        return COLLISION_BODY_HELPER<TV>::Adjust_Nodes_For_Collisions(body, X_old, collision_particles,soft_bindings,nodes_to_check,particle_on_surface,collision_tolerance,
            collision_particle_state,particle_to_collision_geometry_id,max_relative_velocity,dt,friction_table,thickness_table);
    }
};

template<class T>
class COLLISION_BODY<VECTOR<T,3> >
{
public:
    typedef VECTOR<T,3> TV;
    static int Adjust_Nodes_For_Collisions(COLLISION_GEOMETRY<TV>& body, ARRAY_VIEW<const TV> X_old,PARTICLES<TV>& collision_particles,SOFT_BINDINGS<TV>& soft_bindings,
        const ARRAY<int>& nodes_to_check,const ARRAY<bool>& particle_on_surface,const T collision_tolerance,ARRAY<COLLISION_PARTICLE_STATE<TV> >& collision_particle_state,
        ARRAY<COLLISION_GEOMETRY_ID>& particle_to_collision_geometry_id,const T max_relative_velocity,const T dt,const HASHTABLE<int,T> *friction_table,
        const HASHTABLE<int,T> *thickness_table)
    {
            if(TETRAHEDRON_COLLISION_BODY<T>* tetrahedron_collision_body=dynamic_cast<TETRAHEDRON_COLLISION_BODY<T>*>(&body))
                return tetrahedron_collision_body->Adjust_Nodes_For_Collisions(X_old,collision_particles,soft_bindings,nodes_to_check,particle_on_surface,
                    collision_tolerance,collision_particle_state,particle_to_collision_geometry_id,max_relative_velocity,dt,friction_table,thickness_table);
            else
                return COLLISION_BODY_HELPER<TV>::Adjust_Nodes_For_Collisions(body, X_old, collision_particles,soft_bindings,nodes_to_check,particle_on_surface,collision_tolerance,
                    collision_particle_state,particle_to_collision_geometry_id,max_relative_velocity,dt,friction_table,thickness_table);
    }
};

}

#endif
