//#####################################################################
// Copyright 2008, Craig Schroeder, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FRACTURE_REGION
//##################################################################### 
#ifndef __FRACTURE_REGION__
#define __FRACTURE_REGION__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/PARTICLE_PARTITION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_POLICY.h>

namespace PhysBAM{

template<class T>
class FRACTURE_REGION
{
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,3> TV_INT;
    typedef typename GRID<TV>::NODE_ITERATOR NODE_ITERATOR;
    typedef typename RIGID_BODY_POLICY<TV>::WORLD_SPACE_INERTIA_TENSOR T_WORLD_SPACE_INERTIA_TENSOR;
public:
    typedef int HAS_TYPED_READ_WRITE;

    TRIANGULATED_SURFACE<T>* triangulated_surface;
    LEVELSET_IMPLICIT_OBJECT<TV>* implicit_object;
    MATRIX<T,3> levelset_RS;TV levelset_T;
    MATRIX<T,3> object_RS;TV object_T;
    TV_INT fracture_offset;
    T particle_intersection_thickness;
    bool need_destroy_data;
    FRAME<TV> extra_levelset_frame;
    PARTICLE_PARTITION<TV>* particle_partition;
    
    FRACTURE_REGION(TRIANGULATED_SURFACE<T>* triangulated_surface_input,LEVELSET_IMPLICIT_OBJECT<TV>* implicit_object_input,const bool create_particle_partition);
    virtual ~FRACTURE_REGION();

//#####################################################################
    ARRAY<FRACTURE_REGION<T>*> Intersect_With_Rigid_Body(const FRACTURE_REGION<T>& body,const bool use_particle_optimization,const bool tessellate_region=false);
    T Compute_Volume() const;
    void Compute_Inertial_Properties(const T density,TV& com,T& mass,T_WORLD_SPACE_INERTIA_TENSOR& inertia) const; 
    void Read(TYPED_ISTREAM& input);
    void Write(TYPED_OSTREAM& output) const;
    void Initialize_Particle_Partition();
//#####################################################################
};
}
#endif
