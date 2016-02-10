//#####################################################################
// Copyright 2007, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_FRACTURE_QUASISTATICS_FORCES
//#####################################################################
#ifndef __RIGID_FRACTURE_QUASISTATICS_FORCES__
#define __RIGID_FRACTURE_QUASISTATICS_FORCES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
namespace PhysBAM{

template<class T> class RIGID_BODY_FRACTURE_OBJECT_3D;

template<class T_input>
class RIGID_FRACTURE_QUASISTATICS_FORCES:public EXAMPLE_FORCES_AND_VELOCITIES<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
    typedef EXAMPLE_FORCES_AND_VELOCITIES<TV> BASE;
    using BASE::Add_External_Forces; // silence -Woverloaded-virtual
public:
    const RIGID_BODY_FRACTURE_OBJECT_3D<T>* fracture_object;
    bool initialized;
    ARRAY<TV>* impulses;

    VECTOR<int,3> null_space_nodes;
    TV null_space_position;

    TV direction1,direction2;
    MATRIX<T,3> dx3_transform,dx2_transform;
    ARRAY<VECTOR<double,3> > rotation_null_space_vector_1,rotation_null_space_vector_2,rotation_null_space_vector_3;

    RIGID_FRACTURE_QUASISTATICS_FORCES():initialized(false)
    {}

//#####################################################################
    void Initialize(const RIGID_BODY_FRACTURE_OBJECT_3D<T>& fracture_object_input,const VECTOR<int,3>& null_space_nodes_input);
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time); // zero out entries corresponding to external positions
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif

