//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_RIGIDS_PARTICLES
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_0X0.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_MATRIX_0X0.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_MATRIX_1X1.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Particles/RIGIDS_PARTICLES_FORWARD.h>
namespace PhysBAM{

void Initialize_Rigids_Particles()
{
    static bool done=false;if(done) return;done=true;
    Register_Attribute_Name(ATTRIBUTE_ID_RIGID_MASS,"rigid_mass");
    Register_Attribute_Name(ATTRIBUTE_ID_RIGID_INERTIA_TENSOR,"rigid_inertia_tensor");
    Register_Attribute_Name(ATTRIBUTE_ID_ANGULAR_MOMENTUM,"angular_momentum");
    Register_Attribute_Name(ATTRIBUTE_ID_KINEMATIC,"kinematic");

    #define READ_WRITE_SCALAR_HELPER(T,RW) \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<DIAGONAL_MATRIX<T,3> >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<MATRIX<T,1,1> >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<MATRIX<T,0,0> >(); \

    #define READ_WRITE_TYPE_HELPER(RW) \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<bool>();

    READ_WRITE_SCALAR_HELPER(float,float);
    READ_WRITE_TYPE_HELPER(float);
    #ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    READ_WRITE_SCALAR_HELPER(float,double);
    READ_WRITE_SCALAR_HELPER(double,float);
    READ_WRITE_SCALAR_HELPER(double,double);
    READ_WRITE_TYPE_HELPER(double);
    #endif

    Initialize_Geometry_Particle();
}
}
