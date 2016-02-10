//#####################################################################
// Copyright 2002, Ronald Fedkiw.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HAMILTONIAN_3D
//#####################################################################
#ifndef __HAMILTONIAN_3D__
#define __HAMILTONIAN_3D__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class TV> class GRID;

template<class T>
class HAMILTONIAN_3D
{
    typedef VECTOR<T,3> TV;
protected:
    GRID<TV>& grid;

    HAMILTONIAN_3D(GRID<TV>& grid_input)
        :grid(grid_input)
    {}

public:
//#####################################################################
    virtual ~HAMILTONIAN_3D(){}
    virtual T H(const T phi_x,const T phi_y,const T phi_z,const int i=0,const int j=0,const int k=0,const T t=0){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T Maxabs_H1(const T phi_x_1,const T phi_x_2,const T phi_y_1,const T phi_y_2,const T phi_z_1,const T phi_z_2,const int i=0,const int j=0,const int k=0,const T t=0){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T Maxabs_H2(const T phi_x_1,const T phi_x_2,const T phi_y_1,const T phi_y_2,const T phi_z_1,const T phi_z_2,const int i=0,const int j=0,const int k=0,const T t=0){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T Maxabs_H3(const T phi_x_1,const T phi_x_2,const T phi_y_1,const T phi_y_2,const T phi_z_1,const T phi_z_2,const int i=0,const int j=0,const int k=0,const T t=0){PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};
}
#endif

