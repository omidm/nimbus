//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Nipun Kwatra, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_3D_EIGENSYSTEM_H  
//##################################################################### 
#ifndef __EULER_3D_EIGENSYSTEM_H__
#define __EULER_3D_EIGENSYSTEM_H__   

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Euler_Equations/EULER_EIGENSYSTEM.h>
namespace PhysBAM{

template<class T_input>
class EULER_3D_EIGENSYSTEM_H:public EULER_EIGENSYSTEM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;typedef VECTOR<T,5> TV_DIMENSION;
    enum WORKAROUND1 {d=TV_DIMENSION::m};
    typedef EULER_EIGENSYSTEM<GRID<VECTOR<T_input,3> > > BASE;
public:
    using BASE::eos;

    bool only_pressure_flux;

    EULER_3D_EIGENSYSTEM_H(bool only_pressure_flux_input=false):only_pressure_flux(only_pressure_flux_input)
    {}

//#####################################################################
    void Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& H,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) PHYSBAM_OVERRIDE;   
    void Flux_Using_Face_Velocity(VECTOR<int,2> range,const int face_index,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& H,const bool use_standard_average,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) PHYSBAM_OVERRIDE;
    T Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell) PHYSBAM_OVERRIDE;
    bool Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int ij,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right) PHYSBAM_OVERRIDE;   
    void Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int ij,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R) PHYSBAM_OVERRIDE;   
//#####################################################################
};
}    
#endif
