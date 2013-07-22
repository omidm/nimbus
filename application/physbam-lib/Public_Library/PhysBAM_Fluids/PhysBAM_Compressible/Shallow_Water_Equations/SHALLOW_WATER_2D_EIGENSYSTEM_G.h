//#####################################################################
// Copyright 2003-2007, Ronald Fedkiw, Eran Guendelman, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_2D_EIGENSYSTEM_G  
//#####################################################################
#ifndef __SHALLOW_WATER_2D_EIGENSYSTEM_G__
#define __SHALLOW_WATER_2D_EIGENSYSTEM_G__   

#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
namespace PhysBAM{

template<class T>
class SHALLOW_WATER_2D_EIGENSYSTEM_G:public EIGENSYSTEM<T,VECTOR<T,3> >
{
    typedef VECTOR<T,3> TV_DIMENSION;
    enum WORKAROUND1 {d=TV_DIMENSION::m};
public:
    T gravity;

    SHALLOW_WATER_2D_EIGENSYSTEM_G(const T gravity_input)
        :gravity(gravity_input)
    {}
    
//#####################################################################
    void Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& G,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) PHYSBAM_OVERRIDE;        
    T Maximum_Magnitude_Eigenvalue(const VECTOR<T,3>& U_cell) PHYSBAM_OVERRIDE;
    bool Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right) PHYSBAM_OVERRIDE; 
    void Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R) PHYSBAM_OVERRIDE; 
//#####################################################################
};
}    
#endif
