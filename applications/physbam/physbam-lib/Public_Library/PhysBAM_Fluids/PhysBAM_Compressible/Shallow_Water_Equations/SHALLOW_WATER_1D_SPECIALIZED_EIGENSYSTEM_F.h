//#####################################################################
// Copyright 2003-2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_1D_SPECIALIZED_EIGENSYSTEM_F  
//##################################################################### 
#ifndef __SHALLOW_WATER_1D_SPECIALIZED_EIGENSYSTEM_F__
#define __SHALLOW_WATER_1D_SPECIALIZED_EIGENSYSTEM_F__   

#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
namespace PhysBAM{

template<class T>
class SHALLOW_WATER_1D_SPECIALIZED_EIGENSYSTEM_F:public EIGENSYSTEM<T,VECTOR<T,2> >
{
    enum WORKAROUND1 {d=2};
public:
    T gravity;
    ARRAY<T,VECTOR<int,1> >& eta_ghost;
    T min_height;

    SHALLOW_WATER_1D_SPECIALIZED_EIGENSYSTEM_F(const T gravity_input,ARRAY<T,VECTOR<int,1> >& eta_ghost_input,const T min_height_input)
        :gravity(gravity_input),eta_ghost(eta_ghost_input),min_height(min_height_input)
    {}
   
//#####################################################################
    void Flux(const int m,const ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& U,ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& F,ARRAY<VECTOR<T,2> ,VECTOR<int,1> >* U_clamped=0) PHYSBAM_OVERRIDE;        
    bool Eigenvalues(const ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& U,const int i,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right) PHYSBAM_OVERRIDE; 
    void Eigenvectors(const ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R) PHYSBAM_OVERRIDE; 
//#####################################################################
};
}    
#endif
