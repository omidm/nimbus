//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REACTIVE_EULER_2D_EIGENSYSTEM_F  
//#################################################################### 
#ifndef __REACTIVE_EULER_2D_EIGENSYSTEM_F__
#define __REACTIVE_EULER_2D_EIGENSYSTEM_F__   

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Reactive_Euler_Equations/REACTIVE_EULER.h>
namespace PhysBAM{

template<class T_input>
class REACTIVE_EULER_2D_EIGENSYSTEM_F:public EIGENSYSTEM<T_input,VECTOR<T_input,5> >,public REACTIVE_EULER<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;typedef VECTOR<T,2> TV;typedef VECTOR<T,5> TV_DIMENSION;
    enum WORKAROUND1 {d=TV_DIMENSION::m};
public:
    using REACTIVE_EULER<GRID<TV> >::eos;

    REACTIVE_EULER_2D_EIGENSYSTEM_F(REACTIVE_EOS<T>& eos_input)
        :REACTIVE_EULER<GRID<TV> >(eos_input)
    {}
    
//#####################################################################
    void Flux(const int m,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& F,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_clamped=0) PHYSBAM_OVERRIDE;        
    T Maximum_Magnitude_Eigenvalue(const TV_DIMENSION& U_cell) PHYSBAM_OVERRIDE;
    bool Eigenvalues(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,ARRAY<T,VECTOR<int,1> >& lambda,ARRAY<T,VECTOR<int,1> >& lambda_left,ARRAY<T,VECTOR<int,1> >& lambda_right) PHYSBAM_OVERRIDE;   
    void Eigenvectors(const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,const int i,MATRIX<T,d,d>& L,MATRIX<T,d,d>& R) PHYSBAM_OVERRIDE;   
//#####################################################################
};   
}
#endif
