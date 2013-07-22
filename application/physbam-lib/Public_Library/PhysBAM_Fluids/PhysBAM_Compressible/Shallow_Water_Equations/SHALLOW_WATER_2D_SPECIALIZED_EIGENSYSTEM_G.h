//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_G  
//##################################################################### 
//
// only the flux of G(U) is different from the F(U) term and only because of the slice_index variable - otherwise,  the eigensystem is the same
//
//#####################################################################
#ifndef __SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_G__
#define __SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_G__   

#include <PhysBAM_Fluids/PhysBAM_Compressible/Shallow_Water_Equations/SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F.h>
namespace PhysBAM{

template<class T>
class SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_G:public SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F<T>
{
public:
    using SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F<T>::slice_index;using SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F<T>::gravity;

    SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_G(const T gravity_input,ARRAY<T,VECTOR<int,2> > &eta_ghost_input,const T min_height_input)
        :SHALLOW_WATER_2D_SPECIALIZED_EIGENSYSTEM_F<T>(gravity_input,eta_ghost_input,min_height_input)
    {}
    
//#####################################################################
    void Flux(const int m,const ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& U,ARRAY<VECTOR<T,2> ,VECTOR<int,1> >& F,ARRAY<VECTOR<T,2> ,VECTOR<int,1> >* U_clamped=0) PHYSBAM_OVERRIDE;
//#####################################################################
};
}    
#endif
