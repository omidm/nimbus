//#####################################################################
// Copyright 2002, Ronald Fedkiw
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Lagrange_Equations/ARTIFICIAL_VISCOSITY_VNR_1D.h>
using namespace PhysBAM;
//#####################################################################
// Function Artificial_Viscosity
//#####################################################################
// Q is (1,m-1) at the cell centers with the pressure
template<class T> void ARTIFICIAL_VISCOSITY_VNR_1D<T>::
Get_Artificial_Viscosity(EOS<T>& eos,GRID_LAGRANGE_1D<T>& grid,const ARRAY<T,VECTOR<int,1> >& mass,const ARRAY<T,VECTOR<int,1> >& velocity,const ARRAY<T,VECTOR<int,1> >& energy,ARRAY<T,VECTOR<int,1> >& Q)
{
    int i;
    int m=grid.m;   
    
    ARRAY<T,VECTOR<int,1> > length(1,m-1);grid.Get_Lengths(length);  
    ARRAY<T,VECTOR<int,1> > u_jump(1,m-1);for(i=1;i<=m-1;i++) u_jump(i)=velocity(i+1)-velocity(i);
    
    for(i=1;i<=m-1;i++){
        if(u_jump(i) >= 0) Q(i)=0;
        else{
            Q(i)=constant*mass(i)/length(i)*sqr(u_jump(i));
            if(limiter){
                T ux_center=u_jump(i)/length(i);
                T r_left=1,r_right=1;
                if(ux_center != 0){
                    if(i != 1) r_left=(u_jump(i-1)/length(i-1))/ux_center;
                    if(i != m-1) r_right=(u_jump(i+1)/length(i+1))/ux_center;}
                T psi=max((T)0,min((r_left+r_right)/2,2*r_left,2*r_right,(T)1));
                Q(i)=(1-psi)*Q(i);}}}
}
//#####################################################################
template class ARTIFICIAL_VISCOSITY_VNR_1D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class ARTIFICIAL_VISCOSITY_VNR_1D<double>;
#endif
