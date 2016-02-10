//#####################################################################
// Copyright 2002-2004, Doug Enright, Ronald Fedkiw, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSERVATION_ENO_LLF_AND_CENTRAL  
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEPARABLE_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/Minmod.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF_AND_CENTRAL.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
using namespace PhysBAM;
//#####################################################################
// Function Conservation_Solver
//#####################################################################
// psi is size (1,m) - U is size 3 by (-2,m+3) with 3 ghost cells - Fx is size 3 by (1,m)
template<class T_GRID,int d> void CONSERVATION_ENO_LLF_AND_CENTRAL<T_GRID,d>::
Conservation_Solver(const int m,const T dx,const ARRAY<bool,VECTOR<int,1> >& psi,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& Fx,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem_explicit,
    const VECTOR<bool,2>& outflow_boundaries,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_flux)
{
    switch(order){
        case 1:Conservation_Solver_Helper<1>(m,dx,psi,U,Fx,eigensystem,eigensystem_explicit,outflow_boundaries);break;
        case 2:Conservation_Solver_Helper<2>(m,dx,psi,U,Fx,eigensystem,eigensystem_explicit,outflow_boundaries);break;
        case 3:Conservation_Solver_Helper<3>(m,dx,psi,U,Fx,eigensystem,eigensystem_explicit,outflow_boundaries);break;
        default: PHYSBAM_FATAL_ERROR();}
}
//#####################################################################
// Function Conservation_Solver_Helper
//#####################################################################
template<class T_GRID,int d> template<int eno_order> void CONSERVATION_ENO_LLF_AND_CENTRAL<T_GRID,d>::
Conservation_Solver_Helper(const int m,const T dx,const ARRAY<bool,VECTOR<int,1> >& psi,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& Fx,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem,
    EIGENSYSTEM<T,TV_DIMENSION>& eigensystem_explicit,const VECTOR<bool,2>& outflow_boundaries)
{
    int k,i,j;

    // divided differences    
    ARRAY<VECTOR<T,eno_order> ,VECTOR<int,2> > DU(1,d,-2,m+3),DF(1,d,-2,m+3);
    ARRAY<TV_DIMENSION,VECTOR<int,1> > F(-2,m+3);eigensystem.Flux(m,U,F); 
    for(i=-2;i<=m+3;i++) for(k=1;k<=d;k++){DU(k,i)(1)=U(i)(k);DF(k,i)(1)=F(i)(k);}
    for(j=2;j<=eno_order;j++) for(k=1;k<=d;k++) for(i=-2;i<=m+4-j;i++){DU(k,i)(j)=(DU(k,i+1)(j-1)-DU(k,i)(j-1))/(j*dx);DF(k,i)(j)=(DF(k,i+1)(j-1)-DF(k,i)(j-1))/(j*dx);}
                
    // calculate the fluxes 
    ARRAY<bool,VECTOR<int,1> > psi_ghost(0,m+1);ARRAY<bool,VECTOR<int,1> >::Put(psi,psi_ghost); // ghost points for the if statement below  
    ARRAY<TV_DIMENSION,VECTOR<int,1> > flux(0,m); // fluxes to the right of each point
    ARRAY<T,VECTOR<int,1> > lambda(1,d),lambda_left(1,d),lambda_right(1,d);
    MATRIX<T,d,d> L,R;
    ARRAY<VECTOR<T,eno_order> ,VECTOR<int,2> > LDU(1,d,-2,m+3),LDF(1,d,-2,m+3);
    for(i=0;i<=m;i++) if(psi_ghost(i) || psi_ghost(i+1)){ // compute flux
        // eigensystem
        if(eigensystem.Eigenvalues(U,i,lambda,lambda_left,lambda_right)){
            eigensystem.Eigenvectors(U,i,L,R);
            // transfer the divided differences into the characteristic fields
            for(j=1;j<=eno_order;j++) for(int ii=i+1-j;ii<=i+1;ii++) if(ii >= -2 && ii <= m+4-j) for(k=1;k<=d;k++){
                LDU(k,ii)(j)=LDF(k,ii)(j)=0;
                for(int kk=1;kk<=d;kk++){LDU(k,ii)(j)+=L(k,kk)*DU(kk,ii)(j);LDF(k,ii)(j)+=L(k,kk)*DF(kk,ii)(j);}}
            // find a flux in each characteristic field
            if(eno_order == 1) for(k=1;k<=d;k++){
                T alpha=Alpha(lambda_left,lambda_right,k,d);
                T flux_left=LDF(k,i)(1)+alpha*LDU(k,i)(1);
                T flux_right=LDF(k,i+1)(1)-alpha*LDU(k,i+1)(1);
                for(int kk=1;kk<=d;kk++) flux(i)(kk)+=(T).5*(flux_left+flux_right)*R(k,kk);}
            else if(eno_order == 2) for(k=1;k<=d;k++){
                T alpha=Alpha(lambda_left,lambda_right,k,d);
                T flux_left=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,LDF(k,i)(1)+alpha*LDU(k,i)(1),LDF(k,i-1)(2)+alpha*LDU(k,i-1)(2),LDF(k,i)(2)+alpha*LDU(k,i)(2));
                T flux_right=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,LDF(k,i+1)(1)-alpha*LDU(k,i+1)(1),-(LDF(k,i+1)(2)-alpha*LDU(k,i+1)(2)),-(LDF(k,i)(2)-alpha*LDU(k,i)(2)));
                for(int kk=1;kk<=d;kk++) flux(i)(kk)+=(T).5*(flux_left+flux_right)*R(k,kk);}
            else if(eno_order == 3) for(k=1;k<=d;k++){
                T alpha=Alpha(lambda_left,lambda_right,k,d);
                T flux_left=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,LDF(k,i)(1)+alpha*LDU(k,i)(1),LDF(k,i-1)(2)+alpha*LDU(k,i-1)(2),LDF(k,i)(2)+alpha*LDU(k,i)(2),
                                                                       LDF(k,i-2)(3)+alpha*LDU(k,i-2)(3),LDF(k,i-1)(3)+alpha*LDU(k,i-1)(3),LDF(k,i)(3)+alpha*LDU(k,i)(3));
                T flux_right=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,LDF(k,i+1)(1)-alpha*LDU(k,i+1)(1),-(LDF(k,i+1)(2)-alpha*LDU(k,i+1)(2)),-(LDF(k,i)(2)-alpha*LDU(k,i)(2)),
                                                                         LDF(k,i+1)(3)-alpha*LDU(k,i+1)(3),LDF(k,i)(3)-alpha*LDU(k,i)(3),LDF(k,i-1)(3)-alpha*LDU(k,i-1)(3));
                for(int kk=1;kk<=d;kk++) flux(i)(kk)+=(T).5*(flux_left+flux_right)*R(k,kk);}}
        else{ // use the central scheme
            // change parameters
            int save_alpha=field_by_field_alpha;field_by_field_alpha=0; // use maximum alpha
            T save_amplification_factor=amplification_factor;amplification_factor=central_amplification_factor;
            // find a flux in each component
            if(central_order == 1) for(k=1;k<=d;k++){
                T alpha=Alpha(lambda_left,lambda_right,k,d);
                T flux_left=DF(k,i)(1)+alpha*DU(k,i)(1);
                T flux_right=DF(k,i+1)(1)-alpha*DU(k,i+1)(1);
                flux(i)(k)=(T).5*(flux_left+flux_right);}
            else if(central_order == 2) for(k=1;k<=d;k++){
                T alpha=Alpha(lambda_left,lambda_right,k,d);
                T flux_left=Minmod(dx,DF(k,i)(1)+alpha*DU(k,i)(1),DF(k,i-1)(2)+alpha*DU(k,i-1)(2),DF(k,i)(2)+alpha*DU(k,i)(2));
                T flux_right=Minmod(dx,DF(k,i+1)(1)-alpha*DU(k,i+1)(1),-(DF(k,i+1)(2)-alpha*DU(k,i+1)(2)),-(DF(k,i)(2)-alpha*DU(k,i)(2)));
                flux(i)(k)=(T).5*(flux_left+flux_right);}
            // change parameters back to the usaul ones
            field_by_field_alpha=save_alpha; 
            amplification_factor=save_amplification_factor;}}

    // difference the fluxes
    T one_over_dx=1/dx;
    for(i=1;i<=m;i++) if(psi_ghost(i)) Fx(i)=(flux(i)-flux(i-1))*one_over_dx;
}
//#####################################################################
#define INSTANTIATION_HELPER(T_GRID) \
    template class CONSERVATION_ENO_LLF_AND_CENTRAL<T_GRID,1>; \
    template class CONSERVATION_ENO_LLF_AND_CENTRAL<T_GRID,2>; \
    template class CONSERVATION_ENO_LLF_AND_CENTRAL<T_GRID,3>;
#define P(...) __VA_ARGS__
INSTANTIATION_HELPER(P(GRID<VECTOR<float,1> >))
INSTANTIATION_HELPER(P(GRID<VECTOR<float,2> >))
INSTANTIATION_HELPER(P(GRID<VECTOR<float,3> >))
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(P(GRID<VECTOR<double,1> >))
INSTANTIATION_HELPER(P(GRID<VECTOR<double,2> >))
INSTANTIATION_HELPER(P(GRID<VECTOR<double,3> >))
#endif
