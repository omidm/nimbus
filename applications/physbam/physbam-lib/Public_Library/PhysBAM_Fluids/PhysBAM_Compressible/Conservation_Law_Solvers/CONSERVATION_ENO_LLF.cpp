//#####################################################################
// Copyright 2002-2009, Doug Enright, Ronald Fedkiw, Jon Gretarsson, Nipun Kwatra, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONSERVATION_ENO_LLF  
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_SEPARABLE_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/CONSERVATION_ENO_LLF.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/EIGENSYSTEM.h>
using namespace PhysBAM;
//#####################################################################
// Function Conservation_Solver
//#####################################################################
// psi is size (1,m) - U is size 3 by (-2,m+3) with 3 ghost cells - Fx is size 3 by (1,m)
template<class T_GRID,int d> void CONSERVATION_ENO_LLF<T_GRID,d>::
Conservation_Solver(const int m,const T dx,const ARRAY<bool,VECTOR<int,1> >& psi,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,
    ARRAY<TV_DIMENSION,VECTOR<int,1> >& Fx,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem_explicit,
    const VECTOR<bool,2>& outflow_boundaries,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_flux)
{
    switch(order){
        case 1:Conservation_Solver_Helper<1>(m,dx,psi,U,Fx,eigensystem,eigensystem_explicit,outflow_boundaries,U_flux);break;
        case 2:Conservation_Solver_Helper<2>(m,dx,psi,U,Fx,eigensystem,eigensystem_explicit,outflow_boundaries,U_flux);break;
        case 3:Conservation_Solver_Helper<3>(m,dx,psi,U,Fx,eigensystem,eigensystem_explicit,outflow_boundaries,U_flux);break;
        default:PHYSBAM_FATAL_ERROR();}
}
//#####################################################################
// Function Conservation_Solver_Helper
//#####################################################################
// Does not implement MENO for fully explicit case. Also does not implement the outflow_boundaries fix. Use Experimental version instead for them.
// TODO(kwatra): can try implementing Flux_Divided_By_Velocity for fully explicit eigensystems and see if multiplying by face_velocity later works.
template<class T_GRID,int d> template<int eno_order> void CONSERVATION_ENO_LLF<T_GRID,d>::
Conservation_Solver_Helper(const int m,const T dx,const ARRAY<bool,VECTOR<int,1> >& psi,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,
    ARRAY<TV_DIMENSION,VECTOR<int,1> >& Fx,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem_explicit,
    const VECTOR<bool,2>& outflow_boundaries,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_flux)
{
    int number_ghost_cells=3;
    int start_index_ghost=U.domain.min_corner.x,end_index_ghost=U.domain.max_corner.x;
    // divided differences    
    ARRAY<TV_DIMENSION,VECTOR<int,2> > DU(start_index_ghost,end_index_ghost,1,eno_order),DF(start_index_ghost,end_index_ghost,1,eno_order);
    ARRAY<TV_DIMENSION,VECTOR<int,1> > F(start_index_ghost,end_index_ghost);
    if(use_face_velocity_for_fluxes) eigensystem.Flux_Divided_By_Velocity(m,U,F,U_flux);
    else eigensystem.Flux(m,U,F,U_flux);
    for(int i=start_index_ghost;i<=end_index_ghost;i++){DU(i,1)=U(i);DF(i,1)=F(i);}
    for(int j=2;j<=order;j++) for(int i=start_index_ghost;i<=end_index_ghost+1-j;i++){
        DU(i,j)=(DU(i+1,j-1)-DU(i,j-1))/(j*dx);DF(i,j)=(DF(i+1,j-1)-DF(i,j-1))/(j*dx);}

    int start_index=U.domain.min_corner.x+number_ghost_cells,end_index=U.domain.max_corner.x-number_ghost_cells;
    ARRAY<bool,VECTOR<int,1> > psi_ghost(start_index-1,end_index+1);ARRAY<bool,VECTOR<int,1> >::Put(psi,psi_ghost);
    ARRAY<T,VECTOR<int,1> > lambda(1,d),lambda_left(1,d),lambda_right(1,d);

    // get globally max alpha
    TV_DIMENSION max_alpha;
    if(use_global_llf) for(int i=start_index;i<=end_index;i++) if(psi_ghost(i) || psi_ghost(i+1)){
        if(use_explicit_eigensystem_for_alphas) eigensystem_explicit.Eigenvalues(U,i,lambda,lambda_left,lambda_right);
        else eigensystem.Eigenvalues(U,i,lambda,lambda_left,lambda_right);
        for(int k=1;k<=d;k++){
            T alpha=Alpha(lambda_left,lambda_right,k,d);
            max_alpha[k]=max(max_alpha[k],alpha);}}

    ARRAY<TV_DIMENSION,VECTOR<int,1> > flux(U.domain.min_corner.x+2,U.domain.max_corner.x-3); // fluxes to the right of each point
    MATRIX<T,d,d> L,R;
    ARRAY<TV_DIMENSION,VECTOR<int,2> > LDU(start_index_ghost,end_index_ghost,1,eno_order),LDF(start_index_ghost,end_index_ghost,1,eno_order);
    const bool all_eigenvalues_same=eigensystem.All_Eigenvalues_Same();
    for(int i=start_index-1;i<=end_index;i++) if(psi_ghost(i) || psi_ghost(i+1)){ // compute flux
        // transfer the divided differences into the characteristic fields
        T velocity_multiplier=(T)1;
        if(use_face_velocity_for_fluxes) velocity_multiplier=eigensystem.Get_Face_Velocity_Component(i,use_standard_average,U);
        if(!all_eigenvalues_same){
            eigensystem.Eigenvectors(U,i,L,R);
            for(int j=1;j<=order;j++) for(int ii=i+1-j;ii<=i+1;ii++) if(ii >= start_index_ghost && ii <= end_index_ghost+1-j){
                LDU(ii,j)=L*DU(ii,j);LDF(ii,j)=L*velocity_multiplier*DF(ii,j);}}
        else{
            for(int j=1;j<=order;j++) for(int ii=i+1-j;ii<=i+1;ii++) if(ii >= start_index_ghost && ii <= end_index_ghost+1-j){
                LDU(ii,j)=DU(ii,j);LDF(ii,j)=velocity_multiplier*DF(ii,j);}}

        // compute max_alpha for use_global_llf=false case
        if(!use_global_llf){
            eigensystem.Eigenvalues(U,i,lambda,lambda_left,lambda_right);
            for(int k=1;k<=d;k++) max_alpha[k]=Alpha(lambda_left,lambda_right,k,d);}

        // find a flux in each characteristic field
        TV_DIMENSION Lflux; // flux without multiplication by R^T
        if(eno_order == 1) for(int k=1;k<=d;k++){
            T alpha=max_alpha[k];
            T flux_left=LDF(i,1)(k)+alpha*LDU(i,1)(k);
            T flux_right=LDF(i+1,1)(k)-alpha*LDU(i+1,1)(k);
            Lflux(k)=(T).5*(flux_left+flux_right);}
        else if(eno_order == 2) for(int k=1;k<=d;k++){
            T alpha=max_alpha[k];
            T flux_left=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,LDF(i,1)(k)+alpha*LDU(i,1)(k),LDF(i-1,2)(k)+alpha*LDU(i-1,2)(k),
                LDF(i,2)(k)+alpha*LDU(i,2)(k));
            T flux_right=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,LDF(i+1,1)(k)-alpha*LDU(i+1,1)(k),-(LDF(i+1,2)(k)-alpha*LDU(i+1,2)(k)),
                -(LDF(i,2)(k)-alpha*LDU(i,2)(k)));
            Lflux(k)=(T).5*(flux_left+flux_right);}
        else if(eno_order == 3) for(int k=1;k<=d;k++){
            T alpha=max_alpha[k];
            T flux_left=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,LDF(i,1)(k)+alpha*LDU(i,1)(k),LDF(i-1,2)(k)+alpha*LDU(i-1,2)(k),
                LDF(i,2)(k)+alpha*LDU(i,2)(k),LDF(i-2,3)(k)+alpha*LDU(i-2,3)(k),LDF(i-1,3)(k)+alpha*LDU(i-1,3)(k),LDF(i,3)(k)+alpha*LDU(i,3)(k));
            T flux_right=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,LDF(i+1,1)(k)-alpha*LDU(i+1,1)(k),-(LDF(i+1,2)(k)-alpha*LDU(i+1,2)(k)),
                -(LDF(i,2)(k)-alpha*LDU(i,2)(k)),LDF(i+1,3)(k)-alpha*LDU(i+1,3)(k),LDF(i,3)(k)-alpha*LDU(i,3)(k),LDF(i-1,3)(k)-alpha*LDU(i-1,3)(k));
            Lflux(k)=(T).5*(flux_left+flux_right);}
        else for(int k=1;k<=d;++k) Lflux(k)=(T)0;

        if(all_eigenvalues_same) flux(i)=Lflux;
        else flux(i)=R.Transpose_Times(Lflux);}

    // difference the fluxes
    T one_over_dx=1/dx;
    for(int i=U.domain.min_corner.x+3;i<=U.domain.max_corner.x-3;i++) if(psi_ghost(i)) Fx(i)=(flux(i)-flux(i-1))*one_over_dx;
    if(save_fluxes) for(int i=U.domain.min_corner.x+2;i<=U.domain.max_corner.x-3;i++) if(psi_ghost(i) || psi_ghost(i+1)) flux_temp(i)=flux(i);
}
//#####################################################################
// Function Conservation_Solver_Helper
//#####################################################################
template<class T_GRID,int d> template<int eno_order> void CONSERVATION_ENO_LLF<T_GRID,d>::
Conservation_Solver_Helper_Experimental(const int m,const T dx,const ARRAY<bool,VECTOR<int,1> >& psi,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U,ARRAY<TV_DIMENSION,VECTOR<int,1> >& Fx,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem,EIGENSYSTEM<T,TV_DIMENSION>& eigensystem_explicit,
    const VECTOR<bool,2>& outflow_boundaries,ARRAY<TV_DIMENSION,VECTOR<int,1> >* U_flux)
{
    ARRAY<TV_DIMENSION,VECTOR<int,1> > F(U.domain.min_corner.x,U.domain.max_corner.x);
    if(!use_face_velocity_for_fluxes){
        if(U_flux) eigensystem.Flux(m,U,F,U_flux);
        else eigensystem.Flux(m,U,F);}
    // TODO: Currently using U_extrapolated for Jacobian always. We might want to move this decision to the boundary, as in some boundaries like MPI we wont want to use extrapolated values but real
    // values.
    ARRAY<TV_DIMENSION,VECTOR<int,1> > U_extrapolated(U);
    if(outflow_boundaries(1)) U_extrapolated(0)=U(1);if(outflow_boundaries(2)) U_extrapolated(m+1)=U(m);

    ARRAY<bool,VECTOR<int,1> > psi_ghost(U.domain.min_corner.x+2,U.domain.max_corner.x-2);ARRAY<bool,VECTOR<int,1> >::Put(psi,psi_ghost); // ghost points for the if statement below  
    ARRAY<TV_DIMENSION,VECTOR<int,1> > flux(U.domain.min_corner.x+2,U.domain.max_corner.x-3); // fluxes to the right 0 of each point
    ARRAY<T,VECTOR<int,1> > lambda(1,d),lambda_left(1,d),lambda_right(1,d);
    MATRIX<T,d,d> L,R;
    ARRAY<VECTOR<T,eno_order> ,VECTOR<int,2> > DLU(1,d,U.domain.min_corner.x,U.domain.max_corner.x),DLF(1,d,U.domain.min_corner.x,U.domain.max_corner.x);
    ARRAY<VECTOR<T,eno_order> ,VECTOR<int,2> > DU(1,d,U.domain.min_corner.x,U.domain.max_corner.x),DF(1,d,U.domain.min_corner.x,U.domain.max_corner.x);
    ARRAY<VECTOR<T,eno_order> ,VECTOR<int,2> > *Dstate_ptr,*Dflux_ptr;

    TV_DIMENSION max_alpha;
    if(use_global_llf) for(int i=U.domain.min_corner.x+2;i<=U.domain.max_corner.x-3;i++) if(psi_ghost(i) || psi_ghost(i+1)){
        if(use_explicit_eigensystem_for_alphas) eigensystem_explicit.Eigenvalues(U_extrapolated,i,lambda,lambda_left,lambda_right);
        else eigensystem.Eigenvalues(U_extrapolated,i,lambda,lambda_left,lambda_right);
        for(int k=1;k<=d;k++){
            T alpha=Alpha(lambda_left,lambda_right,k,d);
            max_alpha[k]=max(max_alpha[k],alpha);}}
    for(int i=U.domain.min_corner.x+2;i<=U.domain.max_corner.x-3;i++) if(psi_ghost(i) || psi_ghost(i+1)){ // compute flux
        eigensystem.Eigenvalues(U_extrapolated,i,lambda,lambda_left,lambda_right);
        if(use_face_velocity_for_fluxes){
            if(U_flux) eigensystem.Flux_Using_Face_Velocity(VECTOR<int,2>(i+1-eno_order,i+eno_order),i,U,F,use_standard_average,U_flux);
            else eigensystem.Flux_Using_Face_Velocity(VECTOR<int,2>(i+1-eno_order,i+eno_order),i,U,F,use_standard_average);}
        if(!eigensystem.All_Eigenvalues_Same()){
            eigensystem.Eigenvectors(U_extrapolated,i,L,R);
            // transfer the fluxes into the characteristic fields
            for(int ii=i+1-eno_order;ii<=i+eno_order;ii++){
                for(int k=1;k<=d;k++){
                    DLU(k,ii)(1)=0;DLF(k,ii)(1)=0;
                    // TODO need to do something with psi-ghost here instead of ii<1 and ii>m
                    if(ii<1 && outflow_boundaries(1) && lambda(k)<0) for(int kk=1;kk<=d;kk++){DLU(k,ii)(1)+=L(k,kk)*U(1)(kk);DLF(k,ii)(1)+=L(k,kk)*F(1)(kk);}
                    else if(ii>m && outflow_boundaries(2) && lambda(k)>0) for(int kk=1;kk<=d;kk++){DLU(k,ii)(1)+=L(k,kk)*U(m)(kk);DLF(k,ii)(1)+=L(k,kk)*F(m)(kk);}
                    else for(int kk=1;kk<=d;kk++){DLU(k,ii)(1)+=L(k,kk)*U(ii)(kk);DLF(k,ii)(1)+=L(k,kk)*F(ii)(kk);}}}
            // compute the divided differences
            for(int j=2;j<=eno_order;j++) for(int k=1;k<=d;k++) for(int ii=i+1-eno_order;ii<=i+eno_order-j+1;ii++){
                DLU(k,ii)(j)=(DLU(k,ii+1)(j-1)-DLU(k,ii)(j-1))/(j*dx);DLF(k,ii)(j)=(DLF(k,ii+1)(j-1)-DLF(k,ii)(j-1))/(j*dx);}
            Dstate_ptr=&DLU;Dflux_ptr=&DLF;}
        else{
            for(int ii=i+1-eno_order;ii<=i+eno_order;ii++){
                for(int k=1;k<=d;k++){
                    DU(k,ii)(1)=0;DF(k,ii)(1)=0;
                    if(ii<1 && outflow_boundaries(1) && lambda(k)<0){DU(k,ii)(1)=U(1)(k);DF(k,ii)(1)=F(1)(k);}
                    else if(ii>m && outflow_boundaries(2) && lambda(k)>0){DU(k,ii)(1)=U(m)(k);DF(k,ii)(1)=F(m)(k);}
                    else{DU(k,ii)(1)=U(ii)(k);DF(k,ii)(1)=F(ii)(k);}}}
            // compute the divided differences
            for(int j=2;j<=eno_order;j++) for(int k=1;k<=d;k++) for(int ii=i+1-eno_order;ii<=i+eno_order-j+1;ii++){
                DU(k,ii)(j)=(DU(k,ii+1)(j-1)-DU(k,ii)(j-1))/(j*dx);DF(k,ii)(j)=(DF(k,ii+1)(j-1)-DF(k,ii)(j-1))/(j*dx);}
            Dstate_ptr=&DU;Dflux_ptr=&DF;}
        ARRAY<VECTOR<T,eno_order> ,VECTOR<int,2> > &Dstate=*Dstate_ptr,&Dflux=*Dflux_ptr;
        // find a flux in each characteristic field
        if(eno_order == 1) for(int k=1;k<=d;k++){
            T alpha;
            if(use_global_llf) alpha=max_alpha[k];
            else alpha=Alpha(lambda_left,lambda_right,k,d);
            T flux_left=Dflux(k,i)(1)+alpha*Dstate(k,i)(1);
            T flux_right=Dflux(k,i+1)(1)-alpha*Dstate(k,i+1)(1);
            if(!eigensystem.All_Eigenvalues_Same()) for(int kk=1;kk<=d;kk++) flux(i)(kk)+=(T).5*(flux_left+flux_right)*R(k,kk);
            else flux(i)(k)=(T).5*(flux_left+flux_right);}
        else if(eno_order == 2) for(int k=1;k<=d;k++){
            T alpha;
            if(use_global_llf) alpha=max_alpha[k];
            else alpha=Alpha(lambda_left,lambda_right,k,d);
            T flux_left=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,Dflux(k,i)(1)+alpha*Dstate(k,i)(1),Dflux(k,i-1)(2)+alpha*Dstate(k,i-1)(2),Dflux(k,i)(2)+alpha*Dstate(k,i)(2));
            T flux_right=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,Dflux(k,i+1)(1)-alpha*Dstate(k,i+1)(1),-(Dflux(k,i+1)(2)-alpha*Dstate(k,i+1)(2)),-(Dflux(k,i)(2)-alpha*Dstate(k,i)(2)));
            if(!eigensystem.All_Eigenvalues_Same()) for(int kk=1;kk<=d;kk++) flux(i)(kk)+=(T).5*(flux_left+flux_right)*R(k,kk);
            else flux(i)(k)=(T).5*(flux_left+flux_right);}
        else if(eno_order == 3) for(int k=1;k<=d;k++){
            T alpha;
            if(use_global_llf) alpha=max_alpha[k];
            else alpha=Alpha(lambda_left,lambda_right,k,d);
            T flux_left=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,Dflux(k,i)(1)+alpha*Dstate(k,i)(1),Dflux(k,i-1)(2)+alpha*Dstate(k,i-1)(2),Dflux(k,i)(2)+alpha*Dstate(k,i)(2),
                    Dflux(k,i-2)(3)+alpha*Dstate(k,i-2)(3),Dflux(k,i-1)(3)+alpha*Dstate(k,i-1)(3),Dflux(k,i)(3)+alpha*Dstate(k,i)(3));
            T flux_right=ADVECTION_SEPARABLE_UNIFORM<GRID<TV>,T>::ENO(dx,Dflux(k,i+1)(1)-alpha*Dstate(k,i+1)(1),-(Dflux(k,i+1)(2)-alpha*Dstate(k,i+1)(2)),-(Dflux(k,i)(2)-alpha*Dstate(k,i)(2)),
                    Dflux(k,i+1)(3)-alpha*Dstate(k,i+1)(3),Dflux(k,i)(3)-alpha*Dstate(k,i)(3),Dflux(k,i-1)(3)-alpha*Dstate(k,i-1)(3));
            if(!eigensystem.All_Eigenvalues_Same()) for(int kk=1;kk<=d;kk++) flux(i)(kk)+=(T).5*(flux_left+flux_right)*R(k,kk);
            else flux(i)(k)=(T).5*(flux_left+flux_right);}}
    
    // difference the fluxes
    T one_over_dx=1/dx;
    for(int i=U.domain.min_corner.x+3;i<=U.domain.max_corner.x-3;i++) if(psi_ghost(i)) Fx(i)=(flux(i)-flux(i-1))*one_over_dx;
    if(save_fluxes) for(int i=U.domain.min_corner.x+2;i<=U.domain.max_corner.x-3;i++) if(psi_ghost(i) || psi_ghost(i+1)) flux_temp(i)=flux(i);
}
//#####################################################################
#define INSTANTIATION_HELPER(T_GRID) \
    template class CONSERVATION_ENO_LLF<T_GRID,1>; \
    template class CONSERVATION_ENO_LLF<T_GRID,2>; \
    template class CONSERVATION_ENO_LLF<T_GRID,3>; \
    template class CONSERVATION_ENO_LLF<T_GRID,4>; \
    template class CONSERVATION_ENO_LLF<T_GRID,5>; \
    template class CONSERVATION_ENO_LLF<T_GRID,6>;
#define P(...) __VA_ARGS__
INSTANTIATION_HELPER(P(GRID<VECTOR<float,1> >))
INSTANTIATION_HELPER(P(GRID<VECTOR<float,2> >))
INSTANTIATION_HELPER(P(GRID<VECTOR<float,3> >))
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(P(GRID<VECTOR<double,1> >))
INSTANTIATION_HELPER(P(GRID<VECTOR<double,2> >))
INSTANTIATION_HELPER(P(GRID<VECTOR<double,3> >))
#endif
