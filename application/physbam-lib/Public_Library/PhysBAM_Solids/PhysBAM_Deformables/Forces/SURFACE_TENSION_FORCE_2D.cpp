//#####################################################################
// Copyright 2010-2011, Craig Shroeder, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SURFACE_TENSION_FORCE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> SURFACE_TENSION_FORCE<VECTOR<T,2> >::
SURFACE_TENSION_FORCE(SEGMENTED_CURVE_2D<T>& surface_input,T surface_tension_coefficient_input)
    :BASE(dynamic_cast<PARTICLES<VECTOR<T,2> >&>(surface_input.particles)),surface(surface_input),particles(surface_input.particles),surface_tension_coefficient(surface_tension_coefficient_input),dt(0),apply_explicit_forces(true),apply_implicit_forces(true),apply_tangential_implicit(true)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> SURFACE_TENSION_FORCE<VECTOR<T,2> >::
~SURFACE_TENSION_FORCE()
{
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,2> >::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    if(apply_explicit_forces)
        for(int i=1;i<=surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
            TV f=(particles.X(k.x)-particles.X(k.y))*coefficients(i);
            F(k.x)-=f;
            F(k.y)+=f;}
}
//#####################################################################
// Function Tangential_Helper
//#####################################################################
template<class T> static void Tangential_Helper(MATRIX<T,2,1>& tangential,const VECTOR<T,2>& normal)
{
    tangential.Column(1)=normal.Orthogonal_Vector();
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,2> >::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    coefficients.Resize(surface.mesh.elements.m);
    normal.Resize(surface.mesh.elements.m);
    for(int i=1;i<=surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
        normal(i)=(particles.X(k.x)-particles.X(k.y)).Orthogonal_Vector();
        coefficients(i)=surface_tension_coefficient/normal(i).Normalize();}
    if(apply_implicit_forces){
        sqrt_coefficients.Resize(surface.mesh.elements.m);
        for(int i=1;i<=surface.mesh.elements.m;i++)
            sqrt_coefficients(i)=dt*sqrt(coefficients(i));}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,2> >::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    if(apply_implicit_forces)
        for(int i=1;i<=surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
            TV f;
            if(apply_tangential_implicit) f=(V(k.x)-V(k.y))*coefficients(i)*dt;
            else f=(TV::Dot_Product((V(k.x)-V(k.y)),normal(i))*coefficients(i)*dt)*normal(i);
            F(k.x)-=f;
            F(k.y)+=f;}
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,2> >::
Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const
{
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,2> >::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,2> >::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class T> int SURFACE_TENSION_FORCE<VECTOR<T,2> >::
Velocity_Dependent_Forces_Size() const
{
    int size=0;
    if(apply_implicit_forces){
        if(!apply_tangential_implicit) size+=surface.mesh.elements.m;
        else size+=2*surface.mesh.elements.m;}
    return size;
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,2> >::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
    if(apply_implicit_forces){
        if(!apply_tangential_implicit) 
            for(int i=1;i<=surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
                aggregate(i)=TV::Dot_Product((V(k.x)-V(k.y)),normal(i))*sqrt_coefficients(i);}
        else for(int i=1;i<=surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
            TV value=(V(k.x)-V(k.y))*sqrt_coefficients(i);
            aggregate(i)=value(1);aggregate(surface.mesh.elements.m+i)=value(2);}}

}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,2> >::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const
{
    if(apply_implicit_forces){
        if(!apply_tangential_implicit) 
            for(int i=1;i<=surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
                TV f=(aggregate(i)*sqrt_coefficients(i))*normal(i);
                F(k.x)+=f;
                F(k.y)-=f;}
        else for(int i=1;i<=surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
            TV f;f(1)=(aggregate(i)*sqrt_coefficients(i));f(2)=(aggregate(surface.mesh.elements.m+i)*sqrt_coefficients(i));
            F(k.x)+=f;
            F(k.y)-=f;}}
}
//#####################################################################
// Function Add_Raw_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,2> >::
Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const
{
    if(apply_implicit_forces){
        if(!apply_tangential_implicit) 
            for(int i=1;i<=surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
                int off_kx=(k.x-1)*TV::m,off_ky=(k.y-1)*TV::m;
                TV sn=sqrt_coefficients(i)*normal(i);
                for(int j=1;j<=TV::m;j++){
                    data.Append(TRIPLE<int,int,T>(i,off_kx+j,sn(j)));
                    data.Append(TRIPLE<int,int,T>(i,off_ky+j,-sn(j)));}}
        else for(int i=1;i<=surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
                int off_kx=(k.x-1)*TV::m,off_ky=(k.y-1)*TV::m;
                T s=sqrt_coefficients(i);
                for(int j=1;j<=TV::m;j++){
                    data.Append(TRIPLE<int,int,T>(i,off_kx+j,s));
                    data.Append(TRIPLE<int,int,T>(i,off_ky+j,-s));}}}
}
//#####################################################################
// Function Construct_Implicit_Matrix_First_Half
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,2> >::
Construct_Implicit_Matrix_First_Half(SPARSE_MATRIX_FLAT_MXN<T>& matrix,int n)
{
    int pm=(n==0)?particles.X.m:n;
    matrix.Reset(pm*TV::m);
    if(1){//segment
        if(apply_tangential_implicit){
            for(int i=1;i<=surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
                T sn=sqrt_coefficients(i);
                for(int axis=1;axis<=TV::m;axis++){
                    matrix.Append_Entry_To_Current_Row(pm*(axis-1)+k.x,sn);
                    matrix.Append_Entry_To_Current_Row(pm*(axis-1)+k.y,-sn);
                    matrix.Finish_Row();}}}
        else{
            for(int i=1;i<=surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
                TV sn=sqrt_coefficients(i)*normal(i);
                for(int axis=1;axis<=TV::m;axis++){
                    matrix.Append_Entry_To_Current_Row(pm*(axis-1)+k.x,sn(axis));
                    matrix.Append_Entry_To_Current_Row(pm*(axis-1)+k.y,-sn(axis));
                    matrix.Finish_Row();}}}}
    else{//vertex
        surface.mesh.Initialize_Incident_Elements();
        for(int i=1;i<=surface.mesh.incident_elements->m;i++){
            ARRAY<int>& list=(*surface.mesh.incident_elements)(i);
            if(list.m==0) continue;
            for(int axis=1;axis<=TV::m;axis++){
                for(int j=1;j<=list.m;j++){
                    VECTOR<int,2> k=surface.mesh.elements(list(j));
                    T sn=sqrt_coefficients(list(j));
                    if(k.y==i) matrix.Append_Entry_To_Current_Row(pm*(axis-1)+k.x,sn);
                    else matrix.Append_Entry_To_Current_Row(pm*(axis-1)+k.y,sn);}
                matrix.Finish_Row();}}
        surface.mesh.Delete_Auxiliary_Structures();}
    matrix.Sort_Entries();
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,2> >::
Initialize_CFL(ARRAY_VIEW<typename DEFORMABLES_FORCES<VECTOR<T,2> >::FREQUENCY_DATA> frequency)
{
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class T> T SURFACE_TENSION_FORCE<VECTOR<T,2> >::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,2> >::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    surface.mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,2> >::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class T> T SURFACE_TENSION_FORCE<VECTOR<T,2> >::
Potential_Energy(const T time) const
{
    T pe=0;
    if(apply_explicit_forces)
        for(int i=1;i<=surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
            pe+=surface_tension_coefficient*(particles.X(k.x)-particles.X(k.y)).Magnitude();}
    return pe;
}
//#####################################################################
// Function Compute_Curvatures
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,2> >::
Compute_Curvatures(ARRAY_VIEW<T> curvature) const
{
    ARRAY<VECTOR<T,2> > F(surface.mesh.elements.m);
    Add_Velocity_Independent_Forces(F,0);
    ARRAY<T> A(coefficients.m);
    for(int i=1;i<=surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
        T len=(particles.X(k.x)-particles.X(k.y)).Magnitude()/2;
        A(k.x)+=len;
        A(k.y)+=len;}
    for(int i=1;i<=F.m;i++) curvature(i)=abs(F(i).Magnitude()/A(i)/surface_tension_coefficient);
}
//#####################################################################
// Function Dump_Curvatures
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,2> >::
Dump_Curvatures() const
{
    T mn=FLT_MAX,mx=-mn,av=0;
    int n=0;
    ARRAY<VECTOR<T,2> > F(surface.mesh.elements.m);
    Add_Velocity_Independent_Forces(F,0);
    ARRAY<T> A(coefficients.m);
    for(int i=1;i<=surface.mesh.elements.m;i++){VECTOR<int,2> k=surface.mesh.elements(i);
        T len=(particles.X(k.x)-particles.X(k.y)).Magnitude()/2;
        A(k.x)+=len;
        A(k.y)+=len;}
    for(int i=1;i<=F.m;i++){
        T K=abs(F(i).Magnitude()/A(i)/surface_tension_coefficient);
        TV dx=particles.X(i)-(T).02;
        av+=K;
        n++;
        if(K>mx) mx=K;
        if(K<mn) mn=K;}
    if(n) {std::stringstream ss;ss<<"cstats "<<mn<<"   "<<mx<<std::endl;LOG::filecout(ss.str());}
    if(coefficients.m) {std::stringstream ss;ss<<"length estimates  "<<ARRAYS_COMPUTATIONS::Max(coefficients)/ARRAYS_COMPUTATIONS::Min(coefficients)<<std::endl;LOG::filecout(ss.str());}
}
template class SURFACE_TENSION_FORCE<VECTOR<float,2> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SURFACE_TENSION_FORCE<VECTOR<double,2> >;
#endif
