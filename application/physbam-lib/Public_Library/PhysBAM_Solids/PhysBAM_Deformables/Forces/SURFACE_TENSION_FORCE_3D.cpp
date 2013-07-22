//#####################################################################
// Copyright 2011, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SURFACE_TENSION_FORCE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> SURFACE_TENSION_FORCE<VECTOR<T,3> >::
SURFACE_TENSION_FORCE(TRIANGULATED_SURFACE<T>& surface_input,T surface_tension_coefficient_input)
    :BASE(dynamic_cast<PARTICLES<VECTOR<T,3> >&>(surface_input.particles)),surface(surface_input),particles(surface_input.particles),surface_tension_coefficient(surface_tension_coefficient_input),dt(0),apply_explicit_forces(true),apply_implicit_forces(true),use_sph_laplace(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> SURFACE_TENSION_FORCE<VECTOR<T,3> >::
~SURFACE_TENSION_FORCE()
{
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,3> >::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    if(apply_explicit_forces){
        ARRAY_VIEW<TV> X(particles.X);
        if(use_sph_laplace) for(int i=1;i<=sph_coefficients.m;i++){
            const SPH_COEFFICIENT_ENTRY& entry=sph_coefficients(i);
            TV f=(X(entry.i)-X(entry.j))*entry.weight;
            F(entry.i)-=f;
            F(entry.j)+=f;}
        else for(int i=1;i<=fvm_coefficients.m;i++){
            const VECTOR<int,2>& e=segments(i);
            const T& c=fvm_coefficients(i);
            TV f=c*(X(e(1))-X(e(2)));
            F(e(1))+=f;
            F(e(2))-=f;}}
}
//#####################################################################
// Function Tangential_Helper
//#####################################################################
template<class T> static void Tangential_Helper(MATRIX<T,3,2>& tangential,const VECTOR<T,3>& normal)
{
    tangential.Column(1)=normal.Unit_Orthogonal_Vector();
    tangential.Column(2)=normal.Cross_Product(normal,tangential.Column(1));
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,3> >::
Update_Position_Based_State(const T time,const bool is_position_update)
{
    if(use_sph_laplace){
        vertex_areas.Resize(particles.X.m);
        ARRAYS_COMPUTATIONS::Fill(vertex_areas,(T)0);
        for(int i=1;i<=surface.mesh.elements.m;i++){VECTOR<int,3> k=surface.mesh.elements(i);
            T area=surface.Area(i)/3;
            for(int axis=1;axis<=3;axis++) vertex_areas(k(axis))+=area;}
        ARRAY<ARRAY<int>,TV_INT> index(grid.Domain_Indices(1));
        for(int i=1;i<=particles.X.m;i++){
            const TV& p=particles.X(i);
            TV_INT cell=grid.Cell(p,1);
            index(cell).Append(i);}
        T h=(grid.dX*2).Product();
        T constant=1/(4*(T)pi*h*h);
        int int_radius=2;
        sph_coefficients.Resize(0);
        for(int i=1;i<=particles.X.m;i++){
            const TV& p=particles.X(i);
            TV_INT cell=grid.Cell(p,1);
            RANGE<TV_INT> cell_range(clamp_min(cell-int_radius,grid.Domain_Indices(1).min_corner),clamp_max(cell+int_radius,grid.Domain_Indices(1).max_corner));
            for(UNIFORM_GRID_ITERATOR_CELL<TV> iter(grid,cell_range);iter.Valid();iter.Next()){
                TV_INT cell1=iter.Cell_Index();
                if(index(cell1).Size()==0) continue;
                for(int block_i=1;block_i<=index(cell1).Size();block_i++){int j=index(cell1)(block_i);
                    if(i>=j) continue;
                    T weight=surface_tension_coefficient*vertex_areas(i)*vertex_areas(j)*exp(-(p-particles.X(j)).Magnitude_Squared()/(4*h))*constant;
                    if(weight<=1e-6) continue;
                    SPH_COEFFICIENT_ENTRY entry;
                    entry.i=i;entry.j=j;entry.weight=weight;
                    sph_coefficients.Append(entry);}}}}
    else{
        surface.mesh.Initialize_Segment_Mesh();
        surface.mesh.Initialize_Edge_Triangles();
        segments=surface.mesh.segment_mesh->elements;
        fvm_coefficients.Resize(segments.m);
        for(int i=1;i<=segments.m;i++){
            VECTOR<int,2> e=segments(i);
            ARRAY<int>& list=(*surface.mesh.edge_triangles)(i);
            T coefficient=0;
            for(int j=1;j<=list.m;j++) coefficient+=Cotangent(list(j),e(1),e(2));
            fvm_coefficients(i)=min((T)0,-surface_tension_coefficient*coefficient*(T)0.5);}}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,3> >::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    if(apply_implicit_forces){};
}
//#####################################################################
// Function Construct_Implicit_Matrix_First_Half
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,3> >::
Construct_Implicit_Matrix_First_Half(SPARSE_MATRIX_FLAT_MXN<T>& matrix,int n)
{
    if(use_sph_laplace){
        int pm=(n==0)?particles.X.m:n;
        matrix.Reset(pm*3);
        for(int i=1;i<=sph_coefficients.m;i++){
            const SPH_COEFFICIENT_ENTRY& entry=sph_coefficients(i);
            for(int axis=1;axis<=3;axis++){
                matrix.Append_Entry_To_Current_Row(pm*(axis-1)+entry.i,dt*sqrt(entry.weight));
                matrix.Append_Entry_To_Current_Row(pm*(axis-1)+entry.j,-dt*sqrt(entry.weight));
                matrix.Finish_Row();}}
        matrix.Sort_Entries();}
    else{
        int pm=(n==0)?particles.X.m:n;
        matrix.Reset(pm*3);
        if(1){// segment aggregation
            for(int i=1;i<=fvm_coefficients.m;i++){
                const VECTOR<int,2>& e=segments(i);
                const T& c=fvm_coefficients(i);
                if(abs(c)<1e-16) continue;
                for(int axis=1;axis<=3;axis++){
                    matrix.Append_Entry_To_Current_Row(pm*(axis-1)+e(1),dt*sqrt(-c));
                    matrix.Append_Entry_To_Current_Row(pm*(axis-1)+e(2),-dt*sqrt(-c));
                    matrix.Finish_Row();}}}
        else{// vertex aggregation
            surface.mesh.segment_mesh->Initialize_Incident_Elements();
            for(int i=1;i<=surface.mesh.segment_mesh->incident_elements->m;i++){
                ARRAY<int>& list=(*surface.mesh.segment_mesh->incident_elements)(i);
                if(list.m==0) continue;
                for(int axis=1;axis<=TV::m;axis++){
                    for(int j=1;j<=list.m;j++){
                        VECTOR<int,2> k=segments(list(j));
                        T sn=dt*sqrt(-fvm_coefficients(list(j)));
                        if(k.y==i) matrix.Append_Entry_To_Current_Row(pm*(axis-1)+k.x,sn);
                        else matrix.Append_Entry_To_Current_Row(pm*(axis-1)+k.y,sn);}
                    matrix.Finish_Row();}}
            surface.mesh.segment_mesh->Delete_Auxiliary_Structures();}
        surface.mesh.Delete_Auxiliary_Structures();
        matrix.Sort_Entries();}
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,3> >::
Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const
{
}
//#####################################################################
// Function Enforce_Definiteness
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,3> >::
Enforce_Definiteness(const bool enforce_definiteness_input)
{
}
//#####################################################################
// Function Add_Implicit_Velocity_Independent_Forces
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,3> >::
Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
}
//#####################################################################
// Function Velocity_Dependent_Forces_Size
//#####################################################################
template<class T> int SURFACE_TENSION_FORCE<VECTOR<T,3> >::
Velocity_Dependent_Forces_Size() const
{
    int size=0;
    //if(apply_implicit_forces) size+=surface.mesh.elements.m;
    return size;
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,3> >::
Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const
{
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces_Second_Half
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,3> >::
Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV > F,const T time) const
{
}
//#####################################################################
// Function Add_Raw_Velocity_Dependent_Forces_First_Half
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,3> >::
Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const
{
}
//#####################################################################
// Function Initialize_CFL
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,3> >::
Initialize_CFL(ARRAY_VIEW<typename DEFORMABLES_FORCES<VECTOR<T,3> >::FREQUENCY_DATA> frequency)
{
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class T> T SURFACE_TENSION_FORCE<VECTOR<T,3> >::
CFL_Strain_Rate() const
{
    return FLT_MAX;
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,3> >::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    surface.mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,3> >::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV >* mpi_solids)
{
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class T> T SURFACE_TENSION_FORCE<VECTOR<T,3> >::
Potential_Energy(const T time) const
{
    T pe=0;
    return pe;
}
//#####################################################################
// Function Dump_Curvatures
//#####################################################################
template<class T> void SURFACE_TENSION_FORCE<VECTOR<T,3> >::
Dump_Curvatures() const
{
}
template class SURFACE_TENSION_FORCE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SURFACE_TENSION_FORCE<VECTOR<double,3> >;
#endif
