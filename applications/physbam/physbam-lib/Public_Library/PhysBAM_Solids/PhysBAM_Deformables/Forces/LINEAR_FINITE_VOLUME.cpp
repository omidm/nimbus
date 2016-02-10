//#####################################################################
// Copyright 2004-2007, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Igor Neverov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_2X2.h>
#include <PhysBAM_Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> LINEAR_FINITE_VOLUME<TV,d>::
LINEAR_FINITE_VOLUME(T_OBJECT& object,const T youngs_modulus,const T poissons_ratio,const T Rayleigh_coefficient)
    :DEFORMABLES_FORCES<TV>(dynamic_cast<PARTICLES<TV>&>(object.particles)),object(object),mesh(object.mesh),use_uniform_density(false),density_list(0)
{
    assert(poissons_ratio!=-1 && poissons_ratio!=.5);
    lambda=youngs_modulus*poissons_ratio/((1+poissons_ratio)*(1-2*poissons_ratio));
    mu=youngs_modulus/(2*(1+poissons_ratio));
    alpha=Rayleigh_coefficient*lambda;beta=Rayleigh_coefficient*mu;
    Initialize_Material_State(particles.X);
    mesh.Initialize_Incident_Elements();
    if(use_uniform_density){
        T total_mass=0;
        ARRAY<int> mesh_particles;mesh.elements.Flattened().Get_Unique(mesh_particles);
        for(int i=1;i<=mesh_particles.m;i++) total_mass+=particles.mass(mesh_particles(i));
        density=total_mass/object.Total_Size();
        if(density==0) density=TV::dimension==1?1:TV::dimension==2?100:1000;}
    else{
        density_list=new ARRAY<T>(mesh.elements.m);
        for(int i=1;i<=mesh.elements.m;i++){
            const VECTOR<int,d+1>& nodes=mesh.elements(i);
            T volume=object.Signed_Size(i);
            for(int j=1;j<=nodes.m;j++) (*density_list)(i)+=particles.mass(nodes(j))/(*mesh.incident_elements)(nodes(j)).m/volume;
            if((*density_list)(i)==0) (*density_list)(i)=TV::dimension==1?1:TV::dimension==2?100:1000;}}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> LINEAR_FINITE_VOLUME<TV,d>::
~LINEAR_FINITE_VOLUME()
{}
//#####################################################################
// Function Initialize_Material_State
//#####################################################################
namespace{
template<class T,int d> inline VECTOR<T,d> Normal(const MATRIX<T,d>& A)
{
    return VECTOR<T,d>();
}
template<class T> inline VECTOR<T,3> Normal(const MATRIX<T,3,2>& A)
{
    return A.Weighted_Normal().Normalized();
}
template<class T,int d> inline MATRIX<T,d> Pseudoinverse(const MATRIX<T,d>& A)
{
    return A.Inverse();
}
template<class T> inline MATRIX<T,2,3> Pseudoinverse(const MATRIX<T,3,2>& A)
{
    T scale=A.Parallelepiped_Measure();assert(scale);
    return (T)1/scale*A.Cofactor_Matrix().Transposed();
}
}
template<class TV,int d> void LINEAR_FINITE_VOLUME<TV,d>::
Initialize_Material_State(ARRAY_VIEW<const TV> X)
{
    Dm_inverse.Resize(mesh.elements.m,false,false);Bm.Resize(mesh.elements.m,false,false);
    if(TV::m>d) normals.Resize(mesh.elements.m,false,false);
    for(int t=1;t<=mesh.elements.m;t++){
        MATRIX<T,TV::m,d> Dm=Ds(X,t);
        if(TV::m>d) normals(t)=Normal(Dm);
        Dm_inverse(t)=Pseudoinverse(Dm);
        Bm(t)=-(T)1/FACTORIAL<d>::value*Dm.Cofactor_Matrix();}
}
//#####################################################################
// Function Add_Dependencies
//#####################################################################
template<class TV,int d> void LINEAR_FINITE_VOLUME<TV,d>::
Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
{
    mesh.Add_Dependencies(dependency_mesh);
}
//#####################################################################
// Function Update_Mpi
//#####################################################################
template<class TV,int d> void LINEAR_FINITE_VOLUME<TV,d>::
Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids)
{
    force_elements.Update(mesh.elements,particle_is_simulated);
    force_particles.Update(mesh.elements.Flattened(),particle_is_simulated);
}
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV,int d> void LINEAR_FINITE_VOLUME<TV,d>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        MATRIX<T,TV::m,d> G=Stress(t)*Bm(t);
        STRAIN_MEASURE<TV,d>::Distribute_Force(F,mesh.elements(t),G);}
}
//#####################################################################
// Function Add_Velocity_Dependent_Forces
//#####################################################################
template<class TV,int d> void LINEAR_FINITE_VOLUME<TV,d>::
Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        SYMMETRIC_MATRIX<T,TV::m> cauchy_strain_rate=(Ds(V,t)*Dm_inverse(t)).Symmetric_Part();
        MATRIX<T,TV::m,d> G=(2*beta*cauchy_strain_rate+alpha*cauchy_strain_rate.Trace())*Bm(t);
        STRAIN_MEASURE<TV,d>::Distribute_Force(F,mesh.elements(t),G);}
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class TV,int d> void LINEAR_FINITE_VOLUME<TV,d>::
Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const
{
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        MATRIX<T,TV::m,d> dG=Stress_Differential(dX,t)*Bm(t);
        STRAIN_MEASURE<TV,d>::Distribute_Force(dF,mesh.elements(t),dG);}
}
//#####################################################################
// Function Add_Force_Differential
//#####################################################################
template<class TV,int d> void LINEAR_FINITE_VOLUME<TV,d>::
Add_Force_Differential(const ARRAY<SYMMETRIC_MATRIX<T,TV::m> >& stress_differential,ARRAY_VIEW<TV> dF) const
{
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        MATRIX<T,TV::m,d> dG=stress_differential(t)*Bm(t);
        STRAIN_MEASURE<TV,d>::Distribute_Force(dF,mesh.elements(t),dG);}
}
//#####################################################################
// Function Intialize_CFL
//#####################################################################
namespace{
template<class T,int d> inline T Simplex_Minimum_Altitude(const MATRIX<T,d>& Dm_inverse)
{
    return Dm_inverse.Inverse().Simplex_Minimum_Altitude();
}
template<class T> inline T Simplex_Minimum_Altitude(const MATRIX<T,2,3>& Dm_inverse)
{
    return Dm_inverse.transpose.R_From_QR_Factorization().Inverse().Simplex_Minimum_Altitude();
}
}
template<class TV,int d> void LINEAR_FINITE_VOLUME<TV,d>::
Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency)
{
    // TODO: MPI
    T one_over_cfl_number=1/cfl_number,one_over_cfl_number_squared=sqr(one_over_cfl_number);

    ARRAY<FREQUENCY_DATA> fragment_particle_frequency(frequency.Size());
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        T local_density=use_uniform_density?density:(*density_list)(t);
        T altitude_squared=sqr(Simplex_Minimum_Altitude(Dm_inverse(t)));
        T elastic_squared=(lambda+2*mu)/(local_density*altitude_squared)*one_over_cfl_number_squared;
        T damping=(alpha+2*beta)/(local_density*altitude_squared)*one_over_cfl_number;
        const VECTOR<int,d+1>& nodes=mesh.elements(t);
        for(int j=1;j<=nodes.m;j++){FREQUENCY_DATA& data=fragment_particle_frequency(nodes[j]);
            data.elastic_squared=max(data.elastic_squared,elastic_squared);data.damping=max(data.damping,damping);}}
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
        frequency(p).elastic_squared+=fragment_particle_frequency(p).elastic_squared;
        frequency(p).damping+=fragment_particle_frequency(p).damping;}
}
//#####################################################################
// Function CFL_Strain_Rate
//#####################################################################
template<class TV,int d> typename TV::SCALAR LINEAR_FINITE_VOLUME<TV,d>::
CFL_Strain_Rate() const
{
    T max_strain_rate=0;
    for(ELEMENT_ITERATOR iterator(force_elements);iterator.Valid();iterator.Next()){int t=iterator.Data();
        max_strain_rate=max(max_strain_rate,(Ds(particles.V,t)*Dm_inverse(t)).Max_Abs());}
    return Robust_Divide(max_strain_per_time_step,max_strain_rate);
}
//#####################################################################
template class LINEAR_FINITE_VOLUME<VECTOR<float,2>,2>;
template class LINEAR_FINITE_VOLUME<VECTOR<float,3>,2>;
template class LINEAR_FINITE_VOLUME<VECTOR<float,3>,3>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LINEAR_FINITE_VOLUME<VECTOR<double,2>,2>;
template class LINEAR_FINITE_VOLUME<VECTOR<double,3>,2>;
template class LINEAR_FINITE_VOLUME<VECTOR<double,3>,3>;
#endif
