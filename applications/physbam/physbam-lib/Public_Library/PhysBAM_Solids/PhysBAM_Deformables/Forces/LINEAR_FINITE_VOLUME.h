//#####################################################################
// Copyright 2004-2007, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Igor Neverov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_FINITE_VOLUME
//#####################################################################
#ifndef __LINEAR_FINITE_VOLUME__
#define __LINEAR_FINITE_VOLUME__

#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <PhysBAM_Tools/Math_Tools/FACTORIAL.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/CONSTITUTIVE_MODELS_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV,int d>
class LINEAR_FINITE_VOLUME:public DEFORMABLES_FORCES<TV>,public FINITE_VOLUME_TAG
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d>::OBJECT T_OBJECT;
    typedef typename MESH_POLICY<d>::MESH T_MESH;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;using BASE::max_strain_per_time_step;using BASE::cfl_number;
    typedef typename FORCE_ELEMENTS::ITERATOR ELEMENT_ITERATOR;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

    T_OBJECT& object;
    T_MESH& mesh;
    T lambda,mu; // Lame coefficients
    T alpha,beta;
    bool use_uniform_density;
private:
    ARRAY<TV> normals;
    ARRAY<MATRIX<T,d,TV::m> > Dm_inverse;
    ARRAY<MATRIX<T,TV::m,d> > Bm;
    ARRAY<MATRIX<T,d,TV::m> >* Dm_inverse_save;
    ARRAY<MATRIX<T,TV::m,d> >* Bm_save;
    FORCE_ELEMENTS force_elements;
    FORCE_ELEMENTS force_particles;
    ARRAY<T>* density_list;
    T density;
public:

    LINEAR_FINITE_VOLUME(T_OBJECT& object,const T youngs_modulus,const T poissons_ratio,const T Rayleigh_coefficient);
    virtual ~LINEAR_FINITE_VOLUME();

    void Initialize_Dm_Inverse_Save()
    {Dm_inverse_save=new ARRAY<MATRIX<T,d,TV::m> >(Dm_inverse);}

    void Copy_Dm_Inverse_Save_Into_Dm_Inverse(const ARRAY<int>& map)
    {Dm_inverse=Dm_inverse_save->Subset(map);}

    void Initialize_Bm_Save()
    {Bm_save=new ARRAY<MATRIX<T,TV::m,d> >(Bm);}

    void Copy_Bm_Save_Into_Bm(const ARRAY<int>& map)
    {Bm=Bm_save->Subset(map);}

    MATRIX<T,TV::m,d> Ds(ARRAY_VIEW<const TV> X,const int simplex) const
    {return STRAIN_MEASURE<TV,d>::Ds(X,mesh.elements(simplex));}

    SYMMETRIC_MATRIX<T,TV::m> Stress(const int simplex) const
    {SYMMETRIC_MATRIX<T,TV::m> cauchy_strain=(Ds(particles.X,simplex)*Dm_inverse(simplex)).Symmetric_Part()-1;
    if(TV::m>d) cauchy_strain+=SYMMETRIC_MATRIX<T,TV::m>::Outer_Product(normals(simplex));
    return 2*mu*cauchy_strain+lambda*cauchy_strain.Trace();}

    SYMMETRIC_MATRIX<T,TV::m> Stress_Differential(ARRAY_VIEW<const TV> dX,const int simplex) const
    {SYMMETRIC_MATRIX<T,TV::m> cauchy_strain_differential=(Ds(dX,simplex)*Dm_inverse(simplex)).Symmetric_Part();
    return 2*mu*cauchy_strain_differential+lambda*cauchy_strain_differential.Trace();}

//#####################################################################
    void Initialize_Material_State(ARRAY_VIEW<const TV> X);
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const PHYSBAM_OVERRIDE;
    void Add_Force_Differential(const ARRAY<SYMMETRIC_MATRIX<T,TV::m> >& stress_differential,ARRAY_VIEW<TV> dF) const;
    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
//#####################################################################
};

template<class T,class T_OBJECT> LINEAR_FINITE_VOLUME<typename T_OBJECT::VECTOR_T,T_OBJECT::MESH::dimension>*
Create_Linear_Finite_Volume(T_OBJECT& object,const T youngs_modulus=3e6,const T poissons_ratio=.475,const T Rayleigh_coefficient=.05,const bool limit_time_step_by_strain_rate=true,
    const T max_strain_per_time_step=.1,const bool use_rest_state_for_strain_rate=true)
{
    typedef typename T_OBJECT::VECTOR_T TV;static const int d=T_OBJECT::MESH::dimension;
    LINEAR_FINITE_VOLUME<TV,d>* fvm=new LINEAR_FINITE_VOLUME<TV,d>(object,youngs_modulus,poissons_ratio,Rayleigh_coefficient);
    fvm->Limit_Time_Step_By_Strain_Rate(limit_time_step_by_strain_rate,max_strain_per_time_step);
    fvm->Use_Rest_State_For_Strain_Rate(use_rest_state_for_strain_rate);
    return fvm;
}

}
#endif
