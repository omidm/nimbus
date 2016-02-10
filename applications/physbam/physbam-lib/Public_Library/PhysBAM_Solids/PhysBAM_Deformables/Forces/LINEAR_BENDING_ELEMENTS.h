//#####################################################################
// Copyright 2007, Geoffrey Irving, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_BENDING_ELEMENTS
//#####################################################################
#ifndef __LINEAR_BENDING_ELEMENTS__
#define __LINEAR_BENDING_ELEMENTS__

#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV>
class LINEAR_BENDING_ELEMENTS:public DEFORMABLES_FORCES<TV>
{
    typedef typename TV::SCALAR T;
    enum WORKAROUND {d=TV::m};
    typedef typename MESH_POLICY<d-1>::MESH T_MESH;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

    T_MESH& mesh;
    T stiffness,damping;
private:
    ARRAY<T> stiffness_matrix_diagonal;
    SPARSE_MATRIX_FLAT_NXN<T> stiffness_matrix_upper;
public:

    LINEAR_BENDING_ELEMENTS(PARTICLES<TV>& particles,T_MESH& mesh);
    ~LINEAR_BENDING_ELEMENTS();

    void Set_Stiffness(const T stiffness_input)
    {stiffness=stiffness_input;}

    void Set_Damping(const T damping_input)
    {damping=damping_input;}

    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE
    {}

    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE// linear bending elements require large edge spring stiffness to work, so we assume those dominate the CFL
    {return FLT_MAX;}

    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE
    {}

//#####################################################################
    void Compute_Stiffness_Matrix(ARRAY_VIEW<const TV> X);
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    T Compute_Energy() const;
private:
    void Add_Force(const T scale,ARRAY_VIEW<const TV> X,ARRAY_VIEW<TV> F) const;
//#####################################################################
};

template<class TV> LINEAR_BENDING_ELEMENTS<TV>*
Create_Linear_Bending_Elements(PARTICLES<TV>& particles,typename MESH_POLICY<TV::m-1>::MESH& mesh,const typename TV::SCALAR stiffness,
    const typename TV::SCALAR damping)
{
    LINEAR_BENDING_ELEMENTS<TV>* bend=new LINEAR_BENDING_ELEMENTS<TV>(particles,mesh);
    bend->Compute_Stiffness_Matrix(particles.X);
    bend->Set_Stiffness(stiffness);
    bend->Set_Damping(damping);
    return bend;
}

template<class T_OBJECT> LINEAR_BENDING_ELEMENTS<typename T_OBJECT::VECTOR_T>*
Create_Linear_Bending_Elements(T_OBJECT& object,const typename T_OBJECT::VECTOR_T::SCALAR stiffness,const typename T_OBJECT::VECTOR_T::SCALAR damping)
{
    return Create_Linear_Bending_Elements(object.particles,object.mesh,stiffness,damping);
}

}
#endif
