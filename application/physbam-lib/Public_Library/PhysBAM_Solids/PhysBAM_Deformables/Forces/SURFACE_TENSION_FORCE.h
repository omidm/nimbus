//#####################################################################
// Copyright 2011, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SURFACE_TENSION_FORCE
//#####################################################################
#ifndef __SURFACE_TENSION_FORCE__
#define __SURFACE_TENSION_FORCE__

#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SURFACE_TENSION_FORCE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/POINTWISE_DEFORMABLE_FORCE.h>
namespace PhysBAM{

template<class TV>
class SURFACE_TENSION_FORCE:public DEFORMABLES_FORCES<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    struct UNUSABLE{};

    UNUSABLE& surface;
    GEOMETRY_PARTICLES<TV>& particles;
    bool apply_implicit_forces;
    bool apply_tangential_implicit,use_sph_laplace;
    GRID<TV> grid;
    T dt;

    SURFACE_TENSION_FORCE(MESH_OBJECT<TV,SIMPLEX_MESH<TV::m> >& surface_input,T surface_tension_coefficient_input);
    virtual ~SURFACE_TENSION_FORCE();

//#####################################################################
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE{}
    void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const PHYSBAM_OVERRIDE{}
    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE{}
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE{}
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE{}
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE{}
    int Velocity_Dependent_Forces_Size() const PHYSBAM_OVERRIDE{return 0;}
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const PHYSBAM_OVERRIDE{}
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE{}
    void Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const PHYSBAM_OVERRIDE{}
    void Initialize_CFL(ARRAY_VIEW<typename BASE::FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE{}
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE{return FLT_MAX;}
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE{}
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE{}
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE{return 0;}
    void Dump_Curvatures() const{}
    void Compute_Curvatures(ARRAY_VIEW<T> curvature) const{}
    void Construct_Implicit_Matrix_First_Half(SPARSE_MATRIX_FLAT_MXN<T>& matrix,int n=0){}
//#####################################################################
};

template<class T_input>
class SURFACE_TENSION_FORCE<VECTOR<T_input,2> >:public DEFORMABLES_FORCES<VECTOR<T_input,2> >
{
public:
    typedef T_input T;
    typedef VECTOR<T,2> TV;
    typedef DEFORMABLES_FORCES<TV> BASE;
    SEGMENTED_CURVE_2D<T>& surface;
    GEOMETRY_PARTICLES<TV>& particles;
    T surface_tension_coefficient;
    ARRAY<T> coefficients;
    ARRAY<TV> normal;
    ARRAY<MATRIX<T,TV::m,TV::m-1> > tangential;
    ARRAY<T> sqrt_coefficients;
    GRID<TV> grid;
    T dt;
    bool apply_explicit_forces;
    bool apply_implicit_forces;
    bool apply_tangential_implicit,use_sph_laplace;

    SURFACE_TENSION_FORCE(SEGMENTED_CURVE_2D<T>& surface_input,T surface_tension_coefficient_input);
    virtual ~SURFACE_TENSION_FORCE();

//#####################################################################
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const PHYSBAM_OVERRIDE;
    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    int Velocity_Dependent_Forces_Size() const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const PHYSBAM_OVERRIDE;
    void Initialize_CFL(ARRAY_VIEW<typename DEFORMABLES_FORCES<VECTOR<T,2> >::FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
    void Dump_Curvatures() const;
    void Compute_Curvatures(ARRAY_VIEW<T> curvature) const;
    void Construct_Implicit_Matrix_First_Half(SPARSE_MATRIX_FLAT_MXN<T>& matrix,int n=0);
//#####################################################################
};

template<class T_input>
class SURFACE_TENSION_FORCE<VECTOR<T_input,3> >:public DEFORMABLES_FORCES<VECTOR<T_input,3> >
{
public:
    typedef T_input T;
    typedef VECTOR<T,3> TV;
    typedef VECTOR<int,3> TV_INT;
    struct SPH_COEFFICIENT_ENTRY{
        int i,j;
        T weight;
    };
    typedef DEFORMABLES_FORCES<TV> BASE;
    TRIANGULATED_SURFACE<T>& surface;
    GEOMETRY_PARTICLES<TV>& particles;
    T surface_tension_coefficient;
    ARRAY<SPH_COEFFICIENT_ENTRY> sph_coefficients;
    ARRAY<T> fvm_coefficients;
    ARRAY<VECTOR<int,2> > segments;
    ARRAY<T> vertex_areas;
    T dt;
    GRID<TV> grid;
    bool apply_explicit_forces,apply_implicit_forces,apply_tangential_implicit,use_sph_laplace;

    SURFACE_TENSION_FORCE(TRIANGULATED_SURFACE<T>& surface_input,T surface_tension_coefficient_input);
    virtual ~SURFACE_TENSION_FORCE();

//#####################################################################
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const PHYSBAM_OVERRIDE;
    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    int Velocity_Dependent_Forces_Size() const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Raw_Velocity_Dependent_Forces_First_Half(ARRAY<TRIPLE<int,int,T> >& data) const PHYSBAM_OVERRIDE;
    void Initialize_CFL(ARRAY_VIEW<typename DEFORMABLES_FORCES<VECTOR<T,3> >::FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
    void Dump_Curvatures() const;
    void Compute_Curvatures(ARRAY_VIEW<T> curvature) const{}
    void Construct_Implicit_Matrix_First_Half(SPARSE_MATRIX_FLAT_MXN<T>& matrix,int n=0);
//#####################################################################
protected:
    T Cotangent(int t,int a,int b){
        int c;
        for(int i=1;i<=TV::m;i++){
            c=surface.mesh.elements(t)(i);
            if(c!=a && c!=b) break;}
        T A2,B2,C2,AB;
        A2=(particles.X(b)-particles.X(c)).Magnitude_Squared();
        B2=(particles.X(a)-particles.X(c)).Magnitude_Squared();
        C2=(particles.X(a)-particles.X(b)).Magnitude_Squared();
        AB=sqrt(A2*B2);
        return (A2+B2-C2)*sqrt(C2)*(T)0.5/AB;
    }
};

}
#endif
