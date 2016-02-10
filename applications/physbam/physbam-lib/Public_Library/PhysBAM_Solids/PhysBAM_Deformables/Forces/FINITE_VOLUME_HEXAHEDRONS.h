//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FINITE_VOLUME_HEXAHEDRONS
//#####################################################################
#ifndef __FINITE_VOLUME_HEXAHEDRONS__
#define __FINITE_VOLUME_HEXAHEDRONS__

#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class T> class STRAIN_MEASURE_HEXAHEDRONS;
template<class T, int d> class CONSTITUTIVE_MODEL;
template<class T, int d> class DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE;
template<class T, int d> class DIAGONALIZED_STRESS_DERIVATIVE;
template<class T, int d> class ISOTROPIC_CONSTITUTIVE_MODEL;
template<class T, int d> class ANISOTROPIC_CONSTITUTIVE_MODEL;

template<class T_input>
class FINITE_VOLUME_HEXAHEDRONS:public DEFORMABLES_FORCES<VECTOR<T_input,3> >,public FINITE_VOLUME_TAG
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::particles;using BASE::cfl_number;using BASE::max_strain_per_time_step;
    typedef typename FORCE_ELEMENTS::ITERATOR ELEMENT_ITERATOR;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

    STRAIN_MEASURE_HEXAHEDRONS<T>& strain_measure;
    CONSTITUTIVE_MODEL<T,3>& constitutive_model;
    ARRAY<VECTOR<T,8> > Be_scales;
    ARRAY<VECTOR<MATRIX<T,3>,8> > U;
    ARRAY<VECTOR<ARRAY<TV>,8> > De_inverse_hat; // actually H (Dm H)^-1 Fp^-1
    ARRAY<VECTOR<DIAGONAL_MATRIX<T,3>,8> > Fe_hat;
    ARRAY<VECTOR<DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3>,8> >* dPi_dFe;
    ARRAY<VECTOR<DIAGONALIZED_STRESS_DERIVATIVE<T,3>,8> >* dP_dFe;
    ARRAY<VECTOR<T,8> >* Be_scales_save;
    ARRAY<VECTOR<MATRIX<T,3>,8> >* V;
    ARRAY<VECTOR<SYMMETRIC_MATRIX<T,3>,8> >* node_stiffness;
    ARRAY<VECTOR<MATRIX<T,3>,8> >* edge_stiffness;
    bool use_uniform_density;
    bool destroy_data;
private:
    ISOTROPIC_CONSTITUTIVE_MODEL<T,3>* isotropic_model;
    ANISOTROPIC_CONSTITUTIVE_MODEL<T,3>* anisotropic_model;
    FORCE_ELEMENTS force_elements;
    FORCE_ELEMENTS force_particles;
    ARRAY<T>* density_list;
    T density;
public:

    FINITE_VOLUME_HEXAHEDRONS(STRAIN_MEASURE_HEXAHEDRONS<T>& strain_measure,CONSTITUTIVE_MODEL<T,3>& constitutive_model);
    ~FINITE_VOLUME_HEXAHEDRONS();

    void Save_V()
    {if(!V) V=new ARRAY<VECTOR<MATRIX<T,3>,8> >();}

    void Use_Quasistatics()
    {Save_V();Save_Stress_Derivative();}

    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE
    {constitutive_model.enforce_definiteness=enforce_definiteness_input;}

    void Update_Be_Scales()
    {Be_scales=strain_measure.DmH_determinant;Be_scales*=(T)-1;}

    void Initialize_Be_Scales_Save()
    {Be_scales_save=new ARRAY<VECTOR<T,8> >(Be_scales);}

    void Copy_Be_Scales_Save_Into_Be_Scales(const ARRAY<int>& map)
    {Be_scales=Be_scales_save->Subset(map);}

//#####################################################################
    void Save_Stress_Derivative();
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const PHYSBAM_OVERRIDE;
    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
//#####################################################################
};

template<class T> FINITE_VOLUME_HEXAHEDRONS<T>*
Create_Finite_Volume(HEXAHEDRALIZED_VOLUME<T>& hexahedralized_volume,CONSTITUTIVE_MODEL<T,3>* constitutive_model,const bool limit_time_step_by_strain_rate=true,
    const T max_strain_per_time_step=.1,const bool use_rest_state_for_strain_rate=true);

template<class T> FINITE_VOLUME_HEXAHEDRONS<T>*
Create_Quasistatic_Finite_Volume(HEXAHEDRALIZED_VOLUME<T>& hexahedralized_volume,CONSTITUTIVE_MODEL<T,3>* constitutive_model);
}
#endif
