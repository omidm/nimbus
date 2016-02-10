//#####################################################################
// Copyright 2007, Craig Schroeder, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __LINEAR_TET_SPRINGS__
#define __LINEAR_TET_SPRINGS__

#include <PhysBAM_Tools/Data_Structures/FORCE_ELEMENTS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{
template<class T> class TETRAHEDRALIZED_VOLUME;

class TETRAHEDRON_MESH;
template<class T_input>
class LINEAR_TET_SPRINGS:public DEFORMABLES_FORCES<VECTOR<T_input,3> >,public SPRINGS_TAG
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef DEFORMABLES_FORCES<TV> BASE;
    using BASE::Invalidate_CFL;using BASE::particles;using BASE::use_rest_state_for_strain_rate;using BASE::max_strain_per_time_step;using BASE::use_implicit_velocity_independent_forces;
protected:
    using BASE::cfl_number;
    enum WORKAROUND{spring_count=7};
    typedef typename FORCE_ELEMENTS::ITERATOR ELEMENT_ITERATOR;
    typedef typename BASE::FREQUENCY_DATA FREQUENCY_DATA;

public:
    TETRAHEDRON_MESH& mesh;

    struct SPRING_PARAMETER
    {
        T youngs_modulus;
        T restlength,visual_restlength;
        T damping;
    };
    struct SPRING_STATE
    {
        SPRING_STATE()
           :id(0),coefficient(0),current_length(0)
        {}
        unsigned char id; // point faces 1,2,3,4, edge edges 5,6,7, colinear spring 8, no spring 9
        VECTOR<int,4> spring_nodes;
        T coefficient;
        T current_length;
        TV direction;
        TV weights;
    };

    ARRAY<VECTOR<SPRING_PARAMETER,spring_count> > spring_parameters;

protected:
    FORCE_ELEMENTS force_elements;
    
    ARRAY<SPRING_STATE> spring_states;
    ARRAY<VECTOR<T,6> > edge_restlength_squared; // (1,2),(1,3),(1,4),(2,3),(2,4),(3,4)
    bool use_springs_compressed_beyond_threshold; // only use the springs compressed enough
    T spring_compression_fraction_threshold;
    T minimum_edge_compression_squared; // threshold for judging cross product robust
    T minimum_sin; // threshold for judging cross product robust

public:
    static VECTOR<int,4> Spring_Nodes(unsigned char pair_id,const VECTOR<int,4>& n);
    static VECTOR<int,2> Edge_Indices(unsigned char pair_id);

    void Set_Stiffness(const T youngs_modulus_input)
    {Invalidate_CFL();
    for(int i=1;i<=spring_parameters.m;i++) for(int s=1;s<=spring_count;s++) spring_parameters(i)(s).youngs_modulus=youngs_modulus_input;}

    void Set_Stiffness(const ARRAY<VECTOR<T,spring_count> >& youngs_modulus_input)
    {Invalidate_CFL();
    for(int i=1;i<=spring_parameters.m;i++) for(int s=1;s<=spring_count;s++) spring_parameters(i)(s).youngs_modulus=youngs_modulus_input(i)(s);}

    void Set_Restlength(const ARRAY<VECTOR<T,spring_count> >& restlength_input)
    {Invalidate_CFL();
    for(int i=1;i<=spring_parameters.m;i++) for(int s=1;s<=spring_count;s++) spring_parameters(i)(s).restlength=restlength_input(i)(s);}

    virtual void Clamp_Restlength(const T clamped_restlength)
    {Invalidate_CFL();
    for(int i=1;i<=spring_parameters.m;i++) for(int s=1;s<=spring_count;s++) spring_parameters(i)(s).restlength=max(spring_parameters(i)(s).visual_restlength,clamped_restlength);}

//#####################################################################
    LINEAR_TET_SPRINGS(PARTICLES<TV>& particles,TETRAHEDRON_MESH& mesh,const bool implicit);
    virtual ~LINEAR_TET_SPRINGS(){}
    void Set_Restlength_From_Particles();
    void Set_Restlength_From_Material_Coordinates(ARRAY_VIEW<const TV> material_coordinates);
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
    void Initialize_CFL(ARRAY_VIEW<FREQUENCY_DATA> frequency) PHYSBAM_OVERRIDE;
    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE;
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const;
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Use_Springs_Compressed_Beyond_Threshold_Only(const bool use_springs_compressed_beyond_threshold_input=true,const T threshold_fraction=.25);
    void Set_Overdamping_Fraction(const T overdamping_fraction); // 1 is critically damped
    void Clamp_Restlength_With_Fraction_Of_Springs(const T fraction=.01);
    void Print_Restlength_Statistics(); // 1 is critically damped
    int Find_Shortest_Spring(const int tet,const VECTOR<int,4> element_nodes,VECTOR<int,4>& spring_nodes,T& minimum_signed_distance,TV& minimum_normal,TV& weights) const;
//#####################################################################
};

template<class T> LINEAR_TET_SPRINGS<T>*
Create_Tet_Springs(PARTICLES<VECTOR<T,3> >& particles,TETRAHEDRON_MESH& mesh,const T stiffness,const T overdamping_fraction,
    const bool use_compressed_by_threshold_only=true,const T fraction_compression=.1,const bool limit_time_step_by_strain_rate=true,const T max_strain_per_time_step=.1,
    const bool use_rest_state_for_strain_rate=true,const T restlength_enlargement_fraction=0,const bool verbose=true,const bool implicit=false);

template<class T> LINEAR_TET_SPRINGS<T>*
Create_Tet_Springs(TETRAHEDRALIZED_VOLUME<T>& volume,const T stiffness,
    const T overdamping_fraction,const bool use_compressed_by_threshold_only=true,const T fraction_compression=.1,const bool limit_time_step_by_strain_rate=true,
    const T max_strain_per_time_step=.1,const bool use_rest_state_for_strain_rate=true,const T restlength_enlargement_fraction=0,const bool verbose=true,const bool implicit=false);
}
#endif
