//#####################################################################
// Copyright 2003-2009, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonthan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_BODY_COLLECTION
//#####################################################################
#ifndef __DEFORMABLE_BODY_COLLECTION__
#define __DEFORMABLE_BODY_COLLECTION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#include <PhysBAM_Tools/Utilities/Find_Type.h>
#endif
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLES_FORCES.h>
namespace PhysBAM{

template<class TV> class BINDING_LIST;
template<class TV> class DEFORMABLE_GEOMETRY_COLLECTION;
template<class TV> class COLLISION_GEOMETRY_COLLECTION;
template<class TV> class DEFORMABLE_OBJECT_COLLISIONS;
template<class TV> class PARTICLES;
template<class TV> class SOFT_BINDINGS;
template<class TV> class TRIANGLE_REPULSIONS;
template<class TV> class TRIANGLE_COLLISIONS;
template<class TV> class TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY;
template<class TV> class COLLISION_PENALTY_FORCES;
template<class TV> class DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES;
template<class TV> class MPI_SOLIDS;
template<class TV> class DEFORMABLES_FORCES;
template<class TV> class TRIANGLE_COLLISION_PARAMETERS;
class ARRAY_COLLECTION;

template<class TV>
class DEFORMABLE_BODY_COLLECTION:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename DEFORMABLES_FORCES<TV>::FREQUENCY_DATA T_FREQUENCY_DEFORMABLE;
public:
    PARTICLES<TV>& particles;
    DEFORMABLE_GEOMETRY_COLLECTION<TV>& deformable_geometry;
    bool simulate; // TODO: use one of those per fragment
    bool owns_data;

    // These are pruned for MPI.
    ARRAY<int> simulated_particles;
    ARRAY<int> dynamic_particles;

    BINDING_LIST<TV>& binding_list;
    SOFT_BINDINGS<TV>& soft_bindings;
    ARRAY<DEFORMABLES_FORCES<TV>*> deformables_forces;

    MPI_SOLIDS<TV>* mpi_solids;
    ARRAY<T_FREQUENCY_DEFORMABLE> frequency; // hertz for deformable CFL
    T cfl_number;
    T cfl_elastic,cfl_damping;
    bool implicit_damping;

    bool print_diagnostics;
    bool print_residuals;
    bool print_energy;
    int iterations_used_diagnostic;

    DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES<TV>* deformables_example_forces_and_velocities;
    TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY<TV>& triangle_repulsions_and_collisions_geometry;
    TRIANGLE_REPULSIONS<TV>& triangle_repulsions;
    TRIANGLE_COLLISIONS<TV>& triangle_collisions;
    DEFORMABLE_OBJECT_COLLISIONS<TV>& collisions;
    ARRAY<COLLISION_PENALTY_FORCES<TV>*> collision_penalty_forces;
    bool use_embedded_collisions;
    bool use_nonembedded_self_collision; // TODO: have one of these per fragment
    bool check_stale;

    DEFORMABLE_BODY_COLLECTION(DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES<TV>* deformables_example_forces_and_velocities_input,COLLISION_GEOMETRY_COLLECTION<TV>& collision_body_list,
        ARRAY_COLLECTION* array_collection=0);
    virtual ~DEFORMABLE_BODY_COLLECTION();

    template<class T_FORCE> T_FORCE
    Find_Force(const int index=1)
    {return Find_Type<T_FORCE>(deformables_forces,index);}

    template<class T_FORCE> const T_FORCE
    Find_Force(const int index=1) const
    {return Find_Type<T_FORCE>(deformables_forces,index);}

//#####################################################################
    void Set_CFL_Number(const T cfl_number_input=.5);
    void Read(const STREAM_TYPE stream_type,const std::string& prefix,const int frame,const int static_frame,const bool include_static_variables,const bool read_from_every_process);
    void Write(const STREAM_TYPE stream_type,const std::string& prefix,const int frame,const int static_frame,const bool include_static_variables,const bool write_from_every_process) const;
    int Add_Force(DEFORMABLES_FORCES<TV>* force);
    void Initialize(TRIANGLE_COLLISION_PARAMETERS<TV>& triangle_collisions_parameters);
    void Adjust_Mesh_For_Self_Collision();
    void Adjust_Mesh_For_Self_Repulsion(); // same as Adjust_Mesh_For_Self_Collision, but assumes only velocities change
    void Update_Collision_Penalty_Forces_And_Derivatives();
    void Update_Stored_Repulsion_Connectivity(const TRIANGLE_REPULSIONS<TV>& repulsions);
    void Read_Dynamic_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame);
    void Write_Dynamic_Variables(const STREAM_TYPE stream_type,const std::string& prefix,const int frame) const;
    void Update_Simulated_Particles();
    void Update_Simulated_Particles(DEFORMABLES_EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities);
    void Set_Mpi_Solids(MPI_SOLIDS<TV>* mpi_solids);
    void Update_CFL();
    T CFL(const bool verbose=false);
    T CFL_Elastic_And_Damping();
    T CFL_Elastic();
    T CFL_Damping();
    T CFL_Strain_Rate();

    void Update_Position_Based_State(const T time,const bool is_position_update);
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F_full,const T time) const;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V_full,ARRAY_VIEW<TV> F_full,const T time) const;
    void Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V_full,ARRAY_VIEW<TV> F_full,const T scale,const T time) const;

    void Save_Potential_Energy(const T time);
    void Compute_Energy_Error(ARRAY_VIEW<const TV> velocity_save,const T time,const T dt);
    void Add_Energy_Correction_Force(ARRAY_VIEW<const TV> velocity_save,const int energy_correction_iterations,const T time,const T dt);
    void Test_Energy(const T time);
    void Test_Force_Derivatives(const T time);
    void Compute_Previously_Applied_Forces();
    void Setup_Set_Velocity_From_Positions(const T time,const bool is_position_update,const bool reset_alphas);
    void Store_Velocities();
    void Read(const STREAM_TYPE,const std::string& prefix,const std::string& static_prefix,const int frame,const int static_frame,const bool include_static_variables,
        const bool read_from_every_process);
    void Write(const STREAM_TYPE,const std::string& prefix,const std::string& static_prefix,const int frame,const int static_frame,const bool include_static_variables,
        const bool write_from_every_process) const;
//#####################################################################
};
}
#endif
