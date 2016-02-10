//#####################################################################
// Copyright 2003-2007, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_ALTITUDE_SPRINGS_2D
//#####################################################################
#ifndef __LINEAR_ALTITUDE_SPRINGS_2D__
#define __LINEAR_ALTITUDE_SPRINGS_2D__

#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS.h>
namespace PhysBAM{

template<class T_input>
class LINEAR_ALTITUDE_SPRINGS_2D:public LINEAR_ALTITUDE_SPRINGS<VECTOR<T_input,2>,2>
{
    typedef T_input T;
private:
    typedef VECTOR<T,2> TV;
    typedef LINEAR_ALTITUDE_SPRINGS<TV,2> BASE;typedef typename BASE::SPRING_STATE SPRING_STATE;typedef typename BASE::SPRING_PARAMETER SPRING_PARAMETER;
    using BASE::use_springs_compressed_beyond_threshold;using BASE::spring_compression_fraction_threshold;using BASE::print_number_used;
    using BASE::force_elements;
public:
    using BASE::particles;using BASE::Invalidate_CFL;
    using BASE::max_strain_per_time_step;using BASE::use_rest_state_for_strain_rate;using BASE::spring_states;
    using BASE::parameters;using BASE::plastic_parameters;using BASE::use_plasticity;using BASE::compute_half_forces;
    using BASE::mesh;typedef typename BASE::ELEMENT_ITERATOR ELEMENT_ITERATOR;

    LINEAR_ALTITUDE_SPRINGS_2D(PARTICLES<TV>& particles,TRIANGLE_MESH& mesh)
        :LINEAR_ALTITUDE_SPRINGS<TV,2>(particles,mesh)
    {}

    virtual ~LINEAR_ALTITUDE_SPRINGS_2D()
    {}

//#####################################################################
    void Set_Stiffness_Based_On_Reduced_Mass(const T scaling_coefficient); // assumes mass and restlength are already defined
    void Set_Restlength_From_Particles();
    void Set_Restlength_From_Material_Coordinates(ARRAY_VIEW<TV> material_coordinates);
    void Set_Overdamping_Fraction(const T overdamping_fraction); // 1 is critically damped
    void Ensure_Minimum_Overdamping_Fraction(const T overdamping_fraction=1); // 1 is critically damped
    void Update_Position_Based_State(const T time,const bool is_position_update) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    int Velocity_Dependent_Forces_Size() const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE;
//#####################################################################
};

template<class T> LINEAR_ALTITUDE_SPRINGS_2D<T>*
Create_Altitude_Springs(PARTICLES<VECTOR<T,2> >& particles,TRIANGLE_MESH& mesh,const T stiffness=200,const T overdamping_fraction=2,
    const bool use_compressed_by_threshold_only=true,const T fraction_compression=.1,const bool limit_time_step_by_strain_rate=true,
    const T max_strain_per_time_step=.1,const bool use_rest_state_for_strain_rate=true,const T restlength_enlargement_fraction=0,const bool verbose=true)
{
    return Create_Altitude_Springs_Base(particles,mesh,stiffness,overdamping_fraction,use_compressed_by_threshold_only,fraction_compression,
        limit_time_step_by_strain_rate,max_strain_per_time_step,use_rest_state_for_strain_rate,restlength_enlargement_fraction,verbose);
}

template<class T> LINEAR_ALTITUDE_SPRINGS_2D<T>*
Create_Altitude_Springs(TRIANGULATED_AREA<T>& object,const T stiffness=200,const T overdamping_fraction=2,
    const bool use_compressed_by_threshold_only=true,const T fraction_compression=.1,const bool limit_time_step_by_strain_rate=true,
    const T max_strain_per_time_step=.1,const bool use_rest_state_for_strain_rate=true,const T restlength_enlargement_fraction=0,const bool verbose=true)
{
    return Create_Altitude_Springs(dynamic_cast<PARTICLES<VECTOR<T,2> >&>(object.particles),object.mesh,stiffness,overdamping_fraction,use_compressed_by_threshold_only,fraction_compression,
        limit_time_step_by_strain_rate,max_strain_per_time_step,use_rest_state_for_strain_rate,restlength_enlargement_fraction,verbose);
}

}
#endif
