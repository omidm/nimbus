#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2006, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MELTING_EXAMPLE
//#####################################################################
#ifndef __MELTING_EXAMPLE__
#define __MELTING_EXAMPLE__

#include <PhysBAM_Tools/Data_Structures/DATA_STRUCTURES_FORWARD.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_POLICY.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/MELTING_PARAMETERS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class TV,int d> class EMBEDDED_MATERIAL_SURFACE;

template<class TV_input,int d>
class MELTING_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV_input> >
{
    typedef TV_input TV;typedef typename TV::SCALAR T;
    typedef typename RED_GREEN_POLICY<VECTOR<T,d> >::GRID_T T_RED_GREEN_GRID;
    typedef typename RED_GREEN_POLICY<VECTOR<T,d> >::RED_SIMPLEX T_RED_SIMPLEX;
    typedef typename MATRIX_POLICY<VECTOR<T,d> >::DIAGONAL_MATRIX T_DIAGONAL_MATRIX;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::output_directory;using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::use_melting;using BASE::stream_type;

    MELTING_PARAMETERS<TV,d> melting_parameters;
    T maximum_velocity;
private:
    T near_interface_threshold;
    ARRAY<T_DIAGONAL_MATRIX>* refinement_Fe_hat;
public:
    bool initialized;

    MELTING_EXAMPLE(const STREAM_TYPE stream_type,const int number_of_regions,const typename FLUIDS_PARAMETERS<GRID<TV> >::TYPE type);
    virtual ~MELTING_EXAMPLE();

    void Initialize_Phi() PHYSBAM_OVERRIDE
    {PHYSBAM_FATAL_ERROR("wrong version of Initialize_Phi called");}

    void Setup_Initial_Refinement() PHYSBAM_OVERRIDE
    {PHYSBAM_FATAL_ERROR("wrong version of Setup_Initial_Refinement called");}

//#####################################################################
    virtual void Initialize_Deformable_And_Rigid_Bodies(){}
    virtual void Initialize_Forces(){}
    virtual void Initialize_Grid(const int object,T_RED_GREEN_GRID& grid)=0;
    virtual void Initialize_Levelset_Velocity(const int object,ARRAY<VECTOR<T,d> >& V)=0;
    virtual void Initialize_Phi(const int object,ARRAY<T>& phi)=0;
    virtual void Initialize_Particle_Positions_And_Velocities(const int object)=0;
    virtual void Refinement_Criteria_Precomputation(const int object,const T dt);
    virtual bool Refinement_Criteria(const int object,const T_RED_SIMPLEX* triangle);
    virtual void Melting_Levelset_Substep(const int object,const T dt,const T time);
    void Initialize_Bodies() PHYSBAM_OVERRIDE;
    void Add_Deformable_Melting_Object(const int index);
    void Add_Rigid_Melting_Object(const int index);
private:
    int Add_Melting_Object(const typename MELTING_PARAMETERS<TV,d>::BODY_TYPE type,const int index,EMBEDDED_MATERIAL_SURFACE<TV,d>& embedding);
public:
    void Mark_Nodes_Inside_Objects(ARRAY<bool,VECTOR<int,2> >& inside_objects); // only in 2d
    void Update_Solids_Topology_For_Melting(const T dt,const T time,const bool reinitialize);
    void Update_Topology(const T dt);
    bool Interface_Refinement_Criteria(const int object,const T_RED_SIMPLEX* simplex);
    bool Deformation_Refinement_Criteria(const int object,const T_RED_SIMPLEX* simplex);
    void Setup_Initial_Refinement(const int object); // unrelated to base class version
    void Read_Output_Files_Solids(const int frame) PHYSBAM_OVERRIDE;
    void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE;
private:
    struct INITIAL_REFINEMENT_CRITERIA_HELPER{ARRAY<T>* phi;int* maximum_depth;};
    static bool Initial_Refinement_Criteria_Helper(INITIAL_REFINEMENT_CRITERIA_HELPER* helper,const typename RED_GREEN_POLICY<VECTOR<T,d> >::RED_SIMPLEX* simplex);
    static bool Refinement_Criteria_Helper(PAIR<MELTING_EXAMPLE<TV,d>*,int>* helper,const typename RED_GREEN_POLICY<VECTOR<T,d> >::RED_SIMPLEX* simplex);
//#####################################################################
};
}
#endif
#endif
