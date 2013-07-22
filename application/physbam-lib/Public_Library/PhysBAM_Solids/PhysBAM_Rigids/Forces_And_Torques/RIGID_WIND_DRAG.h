//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_WIND_DRAG
//#####################################################################
#ifndef __RIGID_WIND_DRAG__
#define __RIGID_WIND_DRAG__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_SIMPLEX_POLICY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/INTERPOLATION_COLLIDABLE_UNIFORM_FORWARD.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGIDS_FORCES.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY;
template<class TV> class GRID;

template<class TV>
class RIGID_WIND_DRAG:public RIGIDS_FORCES<TV>
{
    typedef VECTOR<int,TV::m> TV_INT;
    typedef typename TV::SCALAR T;
public:
    typedef RIGIDS_FORCES<TV> BASE;
    typedef typename BASE::FREQUENCY_DATA RIGID_FREQUENCY_DATA;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::dimension-1>::OBJECT T_SIMPLICIAL_OBJECT;
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::dimension-1>::SIMPLEX T_SIMPLEX;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;

    RIGID_BODY<TV>* rigid_body;
    bool use_constant_wind;
    T constant_wind_viscosity;
    TV constant_wind;
    bool use_spatially_varying_wind;
    T spatially_varying_wind_viscosity;
    GRID<TV> V_grid;
    T_ARRAYS_VECTOR *spatially_varying_wind;
    T wind_density;
    T_ARRAYS_SCALAR *spatially_varying_wind_density,*spatially_varying_wind_pressure;
    T linear_normal_viscosity; // uses vertex normals
    mutable T surface_area;
private:
    static const LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM<GRID<TV>,T> interpolation;
    static const LINEAR_INSIDE_CONSTANT_OUTSIDE_INTERPOLATION_UNIFORM<GRID<TV>,TV> vector_interpolation;
    struct OPTIMIZATION{
        OPTIMIZATION()
            :area_over_m((T)0)
        {}

        T area_over_m;TV inward_normal,center,wind_velocity;
    };
    ARRAY<OPTIMIZATION> optimization;
    ARRAY<TV> vertex_normals; // TODO: use triangulated_surface.vertex_normals
public:

    RIGID_WIND_DRAG(RIGID_BODY<TV>& rigid_body_input);
    virtual ~RIGID_WIND_DRAG();

    void Use_Constant_Wind(const T viscosity_input,const TV& wind_input=TV())
    {use_constant_wind=true;constant_wind_viscosity=viscosity_input;constant_wind=wind_input;}

    void Use_Spatially_Varying_Wind(const T viscosity_input,const GRID<TV>& grid_input,T_ARRAYS_VECTOR& V_input)
    {use_spatially_varying_wind=true;spatially_varying_wind_viscosity=viscosity_input;
    V_grid=grid_input;spatially_varying_wind=&V_input;}

    void Set_Wind_Density(const T wind_density_input)
    {wind_density=wind_density_input;}

    void Set_Wind_Density(T_ARRAYS_SCALAR& density_input)
    {spatially_varying_wind_density=&density_input;}

    void Set_Wind_Pressure(T_ARRAYS_SCALAR& pressure_input) // only valid for volumetric objects
    {spatially_varying_wind_pressure=&pressure_input;}

    void Use_Linear_Normal_Viscosity(const T viscosity_input)
    {linear_normal_viscosity=viscosity_input;}

private:
    TV Spatially_Varying_Wind_Velocity(const TV& X) const
    {return vector_interpolation.Clamped_To_Array(V_grid,*spatially_varying_wind,X);}

    T Spatially_Varying_Wind_Density(const TV& X) const
    {return interpolation.Clamped_To_Array(V_grid,*spatially_varying_wind_density,X);}

    T Spatially_Varying_Wind_Pressure(const TV& X) const
    {return interpolation.Clamped_To_Array(V_grid,*spatially_varying_wind_pressure,X);}

    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE
    {}

    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE
    {return FLT_MAX;}

    void Initialize_CFL(ARRAY_VIEW<RIGID_FREQUENCY_DATA> rigid_frequency) PHYSBAM_OVERRIDE
    {}

    void Update_Mpi(const ARRAY<bool>& particle_is_simulated)
    {}

    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const
    {}

    TV Add_Velocity_Independent_Forces_Helper(TV relative_velocity,int t) const;

public:

//#####################################################################
    void Update_Position_Based_State(const T time) PHYSBAM_OVERRIDE;
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
