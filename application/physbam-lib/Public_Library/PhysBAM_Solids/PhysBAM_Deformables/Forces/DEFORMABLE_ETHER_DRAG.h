//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_ETHER_DRAG
//#####################################################################
#ifndef __DEFORMABLE_ETHER_DRAG__
#define __DEFORMABLE_ETHER_DRAG__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/POINTWISE_DEFORMABLE_FORCE.h>
namespace PhysBAM{

template<class TV>
class DEFORMABLE_ETHER_DRAG:public POINTWISE_DEFORMABLE_FORCE<TV>
{
private:
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef POINTWISE_DEFORMABLE_FORCE<TV> BASE;
    using BASE::particles;using BASE::mpi_solids;using BASE::force_particles;
    typedef typename FORCE_ELEMENTS::ITERATOR ELEMENT_ITERATOR;

    bool use_constant_wind;
    T constant_wind_viscosity,constant_wind_angular_viscosity;
    TV constant_wind;
    bool use_spatially_varying_wind;
    T spatially_varying_wind_viscosity;
    RANGE<TV> spatially_varying_wind_domain;
    GRID<TV> V_grid;
    ARRAY<TV,TV_INT>* spatially_varying_wind;
    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,TV> vector_interpolation;

    DEFORMABLE_ETHER_DRAG(PARTICLES<TV>& particles_input,ARRAY<int>* influenced_particles_input,T dynamic_ether_viscosity=0,T angular_viscosity=0)
        :POINTWISE_DEFORMABLE_FORCE<TV>(particles_input,influenced_particles_input),use_constant_wind(dynamic_ether_viscosity!=0),
        constant_wind_viscosity(dynamic_ether_viscosity),constant_wind_angular_viscosity(angular_viscosity),use_spatially_varying_wind(false),spatially_varying_wind_viscosity(0)
    {}

    DEFORMABLE_ETHER_DRAG(PARTICLES<TV>& particles_input,const bool influence_all_particles_input,T dynamic_ether_viscosity=0,T angular_viscosity=0)
        :POINTWISE_DEFORMABLE_FORCE<TV>(particles_input,influence_all_particles_input),use_constant_wind(dynamic_ether_viscosity!=0),
        constant_wind_viscosity(dynamic_ether_viscosity),constant_wind_angular_viscosity(angular_viscosity),use_spatially_varying_wind(false),spatially_varying_wind_viscosity(0)
    {}

    virtual ~DEFORMABLE_ETHER_DRAG()
    {}

    void Use_No_Drag()
    {use_constant_wind=use_spatially_varying_wind=false;}

    void Use_Constant_Wind(const T viscosity_input,const TV& wind_input=TV())
    {if(!viscosity_input && wind_input==TV()) Use_No_Drag();
    else{use_constant_wind=true;constant_wind_viscosity=viscosity_input;constant_wind=wind_input;}}

    void Use_Spatially_Varying_Wind(const T viscosity_input,const RANGE<TV>& domain_input,const GRID<TV>& grid_input,ARRAY<TV,TV_INT>& V_input)
    {use_spatially_varying_wind=true;spatially_varying_wind_viscosity=viscosity_input;
    spatially_varying_wind_domain=domain_input;V_grid=grid_input;spatially_varying_wind=&V_input;}

    TV Spatially_Varying_Wind_Velocity(const TV& X) const
    {return vector_interpolation.Clamped_To_Array(V_grid,*spatially_varying_wind,X);}

    void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const PHYSBAM_OVERRIDE
    {}

    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE
    {}

    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE
    {}

//#####################################################################
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
