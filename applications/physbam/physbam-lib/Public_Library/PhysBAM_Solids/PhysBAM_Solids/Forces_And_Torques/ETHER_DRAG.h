//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ETHER_DRAG
//#####################################################################
#ifndef __ETHER_DRAG__
#define __ETHER_DRAG__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/POINTWISE_FORCE.h>
namespace PhysBAM{

template<class T_GRID>
class ETHER_DRAG:public POINTWISE_FORCE<typename T_GRID::VECTOR_T>
{
private:
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_TV;
public:
    typedef POINTWISE_FORCE<TV> BASE;
    using BASE::particles;using BASE::rigid_body_collection;using BASE::mpi_solids;using BASE::force_particles;using BASE::force_rigid_body_particles;
    typedef typename FORCE_ELEMENTS::ITERATOR ELEMENT_ITERATOR;

    bool use_constant_wind;
    T constant_wind_viscosity,constant_wind_angular_viscosity;
    TV constant_wind;
    bool use_spatially_varying_wind;
    T spatially_varying_wind_viscosity;
    RANGE<TV> spatially_varying_wind_domain;
    T_GRID V_grid;
    T_ARRAYS_TV* spatially_varying_wind;
    LINEAR_INTERPOLATION_UNIFORM<T_GRID,TV> vector_interpolation;

    ETHER_DRAG(PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,ARRAY<int>* influenced_particles_input,
        ARRAY<int>* influenced_rigid_body_particles_input,T dynamic_ether_viscosity=0,T angular_viscosity=0);
    ETHER_DRAG(PARTICLES<TV>& particles_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const bool influence_all_particles_input,
        const bool influence_all_rigid_body_particles_input,T dynamic_ether_viscosity=0,T angular_viscosity=0);
    virtual ~ETHER_DRAG();

    void Use_No_Drag()
    {use_constant_wind=use_spatially_varying_wind=false;}

    void Use_Constant_Wind(const T viscosity_input,const TV& wind_input=TV())
    {if(!viscosity_input && wind_input==TV()) Use_No_Drag();
    else{use_constant_wind=true;constant_wind_viscosity=viscosity_input;constant_wind=wind_input;}}

    void Use_Spatially_Varying_Wind(const T viscosity_input,const RANGE<TV>& domain_input,const T_GRID& grid_input,T_ARRAYS_TV& V_input)
    {use_spatially_varying_wind=true;spatially_varying_wind_viscosity=viscosity_input;
    spatially_varying_wind_domain=domain_input;V_grid=grid_input;spatially_varying_wind=&V_input;}

    TV Spatially_Varying_Wind_Velocity(const TV& X) const
    {return vector_interpolation.Clamped_To_Array(V_grid,*spatially_varying_wind,X);}

//#####################################################################
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE;
    void Add_Force_Differential(ARRAY_VIEW<const TV> dX,ARRAY_VIEW<TV> dF,const T time) const PHYSBAM_OVERRIDE;
    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE;
    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
