//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __WATER_DRIVER__
#define __WATER_DRIVER__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Advection/ADVECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
#include "application/water_multiple/app_utils.h"
#include "application/water_multiple/options.h"
#include "shared/nimbus.h"

namespace PhysBAM{


template<class TV> class WATER_EXAMPLE;

template<class TV>
class WATER_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;
    typedef typename ADVECTION_POLICY<GRID<TV> >::ADVECTION_SEMI_LAGRANGIAN_SCALAR T_ADVECTION_SEMI_LAGRANGIAN_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename GRID_ARRAYS_POLICY<GRID<TV> >::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;

protected:
    WATER_EXAMPLE<TV>& example;

public:
    bool init_phase;
    int current_frame;
    T time;
    int output_number;

    WATER_DRIVER(WATER_EXAMPLE<TV>& example);
    virtual ~WATER_DRIVER();

    void Initialize(const nimbus::Job *job,
                    const nimbus::DataArray &da,
                    bool set_boundary_conditions = true);

    bool InitializeIncompressibleProjectionHelper(
        const application::DataConfig& data_config,
        const GRID<TV>& grid_input,
        INCOMPRESSIBLE_UNIFORM<GRID<TV> >* incompressible,
        PROJECTION_DYNAMICS_UNIFORM<GRID<TV> >* projection);

    bool InitializeParticleLevelsetEvolutionHelper(
        const application::DataConfig& data_config,
        const GRID<TV>& grid_input,
        PARTICLE_LEVELSET_EVOLUTION_UNIFORM<GRID<TV> >*
        particle_levelset_evolution);

    void CalculateFrameImpl(const nimbus::Job *job,
                            const nimbus::DataArray &da,
                            const bool set_boundary_conditions,
                            const T dt);

    void WriteFrameImpl(const nimbus::Job *job,
                        const nimbus::DataArray &da,
                        const bool set_boundary_conditions,
                        const T dt);

    bool ComputeOccupiedBlocksImpl(const nimbus::Job *job,
                                   const nimbus::DataArray &da,
                                   T dt);

    bool AdjustPhiWithObjectsImpl(const nimbus::Job *job,
                                  const nimbus::DataArray &da,
                                  T dt);

    bool AdvectPhiImpl(const nimbus::Job *job,
                       const nimbus::DataArray &da,
                       T dt);

    bool StepParticlesImpl(const nimbus::Job *job,
                           const nimbus::DataArray &da,
                           T dt);

    bool AdvectRemovedParticlesImpl(const nimbus::Job *job,
                           const nimbus::DataArray &da,
                           T dt);

    bool AdvectVImpl(const nimbus::Job *job,
                     const nimbus::DataArray &da,
                     T dt);

    bool ApplyForcesImpl(const nimbus::Job *job,
                    const nimbus::DataArray &da,
                    T dt);

    bool ModifyLevelSetImpl(const nimbus::Job *job,
                            const nimbus::DataArray &da,
                            T dt);

    bool ExtrapolatePhiImpl(const nimbus::Job *job,
                            const nimbus::DataArray &da,
                            T dt);

    bool AdjustPhiImpl(const nimbus::Job *job,
                       const nimbus::DataArray &da,
                       T dt);

    bool DeleteParticlesImpl(const nimbus::Job *job,
                             const nimbus::DataArray &da,
                             T dt);

    bool ReincorporateParticlesImpl(const nimbus::Job *job,
                                    const nimbus::DataArray &da,
                                    T dt);

    bool ProjectionImpl (const nimbus::Job *job,
                         const nimbus::DataArray &da,
                         T dt);

    bool ExtrapolationImpl (const nimbus::Job *job,
                          const nimbus::DataArray &da,
                          T dt);

    void Write_Output_Files(const int frame);

    void Write_Substep(const std::string& title,const int substep,const int level=0);
    void Run(RANGE<TV_INT>& domain,const T dt,const T time);

//#####################################################################
};
}
#endif
