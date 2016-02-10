//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __SMOKE_DRIVER__
#define __SMOKE_DRIVER__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include "application/smoke/app_utils.h"
#include "application/smoke/options.h"
#include "shared/nimbus.h"

namespace PhysBAM{


template<class TV> class SMOKE_EXAMPLE;

template<class TV>
class SMOKE_DRIVER
{
    typedef typename TV::SCALAR T;
    typedef typename TV::template REBIND<int>::TYPE TV_INT;

protected:
    SMOKE_EXAMPLE<TV>& example;

public:
    bool init_phase;
    int current_frame;
    T time;
    int output_number;

    SMOKE_DRIVER(SMOKE_EXAMPLE<TV>& example);
    virtual ~SMOKE_DRIVER();

    void InitializeFirstDistributed(
        const nimbus::Job *job,
        const nimbus::DataArray &da);
    void Initialize(const nimbus::Job *job,
                    const nimbus::DataArray &da);
    void InitializeUseCache(const nimbus::Job *job,
                    const nimbus::DataArray &da);

    bool InitializeProjectionHelper(
	const application::DataConfig& data_config,
	const GRID<TV>& grid_input,
	PROJECTION_UNIFORM<GRID<TV> >* projection);

    void CalculateFrameImpl(const nimbus::Job *job,
                            const nimbus::DataArray &da,
                            const bool set_boundary_conditions,
                            const T dt);

    void WriteOutputSplitImpl(const nimbus::Job *job,
                              const nimbus::DataArray &da,
                              const bool set_boundary_conditions,
                              const T dt, const int rank);

    bool ScalarAdvanceImpl(const nimbus::Job *job,
			   const nimbus::DataArray &da,
			   const T dt);

    bool ConvectImpl(const nimbus::Job *job,
		     const nimbus::DataArray &da,
		     const T dt);

    bool UpdateGhostVelocitiesImpl(const nimbus::Job *job,
				   const nimbus::DataArray &da,
				   T dt);

    bool UpdateGhostDensitiesImpl(const nimbus::Job *job,
				   const nimbus::DataArray &da,
				   T dt);

    bool ProjectionCalculateBoundaryConditionPartOneImpl(
							 const nimbus::Job *job,
							 const nimbus::DataArray &da,
							 T dt);

    bool ProjectionCalculateBoundaryConditionPartTwoImpl(
							 const nimbus::Job *job,
							 const nimbus::DataArray &da,
							 T dt);

    bool ProjectionConstructMatrixImpl(const nimbus::Job *job,
                                       const nimbus::DataArray &da,
                                       T dt);

    bool ProjectionWrapupImpl(const nimbus::Job *job,
			      const nimbus::DataArray &da,
			      T dt);

    void Write_Output_Files(const int frame, int rank = -1);

    void Write_Substep(const std::string& title,const int substep,const int level=0);

//#####################################################################
};
}
#endif
