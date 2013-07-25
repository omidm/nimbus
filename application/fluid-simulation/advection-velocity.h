#ifndef ADVECTION_VELOCITY_H
#include "WATER_DRIVER.h"
#include "WATER_EXAMPLE.h"

namespace ADVECT_VELOCITY_NS {
// Import auxiliary class.
typedef typename ::PhysBAM::WATER_DRIVER::TV TV;
typedef typename ::PhysBAM::WATER_DRIVER::TV_INT TV_INT;

// Import function class.
typedef typename ::PhysBAM::WATER_DRIVER::ADVECT_VELOCITY_WORKER_T::AVERAGING_TYPE AVERAGING_TYPE;
typedef typename ::PhysBAM::WATER_DRIVER::ADVECT_VELOCITY_WORKER_T::INTERPOLATION_TYPE INTERPOLATION_TYPE;

void* advect_velocity_worker(void *arg);
void* advect_velocity_fetcher(void *arg);
void* advect_velocity_fetcher_network(void *arg);
void run(::PhysBAM::WATER_DRIVER &driver, ::PhysBAM::WATER_EXAMPLE &example,
    ::PhysBAM::WATER_DRIVER::ADVECT_VELOCITY_WORKER_T &INFO, AVERAGING_TYPE &averaging,
    INTERPOLATION_TYPE &interpolation, const TV_INT &segment_start);
}
#endif
