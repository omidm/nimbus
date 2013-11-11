#ifndef __INIT_MAIN__
#define __INIT_MAIN__

#include <PhysBAM_Tools/Parallel_Computation/MPI_WORLD.h>
#include "projection_driver.h"

namespace PhysBAM {
// driver and mpi_world are set to be global for fast implementation.
// [TODO] Change it.
extern PROJECTION_DRIVER< VECTOR<float,2> >* driver;
extern MPI_WORLD* mpi_world; 
}  // namespace PhysBAM

void InitMain(int argc, char* argv[]);
void FinishMain();

#endif
