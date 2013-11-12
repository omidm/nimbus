#ifndef __PROJECTION_MAIN__ 
#define __PROJECTION_MAIN__

#include "shared/nimbus.h"
using nimbus::Data;
using nimbus::Application;

namespace PhysBAM {
template<class TV> class PROJECTION_DRIVER;
template<class T, int d> class VECTOR;
class MPI_WORLD;
}  // namespace PhysBAM

class App : public Application {
 public:
  App() {}
  void Load();
  // Initializes projection driver, and mpi world.
  void InitMain(int argc, char* argv[]);
  // Tears down the projection driver and the mpi world.
  void FinishMain();
  // Driver contains all the meta-data and functions to do projection on
  // a grid.
  PhysBAM::PROJECTION_DRIVER< PhysBAM::VECTOR<float,2> >* app_driver;
  // The PhysBAM mpi handler for communication, expected to be removed in 
  // the future.
  PhysBAM::MPI_WORLD* app_mpi_world;
};

// "ProfileData" is a dump data to enforce worker placement for now.
class ProfileData : public Data {
 public:
  ProfileData() {}
  virtual void Create() {}
  virtual Data * Clone() {return new ProfileData;}
};

#endif
