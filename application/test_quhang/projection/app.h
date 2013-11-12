#ifndef __PROJECTION_MAIN__ 
#define __PROJECTION_MAIN__

#include "shared/nimbus.h"
using nimbus::Data;
using nimbus::Job;
using nimbus::Application;

class Main : public Job {
  public:
    Main(Application* app);
    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job* Clone();
};

class Initialization : public Job {
 public:
  Initialization(Application* app); 
  virtual void Execute(Parameter params, const DataArray& da);
  virtual Job* Clone();
};

class SpawnOneIterationIfNeeded : public Job {
 public:
  SpawnOneIterationIfNeeded(Application* app);
  virtual void Execute(Parameter params, const DataArray& da);
  virtual Job* Clone();
};

class OneIteration : public Job {
 public:
  OneIteration(Application* app); 
  virtual void Execute(Parameter params, const DataArray& da);
  virtual Job* Clone();
};

class Finish : public Job {
 public:
  Finish(Application* app); 
  virtual void Execute(Parameter params, const DataArray& da);
  virtual Job* Clone();
};

class App : public Application {
 public:
  App() {}
  void Load();
};

// "ProfileData" is a dump data to enforce worker placement for now.
class ProfileData : public Data {
 public:
  ProfileData() {}
  virtual void Create() {}
  virtual Data * Clone() {return new ProfileData;}
};

#endif
