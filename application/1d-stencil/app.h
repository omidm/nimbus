#ifndef _1DSTENCIL
#define _1DSTENCIL

#include "worker/application.h"
#include "worker/job.h"
#include "worker/data.h"
#include "shared/nimbus_types.h"

using namespace nimbus; // NOLINT

class Main : public Job {
  public:
    Main(Application* app, JobType type);
    virtual void Execute(std::string params, const DataArray& da);
    virtual Job * Clone();
};

class Init : public Job {
  public:
    Init(Application* app, JobType type);
    virtual void Execute(std::string params, const DataArray& da);
    virtual Job * Clone();
};

class Print : public Job {
  public:
    Print(Application* app, JobType type);
    virtual void Execute(std::string params, const DataArray& da);
    virtual Job * Clone();
};

class ApplyLeft : public Job {
  public:
    ApplyLeft(Application* app, JobType type);
    virtual void Execute(std::string params, const DataArray& da);
    virtual Job * Clone();
};

class ApplyRight : public Job {
  public:
    ApplyRight(Application* app, JobType type);
    virtual void Execute(std::string params, const DataArray& da);
    virtual Job * Clone();
};

class UpdateLeft : public Job {
  public:
    UpdateLeft(Application* app, JobType type);
    virtual void Execute(std::string params, const DataArray& da);
    virtual Job * Clone();
};

class UpdateRight : public Job {
  public:
    UpdateRight(Application* app, JobType type);
    virtual void Execute(std::string params, const DataArray& da);
    virtual Job * Clone();
}; 

class ForLoop : public Job {
  public:
    ForLoop(Application* app, JobType type);
    virtual void Execute(std::string params, const DataArray& da);
    virtual Job * Clone();
};

class Vec : public Data {
  public:
    Vec(int);
    int size;
    int *arr;
    virtual void Create();
    virtual Data * Clone();
};

class App : public Application {
  public:
    App();
    virtual void load();
};











#endif
