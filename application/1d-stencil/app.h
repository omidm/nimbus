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
    virtual void execute(std::string params, const DataArray& da);
    virtual Job * clone();
};

class Init : public Job {
  public:
    Init(Application* app, JobType type);
    virtual void execute(std::string params, const DataArray& da);
    virtual Job * clone();
};

class Print : public Job {
  public:
    Print(Application* app, JobType type);
    virtual void execute(std::string params, const DataArray& da);
    virtual Job * clone();
};

class ApplyLeft : public Job {
  public:
    ApplyLeft(Application* app, JobType type);
    virtual void execute(std::string params, const DataArray& da);
    virtual Job * clone();
};

class ApplyRight : public Job {
  public:
    ApplyRight(Application* app, JobType type);
    virtual void execute(std::string params, const DataArray& da);
    virtual Job * clone();
};

class UpdateLeft : public Job {
  public:
    UpdateLeft(Application* app, JobType type);
    virtual void execute(std::string params, const DataArray& da);
    virtual Job * clone();
};

class UpdateRight : public Job {
  public:
    UpdateRight(Application* app, JobType type);
    virtual void execute(std::string params, const DataArray& da);
    virtual Job * clone();
}; 

class ForLoop : public Job {
  public:
    ForLoop(Application* app, JobType type);
    virtual void execute(std::string params, const DataArray& da);
    virtual Job * clone();
};

class Vec : public Data {
  public:
    Vec(int);
    int size;
    int *arr;
    virtual void create();
    virtual Data * clone();
};

class App : public Application {
  public:
    App();
    virtual void load();
};











#endif
