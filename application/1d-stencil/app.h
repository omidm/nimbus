#ifndef _1DSTENCIL
#define _1DSTENCIL

#include "lib/application.h"
#include "lib/job.h"
#include "lib/data.h"

class Main : public Job {
  public:
    Main(Application* app, JobType type);
    virtual void Execute(std::string params, const dataArray& da);
};

class Init : public Job {
  public:
    Init(Application* app, JobType type);
    virtual void Execute(std::string params, const dataArray& da);
};

class Print : public Job {
  public:
    Print(Application* app, JobType type);
    virtual void Execute(std::string params, const dataArray& da);
};

class ApplyLeft : public Job {
  public:
    ApplyLeft(Application* app, JobType type);
    virtual void Execute(std::string params, const dataArray& da);
};

class ApplyRight : public Job {
  public:
    ApplyRight(Application* app, JobType type);
    virtual void Execute(std::string params, const dataArray& da);
};

class UpdateLeft : public Job {
  public:
    UpdateLeft(Application* app, JobType type);
    virtual void Execute(std::string params, const dataArray& da);
};

class UpdateRight : public Job {
  public:
    UpdateRight(Application* app, JobType type);
    virtual void Execute(std::string params, const dataArray& da);
}; 

class ForLoop : public Job {
  public:
    ForLoop(Application* app, JobType type);
    virtual void Execute(std::string params, const dataArray& da);
};

class Vec : public Data {
  public:
    Vec(int);
    int size;
    int *arr;
    void Create();
};

class App : public Application {
  public:
    App();
    virtual void load();
};











#endif
