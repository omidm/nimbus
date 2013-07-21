#ifndef _APPLICATION
#define _APPLICATION

#include <set>
#include "lib/job.h"
#include "lib/data.h"
#include "lib/client.h"

class Scheduler;

class Application {
  public:

    int id;
    int priority;
    DataMap dataMap;
    JobMap jobMap;

    Application();
    virtual void loadApp();
    virtual void run(Client* scheduler);

};

typedef std::map<int, Application*> AppMap;

#endif
