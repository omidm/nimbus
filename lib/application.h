#ifndef _APPLICATION
#define _APPLICATION

#include <set>
#include "job.h"
#include "data.h"


class Application
{
  public:

    int id;
    int priority;
    DataMap dataMap;
    JobMap jobMap;

    Application();
    virtual void loadApp();
    void Run();

};


typedef std::map<int, Application*> AppMap;




#endif
