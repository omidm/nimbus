#ifndef _WORKER
#define _WORKER

#include "cluster.h"
#include "data.h"
#include "job.h"

#define DataSet std::set<Data>
#define JobSet std::set<Job>

class Worker
{
  public:
    Worker();

    int id;
    Computer host;
    unsigned int port;

    DataSet dataSet;
    JobSet jobs;

    void run();

    void addJob(Job *);
    void delJob(Job *);
    

};


typedef std::map<int, Worker*> WorkerMap;



#endif
