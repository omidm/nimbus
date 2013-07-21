#ifndef _SCHEDULER
#define _SCHEDULER

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/thread.hpp>
#include "lib/scheduler_server.h"
#include "lib/application.h"
#include "lib/cluster.h"
#include "lib/worker.h"
#include "lib/parser.h"

#define USER_CM_FILE "cmu.txt"
#define WORKER_CM_FILE "cmw.txt"

using namespace std;

typedef set<string> CmSet;

class Worker;
typedef std::map<int, Worker*> WorkerMap;

class Scheduler {
  public:
    Scheduler(uint listening_port);
    
    Computer host;
    unsigned int port;
    unsigned int appId;

    SchedulerServer* server;
    
    //AppMap appMap;
    WorkerMap workerMap;
    ClusterMap clusterMap;
    
    void run();

    void loadClusterMap(std::string);

    void delWorker(Worker *);
    Worker * addWorker();
    Worker * getWorker(int);

  private:
    void setupUI();
    void setupWI();
    
    boost::thread* worker_interface_thread;
    boost::thread* user_interface_thread;

    void loadUserCommands();
    void loadWorkerCommands();
    
    CmSet userCmSet;
    CmSet workerCmSet;
};

#endif
