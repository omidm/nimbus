#include "scheduler.h"



Scheduler::Scheduler(unsigned int port)
{
  appId = 0;
  this->port = port;
};

void Scheduler::run()
{
  cout << "Running the Scheduler" << endl;

  loadUserCommands();
  loadWorkerCommands();
  
  user_interface_thread = new boost::thread(
      boost::bind(&Scheduler::setupUI,this));

  worker_interface_thread = new boost::thread(
      boost::bind(&Scheduler::setupWI,this));


  user_interface_thread->join();
  worker_interface_thread->join();  

};

void Scheduler::setupWI()
{
  server = new Server(port, this);
  server->run();
};

void Scheduler::setupUI()
{
  while(true)
  {
    cout << "command: ";
    string token ("runapp");
    string str, cm;
    vector<int> args;
    getline(cin, str);
    parseCommand(str, userCmSet, cm, args);
    cout << "you typed: " << cm << endl;
  }

}

void Scheduler::addApp(Application * app)
{
  app->id = appId;
  appMap[appId] = app;
  appId++;
};

void Scheduler::delApp(Application * app)
{
  appMap.erase(app->id);
};

void Scheduler::loadUserCommands()
{
  stringstream cms("loadapp runapp killapp haltapp resumeapp quit");
  while (true) 
  {
    string word;
    cms >> word;
    if(cms.fail()) break;
    userCmSet.insert(word);
  }
};

void Scheduler::loadWorkerCommands()
{
  stringstream cms("runjob killjob haltjob resumejob jobdone createdata copydata deletedata");
  while (true) 
  {
    string word;
    cms >> word;
    if(cms.fail()) break;
    workerCmSet.insert(word);
  }

};


