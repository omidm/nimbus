#include <iostream>
#include <pthread.h>
#include "scheduler.h"
#include "application.h"
#include "../application/1d-stencil/app.h"

using namespace std;

#define LISTENING_PORT 5983


int main (int argc, char *argv[]) {

  std::cout << "Nimbus is up!" << std::endl;
  Scheduler * s = new Scheduler(LISTENING_PORT);
  
  App * app0 = new App();
  s->addApp(app0);
  s->appMap[0]->loadApp();
  
  s->run();
}

