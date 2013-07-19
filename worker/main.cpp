
#include <iostream>
#include <pthread.h>
#include "scheduler.h"
#include "worker.h"
#include "application.h"
#include "../application/1d-stencil/app.h"

using namespace std;

int main (int argc, char *argv[]) {

  std::cout << "Worker is up!" << std::endl;
  Worker * w = new Worker();
  
  App * app0 = new App();
  app0->loadApp();
  
  w->run();

}

