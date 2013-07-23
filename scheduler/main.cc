#include <iostream>
#include <pthread.h>
#include "lib/scheduler.h"
#include "lib/application.h"
#include "../application/1d-stencil/app.h"

using namespace std;

#define LISTENING_PORT 5983


int main (int argc, char *argv[]) {

  std::cout << "Nimbus is up!" << std::endl;
  Scheduler * s = new Scheduler(LISTENING_PORT);
  s->run();
}

