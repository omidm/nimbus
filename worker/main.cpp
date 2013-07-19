#include <iostream>
#include "scheduler.h"


int main (int argc, char *argv[]) {

  std::cout << "Hello World!\n";
  Worker * w = new Worker();
  w->run();

}

