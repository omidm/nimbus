/*
 * Copyright 2013 Stanford University.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the
 *   distribution.
 *
 * - Neither the name of the copyright holders nor the names of
 *   its contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
 * THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

 /*
  * A test application that a worker can send and receive commands.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#include <pthread.h>
#include <iostream>  // NOLINT

#include "lib/nimbus.h"

using ::std::cout;
using ::std::endl;

class TestApplication : public Application {
  public:
  void load() {
    cout << "Loading Nimbus test application." << std::endl;
  }

  void start(SchedulerClient* scheduler) {
    cout << "Starting Nimbus test application." << std::endl;
    const char* commands[] = {
      "no-op",
      "halt 53",
      "run job3 job0,job1,job2 job5,job6 data4,data5 data5 blah",
      "copy   data4          host34   ",
      "copy       data5 192.244.11.2        ",
      "query status job3",
      "run 34",
      NULL
    };

    for (int i = 0; commands[i] != NULL; i++) {
      SchedulerCommand cm(commands[i]);
      cout << "Sending command:  " << cm.ToNetworkData() << std::endl;
      scheduler->sendCommand(&cm);
    }

    while (1) {
      sleep(1);
      SchedulerCommand* comm = scheduler->receiveCommand();
      if (comm->ToNetworkData() != "no-command")
        cout << "Received command: " << comm->ToNetworkData() << std::endl;
    }
  }
};

int main(int argc, char *argv[]) {
  std::cout << "Worker is up!" << std::endl;
  TestApplication * app0 = new TestApplication();
  Worker * w = new Worker(NIMBUS_SCHEDULER_PORT, app0);
  w->run();
}
