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

#include <pthread.h>
#include <iostream>  // NOLINT

#include "simple_worker.h"
#include "shared/nimbus.h"
#include "shared/nimbus_types.h"
#include "worker/application.h"

#include "application/projection_multi_workers/app.h"

int main(int argc, char *argv[]) {
  port_t listening_port;
  int rankID = 0;
  if (argc < 2) {
    std::cout << "ERROR: provide an integer (1 to 4)." <<
      std::endl;
    exit(-1);
  }
  if (*argv[1] == '1') {
    rankID = 1;
    std::cout << "1" << std::endl;
    listening_port = WORKER_PORT_1;
  } else if (*argv[1] == '2') {
    rankID = 2;
    std::cout << "2" << std::endl;
    listening_port = WORKER_PORT_2;
  } else if (*argv[1] == '3') {
    rankID = 3;
    listening_port = WORKER_PORT_3;
  } else if (*argv[1] == '4') {
    rankID = 4;
    listening_port = WORKER_PORT_4;
  } else {
    std::cout << "ERROR: argument should be an integer (1 to 4)." <<
      std::endl;
    exit(-1);
  }
  nimbus_initialize();
  std::cout << "Simple Worker is up!" << std::endl;
  App *app = new App(rankID);
  SimpleWorker * w = new SimpleWorker(NIMBUS_SCHEDULER_IP,
      NIMBUS_SCHEDULER_PORT, listening_port, app);
  w->Run();
}


