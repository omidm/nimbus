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
  * This file has the main function that launches Nimbus scheduler.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#define DEBUG_MODE

#include <iostream> // NOLINT
#include "./scheduler_1d.h"
#include "shared/nimbus.h"
#include "shared/nimbus_types.h"
#include "shared/scheduler_command.h"
#include "shared/parser.h"

void PrintUsage() {
  std::cout << "ERROR: wrong arguments\n";
  std::cout << "Usage:\n";
  std::cout << "./scheduler\n";
  std::cout << "REQUIRED ARGUMENTS:\n";
  std::cout << "\t-port [listening port]\n";
  std::cout << "OPTIONIAL:\n";
  std::cout << "\t-wn [initial worker num, DEFAULT: 2]\n";
}

int main(int argc, char *argv[]) {
  port_t listening_port;
  size_t worker_num;
  bool listening_port_given = false;
  bool worker_num_given = false;

  if (((argc - 1) % 2 != 0) || (argc < 3)) {
    PrintUsage();
    exit(-1);
  }

  for (int i = 1; i < argc; i = i + 2) {
    std::string tag = argv[i];
    std::string val = argv[i+1];
    if (tag == "-port") {
      std::stringstream ss(val);
      ss >> listening_port;
      if (ss.fail()) {
        PrintUsage();
        exit(-1);
      }
      listening_port_given = true;
    } else if (tag == "-wn") {
      std::stringstream ss(val);
      ss >> worker_num;
      if (ss.fail()) {
        PrintUsage();
        exit(-1);
      }
      worker_num_given = true;
    } else {
      PrintUsage();
      exit(-1);
    }
  }

  if (!listening_port_given) {
    PrintUsage();
    exit(-1);
  }

  nimbus::nimbus_initialize();

  Scheduler1D * s = new Scheduler1D(listening_port);

  if (worker_num_given) {
    s->set_min_worker_to_join(worker_num);
    dbg(DBG_SCHED, "Set min initial number of workers to %d.\n", worker_num);
  } else {
    dbg(DBG_SCHED, "Nothig provided for min initial number of workers, using default.\n");
  }

  s->Run();
}

