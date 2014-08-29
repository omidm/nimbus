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
  * A Nimbus worker. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  * Modified: Chinmayee Shah <chinmayee.shah@stanford.edu>
  */

#include <iostream>  // NOLINT
#include <pthread.h>
#include <string>

#include "application/water_multiple/water_app.h"
#include "shared/log.h"
#include "shared/nimbus.h"
#include "shared/nimbus_types.h"
#include "simple_worker.h"
#include "worker/application.h"
#include "worker/worker_manager.h"


void PrintUsage() {
  std::cout << "ERROR: wrong arguments\n";
  std::cout << "Usage:\n";
  std::cout << "./worker\n";
  std::cout << "REQUIRED ARGUMENTS:\n";
  std::cout << "\t-sip [scheduler ip] -sport [scheduler port] -port [listening port]\n";
  std::cout << "OPTIONIAL:\n";
  std::cout << "\t-ip [ip address]\n";
  std::cout << "\t-s [loop counter]\n";
  std::cout << "\t-pn [part num]\n";
  std::cout << "\t-ithread [threading inside a job]\n";
  std::cout << "\t-othread [threading across jobs]\n";
}

int main(int argc, char *argv[]) {
  port_t listening_port, scheduler_port;
  std::string scheduler_ip, ip_address;
  bool ip_address_given = false;
  bool listening_port_given = false;
  bool scheduler_ip_given = false;
  bool scheduler_port_given = false;

  std::string log_id;

  // TODO(omidm): currently not used.
  size_t scale = 40;
  size_t part_num = 64;

  if (((argc - 1) % 2 != 0) || (argc < 3)) {
    PrintUsage();
    exit(-1);
  }

  for (int i = 1; i < argc; i = i + 2) {
    std::string tag = argv[i];
    std::string val = argv[i+1];
    if (tag == "-sip") {
      scheduler_ip = val;
      scheduler_ip_given = true;
    } else if (tag == "-ip") {
      ip_address = val;
      ip_address_given = true;
    } else if (tag == "-sport") {
      std::stringstream ss(val);
      ss >> scheduler_port;
      if (ss.fail()) {
        PrintUsage();
        exit(-1);
      }
      scheduler_port_given = true;
    } else if (tag == "-port") {
      log_id = val;
      std::stringstream ss(val);
      ss >> listening_port;
      if (ss.fail()) {
        PrintUsage();
        exit(-1);
      }
      listening_port_given = true;
    } else if (tag == "-s") {
      std::stringstream ss(val);
      ss >> scale;
      if (ss.fail()) {
        PrintUsage();
        exit(-1);
      }
    } else if (tag == "-pn") {
      std::stringstream ss(val);
      ss >> part_num;
      if (ss.fail()) {
        PrintUsage();
        exit(-1);
      }
    } else if (tag == "-ithread") {
      std::stringstream ss(val);
      ss >> WorkerManager::inside_job_parallism;
      if (ss.fail()) {
        PrintUsage();
        exit(-1);
      }
    } else if (tag == "-othread") {
      std::stringstream ss(val);
      ss >> WorkerManager::across_job_parallism;
      if (ss.fail()) {
        PrintUsage();
        exit(-1);
      }
    } else {
      PrintUsage();
      exit(-1);
    }
  }

  if (!scheduler_ip_given || !scheduler_port_given || !listening_port_given) {
    PrintUsage();
    exit(-1);
  }

  nimbus_initialize();
  std::cout << "Simple Worker is up!" << std::endl;
  application::WaterApp *app = new application::WaterApp();
  SimpleWorker * w = new SimpleWorker(scheduler_ip, scheduler_port, listening_port, app);
  if (ip_address_given) {
    w->set_ip_address(ip_address);
  }

  // TODO: Extra logging information for cache experiments, remove later on
#ifdef CACHE_LOG
  std::string log_file_name = "worker-log";
  log_file_name += "-";
  log_file_name += log_id;
  Log *cache_log = new Log(std::string(log_file_name));
  w->cache_log = cache_log;
  app->translator_log = cache_log;
#endif
  w->Run();
}
