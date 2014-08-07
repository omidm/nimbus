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
  * A Nimbus worker for job spawner application. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include <iostream>  // NOLINT
#include "shared/nimbus.h"
#include "../../application/job_spawner/app.h"

using namespace nimbus; // NOLINT


void PrintUsage() {
  std::cout << "ERROR: wrong arguments\n";
  std::cout << "Usage:\n";
  std::cout << "./worker\n";
  std::cout << "REQUIRED ARGUMENTS:\n";
  std::cout << "\t-sip [scheduler ip] -sport [scheduler port] -port [listening port]\n";
  std::cout << "OPTIONIAL:\n";
  std::cout << "\t-ip [ip address]\n";
  std::cout << "\t-lc [loop counter]\n";
  std::cout << "\t-pn [part num]\n";
  std::cout << "\t-cpp [chunk per part]\n";
  std::cout << "\t-cs [chunk size]\n";
  std::cout << "\t-bw [bandwidth]\n";
  std::cout << "\t-sn [stage num]\n";
  std::cout << "\t-jlu [job length usec]\n";
}


int main(int argc, char *argv[]) {
  port_t listening_port, scheduler_port;
  std::string scheduler_ip, ip_address;
  bool ip_address_given = false;
  bool listening_port_given = false;
  bool scheduler_ip_given = false;
  bool scheduler_port_given = false;

  size_t counter = 30;
  size_t part_num = 100;
  size_t chunk_per_part = 1;
  size_t chunk_size = 50;
  size_t bandwidth = 10;
  size_t stage_num = 10;
  size_t job_length_usec = 0;



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
      std::stringstream ss(val);
      ss >> listening_port;
      if (ss.fail()) {
        PrintUsage();
        exit(-1);
      }
      listening_port_given = true;
    } else if (tag == "-lc") {
      std::stringstream ss(val);
      ss >> counter;
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
    } else if (tag == "-cpp") {
      std::stringstream ss(val);
      ss >> chunk_per_part;
      if (ss.fail()) {
        PrintUsage();
        exit(-1);
      }
    } else if (tag == "-cs") {
      std::stringstream ss(val);
      ss >> chunk_size;
      if (ss.fail()) {
        PrintUsage();
        exit(-1);
      }
    } else if (tag == "-bw") {
      std::stringstream ss(val);
      ss >> bandwidth;
      if (ss.fail()) {
        PrintUsage();
        exit(-1);
      }
    } else if (tag == "-sn") {
      std::stringstream ss(val);
      ss >> stage_num;
      if (ss.fail()) {
        PrintUsage();
        exit(-1);
      }
    } else if (tag == "-jlu") {
      std::stringstream ss(val);
      ss >> job_length_usec;
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
  std::cout << "Job spawner worker is up!" << std::endl;
  JobSpawnerApp * app = new JobSpawnerApp(counter,
                                          part_num,
                                          chunk_per_part,
                                          chunk_size,
                                          bandwidth,
                                          stage_num,
                                          job_length_usec);

  Worker * w = new Worker(scheduler_ip, scheduler_port, listening_port, app);
  if (ip_address_given) {
    w->set_ip_address(ip_address);
  }
  w->Run();
}

