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

#include <boost/program_options.hpp>
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
  namespace po = boost::program_options;


  port_t listening_port;
  port_t controller_port;
  std::string controller_ip;
  std::string ip_address;
  
  WorkerManager::across_job_parallism = 1;

  uint64_t scale;  //

  uint64_t part_num_x;  //
  uint64_t part_num_y;  //
  uint64_t part_num_z;  //

  uint64_t projection_part_num_x;  //
  uint64_t projection_part_num_y;  //
  uint64_t projection_part_num_z;  //

  uint64_t last_frame;  //

  uint64_t max_projection_iterations;  //

  bool global_write;


  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "produce help message")

    // Required arguments
    ("port,p", po::value<port_t>(&listening_port)->required(), "listening port for data exchanger") // NOLINT
    ("cport", po::value<port_t>(&controller_port)->required(), "controller listening port") // NOLINT
    ("cip", po::value<std::string>(&controller_ip)->required(), "controller ip address") // NOLINT

    // Optinal arguments
    ("ip", po::value<std::string>(&ip_address), "forced ip address of the worker, not known by controller") // NOLINT

    ("scale,s", po::value<uint64_t>(&scale), "scale of the simulation") //NOLINT

    ("ithread", po::value<uint64_t>(&WorkerManager::inside_job_parallism), "number of threads within one job") //NOLINT
    ("othread", po::value<uint64_t>(&WorkerManager::across_job_parallism), "number of threads at worker for job execution") //NOLINT

    ("dgw", "deactivate one global write per frame");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 0;
  }

  try {
    po::notify(vm);
  }
  catch(std::exception& e) { // NOLINT
    std::cerr << "ERROR: " << e.what() << "\n";
    return -1;
  }


  nimbus_initialize();
  std::cout << "Simple Worker is up!" << std::endl;
  application::WaterApp *app = new application::WaterApp();
  SimpleWorker * w = new SimpleWorker(controller_ip, controller_port, listening_port, app);

  if (vm.count("ip")) {
    w->set_ip_address(ip_address);
  }

  if (vm.count("dgw")) {
    app->set_global_write(false);
  }

  if (vm.count("scale")) {
    app->set_scale(scale);
  }




  w->Run();
}
