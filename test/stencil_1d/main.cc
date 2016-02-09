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

#include <boost/program_options.hpp>
#include <iostream>  // NOLINT
#include <sstream> // NOLINT
#include "src/shared/nimbus.h"
#include "src/worker/worker_manager.h"
#include "application/stencil_1d/app.h"

using namespace nimbus; // NOLINT

int main(int argc, char *argv[]) {
  namespace po = boost::program_options;

  port_t listening_port;
  port_t controller_port;
  std::string controller_ip;
  std::string ip_address;

  WorkerManager::across_job_parallism = 1;

  size_t iter_num = 150;
  size_t part_num = 4;
  size_t chunk_per_part = 4;
  size_t chunk_size = 5;
  size_t bandwidth = 1;

  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "produce help message")

    // Required arguments
    ("port,p", po::value<port_t>(&listening_port)->required(), "listening port for data exchanger") // NOLINT
    ("cport", po::value<port_t>(&controller_port)->required(), "controller listening port") // NOLINT
    ("cip", po::value<std::string>(&controller_ip)->required(), "controller ip address") // NOLINT

    // Optinal arguments
    ("ip", po::value<std::string>(&ip_address), "forced ip address of the worker, not known by controller") // NOLINT

    ("iter_num,i", po::value<size_t>(&iter_num), "number of iterations [default = 150]") //NOLINT
    ("pn", po::value<size_t>(&part_num), "number of partitions [default = 4]") //NOLINT
    ("cpp", po::value<size_t>(&chunk_per_part), "number of chunks per partition [default = 4]") //NOLINT
    ("cs", po::value<size_t>(&chunk_size), "chunk size in cells [default = 5]") //NOLINT
    ("bw", po::value<size_t>(&bandwidth), "ghost region bandwidth [default = 1]") //NOLINT

    ("othread", po::value<uint64_t>(&WorkerManager::across_job_parallism), "number of threads at worker for job execution") //NOLINT

    ("det", "deactivate execution template");

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
  std::cout << "Stencil1D Worker is up!" << std::endl;
  Stencil1DApp * app = new Stencil1DApp(
      iter_num, part_num, chunk_per_part, chunk_size, bandwidth);
  Worker * w = new Worker(controller_ip, controller_port, listening_port, app);

  if (vm.count("ip")) {
    w->set_ip_address(ip_address);
  }

  if (vm.count("det")) {
    w->set_execution_template_active(false);
  }

  w->Run();
}

