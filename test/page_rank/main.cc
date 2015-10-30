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
  * Nimbus worker for pagerank application.
  * Author: Chinmayee Shah
  */

#include <boost/program_options.hpp>
#include <iostream>     // NOLINT
#include <sstream>      // NOLINT
#include "shared/nimbus.h"
#include "application/page_rank/app.h"
#include "worker/worker_manager.h"

using namespace nimbus; // NOLINT

int main(int argc, char *argv[]) {
  namespace po = boost::program_options;

  // common arguments
  port_t listening_port;
  port_t controller_port;
  std::string controller_ip;
  std::string ip_address;

  // app arguments
  std::string input_dir;
  std::string output_dir;
  size_t num_iterations;

  WorkerManager::across_job_parallism = 1;

  po::options_description desc("Options");
  desc.add_options()
    ("help", "produce help message")

    // Required arguments
    ("port", po::value<port_t>(&listening_port)->required(), "listening port for data exchanger")  // NOLINT
    ("cport", po::value<port_t>(&controller_port)->required(), "controller listening port")        // NOLINT
    ("cip", po::value<std::string>(&controller_ip)->required(), "controller ip address")           // NOLINT

    // Optional arguments
    ("ip", po::value<std::string>(&ip_address), "(optional) forced ip address of the worker, not known by controller")                  // NOLINT
    ("othread", po::value<uint64_t>(&WorkerManager::across_job_parallism), "(optional) number of threads at worker for job execution")  //NOLINT

    // Required app arguments
    ("input", po::value<std::string>(&input_dir)->required(), "input directory containing graph") // NOLINT
    ("output", po::value<std::string>(&output_dir)->required(), "directory to save results") // NOLINT
    ("iterations", po::value<std::size_t>(&num_iterations)->required(), "number of iterations to perform + 1 (outputs saved in final iteration)")   // NOLINT
  ;  // NOLINT

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
    std::cout << desc << "\n";
    return -1;
  }

  nimbus_initialize();
  std::cout << "Simple Worker is up!" << std::endl;
  PageRank *app = new PageRank(input_dir, output_dir,
                               num_iterations);
  Worker *w = new Worker(controller_ip, controller_port, listening_port, app);

  if (vm.count("ip")) {
    w->set_ip_address(ip_address);
  }
  w->Run();
  std::cout << "Worker is done!" << std::endl;
}
