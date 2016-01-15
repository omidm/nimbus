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
  * A Nimbus entry function for a k_means application worker. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include <boost/program_options.hpp>
#include <iostream>  // NOLINT
#include <sstream> // NOLINT
#include "shared/nimbus.h"
#include "../../application/k_means/app.h"
#include "worker/worker_manager.h"

#define DEFAULT_DIMENSION 10
#define DEFAULT_CLUSTER_NUM 2
#define DEFAULT_ITERATION_NUM 10
#define DEFAULT_PARTITION_NUM 10
#define DEFAULT_SAMPLE_NUM_M 1
using namespace nimbus; // NOLINT

int main(int argc, char *argv[]) {
  namespace po = boost::program_options;

  port_t listening_port;
  port_t controller_port;
  std::string controller_ip;
  std::string ip_address;

  size_t dimension;
  size_t cluster_num;
  size_t iteration_num;
  size_t partition_num;
  size_t sample_num_m;

  WorkerManager::across_job_parallism = 1;

  po::options_description desc("Options");
  desc.add_options()
    ("help,h", "produce help message")

    // Required arguments
    ("port,p", po::value<port_t>(&listening_port)->required(), "listening port for data exchanger") // NOLINT
    ("cport", po::value<port_t>(&controller_port)->required(), "controller listening port") // NOLINT
    ("cip", po::value<std::string>(&controller_ip)->required(), "controller ip address") // NOLINT

    // Optinal arguments
    ("ip", po::value<std::string>(&ip_address), "forced ip address of the worker, not known by controller") // NOLINT
    ("dimension,d", po::value<std::size_t>(&dimension)->default_value(DEFAULT_DIMENSION), "dimension of the sample vectors") // NOLINT
    ("iteration,i", po::value<std::size_t>(&iteration_num)->default_value(DEFAULT_ITERATION_NUM), "number of iterations") // NOLINT
    ("cn",
     po::value<std::size_t>(&cluster_num)->default_value(DEFAULT_CLUSTER_NUM), "number of clusters") // NOLINT
    ("pn", po::value<std::size_t>(&partition_num)->default_value(DEFAULT_PARTITION_NUM), "number of partitions") // NOLINT
    ("sn", po::value<std::size_t>(&sample_num_m)->default_value(DEFAULT_SAMPLE_NUM_M), "number of samples in Million") // NOLINT
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
  std::cout << "Simple Worker is up!" << std::endl;
  KMeans *app = new KMeans(dimension,
                           cluster_num,
                           iteration_num,
                           partition_num,
                           sample_num_m);

  Worker * w = new Worker(controller_ip, controller_port, listening_port, app);

  if (vm.count("ip")) {
    w->set_ip_address(ip_address);
  }

  if (vm.count("det")) {
    w->set_execution_template_active(false);
  }

  w->Run();
}

