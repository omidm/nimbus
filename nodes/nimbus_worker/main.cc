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
#include <dlfcn.h>
#include <boost/program_options.hpp>
#include <iostream>  // NOLINT
#include <sstream> // NOLINT
#include "src/shared/nimbus.h"
#include "src/worker/worker_manager.h"

using namespace nimbus; // NOLINT

int main(int argc, char *argv[]) {
  nimbus_initialize();

  namespace po = boost::program_options;

  port_t listening_port;
  port_t controller_port;
  std::string controller_ip;
  std::string application_lib;
  std::string ip_address;

  WorkerManager::across_job_parallism = 1;

  po::options_description desc("Nimbus Worker Options");
  desc.add_options()
    ("help,h", "produce help message")

    // Required arguments
    ("port,p", po::value<port_t>(&listening_port)->required(), "listening port for data exchanger") // NOLINT
    ("cport", po::value<port_t>(&controller_port)->required(), "controller listening port") // NOLINT
    ("cip", po::value<std::string>(&controller_ip)->required(), "controller ip address") // NOLINT
    ("app_lib,l", po::value<std::string>(&application_lib)->required(), "application library") // NOLINT

    // Optinal arguments
    ("ip", po::value<std::string>(&ip_address), "forced ip address of the worker, not known by controller") // NOLINT
    ("othread", po::value<uint64_t>(&WorkerManager::across_job_parallism), "number of threads at worker for job execution") //NOLINT
    ("det", "deactivate execution template");

  int argc_worker = 0;
  while (argc_worker < argc) {
    if (strcmp(argv[argc_worker], "-l") == 0) {
      break;
    }
    argc_worker++;
  }
  argc_worker += 2;
  if (argc_worker > argc) {
    argc_worker = argc;
  }

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc_worker, argv, desc), vm);
  }
  catch(std::exception& e) { // NOLINT
    std::cerr << "ERROR: " << e.what() << "\n";
    return -1;
  }


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


  void * handle = dlopen(application_lib.c_str(), RTLD_LAZY);
  if (handle == NULL) {
    dbg(DBG_ERROR, "ERROR, dlerror: %s\n", dlerror());
    return -1;
  }

  app_builder_t app_builder = (app_builder_t) dlsym(handle, "ApplicationBuilder"); // NOLINT
  if (app_builder == NULL) {
    dbg(DBG_ERROR, "ERROR, dlsym: %s\n", dlerror());
    return -1;
  }

  assert(argc_worker > 1);
  int argc_app = argc - argc_worker + 1;
  Application *app = (*app_builder)(argc_app, argv + argc_worker - 1);
  if (app == NULL) {
    dbg(DBG_ERROR, "ERROR, ApplicationBuilder returned NULL\n");
    return -1;
  }

  Worker * w = new Worker(controller_ip, controller_port, listening_port, app);

  if (vm.count("ip")) {
    w->set_ip_address(ip_address);
  }

  if (vm.count("det")) {
    w->set_execution_template_active(false);
  }

  w->Run();
}

