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


int main(int argc, char *argv[]) {
  namespace po = boost::program_options;


  port_t listening_port;
  port_t controller_port;
  std::string controller_ip;
  std::string ip_address;
  
  WorkerManager::across_job_parallism = 1;

  uint64_t scale;
  uint64_t part_num_x;
  uint64_t part_num_y;
  uint64_t part_num_z;
  uint64_t projection_part_num_x;
  uint64_t projection_part_num_y;
  uint64_t projection_part_num_z;
  uint64_t last_frame;
  uint64_t max_iterations;
  uint64_t iteration_batch;
  uint64_t smart_projection_level;
  float water_level;


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
    ("pnx", po::value<uint64_t>(&part_num_x), "partition number along x") //NOLINT
    ("pny", po::value<uint64_t>(&part_num_y), "partition number along y") //NOLINT
    ("pnz", po::value<uint64_t>(&part_num_z), "partition number along z") //NOLINT
    ("ppnx", po::value<uint64_t>(&projection_part_num_x), "projection partition number along x") //NOLINT
    ("ppny", po::value<uint64_t>(&projection_part_num_y), "projection partition number along y") //NOLINT
    ("ppnz", po::value<uint64_t>(&projection_part_num_z), "projection partition number along z") //NOLINT
    ("last_frame,e", po::value<uint64_t>(&last_frame), "last frame to compute") //NOLINT
    ("maxi", po::value<uint64_t>(&max_iterations), "maximum projection iterations") //NOLINT
    ("ibatch", po::value<uint64_t>(&iteration_batch), "projection iteration batch") //NOLINT
    ("psl", po::value<uint64_t>(&smart_projection_level), "smart projection level") //NOLINT
    ("wl", po::value<float>(&water_level), "initial water level, float between 0 and 1") //NOLINT

    ("ithread", po::value<uint64_t>(&WorkerManager::inside_job_parallism), "number of threads within one job") //NOLINT
    ("othread", po::value<uint64_t>(&WorkerManager::across_job_parallism), "number of threads at worker for job execution") //NOLINT

    ("dpb", "deactivate projection bottleneck job")
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

  if (vm.count("scale")) {
    app->set_scale(scale);
  }

  if (vm.count("pnx")) {
    app->set_part_num_x(part_num_x);
  }

  if (vm.count("pny")) {
    app->set_part_num_y(part_num_y);
  }

  if (vm.count("pnz")) {
    app->set_part_num_z(part_num_z);
  }

  if (vm.count("ppnx")) {
    app->set_projection_part_num_x(projection_part_num_x);
  }

  if (vm.count("ppny")) {
    app->set_projection_part_num_y(projection_part_num_y);
  }

  if (vm.count("ppnz")) {
    app->set_projection_part_num_z(projection_part_num_z);
  }

  if (vm.count("last_frame")) {
    app->set_last_frame(last_frame);
  }

  if (vm.count("psl")) {
    app->set_smart_projection(smart_projection_level);
  }

  if (vm.count("wl")) {
    app->set_water_level(water_level);
  }

  if (vm.count("maxi")) {
    app->set_max_iterations(max_iterations);
  }

  if (vm.count("ibatch")) {
    app->set_iteration_batch(iteration_batch);
  }

  if (vm.count("dpb")) {
    app->set_spawn_projection_loop_bottleneck(false);
  }

  if (vm.count("dgw")) {
    app->set_global_write(false);
  }




  w->Run();
}
