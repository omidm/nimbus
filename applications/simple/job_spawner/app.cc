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
 * An Example application that spawns a lot of jobs in the 3d space.
 *
 * Author: Omid Mashayekhi<omidm@stanford.edu>
 */

#include <boost/program_options.hpp>
#include <vector>
#include "./app.h"
#include "./job.h"
#include "./data.h"

JobSpawnerApp::JobSpawnerApp(size_t counter,
                             size_t part_num,
                             size_t chunk_per_part,
                             size_t chunk_size,
                             size_t bandwidth,
                             size_t stage_num,
                             size_t job_length_usec) {
  counter_ = counter;
  part_num_ = part_num;
  chunk_per_part_ = chunk_per_part;
  chunk_size_ = chunk_size;
  bandwidth_ = bandwidth;
  stage_num_ = stage_num;
  job_length_usec_ = job_length_usec;
};

JobSpawnerApp::~JobSpawnerApp() {
};

void JobSpawnerApp::Load() {
  RegisterJob(NIMBUS_MAIN_JOB_NAME, new Main(this));
  RegisterJob(INIT_JOB_NAME, new Init());
  RegisterJob(LOOP_JOB_NAME, new ForLoop(this));
  RegisterJob(PRINT_JOB_NAME, new Print());
  RegisterJob(STAGE_JOB_NAME, new Stage(this));
  RegisterJob(CONNECTOR_JOB_NAME, new Connector(this));

  RegisterData(DATA_NAME, new Vec());

  // std::cout << "Finished Creating Data and Job Tables" << std::endl;
};


extern "C" Application * ApplicationBuilder(int argc, char *argv[]) {
  namespace po = boost::program_options;

  size_t iter_num = 30;
  size_t part_num = 100;
  size_t chunk_per_part = 1;
  size_t chunk_size = 50;
  size_t bandwidth = 10;
  size_t stage_num = 10;
  size_t job_length_usec = 0;

  po::options_description desc("Stencil 1D Options");
  desc.add_options()
    ("help,h", "produce help message")

    // Optinal arguments
    ("iter_num,i", po::value<size_t>(&iter_num), "number of iterations [default = 150]") //NOLINT
    ("pn", po::value<size_t>(&part_num), "number of partitions [default = 100]") //NOLINT
    ("cpp", po::value<size_t>(&chunk_per_part), "number of chunks per partition [default = 1]") //NOLINT
    ("cs", po::value<size_t>(&chunk_size), "chunk size in cells [default = 50]") //NOLINT
    ("bw", po::value<size_t>(&bandwidth), "ghost region bandwidth [default = 10]") //NOLINT
    ("sn", po::value<size_t>(&stage_num), "number of stages per each iteration block [default = 10]") //NOLINT
    ("jlu", po::value<size_t>(&bandwidth), "each job length in micro seconds [default = 0]"); //NOLINT

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
  }
  catch(std::exception& e) { // NOLINT
    std::cerr << "ERROR: " << e.what() << "\n";
    return NULL;
  }

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return NULL;
  }

  try {
    po::notify(vm);
  }
  catch(std::exception& e) { // NOLINT
    std::cerr << "ERROR: " << e.what() << "\n";
    return NULL;
  }

  return new JobSpawnerApp(iter_num,
                           part_num,
                           chunk_per_part,
                           chunk_size,
                           bandwidth,
                           stage_num,
                           job_length_usec);
}


