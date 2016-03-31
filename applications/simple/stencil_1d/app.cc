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

#define DEFAULT_ITERATION_NUM 150
#define DEFAULT_PARTITION_NUM 4
#define DEFAULT_CHUNK_PER_PARTITION 4
#define DEFAULT_CHUNK_SIZE 5
#define DEFAULT_BANDWIDTH 1
#define DEFAULT_SPIN_WAIT_US 0


Stencil1DApp::Stencil1DApp(const size_t& iteration_num,
                           const size_t& partition_num,
                           const size_t& chunk_per_partition,
                           const size_t& chunk_size,
                           const size_t& bandwidth,
                           const size_t& spin_wait_us) {
  iteration_num_ = iteration_num;
  partition_num_ = partition_num;
  chunk_per_partition_ = chunk_per_partition;
  chunk_size_ = chunk_size;
  bandwidth_ = bandwidth;
  spin_wait_us_ = spin_wait_us;
};

Stencil1DApp::~Stencil1DApp() {
};

void Stencil1DApp::Load() {
  std::cout << "Start Creating Data and Job Tables" << std::endl;

  RegisterJob(NIMBUS_MAIN_JOB_NAME, new Main(this));
  RegisterJob(INIT_JOB_NAME, new Init());
  RegisterJob(LOOP_JOB_NAME, new ForLoop(this));
  RegisterJob(PRINT_JOB_NAME, new Print());
  RegisterJob(STENCIL_JOB_NAME, new Stencil(this));

  RegisterData(DATA_NAME, new Vec());

  std::cout << "Finished Creating Data and Job Tables" << std::endl;
};


extern "C" Application * ApplicationBuilder(int argc, char *argv[]) {
  namespace po = boost::program_options;

  size_t iteration_num;
  size_t partition_num;
  size_t chunk_per_partition;
  size_t chunk_size;
  size_t bandwidth;
  size_t spin_wait_us;

  po::options_description desc("Stencil 1D Options");
  desc.add_options()
    ("help,h", "produce help message")

    // Optinal arguments
    ("iteration,i", po::value<std::size_t>(&iteration_num)->default_value(DEFAULT_ITERATION_NUM), "number of iterations") // NOLINT
    ("pn", po::value<size_t>(&partition_num)->default_value(DEFAULT_PARTITION_NUM), "number of partitions") //NOLINT
    ("cpp", po::value<size_t>(&chunk_per_partition)->default_value(DEFAULT_CHUNK_PER_PARTITION), "number of chunks per partition") //NOLINT
    ("cs", po::value<size_t>(&chunk_size)->default_value(DEFAULT_CHUNK_SIZE), "chunk size in cells") //NOLINT
    ("bw", po::value<size_t>(&bandwidth)->default_value(DEFAULT_BANDWIDTH), "ghost region bandwidth [default = 1]") //NOLINT
    ("spin_wait,w", po::value<std::size_t>(&spin_wait_us)->default_value(DEFAULT_SPIN_WAIT_US), "spin wait in micro seconds, if non zero,replaces the gradient operation with fixed spin wait."); // NOLINT

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

  return new Stencil1DApp(iteration_num,
                          partition_num,
                          chunk_per_partition,
                          chunk_size,
                          bandwidth,
                          spin_wait_us);
}

