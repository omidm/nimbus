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
 * Application is implemented here.
 *
 */

#include <boost/program_options.hpp>
#include <vector>
#include "./app.h"
#include "./job.h"
#include "./data.h"

nimbus::PartitionHandler ph;
AppDataVec AppDataVecPrototype;

Heat::Heat(size_t iter_num,
           int nx, int ny, int nz,
           int pnx, int pny, int pnz) {
  iter_num_ = iter_num;
  nx_ = nx;
  ny_ = ny;
  nz_ = nz;
  pnx_ = pnx;
  pny_ = pny;
  pnz_ = pnz;
};

Heat::~Heat() {
};

void Heat::Load() {
  RegisterJob(NIMBUS_MAIN_JOB_NAME, new Main(this));
  RegisterJob(LOOP_JOB_NAME, new Loop(this));
  RegisterJob(STENCIL_JOB_NAME, new Stencil(this));

  RegisterData(DATA_NAME, new Vec());

  ph.AddPartitions("kRegions",
                   nx_, ny_, nz_,
                   1, 1, 1, // bandwidth
                   pnx_, pny_, pnz_,
                   nimbus::PartitionHandler::INNER);

  ph.AddPartitions("kRegions",
                   nx_, ny_, nz_,
                   1, 1, 1, // bandwidth
                   pnx_, pny_, pnz_,
                   nimbus::PartitionHandler::OUTER);
};


extern "C" Application * ApplicationBuilder(int argc, char *argv[]) {
  namespace po = boost::program_options;

  size_t iter_num = 30;
  int nx = 10;
  int ny = 10;
  int nz = 10;
  int pnx = 2;
  int pny = 2;
  int pnz = 2;

  po::options_description desc("Heat Options");
  desc.add_options()
    ("help,h", "produce help message")

    // Optinal arguments
    ("iter_num,i", po::value<size_t>(&iter_num), "number of iterations [default = 30]") //NOLINT
    ("nx,x", po::value<int>(&nx), "size of x axis [default = 10]") //NOLINT
    ("ny,y", po::value<int>(&ny), "size of y axis [default = 10]") //NOLINT
    ("nz,z", po::value<int>(&nz), "size of z axis [default = 10]") //NOLINT
    ("pnx", po::value<int>(&pnx), "partition number along x axis [default = 2]") //NOLINT
    ("pny", po::value<int>(&pny), "partition number along y axis [default = 2]") //NOLINT
    ("pnz", po::value<int>(&pnz), "partition number along z axis [default = 2]"); //NOLINT

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

  return new Heat(iter_num,
                  nx, ny, nz,
                  pnx, pny, pnz);
}


