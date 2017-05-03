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

ScaffoldApp::ScaffoldApp(std::string message,
                         size_t part_num) {
  message_ = message;
  part_num_ = part_num;
};

ScaffoldApp::~ScaffoldApp() {
};

void ScaffoldApp::Load() {
  RegisterJob(NIMBUS_MAIN_JOB_NAME, new Main(this));
  RegisterJob(PRINT_JOB_NAME, new Print());

  RegisterData(DATA_NAME, new Vec());
};


extern "C" Application * ApplicationBuilder(int argc, char *argv[]) {
  namespace po = boost::program_options;

  size_t part_num = 10;
  std::string message = "Hello world!";

  po::options_description desc("ScaffoldApp Options");
  desc.add_options()
    ("help,h", "produce help message")

    // Optinal arguments
    ("part_num,n", po::value<size_t>(&part_num), "number of partitions [default = 10]") //NOLINT
    ("message,m", po::value<std::string>(&message), "message printed by main job [default: Hello world!]"); //NOLINT

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

  return new ScaffoldApp(message,
                         part_num);
}


