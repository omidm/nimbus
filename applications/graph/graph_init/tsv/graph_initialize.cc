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
 * This program loads a graph file, constructs partitions and logical regions
 * and saves the information to files. It uses GraphPartitioner for these
 * operations.
 */

#include <boost/program_options.hpp>
#include <iostream>  // NOLINT

#include "application/graph_library/graph_partitioner.h"

int main(int argc, char **argv) {
  nimbus::GraphPartitioner graph;

  size_t num_passes = 0;

  namespace po = boost::program_options;
  po::options_description desc("Options");
  desc.add_options()
    ("help", "produce help message")
    ("nodefile", po::value<std::string>()->required(), "node file name")
    ("edgefile", po::value<std::string>()->required(), "edge file name")
    ("partitions", po::value<int>()->required(), "number of partitions")
    ("passes", po::value<size_t>(&num_passes),
     "number of times to refine partitions")
    ("outdir", po::value<std::string>()->required(),
     "output partitions and logical regions to this directory")
    ;  // NOLINT

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }

  try {
    po::notify(vm);
  }
  catch(std::exception& e) { // NOLINT
    std::cout << "ERROR:" << e.what() << "\n";
    std::cout << desc << "\n";
    return 2;
  }

  graph.LoadFromTSV(vm["nodefile"].as<std::string>().c_str(),
                    vm["edgefile"].as<std::string>().c_str());
  graph.PartitionRandomEdgeCutAndRefine(vm["partitions"].as<int>(), num_passes);
  graph.DetermineLogicalObjects();
  graph.SaveGraphInEachPartition(vm["outdir"].as<std::string>().c_str());
  graph.SaveLogicalObjects(vm["outdir"].as<std::string>().c_str());
}
