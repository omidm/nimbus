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
 * Author: Chinmayee Shah
 */

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "applications/graph/page_rank/app.h"
#include "applications/graph/page_rank/edge_data.h"
#include "applications/graph/page_rank/node_data.h"
#include "applications/graph/page_rank/job.h"
#include "src/shared/nimbus.h"

// NOTE: The directory paths should be either absolute or relative to
// "nimbus_worker" executable located at "<nimbus-root>/nodes/nimbus_worker/"
#define DEFAULT_INPUT_DIR "../../applications/graph/page_rank/input-sample"
#define DEFAULT_OUTPUT_DIR "output"
#define DEFAULT_ITERATION_NUM 10

namespace nimbus {

PageRank::PageRank(std::string input_dir, std::string output_dir,
                   size_t num_iterations)
    : input_dir_(input_dir), output_dir_(output_dir),
      num_iterations_(num_iterations),
      graph_helper_(0), num_nodes_(0) {}

PageRank::~PageRank() {
  if (graph_helper_)
    delete graph_helper_;
}

#define SSTR(x) dynamic_cast< std::ostringstream & >( std::ostringstream() << std::dec << x ).str()  // NOLINT

void PageRank::Load() {
  // Register jobs here
  RegisterJob(NIMBUS_MAIN_JOB_NAME, new Main(this));
  RegisterJob(INIT_JOB, new Init(this));
  RegisterJob(FOR_LOOP_JOB, new ForLoop(this));
  RegisterJob(SCATTER_JOB, new Scatter(this));
  RegisterJob(GATHER_JOB, new Gather(this));
  RegisterJob(DUMP_JOB, new Dump(this));
  // Register data here
  RegisterData(EDGES, new EdgeData(EDGES));
  RegisterData(NODES, new NodeData(NODES));
  // read in graph information
  std::cout << "OMID: " << input_dir() << "\n";
  assert(boost::filesystem::is_directory(input_dir()));
  GraphLOs* graph = new GraphLOs();
  set_graph_helper(graph);
  graph->LoadGraphInfo(input_dir());
  std::ifstream num_node_file(SSTR(input_dir() << "/common/num_nodes").c_str());
  num_node_file >> num_nodes_;
  num_node_file.close();
}

std::string PageRank::input_dir() const {
  return input_dir_;
}

std::string PageRank::output_dir() const {
  return output_dir_;
}

size_t PageRank::num_iterations() const {
  return num_iterations_;
}

GraphLOs* PageRank::graph_helper() const {
  return graph_helper_;
}

size_t PageRank::num_nodes() const {
  return num_nodes_;
}

void PageRank::set_graph_helper(GraphLOs* graph_helper) {
  graph_helper_ = graph_helper;
}


extern "C" Application * ApplicationBuilder(int argc, char *argv[]) {
  namespace po = boost::program_options;

  std::string input_dir;
  std::string output_dir;
  size_t iteration_num;

  po::options_description desc("Page Rank Options");
  desc.add_options()
    ("help,h", "produce help message")
    ("input_dir,g", po::value<std::string>(&input_dir)->default_value(DEFAULT_INPUT_DIR), "input directory that holds partitioned graph (absolute path or relative to nimbus_worker executable located at <nimbus-root>/nodes/nimbus_worker)") // NOLINT
    ("output_dir,o", po::value<std::string>(&output_dir)->default_value(DEFAULT_OUTPUT_DIR), "output directory that saves the ranks (absolute path or relative to nimbus_worker executable located at <nimbus-root>/nodes/nimbus_worker)") // NOLINT
    ("iteration,i", po::value<std::size_t>(&iteration_num)->default_value(DEFAULT_ITERATION_NUM), "number of iterations"); // NOLINT

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

  PageRank *app = new PageRank(input_dir,
                               output_dir,
                               iteration_num);

  return app;
}


}  // namespace nimbus

