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
 * GraphLOs reads in number of partitions, number of edge logical objects, and
 * read and write set information for all partitions, and provides helper
 * methods to register all the logical objects and read/ write sets for each
 * partition to spawn jobs.
 * Author: Chinmayee Shah
 */

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <cassert>
#include <cstdlib>
#include <fstream>  // NOLINT
#include <sstream>
#include <string>
#include <vector>
#include "application/graph_library/graph_los.h"
#include "shared/nimbus.h"
#include "worker/job.h"

namespace nimbus {

#define SSTR(x) dynamic_cast< std::ostringstream & >( std::ostringstream() << std::dec << x ).str()  // NOLINT

GraphLOs::GraphLOs() {
}

GraphLOs::~GraphLOs() {
}

/*
 * Reads number of partitions and logical objects, and a list of logical
 * objects each partition reads/ writes for nodes/ edges.
 */
void GraphLOs::LoadGraphInfo(std::string dir_name) {
  assert(boost::filesystem::is_directory(dir_name));
  std::ifstream num_partitions_file, num_edge_los_file,
                edge_rs_file, edge_ws_file;
  // number of partitions
  num_partitions_file.open(SSTR(dir_name << "/common/num_partitions").c_str());  // NOLINT
  num_partitions_file >> num_partitions_;
  edge_lo_read_  = new std::vector<size_t>* [num_partitions_];
  edge_lo_write_ = new std::vector<size_t>* [num_partitions_];
  num_partitions_file.close();
  // total number of edge logical objects
  num_edge_los_file.open(SSTR(dir_name << "/common/num_edge_los").c_str());
  num_edge_los_file >> num_edge_los_;
  num_edge_los_file.close();
  size_t num_lines;
  std::string line;
  // edge read sets
  edge_rs_file.open(SSTR(dir_name << "/common/edge_read_sets").c_str());
  num_lines = 0;
  while (std::getline(edge_rs_file, line)) {
    std::vector<std::string> tokens;
    boost::algorithm::split(tokens, line,
                            boost::is_any_of(": "),
                            boost::token_compress_on);
    size_t partition = boost::lexical_cast<size_t>(tokens[0]);
    assert(partition == num_lines);
    size_t rs_size = boost::lexical_cast<size_t>(tokens[1]);
    assert(rs_size + 2 == tokens.size());
    std::vector<size_t> *edge_lo = new std::vector<size_t>(rs_size);
    edge_lo_read_[partition] = edge_lo;
    for (size_t i = 0; i < rs_size; ++i) {
      (*edge_lo)[i] = boost::lexical_cast<size_t>(tokens[i+2]);
    }
    num_lines++;
  }
  edge_rs_file.close();
  // edge write sets
  edge_ws_file.open(SSTR(dir_name << "/common/edge_write_sets").c_str());
  num_lines = 0;
  while (std::getline(edge_ws_file, line)) {
    std::vector<std::string> tokens;
    boost::algorithm::split(tokens, line,
                            boost::is_any_of(": "),
                            boost::token_compress_on);
    size_t partition = boost::lexical_cast<size_t>(tokens[0]);
    assert(partition == num_lines);
    size_t ws_size = boost::lexical_cast<size_t>(tokens[1]);
    assert(ws_size + 2 == tokens.size());
    std::vector<size_t> *edge_lo = new std::vector<size_t>(ws_size);
    edge_lo_write_[partition] = edge_lo;
    for (size_t i = 0; i < ws_size; ++i) {
      (*edge_lo)[i] = boost::lexical_cast<size_t>(tokens[i+2]);
    }
    num_lines++;
  }
  edge_ws_file.close();
}

/*
 * Define logical objects named name over edges.
 */
void GraphLOs::DefineEdgeLogicalObjects(Job *job, std::string name) {
  std::vector<logical_data_id_t> ids;
  job->GetNewLogicalDataID(&ids, num_edge_los_);
  lolist **lo_map_read_var  = new lolist* [num_partitions_];
  lolist **lo_map_write_var = new lolist* [num_partitions_];
  lo_map_read_[name]  = lo_map_read_var;
  lo_map_write_[name] = lo_map_write_var;
  IDSet<partition_id_t> neighbor;
  for (size_t p = 0; p < num_partitions_; ++p) {
    std::vector<size_t> *edge_los = edge_lo_write_[p];
    lolist* los = new lolist;
    lo_map_write_var[p] = los;
    for (size_t l = 0; l < edge_los->size(); ++l) {
      size_t eid = edge_los->at(l);
      {
        // hack for identifying data objects
        GeometricRegion region(p, 0, 0, 1, 1, 1);
        ID<partition_id_t> partition(p);
        job->DefinePartition(partition, region);
      }
      job->DefineData(name, ids[eid], p, neighbor);
      los->insert(ids[eid]);
    }
  }
  for (size_t p = 0; p < num_partitions_; ++p) {
    std::vector<size_t> *edge_los = edge_lo_read_[p];
    lolist* los = new lolist;
    lo_map_read_var[p] = los;
    for (size_t l = 0; l < edge_los->size(); ++l) {
      size_t eid = edge_los->at(l);
      los->insert(ids[eid]);
    }
  }
}

/*
 * Define logical objects named name over nodes.
 */
void GraphLOs::DefineNodeLogicalObjects(Job *job, std::string name) {
  std::vector<logical_data_id_t> ids;
  job->GetNewLogicalDataID(&ids, num_partitions_);
  lolist **lo_map_read_var  = new lolist* [num_partitions_];
  lolist **lo_map_write_var = new lolist* [num_partitions_];
  lo_map_read_[name]  = lo_map_read_var;
  lo_map_write_[name] = lo_map_write_var;
  IDSet<partition_id_t> neighbor;
  for (size_t p = 0; p < num_partitions_; ++p) {
    job->DefineData(name, ids[p], p, neighbor);
    lo_map_read_var[p]  = new lolist();
    lo_map_read_var[p]->insert(ids[p]);
    lo_map_write_var[p] = new lolist();
    lo_map_write_var[p]->insert(ids[p]);
  }
}

/*
 * Get read set.
 */
const IDSet<logical_data_id_t>* GraphLOs::
GetReadSet(std::string name, partition_id_t p) const {
  return lo_map_read_.at(name)[p];
}

/*
 * Get write set.
 */
const IDSet<logical_data_id_t>* GraphLOs::
GetWriteSet(std::string name, partition_id_t p) const {
  return lo_map_write_.at(name)[p];
}

/*
 * Get number of partitions.
 */
size_t GraphLOs::num_partitions() const {
  return num_partitions_;
}

}  // namespace nimbus
