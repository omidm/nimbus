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
 * GraphPartitioner constructs a graph in-memory from input files, partitions
 * the graph (edge-cut), dtermines disjoint logical regions for vertices and
 * edges (assuming that for a vertex program,  writes happen only to outgoing
 * edges and reads happen from incoming edges, apart from read/ wrte on
 * vertices).
 * Currently, graph construction from TSV files is supported. However, there is
 * no error handling, or handling of special cases like comments and blank
 * lines.
 * Author: Chinmayee Shah
 */

#include <boost/filesystem.hpp>
#include <cassert>
#include <cstdlib>
#include <fstream>  // NOLINT
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "application/graph_library/graph_partitioner.h"

namespace nimbus {

GraphPartitioner::GraphPartitioner() : num_nodes_(0), node_degree_(0),
                                       num_edges_(0), edge_src_(0), edge_dst_(0),
                                       num_partitions_(0), node_partitions_(0),
                                       num_edge_logical_objects_(0),
                                       node_logical_objects_(0),
                                       edge_logical_objects_(0),
                                       edge_logical_object_exists_(0),
                                       edge_lo_write_(0),
                                       edge_lo_read_(0) { }

GraphPartitioner::~GraphPartitioner() {
  if (edge_src_) delete edge_src_;
  if (edge_dst_) delete edge_dst_;
  if (node_degree_) delete node_degree_;
  num_nodes_ = 0;
  num_edges_ = 0;
  if (node_partitions_) delete node_partitions_;

  if (node_logical_objects_) {
    for (size_t p = 0; p < num_partitions_; ++p) {
      if (node_logical_objects_[p]) {
        delete node_logical_objects_[p];
      }
    }
    delete node_logical_objects_;
  }

  if (edge_logical_objects_) {
    for (size_t pp = 0; pp < num_partitions_ * num_partitions_; ++pp) {
      if (edge_logical_objects_[pp]) {
        delete edge_logical_objects_[pp];
      }
    }
    delete edge_logical_objects_;
  }
  num_partitions_ = 0;
  if (edge_logical_object_exists_) {
    delete edge_logical_object_exists_;
  }

  if (edge_lo_write_) {
    for (size_t p = 0; p < num_partitions_; ++p) {
      if (edge_lo_write_[p])
        delete edge_lo_write_[p];
    }
    delete edge_lo_write_;
  }
  if (edge_lo_read_) {
    for (size_t p = 0; p < num_partitions_; ++p) {
      if (edge_lo_read_[p])
        delete edge_lo_read_[p];
    }
    delete edge_lo_read_;
  }
  num_edge_logical_objects_ = 0;
}

/*
 * Constructs a graph from a node and an edge TSV file.
 * Implementation assumes no comments, and no additional (blank) new lines.
 * Implementation does not include any error/ special cases handling.
 */
void GraphPartitioner::LoadFromTSV(std::string node_file_name,
                                   std::string edge_file_name) {
  std::ifstream node_file, edge_file;
  if (!(boost::filesystem::is_regular_file(node_file_name) &&
        boost::filesystem::is_regular_file(edge_file_name))) {
    std::cout << "Invalid input node file or edge file\n";
    exit(-1);
  }
  node_file.open(node_file_name.c_str());
  edge_file.open(edge_file_name.c_str());

  // determine number of nodes and edges
  std::string line;
  num_nodes_ = 0;
  while (std::getline(node_file, line)) {
    num_nodes_++;
  }
  num_edges_ = 0;
  while (std::getline(edge_file, line)) {
    num_edges_++;
  }

  node_file.close();
  edge_file.close();

  node_file.open(node_file_name.c_str());
  edge_file.open(edge_file_name.c_str());

  // make a map from string names to node ids
  std::map<std::string, size_t> node_map;
  node_degree_ = new size_t[num_nodes_];
  size_t node_id = 0;
  std::string node_name;
  while (node_file >> node_name) {
    node_map[node_name] = node_id;
    node_degree_[node_id] = 0;
    node_id++;
  }
  assert(node_id == num_nodes_);

  // determine src node id and dst node id for each edge
  edge_src_ = new size_t[num_edges_];
  edge_dst_ = new size_t[num_edges_];
  size_t edge_id = 0;
  std::string src_name, dst_name;
  while (edge_file >> src_name && edge_file >> dst_name) {
    edge_src_[edge_id] = node_map[src_name];
    edge_dst_[edge_id] = node_map[dst_name];
    node_degree_[node_map[src_name]]++;
    edge_id++;
  }
  assert(edge_id == num_edges_);
  std::cout << "Loaded graph with " << num_nodes_ << " nodes and " <<
    num_edges_ << " edges\n";
}

/*
 * Creates a random edge-cut partition by placing vertices in random partitions.
 */
void GraphPartitioner::PartitionRandomEdgeCut(size_t num_partitions) {
  num_partitions_ = num_partitions;
  srand(0);
  if (!node_partitions_)
    node_partitions_ = new size_t[num_nodes_];
  for (size_t n = 0; n < num_nodes_; ++n) {
    node_partitions_[n] = rand() % num_partitions;  // NOLINT
  }
}

/*
 * Determines disjoint logical objects on edges and vertices. Also determines
 * the set of logical objects in read and write sets for each partition.
 */
void GraphPartitioner::DetermineLogicalObjects() {
  // determine disjoint logical objects for nodes
  node_logical_objects_ = new std::vector<size_t>* [num_partitions_];
  for (size_t p = 0; p < num_partitions_; ++p)
    node_logical_objects_[p] = new std::vector<size_t>;
  for (size_t n = 0; n < num_nodes_; ++n) {
    size_t id = node_partitions_[n];
    node_logical_objects_[id]->push_back(n);
  }
  // determine disjoint logical objects for edges
  edge_logical_objects_ =
    new std::vector<size_t>* [num_partitions_ * num_partitions_];
  edge_logical_object_exists_ = new std::vector<bool>(
      num_partitions_*num_partitions_, false);
  edge_lo_write_ = new std::set<size_t>* [num_partitions_];
  edge_lo_read_ = new std::set<size_t>* [num_partitions_];
  for (size_t pp = 0; pp < num_partitions_ * num_partitions_; ++pp) {
    edge_logical_objects_[pp] = new std::vector<size_t>;
  }
  for (size_t p = 0; p < num_partitions_; ++p) {
    edge_lo_write_[p] = new std::set<size_t>;
    edge_lo_read_[p] = new std::set<size_t>;
  }
  for (size_t e = 0; e < num_edges_; ++e) {
    size_t spid = node_partitions_[edge_src_[e]];
    size_t dpid = node_partitions_[edge_dst_[e]];
    size_t lid = spid * num_partitions_ + dpid;
    edge_logical_objects_[lid]->push_back(e);
    edge_lo_write_[spid]->insert(lid);
    edge_lo_read_[dpid]->insert(lid);
    if (!(*edge_logical_object_exists_)[lid]) {
      (*edge_logical_object_exists_)[lid] = true;
      num_edge_logical_objects_++;
    }
  }
}


#define SSTR(x) dynamic_cast< std::ostringstream & >( std::ostringstream() << std::dec << x ).str()  // NOLINT

/*
 * Outputs 2 files per partition: a nodes file containing ids of nodes in the
 * partition and there degree (number of edges going out),
 * and an edges file containing ids of edges, and source node and destination
 * node ids.
 * This file has detailed information about graph structure (edge id, src and
 * dst id and node degree).
 */
void GraphPartitioner::SaveGraphInEachPartition(std::string dir_name) const {
  boost::filesystem::create_directories(dir_name);
  assert(boost::filesystem::is_directory(dir_name));
  for (size_t p = 0; p < num_partitions_; ++p) {
    std::cout << "Saving graph in partition " << p << "\n";
    boost::filesystem::create_directories(SSTR(dir_name << "/" << p));
    std::string node_fname  = SSTR(dir_name << "/" << p << "/nodes");
    std::string edge_fname  = SSTR(dir_name << "/" << p << "/edges");
    std::ofstream node_file, edge_file;
    node_file.open(node_fname.c_str());
    edge_file.open(edge_fname.c_str());
    for (size_t n = 0; n < num_nodes_; n++) {
      if (node_partitions_[n] == p)
        node_file << n << node_degree_[n] << "\n";
    }
    for (size_t e = 0; e < num_edges_; e++) {
      if (node_partitions_[edge_src_[e]] == p ||
          node_partitions_[edge_dst_[e]] == p)
        edge_file << e << " " << edge_src_[e] << " " << edge_dst_[e] << "\n";
    }
    node_file.close();
    edge_file.close();
  }
}

/*
 * Outputs number of node partitions and edge logical objects in directory
 * common, and a list of logical objects that each partition reads/ writes.
 * Additionally, outputs 3 files per partition:
 *  - a node_lo file containing partition id, size, and node ids in the
 *  node logical object for the partition
 *  (should be the similar to nodes file in SaveGraphInEachPartition)
 *  - an edge_read file containing ids of edge logical objects the partition
 *  reads, each logical object id followed by ids of edges it contains
 *  - an edge_write file containing ids of edge logical objects the partition
 *  writes, each logical object id followed by ids of edges it contains
 */
void GraphPartitioner::SaveLogicalObjects(std::string dir_name) const {
  boost::filesystem::create_directories(dir_name);
  assert(boost::filesystem::is_directory(dir_name));
  {
    boost::filesystem::create_directories(SSTR(dir_name << "/common"));
    // partitions
    std::string num_partitions_fname = SSTR(dir_name << "/common/num_partitions");
    std::ofstream num_partitions_file;
    num_partitions_file.open(num_partitions_fname.c_str());
    num_partitions_file << num_partitions_;
    num_partitions_file.close();
    // number of edge logical objects
    std::string num_edge_los_fname = SSTR(dir_name << "/common/num_edge_los");
    std::ofstream num_edge_los_file;
    num_edge_los_file.open(num_edge_los_fname.c_str());
    num_edge_los_file << num_edge_logical_objects_;
    num_edge_los_file.close();
    // edge read set
    std::string edge_rs_fname = SSTR(dir_name << "/common/edge_read_sets");
    std::ofstream edge_rs_file;
    edge_rs_file.open(edge_rs_fname.c_str());
    for (size_t p = 0; p < num_partitions_; ++p) {
      std::set<size_t> *edge_read = edge_lo_read_[p];
      edge_rs_file << p << " : " << edge_read->size() << " : ";
      for (std::set<size_t>::iterator e = edge_read->begin();
           e != edge_read->end(); ++e) {
        edge_rs_file << (*e) << " ";
      }
      edge_rs_file << "\n";
    }
    // edge write set
    std::string edge_ws_fname = SSTR(dir_name << "/common/edge_write_sets");
    std::ofstream edge_ws_file;
    edge_ws_file.open(edge_ws_fname.c_str());
    for (size_t p = 0; p < num_partitions_; ++p) {
      std::set<size_t> *edge_write = edge_lo_write_[p];
      edge_ws_file << p << " : " << edge_write->size() << " : ";
      for (std::set<size_t>::iterator e = edge_write->begin();
           e != edge_write->end(); ++e) {
        edge_ws_file << (*e) << " ";
      }
      edge_ws_file << "\n";
    }
  }
  for (size_t p = 0; p < num_partitions_; ++p) {
    std::cout << "Saving logical objects to partition " << p << "\n";
    boost::filesystem::create_directories(SSTR(dir_name << "/" << p));
    std::string node_lo_fname    = SSTR(dir_name << "/" << p << "/node_lo");
    std::string edge_read_fname  = SSTR(dir_name << "/" << p << "/edge_read_sets");
    std::string edge_write_fname = SSTR(dir_name << "/" << p << "/edge_write_sets");
    std::ofstream edge_read_file, edge_write_file, node_lo_file, edge_lo_file;
    // node logical object contents
    node_lo_file.open(node_lo_fname.c_str());
    std::vector<size_t> *node_lo = node_logical_objects_[p];
    node_lo_file << p << " : " << node_lo->size() << " : ";
    for (size_t n = 0; n < node_lo->size(); ++n) {
      node_lo_file << node_lo->at(n) << " ";
    }
    node_lo_file.close();
    // edge read set ids and contents
    edge_read_file.open(edge_read_fname.c_str());
    std::set<size_t> *edge_lor = edge_lo_read_[p];
    for (std::set<size_t>::iterator l = edge_lor->begin();
         l != edge_lor->end(); ++l) {
      std::vector<size_t> *edge_lo = edge_logical_objects_[(*l)];
      edge_read_file << (*l) << " : " << edge_lo->size() << " : ";
      for (size_t e = 0; e < edge_lo->size(); ++e) {
        edge_read_file << edge_lo->at(e) << " ";
      }
      edge_read_file << "\n";
    }
    edge_read_file.close();
    // edge write set ids and contents
    edge_write_file.open(edge_write_fname.c_str());
    std::set<size_t> *edge_low = edge_lo_write_[p];
    for (std::set<size_t>::iterator l = edge_low->begin();
         l != edge_low->end(); ++l) {
      std::vector<size_t> *edge_lo = edge_logical_objects_[(*l)];
      edge_write_file << (*l) << " : " << edge_lo->size() << ":";
      for (size_t e = 0; e < edge_lo->size(); ++e) {
        edge_write_file << edge_lo->at(e) << " ";
      }
      edge_write_file << "\n";
    }
    edge_write_file.close();
  }
}

}  // namespace nimbus
