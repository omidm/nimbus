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

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/unordered_map.hpp>
#include <cassert>
#include <cstdlib>
#include <fstream>  // NOLINT
#include <sstream>
#include <string>
#include <vector>
#include "applications/graph/graph_library/graph_los.h"
#include "applications/graph/page_rank/app.h"
#include "applications/graph/page_rank/edge_data.h"
#include "applications/graph/page_rank/node_data.h"
#include "applications/graph/page_rank/job.h"
#include "applications/graph/page_rank/protobuf_compiled/parameter_msg.pb.h"
#include "src/shared/nimbus.h"

#define ALPHA 0.15

namespace nimbus {

void LoadParameter(Parameter* parameter, uint64_t *value) {
  std::string params_str(parameter->ser_data().data_ptr_raw(),
                         parameter->ser_data().size());
  parameter_msg::ParameterMsg msg;
  msg.ParseFromString(params_str);
  *value = msg.elem();
}

void SerializeParameter(Parameter* parameter, size_t value) {
  std::string params_str;
  parameter_msg::ParameterMsg msg;
  msg.set_elem(value);
  msg.SerializeToString(&params_str);
  parameter->set_ser_data(SerializedData(params_str));
}

void GetReadData(const Job &job,
                 const std::string &name,
                 const DataArray &da,
                 DataArray *read) {
  IDSet<physical_data_id_t> read_set = job.read_set();
  size_t rs = read_set.size();
  for (size_t i = 0; i < rs; ++i) {
    Data *d = da[i];
    if (d->name() == name &&
        read_set.contains(d->physical_id())) {
      read->push_back(d);
    }
  }
}

void GetWriteData(const Job &job,
                  const std::string &name,
                  const DataArray &da,
                  DataArray *write) {
  IDSet<physical_data_id_t> read_set = job.read_set();
  IDSet<physical_data_id_t> write_set = job.write_set();
  size_t rs = read_set.size();
  size_t ws = write_set.size();
  assert(rs+ws == da.size());
  for (size_t i = rs; i < rs + ws; ++i) {
    Data *d = da[i];
    if (d->name() == name &&
        write_set.contains(d->physical_id())) {
      write->push_back(d);
    }
  }
}

Main::Main(Application* app) {
  set_application(app);
}

Job* Main::Clone() {
  return new Main(application());
}

void Main::Execute(Parameter params, const DataArray& da) {
  // application, graph_los to help with data objects and jobs
  PageRank* page_rank = static_cast<PageRank*>(application());
  GraphLOs* graph = page_rank->graph_helper();
  size_t num_partitions = graph->num_partitions();

  // define data objects
  graph->DefineEdgeLogicalObjects(this, EDGES);
  graph->DefineNodeLogicalObjects(this, NODES);

  // initialize data -- need to initialize only ranks
  {
    Parameter init_params;
    std::vector<job_id_t> init_job_ids;
    GetNewJobID(&init_job_ids, num_partitions);
    for (size_t p = 0; p < num_partitions; ++p) {
      IDSet<logical_data_id_t> init_rs;
      const IDSet<logical_data_id_t>* nodes_ws = graph->GetWriteSet(NODES, p);
      const IDSet<logical_data_id_t>* edges_ws = graph->GetWriteSet(EDGES, p);
      IDSet<logical_data_id_t> init_ws(*nodes_ws);
      init_ws.insert(*edges_ws);
      SerializeParameter(&init_params, p);
      IDSet<job_id_t> before, after;
      StageJobAndLoadBeforeSet(&before, INIT_JOB, init_job_ids[p],
                               init_rs, init_ws);
      SpawnComputeJob(INIT_JOB, init_job_ids[p], init_rs, init_ws,
                      before, after,
                      init_params, true, GeometricRegion(p,  0, 0, 1, 1, 1));
    }
    MarkEndOfStage();
  }

  // loop job
  {
    size_t num_iterations = page_rank->num_iterations();
    Parameter loop_params;
    SerializeParameter(&loop_params, num_iterations + 1);
    std::vector<job_id_t> loop_job_id;
    GetNewJobID(&loop_job_id, 1);
    IDSet<logical_data_id_t> loop_rs, loop_ws;
    IDSet<job_id_t> before, after;
    StageJobAndLoadBeforeSet(&before, FOR_LOOP_JOB, loop_job_id[0],
                             loop_rs, loop_ws, true);
    SerializeParameter(&loop_params, num_iterations);
    SpawnComputeJob(FOR_LOOP_JOB, loop_job_id[0], loop_rs, loop_ws,
                    before, after, loop_params);
    MarkEndOfStage();
  }
}

ForLoop::ForLoop(Application* app) {
  set_application(app);
}

Job* ForLoop::Clone() {
  return new ForLoop(application());
}

void ForLoop::Execute(Parameter params, const DataArray& da) {
  // application, graph_los to help with data objects and jobs
  PageRank* page_rank = static_cast<PageRank*>(application());
  GraphLOs* graph = page_rank->graph_helper();
  size_t num_partitions = graph->num_partitions();

  uint64_t loop_counter;
  LoadParameter(&params, &loop_counter);

  // run this loop?
  if (loop_counter > 1) {
    std::cout << "Iterations left : " << loop_counter - 1 << "\n";

    StartTemplate("__MARK_STAT_for_loop");

    // scatter job
    {
      Parameter scatter_params;
      std::vector<job_id_t> scatter_job_ids;
      GetNewJobID(&scatter_job_ids, num_partitions);
      for (size_t p = 0; p < num_partitions; ++p) {
        const IDSet<logical_data_id_t>* nodes_rs = graph->GetReadSet(NODES, p);
        const IDSet<logical_data_id_t>* edges_rs = graph->GetWriteSet(EDGES, p);
        IDSet<logical_data_id_t> scatter_rs(*nodes_rs);
        scatter_rs.insert(*edges_rs);
        const IDSet<logical_data_id_t>& scatter_ws(*graph->GetWriteSet(EDGES, p));
        IDSet<job_id_t> before, after;
        StageJobAndLoadBeforeSet(&before, SCATTER_JOB, scatter_job_ids[p],
                                 scatter_rs, scatter_ws);
        SpawnComputeJob(SCATTER_JOB, scatter_job_ids[p], scatter_rs, scatter_ws,
                        before, after,
                        scatter_params, true, GeometricRegion(p, 0, 0, 1, 1, 1));
      }
      MarkEndOfStage();
    }

    // gather job
    {
      Parameter gather_params;
      std::vector<job_id_t> gather_job_ids;
      GetNewJobID(&gather_job_ids, num_partitions);
      for (size_t p = 0; p < num_partitions; ++p) {
        const IDSet<logical_data_id_t>* nodes_rs = graph->GetReadSet(NODES, p);
        const IDSet<logical_data_id_t>* edges_rs = graph->GetReadSet(EDGES, p);
        IDSet<logical_data_id_t> gather_rs(*nodes_rs);
        gather_rs.insert(*edges_rs);
        const IDSet<logical_data_id_t>& gather_ws(*graph->GetWriteSet(NODES, p));
        IDSet<job_id_t> before, after;
        StageJobAndLoadBeforeSet(&before, GATHER_JOB, gather_job_ids[p],
                                 gather_rs, gather_ws);
        SpawnComputeJob(GATHER_JOB, gather_job_ids[p], gather_rs, gather_ws,
                        before, after,
                        gather_params, true, GeometricRegion(p, 0, 0, 1, 1, 1));
      }
      MarkEndOfStage();
    }

    // loop job
    {
      Parameter loop_params;
      std::vector<job_id_t> loop_job_id;
      GetNewJobID(&loop_job_id, 1);
      IDSet<logical_data_id_t> loop_rs, loop_ws;
      IDSet<job_id_t> before, after;
      StageJobAndLoadBeforeSet(&before, FOR_LOOP_JOB, loop_job_id[0],
                               loop_rs, loop_ws, true);
      SerializeParameter(&loop_params, loop_counter - 1);
      SpawnComputeJob(FOR_LOOP_JOB, loop_job_id[0], loop_rs, loop_ws,
                      before, after, loop_params);
      MarkEndOfStage();

      EndTemplate("__MARK_STAT_for_loop");
    }
  } else {
    if (loop_counter > 0) {
      // dump all output
      dbg(DBG_APP, "Saving result.\n");

      // dump job
      {
        Parameter dump_params;
        std::vector<job_id_t> dump_job_ids;
        GetNewJobID(&dump_job_ids, num_partitions);
        for (size_t p = 0; p < num_partitions; ++p) {
          const IDSet<logical_data_id_t>& dump_rs = *(graph->GetReadSet(NODES, p));
          const IDSet<logical_data_id_t> dump_ws;
          IDSet<job_id_t> before, after;
          StageJobAndLoadBeforeSet(&before, DUMP_JOB, dump_job_ids[p],
                                   dump_rs, dump_ws);
          SpawnComputeJob(DUMP_JOB, dump_job_ids[p], dump_rs, dump_ws,
                          before, after,
                          dump_params, true, GeometricRegion(p, 0, 0, 1, 1, 1));
        }
        MarkEndOfStage();
      }

      // loop job
      {
        Parameter loop_params;
        std::vector<job_id_t> loop_job_id;
        GetNewJobID(&loop_job_id, 1);
        IDSet<logical_data_id_t> loop_rs, loop_ws;
        IDSet<job_id_t> before, after;
        StageJobAndLoadBeforeSet(&before, FOR_LOOP_JOB, loop_job_id[0],
                                 loop_rs, loop_ws, true);
        SerializeParameter(&loop_params, loop_counter - 1);
        SpawnComputeJob(FOR_LOOP_JOB, loop_job_id[0], loop_rs, loop_ws,
                        before, after, loop_params);
        MarkEndOfStage();
      }
    } else {
      // terminate application
      TerminateApplication();
    }
  }
}

Init::Init(Application* app) {
  set_application(app);
}

Job* Init::Clone() {
  return new Init(application());
}

#define SSTR(x) dynamic_cast< std::ostringstream & >( std::ostringstream() << std::dec << x ).str()  // NOLINT

struct EdgeNodes {
  size_t src_id;
  size_t dst_id;
};

void Init::Execute(Parameter params, const DataArray& da) {
  PageRank *page_rank = static_cast<PageRank*>(application());

  // graph_los to help with data objects and jobs
  // get partition id for the job
  partition_id_t partition;
  LoadParameter(&params, &partition);
  std::string dir_name = SSTR(page_rank->input_dir() << "/" << partition);
  assert(boost::filesystem::is_directory(dir_name));

  // get write set logical objects consisting of nodes and edges
  std::vector<Data*> nodes_vector;
  GetWriteData(*this, NODES, da, &nodes_vector);
  std::vector<Data*> edges_vector;
  GetWriteData(*this, EDGES, da, &edges_vector);
  boost::unordered_map<size_t, EdgeNodes> edges_partition;

  {
    // read in graph for the nodes, initialize ranks
    std::ifstream node_file;
    node_file.open(SSTR(dir_name << "/nodes").c_str());
    NodeData& nodes = *(static_cast<NodeData*>(nodes_vector[0]));
    nodes.ResetNodes();
    std::string line;
    double rank_init_val = ALPHA/(double)(page_rank->num_nodes());  // NOLINT
    while (std::getline(node_file, line)) {
      std::vector<std::string> tokens;
      boost::algorithm::split(tokens, line,
                              boost::is_any_of(": \n"),
                              boost::token_compress_on);
      size_t node_id = boost::lexical_cast<size_t>(tokens[0]);
      size_t degree  = boost::lexical_cast<size_t>(tokens[1]);
      NodeEntry &entry = nodes[node_id];
      entry.degree = degree;
      entry.rank   = rank_init_val;
    }
    node_file.close();
  }

  {
    // read in graph for the edges
    std::ifstream edge_file;
    edge_file.open(SSTR(dir_name << "/edges").c_str());
    std::string line;
    while (std::getline(edge_file, line)) {
      std::vector<std::string> tokens;
      boost::algorithm::split(tokens, line,
                              boost::is_any_of(": \n"),
                              boost::token_compress_on);
      size_t edge_id = boost::lexical_cast<size_t>(tokens[0]);
      edges_partition[edge_id].src_id  = boost::lexical_cast<size_t>(tokens[1]);
      edges_partition[edge_id].dst_id  = boost::lexical_cast<size_t>(tokens[2]);
    }
    edge_file.close();
  }

  {
    // initialize edge logical objects, initialize delta
    size_t num_edge_data = edges_vector.size();
    std::ifstream edge_file;
    edge_file.open(SSTR(dir_name << "/edge_write_sets").c_str());
    std::string line;
    size_t num = 0;
    while (std::getline(edge_file, line)) {
      assert(num < num_edge_data);
      size_t d = num_edge_data + 1;
      for (size_t i = 0; i < num_edge_data; ++i) {
        if (edges_vector[i]->region().x() == num)
          d = i;
      }
      assert(d < num_edge_data);
      EdgeData& edges = *(static_cast<EdgeData*>(edges_vector[d]));
      std::vector<std::string> tokens;
      boost::algorithm::split(tokens, line,
                              boost::is_any_of(": \n"),
                              boost::token_compress_on);
      size_t num_edges = tokens.size() - 3;
      edges.ResetEdges(num_edges);
      for (size_t e = 0; e < num_edges; ++e) {
        size_t edge_id = boost::lexical_cast<size_t>(tokens[e + 2]);
        EdgeEntry& entry = edges[e];
        entry.src_id = edges_partition[edge_id].src_id;
        entry.dst_id = edges_partition[edge_id].dst_id;
        entry.delta  = 0;
      }
      num++;
    }
  }
}

Scatter::Scatter(Application* app) {
  set_application(app);
}

Job* Scatter::Clone() {
  return new Scatter(application());
}

void Scatter::Execute(Parameter params, const DataArray& da) {
  dbg(DBG_APP, "Scatter job on %ld .\n", region().x());
  // get read set consisting of nodes
  std::vector<Data*> nodes_vector;
  GetReadData(*this, NODES, da, &nodes_vector);
  NodeData& nodes = *(static_cast<NodeData*>(nodes_vector[0]));
  // set contributions from nodes (optimization)
  boost::unordered_map<size_t, NodeEntry>& node_data = nodes.data();
  boost::unordered_map<size_t, NodeEntry>::iterator iter;
  for (iter = node_data.begin(); iter != node_data.end(); ++iter) {
    NodeEntry& entry = iter->second;
    entry.contribution = (1.0 - ALPHA) * entry.rank / (double)(entry.degree);  // NOLINT
  }
  // get wrie set consisting of edges
  std::vector<Data*> edges_vector;
  GetWriteData(*this, EDGES, da, &edges_vector);
  // update delta
  size_t num_edge_data = edges_vector.size();
  for (size_t d = 0; d < num_edge_data; ++d) {
    EdgeData& edges = *(static_cast<EdgeData*>(edges_vector[d]));
    size_t num_edges = edges.num_edges();
    for (size_t e = 0; e < num_edges; ++e) {
      EdgeEntry& entry = edges[e];
      entry.delta = nodes[entry.src_id].contribution;
    }
  }
}

Gather::Gather(Application* app) {
  set_application(app);
}

Job* Gather::Clone() {
  return new Gather(application());
}

void Gather::Execute(Parameter params, const DataArray& da) {
  dbg(DBG_APP, "Gather job on %ld .\n", region().x());
  PageRank *page_rank = static_cast<PageRank*>(application());
  // get read set consisting of edges
  std::vector<Data*> edges_vector;
  GetReadData(*this, EDGES, da, &edges_vector);
  // get write set consisting of nodes
  std::vector<Data*> nodes_vector;
  GetWriteData(*this, NODES, da, &nodes_vector);
  NodeData& nodes = *(static_cast<NodeData*>(nodes_vector[0]));

  // reset node rank to alpha by N
  size_t num_nodes_total = page_rank->num_nodes();
  double reset_val = ALPHA/(double)(num_nodes_total);  // NOLINT
  boost::unordered_map<size_t, NodeEntry>& node_data = nodes.data();
  boost::unordered_map<size_t, NodeEntry>::iterator iter;
  for (iter = node_data.begin(); iter != node_data.end(); ++iter) {
    iter->second.rank = reset_val;
  }

  // update ranks
  size_t num_edge_data = edges_vector.size();
  for (size_t d = 0; d < num_edge_data; ++d) {
    EdgeData& edges = *(static_cast<EdgeData*>(edges_vector[d]));
    size_t num_edges = edges.num_edges();
    for (size_t e = 0; e < num_edges; ++e) {
      EdgeEntry& entry = edges[e];
      nodes[entry.dst_id].rank += entry.delta;  // (1.0 - ALPHA) already in nodeentry.contribution  NOLINT
    }
  }
}

Dump::Dump(Application* app) {
  set_application(app);
}

Job* Dump::Clone() {
  return new Dump(application());
}

void Dump::Execute(Parameter params, const DataArray& da) {
  PageRank *page_rank = static_cast<PageRank*>(application());

  std::vector<Data*> nodes_vector;
  GetReadData(*this, NODES, da, &nodes_vector);
  NodeData& nodes = *(static_cast<NodeData*>(nodes_vector[0]));
  std::string dir_name = SSTR(page_rank->output_dir() << "/" << region().x());
  boost::filesystem::create_directories(dir_name);

  std::ofstream node_ranks_file(SSTR(dir_name << "/ranks").c_str());
  boost::unordered_map<size_t, NodeEntry> &node_data = nodes.data();
  boost::unordered_map<size_t, NodeEntry>::const_iterator iter;
  for (iter = node_data.begin(); iter != node_data.end(); ++iter) {
    node_ranks_file << iter->first << " : " << iter->second.degree << " : " <<
      iter->second.rank << "\n";
  }
  node_ranks_file.close();
}

}  // namespace nimbus
