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

#ifndef NIMBUS_APPLICATION_GRAPH_LIBRARY_GRAPH_PARTITIONER_H_
#define NIMBUS_APPLICATION_GRAPH_LIBRARY_GRAPH_PARTITIONER_H_

#include <set>
#include <string>
#include <vector>

namespace nimbus {

class GraphPartitioner {
  private:
    size_t num_nodes_;     // number of nodes in graph
    size_t *node_degree_;  // degree of each node
    size_t num_edges_;     // number of edges in graph
    size_t *edge_src_;     // source node index for edges
    size_t *edge_dst_;     // destination node index for egdes
    size_t num_partitions_;        // number of partitions
    size_t *node_partitions_;      // node partition id

    size_t num_edge_logical_objects_;  // number of disjoint edge logical objects
    std::vector<size_t> **node_logical_objects_;  // nodes in each logical object (length = num_partitions)  NOLINT
    std::vector<size_t> **edge_logical_objects_;  // edges in each logical object (length = num_partitions^2)  NOLINT
    std::vector<bool> *edge_logical_object_exists_;  // pi*num_partitions + pj = true if logical object pi->pj exist  NOLINT
    std::set<size_t> **edge_lo_write_;  // pth element is a set of edge logical objects that partition p writes (outgoing edges)  NOLINT
    std::set<size_t> **edge_lo_read_;   // pth element is a set of edge logical objects that partition p reads (incoming edges)  NOLINT

  public:
    GraphPartitioner();
    ~GraphPartitioner();
    // load a graph from tsv nodes and edges files
    void LoadFromTSV(std::string node_file_name,
                     std::string edge_file_name);
    // load a graph by reading in number of nodes, number of edges, and edges
    void LoadGraph(std::string num_nodes_file_name,
                   std::string num_edges_file_name,
                   std::string edges_file_name);
    // construct random edge-cut partitions on graph, with p partitions
    void PartitionRandomEdgeCut(size_t num_partitions);
    // construct random edge-cut partitions on graph, with p partitions
    void PartitionRandomEdgeCutRefine(
        size_t num_partitions, size_t step1_partitions, size_t step2_partitions);
    // construct random edge-cut partitions on graph, with p partitions
    void PartitionRandomEdgeCutCoelesce(
        size_t num_partitions, size_t step1_partitions, size_t step2);
    void PartitionUsingInput(
        size_t num_partitions, std::string file_name);
    // construct random edge-cut partitions on graph, with p partitions, and
    // refine partitions based on number of neighboring vertices in another
    // partition
    void PartitionRandomEdgeCutPasses(size_t num_partitions, size_t num_passes);
    // determine logical objects for an edge-cut partitioned graph, and
    // read & write sets
    void DetermineLogicalObjects();
    // save graph over each partition -- for reconstruction purposes
    void SaveGraphInEachPartition(std::string dir_name) const;
    // save logical objects for each partition
    void SaveLogicalObjects(std::string dir_name) const;
};  // class GraphPartitioner

}  // namespace nimbus

#endif  // NIMBUS_APPLICATION_GRAPH_LIBRARY_GRAPH_PARTITIONER_H_
