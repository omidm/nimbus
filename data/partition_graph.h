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
  * PartitionGraph maintains partitions and relationships between partitions
  * for a simulation domain (uniform grid/ mesh/ tree). A vertex represents a
  * partition. A vertex has associated sets of data ids corresponding to the
  * partition.
  *
  * PartitionGraph is an interface for different type specific implementation.
  *
  * PartitionGraph may be used by the scheduler to make good decisions about
  * data placement -- for instance, to place ghost regions with parent main
  * region, and to place data corresponding to the same partition closeby.
  *
  * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
  */

#ifndef NIMBUS_DATA_PARTITION_GRAPH_H_
#define NIMBUS_DATA_PARTITION_GRAPH_H_

#include <set>
#include <map>
#include "shared/nimbus_types.h"

namespace nimbus {
    class Vertex;
    typedef std::set<Vertex> Vertices;
    typedef std::set<data_id_t> DataSet;
    typedef std::map<Vertex, DataSet*> VertexDataSetMap;
    typedef std::map<Vertex, Vertices*> VertexVerticesMap;

    class PartitionGraph {
        public:
            PartitionGraph() {}
            virtual ~PartitionGraph() {}

            // get main nodes corresponding to main partitions
            virtual Vertices* getMainNodes() = 0;
            // get ghost nodes corresponding to ghost partitions
            virtual Vertices* getGhostNodes() = 0;
            // get the map between vertices and corresponding data ids
            virtual VertexDataSetMap* getVertexDataMap() = 0;
            // get parent-child relations between main nodes and ghost nodes
            virtual VertexVerticesMap* getParentChildMap() = 0;
    };
}  // namespace nimbus

#endif  // NIMBUS_DATA_PARTITION_GRAPH_H_
