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
  * PartitionGraph for uniform grid. Supports only static partitioning.
  *
  * Author: Chinmayee Shah <chinmayee.shah@stanford.edu>
  */

#ifndef NIMBUS_DATA_PARTITION_GRAPH_UNIFORM_GRID_3D_H_
#define NIMBUS_DATA_PARTITION_GRAPH_UNIFORM_GRID_3D_H_

#include "shared/nimbus_types.h"
#include "data/partition_graph_flat.h"

namespace nimbus {

    // TODO(chinmayee): Implement Coord3d and Vertex class
    class Coord3d;

    class PartitionGraphUniformGrid3d : public PartitionGraphFlat {
        public:
            PartitionGraphUniformGrid3d();
            // add main nodes and neighbor information
            bool addMainPartition(Vertex partition, Vertices neighbors);
            // add ghost nodes, parent-child and ghost-neighbor information
            bool addGhostPartition
                (Vertex partition, Vertex parent, Vertices neighbors);
            // finalize the partition graph vertices and relations, programmer
            // cannot make changes to partition graph after this method is
            // called
            void finalize();
            // create uniform partitions and add all nodes and relation
            // information
            bool createUniformPartitions
                (partition_t x, partition_t y, partition_t z);

        private:
            // range of grid, TODO(chinmayee): change these to Coord3d
            int w_min_, w_max_;
            // set of main nodes and ghost nodes
            Vertices *main_nodes_, *ghost_nodes_;
    };
}  // namespace nimbus

#endif  // NIMBUS_DATA_PARTITION_GRAPH_UNIFORM_GRID_3D_H_
