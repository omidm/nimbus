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

#ifndef NIMBUS_APPLICATION_GRAPH_LIBRARY_GRAPH_LOS_H_
#define NIMBUS_APPLICATION_GRAPH_LIBRARY_GRAPH_LOS_H_

#include <map>
#include <string>
#include <vector>
#include "shared/nimbus.h"

namespace nimbus {

class GraphLOs {
  private:
    typedef IDSet<logical_data_id_t> lolist;
    size_t num_partitions_;
    size_t num_edge_los_;
    std::vector<size_t> **edge_lo_write_;
    std::vector<size_t> **edge_lo_read_;
    std::map< std::string, lolist** > lo_map_read_;
    std::map< std::string, lolist** > lo_map_write_;

  public:
    GraphLOs();
    ~GraphLOs();
    // read in saved (common) graph information for registering data/
    // spawning jobs
    void LoadGraphInfo(std::string dir_name);
    // define logical objects for a field over edges
    void DefineEdgeLogicalObjects(Job *job, std::string name);
    // register logical objects for a field over vertices
    void DefineNodeLogicalObjects(Job *job, std::string name);
    // get read set for a variable over a partition
    const IDSet<logical_data_id_t>* GetReadSet(std::string name,
                                               partition_id_t p) const;
    // get write set for a variable over a partition
    const IDSet<logical_data_id_t>* GetWriteSet(std::string name,
                                                partition_id_t p) const;
    // accessors
    size_t num_partitions() const;
};  // class GraphLOs

}  // namespace nimbus

#endif  // NIMBUS_APPLICATION_GRAPH_LIBRARY_GRAPH_LOS_H_
