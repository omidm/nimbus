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

#include <boost/unordered_map.hpp>
#include "application/page_rank/node_data.h"
#include "shared/nimbus.h"

namespace nimbus {

NodeData::NodeData() : num_nodes_(0), nodes_(0) {}

NodeData::~NodeData() {}

void NodeData::Create() {
  if (num_nodes_!= 0) {
    num_nodes_ = 0;
  }
  nodes_.clear();
}

void NodeData::Destroy() {
  if (num_nodes_!= 0) {
    num_nodes_ = 0;
  }
}

Data* NodeData::Clone() {
  return new NodeData();
}

void NodeData::Copy(Data* from) {
  NodeData* data = static_cast<NodeData*>(from);
  num_nodes_ = data->num_nodes_;
  nodes_.clear();
  nodes_ = data->nodes_;
}

bool NodeData::Serialize(SerializedData* ser_data) {
  return false;
}

bool NodeData::DeSerialize(const SerializedData &ser_data, Data** result) {
  return false;
}

}  // namespace nimbus
