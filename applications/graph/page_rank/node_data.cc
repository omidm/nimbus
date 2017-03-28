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
#include <string>
#include "applications/graph/page_rank/node_data.h"
#include "applications/graph/page_rank/protobuf_compiled/data_msgs.pb.h"
#include "src/shared/nimbus.h"

namespace nimbus {

NodeData::NodeData(std::string name) {
  set_name(name);
}

NodeData::~NodeData() {}

void NodeData::Create() {
  nodes_.clear();
}

void NodeData::Destroy() {}

Data* NodeData::Clone() {
  return new NodeData(name());
}

void NodeData::Copy(Data* from) {
  assert(false);
  NodeData* data = static_cast<NodeData*>(from);
  nodes_.clear();
  nodes_ = data->nodes_;
}

bool NodeData::Serialize(SerializedData* ser_data) {
  data_msgs::NodeDataMsg node_data_msg;
  boost::unordered_map<size_t, NodeEntry>::const_iterator iter;
  for (iter = nodes_.begin(); iter != nodes_.end(); ++iter) {
    data_msgs::NodeMsg *node_msg = node_data_msg.add_nodes();
    const NodeEntry &entry = iter->second;
    node_msg->set_id(iter->first);
    node_msg->set_degree(entry.degree);
    node_msg->set_rank(entry.rank);
  }
  std::string str;
  node_data_msg.SerializeToString(&str);
  char* ptr = new char[str.length()];
  memcpy(ptr, str.c_str(), str.length());
  ser_data->set_data_ptr(ptr);
  ser_data->set_size(str.length());
  return true;
}

bool NodeData::DeSerialize(const SerializedData &ser_data, Data** result) {
  NodeData *nd = new NodeData(name());
  *result = nd;
  data_msgs::NodeDataMsg node_data_msg;
  std::string str(ser_data.data_ptr_raw(), ser_data.size());
  node_data_msg.ParseFromString(str);
  nd->ResetNodes();
  size_t num_nodes = node_data_msg.nodes_size();
  for (size_t i = 0; i < num_nodes; i++) {
    data_msgs::NodeMsg node_msg = node_data_msg.nodes(i);
    NodeEntry entry;
    entry.degree = node_msg.degree();
    entry.rank = node_msg.rank();
    (*nd)[node_msg.id()] = entry;
  }
  return true;
}

void NodeData::ResetNodes() {
  nodes_.clear();
}

boost::unordered_map<size_t, NodeEntry> &NodeData::data() {
  return nodes_;
}

}  // namespace nimbus
