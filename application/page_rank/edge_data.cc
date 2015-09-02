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

#include <string>
#include "application/page_rank/edge_data.h"
#include "application/page_rank/protobuf_compiled/data_msgs.pb.h"
#include "shared/nimbus.h"

namespace nimbus {

EdgeData::EdgeData() : num_edges_(0), edges_(0) {}

EdgeData::~EdgeData() {
  if (edges_) {
    delete edges_;
  }
}

void EdgeData::Create() {
  if (num_edges_!= 0) {
    num_edges_ = 0;
  }
  if (edges_) {
    delete edges_;
  }
}

void EdgeData::Destroy() {
  if (num_edges_!= 0) {
    num_edges_ = 0;
  }
  if (edges_) {
    delete edges_;
    edges_ = NULL;
  }
}

Data* EdgeData::Clone() {
  return new EdgeData();
}

void EdgeData::Copy(Data* from) {
  EdgeData* data = static_cast<EdgeData*>(from);
  num_edges_ = data->num_edges_;
  if (edges_) {
    delete edges_;
  }
  edges_ = new EdgeEntry[num_edges_];
  memcpy(edges_, data->edges_, num_edges_ * sizeof(EdgeEntry));
}

bool EdgeData::Serialize(SerializedData* ser_data) {
  data_msgs::EdgeDataMsg edge_data_msg;
  for (size_t i = 0; i < num_edges_; ++i) {
    EdgeEntry &entry = edges_[i];
    data_msgs::EdgeMsg *edge_msg = edge_data_msg.add_edges();
    edge_msg->set_src_id(entry.src_id);
    edge_msg->set_dst_id(entry.dst_id);
    edge_msg->set_delta(entry.delta);
  }
  std::string str;
  edge_data_msg.SerializeToString(&str);
  char* ptr = new char[str.length()];
  memcpy(ptr, str.c_str(), str.length());
  ser_data->set_data_ptr(ptr);
  ser_data->set_size(str.length());
  return true;
}

bool EdgeData::DeSerialize(const SerializedData &ser_data, Data** result) {
  data_msgs::EdgeDataMsg edge_data_msg;
  std::string str(ser_data.data_ptr_raw(), ser_data.size());
  edge_data_msg.ParseFromString(str);
  ResetEdges(edge_data_msg.edges_size());
  for (size_t i = 0; i < num_edges_; i++) {
    data_msgs::EdgeMsg edge_msg = edge_data_msg.edges(i);
    EdgeEntry &entry = edges_[i];
    entry.src_id = edge_msg.src_id();
    entry.dst_id = edge_msg.dst_id();
    entry.delta  = edge_msg.delta();
  }
  return true;
}

void EdgeData::ResetEdges(size_t num_edges) {
  num_edges_ = num_edges;
  if (edges_)
    delete edges_;
  edges_ = new EdgeEntry[num_edges_];
  memset(edges_, 0, num_edges_ * sizeof(EdgeEntry));
}

size_t EdgeData::num_edges() {
  return num_edges_;
}

}  // namespace nimbus
