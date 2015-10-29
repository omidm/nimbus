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

#ifndef NIMBUS_APPLICATION_PAGE_RANK_NODE_DATA_H_
#define NIMBUS_APPLICATION_PAGE_RANK_NODE_DATA_H_

#include <boost/unordered_map.hpp>
#include <string>
#include "shared/nimbus.h"

#define NODES "nodes"

namespace nimbus {

struct NodeEntry {
  size_t degree;
  double rank;
  double contribution;  // must be initialized for every use
};

class NodeData : public Data {
  public:
    explicit NodeData(std::string name);
    ~NodeData();

    virtual void Create();
    virtual void Destroy();
    virtual Data* Clone();
    virtual void Copy(Data* from);
    virtual bool Serialize(SerializedData* ser_data);
    virtual bool DeSerialize(const SerializedData& ser_data, Data** result);

    void ResetNodes();
    inline NodeEntry& operator[](size_t nid) { return nodes_[nid]; }
    boost::unordered_map<size_t, NodeEntry> &data();

  private:
    boost::unordered_map<size_t, NodeEntry> nodes_;
};  // class NodeData

}  // namespace nimbus

#endif  // NIMBUS_APPLICATION_PAGE_RANK_NODE_DATA_H_
