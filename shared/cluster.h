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
  * Nimbus abstraction of computational resources.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#ifndef NIMBUS_SHARED_CLUSTER_H_
#define NIMBUS_SHARED_CLUSTER_H_

#include <unordered_map>
#include <set>
#include <string>
#include "shared/nimbus_types.h"
#include "scheduler/scheduler_worker.h"

namespace nimbus {

typedef uint32_t cluster_map_id_t;

class Node;
class Link {
 public:
  Link(Node* src, Node* dest, uint32_t speed, cluster_map_id_t ID) {
    src_ = src;
    dest_ = dest;
    mbps_ = speed;
    id_ = ID;
  }
  ~Link();

  Node* source() {return src_;}
  Node* destination() {return dest_;}
  uint32_t mbps() {return mbps_;}
  cluster_map_id_t id() {return id_;}

 private:
  Node * src_;
  Node * dest_;
  uint32_t mbps_;
  cluster_map_id_t id_;
};

class Computer;
class Switch;

typedef std::set<Node*> NodeSet;
typedef std::set<Link> LinkSet;
typedef std::set<Link*> LinkPtrSet;
typedef std::unordered_map<worker_id_t, Computer*> ComputerMap;
typedef std::unordered_map<switch_id_t, Switch*> SwitchMap;
typedef std::unordered_map<cluster_map_id_t, Node*> NodeMap;

enum NodeType {CLUSTER_COMPUTER, CLUSTER_SWITCH};

class Node {
 public:
  Node();
  Node(uint32_t ipv4, cluster_map_id_t id);
  virtual ~Node();
  virtual NodeType type() = 0;
  virtual cluster_map_id_t id() {return id_;}

  virtual uint32_t mbps() {return mbps_;}
  virtual uint32_t ipv4() {return ipv4_;}

  virtual void set_mbps(uint32_t m) {mbps_ = m;}
  virtual void set_ipv4(uint32_t v) {ipv4_ = v;}

  virtual LinkPtrSet* in_links() {return &in_links_;}
  virtual LinkPtrSet* out_links() {return &out_links_;}

 private:
  uint32_t ipv4_;
  uint32_t mbps_;
  cluster_map_id_t id_;
  LinkPtrSet out_links_;
  LinkPtrSet in_links_;
};



class Computer : public Node {
 public:
  Computer();
  Computer(uint32_t ipv4, cluster_map_id_t id);
  virtual ~Computer();

  virtual NodeType type() {return CLUSTER_COMPUTER;}
  virtual SchedulerWorker* worker() {return worker_;}
  virtual uint64_t memory_size() {return memory_size_;}
  virtual uint32_t core_count() {return core_count_;}
  virtual uint32_t freq_MHz() {return freq_MHz_;}

  virtual void set_memory(uint64_t s) {memory_size_ = s;}
  virtual void set_cores(uint32_t c) {core_count_ = c;}
  virtual void set_freq_MHz(uint32_t f) {freq_MHz_ = f;}
  virtual void set_worker(SchedulerWorker* s) {worker_ = s;}

 private:
  uint64_t memory_size_;
  uint32_t core_count_;
  uint32_t freq_MHz_;
  SchedulerWorker* worker_;
};

class Switch : public Node {
 public:
  Switch(uint32_t ipv4, cluster_map_id_t id);
  ~Switch();

  virtual NodeType type() {return CLUSTER_SWITCH;}

  virtual uint64_t cross_bandwidth() {return cross_bandwidth_;}
  virtual void set_cross_bandwidth(uint64_t b) {cross_bandwidth_ = b;}

  virtual uint32_t ports() {return port_count_;}
  virtual void set_ports(uint32_t p) {port_count_ = p;}

  virtual uint32_t delay_ns() {return delay_ns_;}
  virtual void set_delay_ns(uint32_t d) {delay_ns_ = d;}

  virtual switch_id_t switch_id() {return switch_id_;}
  virtual void set_switch_id(switch_id_t i) {switch_id_ = i;}

 private:
  uint32_t port_count_;
  uint64_t cross_bandwidth_;
  uint32_t delay_ns_;
  switch_id_t switch_id_;
};


class ClusterMap {
 public:
  ClusterMap();
  ~ClusterMap();

  cluster_map_id_t CreateComputer(SchedulerWorker* worker,
                                  uint64_t memory,
                                  uint32_t cores,
                                  uint32_t MHz,
                                  uint32_t mbps,
                                  uint32_t ipv4);
  cluster_map_id_t CreateSwitch(switch_id_t id,
                                uint32_t ports,
                                uint32_t mbps,
                                uint32_t nsdelay,
                                uint32_t ipv4);
  cluster_map_id_t AddLink(cluster_map_id_t src,
                           cluster_map_id_t dest,
                           uint32_t mbps);

  bool Delete(cluster_map_id_t id);

  cluster_map_id_t LookupWorkerId(worker_id_t id);
  cluster_map_id_t LookupSwitchId(switch_id_t id);

  Computer* LookupComputer(cluster_map_id_t id);
  Switch* LookupSwitch(cluster_map_id_t id);
  Node*   LookupNode(cluster_map_id_t id);

  uint32_t LatencyNs(Node * source, Node * destination) {return 0;}
  uint32_t CapacityMbps(Node * source, Node * destination) {return 0;}

 private:
  LinkPtrSet link_set_;
  NodeSet node_set_;
  ComputerMap computer_map_;
  SwitchMap switch_map_;
  NodeMap node_map_;
  cluster_map_id_t id_;
};

}  // namespace nimbus
#endif  // NIMBUS_SHARED_CLUSTER_H_
