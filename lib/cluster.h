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


#ifndef NIMBUS_LIB_CLUSTER_H_
#define NIMBUS_LIB_CLUSTER_H_

#include <inttypes.h>
#include <iostream> // NOLINT
#include <string>
#include <set>

namespace nimbus {

class Computer;
typedef std::set<Computer*> Hosts;

class Node;
class Link {
  Node * ns;
  Node * ne;
  unsigned int capacity;
};

#define LinkSet std::set<Link>
#define LinkSetP std::set<Link*>

enum NodeType {COMPUTER, SWITCH};

class Node {
 public:
  Node();
  ~Node();
  virtual NodeType getType();

 private:
  std::string ip_address_;
  int id_;
  LinkSet linkSet_;
};

typedef std::set<Node*> NodeSet;

class Computer : public Node {
 public:
  Computer();
  ~Computer();

  virtual NodeType type();
  virtual uint64_t memory_size();
  virtual uint32_t level1_cacheSize();
  virtual uint32_t level2_cacheSize();
  virtual uint32_t level3_cacheSize();
  virtual uint32_t core_count();

 private:
  uint64_t memory_size_;
  uint32_t level1_cache_size_;
  uint32_t level2_cache_size_;
  uint32_t level3_cache_size_;
  uint32_t core_count_;
};

class Switch : public Node {
 public:
  Switch();
  ~Switch();

  virtual NodeType type();
  virtual uint32_t port_count();
  virtual uint64_t cross_section_bandwidth();

 private:
  uint32_t port_count_;
  uint64_t cross_section_bandwidth_;
};


class ClusterMap {
 public:
  ClusterMap();
  ~ClusterMap();

  void addNode(Node * node);
  void deleteNode(Node * node);
  void addLink(Link * link);
  void deleteLink(Link * link);
  uint64_t latencyNs(Node * source, Node * destination);
  uint64_t capacityBps(Node * source, Node * destination);
  void route(Node * source, Node * destination, NodeSet* storage);

 private:
  LinkSet link_set_;
  NodeSet node_set_;
};

}  // namespace nimbus
#endif  // NIMBUS_LIB_CLUSTER_H_


