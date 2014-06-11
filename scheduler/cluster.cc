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

/***********************************************************************
 * AUTHOR: Philip Levis <pal>
 *   FILE: .//cluster.cc
 *   DATE: Tue Nov 12 18:45:45 2013
 *  DESCR:
 ***********************************************************************/
#include "scheduler/cluster.h"
#include "shared/dbg.h"

namespace nimbus {
/**
 * \fn nimbus::Node::Node(uint32_t ipv4,
                   cluster_map_id_t id)
 * \brief Brief description.
 * \param ipv4
 * \param id
 * \return
*/

nimbus::Node::Node() {}
nimbus::Node::Node(uint32_t addr,
                   cluster_map_id_t c) {
  ipv4_ = addr;
  id_ = c;
}


/**
 * \fn nimbus::Node::~Node()
 * \brief Brief description.
 * \return
*/
nimbus::Node::~Node() {
  LinkPtrSet::iterator it;
  for (it = in_links_.begin();
       it != in_links_.end();
       it = in_links_.begin()) {
    Link* link = *it;
    delete link;
  }

  for (it = out_links_.begin();
       it != out_links_.end();
       it = out_links_.begin()) {
    Link* link = *it;
    delete link;
  }
}

nimbus::Computer::Computer() {}

/**
 * \fn nimbus::Computer::Computer()
 * \brief Brief description.
 * \return
*/
nimbus::Computer::Computer(uint32_t ipv4, cluster_map_id_t id)
  : Node(ipv4, id) {}


/**
 * \fn nimbus::Computer::~Computer()
 * \brief Brief description.
 * \return
*/
nimbus::Computer::~Computer() {
}


/**
 * \fn nimbus::Switch::Switch()
 * \brief Brief description.
 * \return
*/
nimbus::Switch::Switch(uint32_t ipv4, cluster_map_id_t id) : Node(ipv4, id) {}



/**
 * \fn nimbus::Switch::~Switch()
 * \brief Brief description.
 * \return
*/
nimbus::Switch::~Switch() {
}


/**
 * \fn nimbus::ClusterMap::ClusterMap()
 * \brief Brief description.
 * \return
*/
nimbus::ClusterMap::ClusterMap() {
  id_ = 0;
}


/**
 * \fn nimbus::ClusterMap::~ClusterMap()
 * \brief Brief description.
 * \return
*/
nimbus::ClusterMap::~ClusterMap() {
  NodeMap::iterator it = node_map_.begin();
  for (; it != node_map_.end(); ++it) {
    Node* n = (*it).second;
    dbg(DBG_MEMORY, "Deleting node %i 0x%p\n", (*it).first, n);
    delete n;
  }
}


/**
 * \fn cluster_map_id_t nimbus::ClusterMap::CreateComputer(SchedulerWorker *worker,
                                 uint64_t memory,
                                 uint32_t cores,
                                 uint32_t MHz,
                                 uint32_t mbps,
                                 uint32_t ipv4)
 * \brief Brief description.
 * \param worker
 * \param memory
 * \param cores
 * \param MHz
 * \param mbps
 * \param ipv4
 * \return
*/
cluster_map_id_t nimbus::ClusterMap::CreateComputer(SchedulerWorker *work,
                                                    uint64_t memory,
                                                    uint32_t cores,
                                                    uint32_t MHz,
                                                    uint32_t mbps,
                                                    uint32_t ipv4) {
  Computer* computer = new Computer(ipv4, id_);
  id_++;
  dbg(DBG_MEMORY, "Allocating computer %i: 0x%p\n", computer->id(), computer);
  computer->set_memory(memory);
  computer->set_cores(cores);
  computer->set_freq_MHz(MHz);
  computer->set_mbps(mbps);
  computer->set_ipv4(ipv4);
  computer->set_worker(work);

  computer_map_[work->worker_id()] = computer;
  node_map_[computer->id()] = computer;
  return computer->id();
}


/**
 * \fn cluster_map_id_t nimbus::ClusterMap::CreateSwitch(switch_id_t id,
                                 uint32_t ports,
                                 uint32_t mbps,
                                 uint32_t nsdelay,
                                 uint32_t ipv4)
 * \brief Brief description.
 * \param id
 * \param ports
 * \param mbps
 * \param nsdelay
 * \param ipv4
 * \return
*/
cluster_map_id_t nimbus::ClusterMap::CreateSwitch(switch_id_t id,
                                                  uint32_t ports,
                                                  uint32_t mbps,
                                                  uint32_t nsdelay,
                                                  uint32_t ipv4) {
  Switch* swich = new Switch(ipv4, id_);
  dbg(DBG_MEMORY, "Allocating switch %i: 0x%p\n", swich->id(), swich);
  id_++;

  swich->set_ports(ports);
  swich->set_mbps(mbps);
  swich->set_delay_ns(nsdelay);
  swich->set_ipv4(ipv4);
  swich->set_switch_id(id);

  switch_map_[id] = swich;
  node_map_[swich->id()] = swich;
  return swich->id();
}


/**
 * \fn cluster_map_id_t nimbus::ClusterMap::AddLink(cluster_map_id_t src,
                            cluster_map_id_t dest,
                            uint32_t mbps)
 * \brief Brief description.
 * \param src
 * \param dest
 * \param mbps
 * \return
*/
cluster_map_id_t nimbus::ClusterMap::AddLink(cluster_map_id_t src,
                                             cluster_map_id_t dest,
                                             uint32_t mbps) {
  Node* s = node_map_[src];
  Node* d = node_map_[dest];

  Link* link = new Link(s, d, mbps, id_);
  id_++;

  s->out_links()->insert(link);
  d->in_links()->insert(link);

  link_set_.insert(link);
  return link->id();
}


/**
 * \fn bool nimbus::ClusterMap::Delete(cluster_map_id_t id)
 * \brief Brief description.
 * \param id
 * \return
*/
bool nimbus::ClusterMap::Delete(cluster_map_id_t id) {
  if (node_map_.find(id) == node_map_.end()) {
    dbg(DBG_ERROR, "Tried to delete non-existed cluster map entry %u\n", id);
    return false;
  } else {
    Node* node = node_map_[id];
    node_map_.erase(id);

    if (node->type() == CLUSTER_COMPUTER) {
      Computer* computer = static_cast<Computer*>(node);
      computer_map_.erase(computer->worker()->worker_id());
    } else if (node->type() == CLUSTER_SWITCH) {
      Switch* swich = static_cast<Switch*>(node);
      switch_map_.erase(swich->switch_id());
     }
    // Deleting the node will delete all of its links.
    // Each link delete will remove it from both nodes.
    dbg(DBG_MEMORY, "Deleting node %i 0x%p\n", id, node);
    delete node;
    return true;
  }
}


/**
 * \fn cluster_map_id_t nimbus::ClusterMap::LookupWorkerId(worker_id_t id)
 * \brief Brief description.
 * \param id
 * \return
*/
cluster_map_id_t nimbus::ClusterMap::LookupWorkerId(worker_id_t id) {
  if (computer_map_.find(id) == computer_map_.end()) {
    return 0;
  } else {
    Computer* comp = computer_map_[id];
    return comp->id();
  }
}


/**
 * \fn cluster_map_id_t nimbus::ClusterMap::LookupSwitchId(switch_id_t id)
 * \brief Brief description.
 * \param id
 * \return
*/
cluster_map_id_t nimbus::ClusterMap::LookupSwitchId(switch_id_t id) {
  if (switch_map_.find(id) == switch_map_.end()) {
    return 0;
  } else {
    Switch* swich = switch_map_[id];
    return swich->id();
  }
}


/**
 * \fn Computer * nimbus::ClusterMap::LookupComputer(cluster_map_id_t id)
 * \brief Brief description.
 * \param id
 * \return
*/
Computer * nimbus::ClusterMap::LookupComputer(cluster_map_id_t id) {
  if (node_map_.find(id) == node_map_.end()) {
    return NULL;
  } else {
    Node* n = node_map_[id];
    return static_cast<Computer*>(n);
  }
}


/**
 * \fn Switch * nimbus::ClusterMap::LookupSwitch(cluster_map_id_t id)
 * \brief Brief description.
 * \param id
 * \return
*/
Switch * nimbus::ClusterMap::LookupSwitch(cluster_map_id_t id) {
  if (node_map_.find(id) == node_map_.end()) {
    return NULL;
  } else {
    Node* n = node_map_[id];
    return static_cast<Switch*>(n);
  }
}

Node* nimbus::ClusterMap::LookupNode(cluster_map_id_t id) {
  if (node_map_.find(id) == node_map_.end()) {
    return NULL;
  } else {
    Node* n = node_map_[id];
    return n;
  }
}

nimbus::Link::~Link() {
  src_->out_links()->erase(this);
  dest_->in_links()->erase(this);
}

}  // namespace nimbus
