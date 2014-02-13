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
 *   FILE: .//worker_graph.cc
 *   DATE: Mon Nov 11 21:33:44 2013
 *  DESCR:
 ***********************************************************************/
#include "scheduler/worker_graph.h"

namespace nimbus {
/**
 * \fn nimbus::WorkerGraph::WorkerGraph()
 * \brief Brief description.
 * \return
*/
nimbus::WorkerGraph::WorkerGraph() {
}


/**
 * \fn nimbus::WorkerGraph::~WorkerGraph()
 * \brief Brief description.
 * \return
*/
nimbus::WorkerGraph::~WorkerGraph() {
}


/**
 * \fn bool nimbus::WorkerGraph::AddWorker(SchedulerWorker *worker)
 * \brief Brief description.
 * \param worker
 * \return
*/
bool nimbus::WorkerGraph::AddWorker(SchedulerWorker *worker) {
  if (worker_table_.find(worker->worker_id()) != worker_table_.end()) {
    dbg(DBG_ERROR,
        "Trying to insert existing worker %lu into worker graph.\n",
        worker->worker_id());
    return false;
  } else {
    worker_table_[worker->worker_id()] = worker;
    return true;
  }
}


/**
 * \fn bool nimbus::WorkerGraph::RemoveWorker(SchedulerWorker *worker)
 * \brief Brief description.
 * \param worker
 * \return
*/
bool nimbus::WorkerGraph::RemoveWorker(SchedulerWorker *worker) {
  if (worker_table_.find(worker->worker_id()) == worker_table_.end()) {
    dbg(DBG_ERROR,
        "Trying to remove nonexistant worker %lu from worker graph.\n",
        worker->worker_id());
    return false;
  } else {
    worker_table_.erase(worker->worker_id());
    return true;
  }
}


/**
 * \fn int nimbus::WorkerGraph::AllWorkers(WorkerIdVector *dest)
 * \brief Brief description.
 * \param dest
 * \return
*/
int nimbus::WorkerGraph::AllWorkers(WorkerIdVector *dest) {
  int count = 0;
  dest->clear();
  SchedulerWorkerTable::iterator it = worker_table_.begin();
  for (; it != worker_table_.end(); ++it) {
    worker_id_t id = (*it).first;
    dest->push_back(id);
    count++;
  }
  return count;
}


/**
 * \fn bool nimbus::WorkerGraph::ProcessMessage(WorkerGraphMessage *message)
 * \brief Brief description.
 * \param message
 * \return
*/
bool nimbus::WorkerGraph::ProcessMessage(WorkerGraphMessage *message) {
  switch (message->type()) {
  case WorkerGraphMessage::WORKER_REGISTER:
    if (!message->has_worker_register()) {
      dbg(DBG_ERROR, "Malformed worker graph message received: WORKER_REGISTER message with no worker_register field.\n");  // NOLINT
      return false;
    } else {
      const WorkerRegisterMessage& msg = message->worker_register();
      return ProcessWorkerRegisterMessage(msg);
    }

  case WorkerGraphMessage::SWITCH_REGISTER:
    if (!message->has_switch_register()) {
      dbg(DBG_ERROR, "Malformed worker graph message received: SWITCH_REGISTER message with no switch_register field.\n");  // NOLINT
      return false;
    } else {
      const SwitchRegisterMessage& msg = message->switch_register();
      return ProcessSwitchRegisterMessage(msg);
    }

  case WorkerGraphMessage::WORKER_LINK:
    if (!message->has_worker_link()) {
      dbg(DBG_ERROR, "Malformed worker graph message received: WORKER_LINK message with no worker_link field.\n");  // NOLINT
      return false;
    } else {
      const WorkerLinkMessage& msg = message->worker_link();
      return ProcessWorkerLinkMessage(msg);
    }

  case WorkerGraphMessage::SWITCH_LINK:
    if (!message->has_switch_link()) {
      dbg(DBG_ERROR, "Malformed worker graph message received: SWITCH_LINK message with no switch_link field.\n");  // NOLINT
      return false;
    } else {
      const SwitchLinkMessage& msg = message->switch_link();
      return ProcessSwitchLinkMessage(msg);
    }

    break;
  case WorkerGraphMessage::UPDATE:
    if (!message->has_update()) {
      dbg(DBG_ERROR, "Malformed worker graph message received: UPDATE message with no update field.\n");  // NOLINT
      return false;
    } else {
      const UpdateMessage& msg = message->update();
      return ProcessUpdateMessage(msg);
    }

  case WorkerGraphMessage::INVALID:
    dbg(DBG_ERROR, "Malformed worker graph message received: INVALID type.\n");
    return false;

  default:
    dbg(DBG_ERROR, "Malformed worker graph message received: unrecognized type %i.\n", (int)message->type());  // NOLINT
    return false;
  }
  return false;
}


/**
 * \fn SchedulerWorker * nimbus::WorkerGraph::WorkerById(worker_id_t id)
 * \brief Brief description.
 * \param id
 * \return
*/
SchedulerWorker * nimbus::WorkerGraph::WorkerById(worker_id_t id) {
  if (worker_table_.find(id) == worker_table_.end()) {
    return NULL;
  } else {
    return worker_table_[id];
  }
}


/**
 * \fn const Computer * nimbus::WorkerGraph::ComputerById(worker_id_t id)
 * \brief Brief description.
 * \param id
 * \return
*/
const Computer * nimbus::WorkerGraph::ComputerById(worker_id_t id) {
  if (computer_table_.find(id) == computer_table_.end()) {
    return NULL;
  } else {
    return computer_table_[id];
  }
}


/**
 * \fn uint32_t nimbus::WorkerGraph::InterComputerMbps(worker_id_t src,
                                       worker_id_t dest)
 * \brief Brief description.
 * \param src
 * \param dest
 * \return
*/
uint32_t nimbus::WorkerGraph::InterComputerMbps(worker_id_t src,
                                                worker_id_t dest) {
  return 0;
}


/**
 * \fn uint32_t nimbus::WorkerGraph::InterComputerMicroSec(worker_id_t src,
                                           worker_id_t dest)
 * \brief Brief description.
 * \param src
 * \param dest
 * \return
*/
uint32_t nimbus::WorkerGraph::InterComputerMicroSec(worker_id_t src,
                                                    worker_id_t dest) {
  return 0;
}
/**
 * \fn bool nimbus::WorkerGraph::ProcessWorkerRegisterMessage(WorkerRegisterMessage &msg)
 * \brief Brief description.
 * \param msg
 * \return
 message WorkerRegisterMessage {
  required fixed32 ipv4      = 1;
  required uint32 worker_id  = 2;
  required uint32 cores      = 3 [default = 0];
  required uint32 core_clock = 4 [default = 0];
  required uint64 memory     = 5 [default = 0];
  required uint32 mbps       = 6 [default = 0];
}
*/
bool nimbus::WorkerGraph::ProcessWorkerRegisterMessage(const WorkerRegisterMessage& msg) {
  SchedulerWorker* worker = WorkerById(msg.worker_id());
  return cluster_map_.CreateComputer(worker,
                                     msg.memory(),
                                     msg.cores(),
                                     msg.core_clock(),
                                     msg.mbps(),
                                     msg.ipv4());
}


/**
 * \fn bool nimbus::WorkerGraph::ProcessSwitchRegisterMessage(SwitchRegisterMessage &msg)
 * \brief Brief description.
 * \param msg
 * \return
*/
bool nimbus::WorkerGraph::ProcessSwitchRegisterMessage(const SwitchRegisterMessage& msg) {
  return cluster_map_.CreateSwitch(msg.switch_id(),
                                   msg.ports(),
                                   msg.mbps(),
                                   msg.nsdelay(),
                                   msg.ipv4());
}



/**
 * \fn bool nimbus::WorkerGraph::ProcessWorkerLinkMessage(WorkerLinkMessage &msg)
 * \brief Brief description.
 * \param msg
 * \return
*/
bool nimbus::WorkerGraph::ProcessWorkerLinkMessage(const WorkerLinkMessage& msg) {
  cluster_map_id_t worker = cluster_map_.LookupWorkerId(msg.worker_id());
  cluster_map_id_t swich = cluster_map_.LookupSwitchId(msg.switch_id());
  return cluster_map_.AddLink(worker, swich, msg.mbps());
}


/**
 * \fn bool nimbus::WorkerGraph::ProcessSwitchLinkMessage(SwitchLinkMessage &msg)
 * \brief Brief description.
 * \param msg
 * \return
*/
bool nimbus::WorkerGraph::ProcessSwitchLinkMessage(const SwitchLinkMessage& msg) {
  cluster_map_id_t swich1 = cluster_map_.LookupSwitchId(msg.switch_id1());
  cluster_map_id_t swich2 = cluster_map_.LookupSwitchId(msg.switch_id2());
  return cluster_map_.AddLink(swich1, swich2, msg.mbps());
}


/**
 * \fn bool nimbus::WorkerGraph::ProcessUpdateMessage(UpdateMessage &msg)
 * \brief Brief description.
 * \param msg
 * \return
*/
bool nimbus::WorkerGraph::ProcessUpdateMessage(const UpdateMessage& msg) {
  // Currently do not process real-time updates
  return true;
}


}  // namespace nimbus
