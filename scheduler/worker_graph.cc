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
  int count;
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
  return true;
}


/**
 * \fn const SchedulerWorker * nimbus::WorkerGraph::WorkerById(worker_id_t id)
 * \brief Brief description.
 * \param id
 * \return
*/
const SchedulerWorker * nimbus::WorkerGraph::WorkerById(worker_id_t id) {
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
