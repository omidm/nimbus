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
  * Simple Nimbus scheduler that is supposed to run the application over a
  * single worker. It is intended to check the command exchange interface, the
  * mapping logics and generally the system abstraction soundness.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "./scheduler_single.h"
#define WORKER_NUM 4

SimpleScheduler::SimpleScheduler(unsigned int p)
: Scheduler(p) {
}

void SimpleScheduler::SchedulerCoreProcessor() {
  Log::dbg_printLine("Simple Scheduler Core Processor");

  while (server_->workers()->begin() == server_->workers()->end()) {
    std::cout << "Waiting for the first worker to connect ..." << std::endl;
    sleep(1);
  }

  SchedulerWorkerList::iterator iter = server_->workers()->begin();
  while (true) {
    sleep(1);
    std::cout << "Waiting for the handshake to complete ..." << std::endl;
    if ((*iter)->handshake_done()) {
      break;
    } else {
      ID<worker_id_t> worker_id;
      worker_id.set_elem((*iter)->worker_id());
      std::string ip("you-know");
      ID<port_t> port(0);
      HandshakeCommand cm(worker_id, ip, port);
      std::cout << "Sending command: " << cm.toStringWTags() << std::endl;
      server_->SendCommand(*iter, &cm);
    }

    SchedulerCommandList storage;
    if (server_->ReceiveCommands(&storage, (uint32_t)10)) {
      SchedulerCommandList::iterator iter = storage.begin();
      for (; iter != storage.end(); iter++) {
        SchedulerCommand* comm = *iter;
        std::cout << "Received command: " << comm->toStringWTags() << std::endl;
        ProcessSchedulerCommand(comm);
        delete comm;
      }
    }
  }


  std::vector<job_id_t> j;
  id_maker_.GetNewJobID(&j, 1);
  ID<job_id_t> id(j[0]);
  IDSet<physical_data_id_t> read, write;
  IDSet<job_id_t> before, after;
  Parameter params;
  ComputeJobCommand cm("main", id, read, write, before, after, params);
  std::cout << "Sending command: " << cm.toStringWTags() << std::endl;
  server_->SendCommand(*(server_->workers()->begin()), &cm);

  while (true) {
    SchedulerCommandList storage;
    if (server_->ReceiveCommands(&storage, (uint32_t)10)) {
      SchedulerCommandList::iterator iter = storage.begin();
      for (; iter != storage.end(); iter++) {
        SchedulerCommand* comm = *iter;
        std::cout << "Received command: " << comm->toStringWTags() << std::endl;
        ProcessSchedulerCommand(comm);
        delete comm;
      }
    }

    SchedulerCommandList::iterator iter;
    for (iter = pending_compute_jobs_.begin();
        iter != pending_compute_jobs_.end();) {
      ComputeJobCommand* comm = reinterpret_cast<ComputeJobCommand*>(*iter);
      bool data_created = true;
      IDSet<physical_data_id_t>::IDSetIter it;
      IDSet<physical_data_id_t> read = comm->read_set();
      for (it = read.begin(); it != read.end(); it++) {
        if (create_data_.count(*it) == 0) {
          data_created = false;
          break;
        } else if (!create_data_[*it]) {
          data_created = false;
          break;
        }
      }
      IDSet<physical_data_id_t> write = comm->write_set();
      for (it = write.begin(); it != write.end(); it++) {
        if (create_data_.count(*it) == 0) {
          data_created = false;
          break;
        } else if (!create_data_[*it]) {
          data_created = false;
          break;
        }
      }
      if (data_created) {
        std::cout << "Sending command: " << comm->toStringWTags() << std::endl;
        server_->SendCommand(*(server_->workers()->begin()), comm);
        pending_compute_jobs_.erase(iter++);
        delete comm;
      } else {
        ++iter;
      }
      // sleep(1);
    }
  }
}

void SimpleScheduler::ProcessSpawnComputeJobCommand(SpawnComputeJobCommand* cm) {
  SchedulerCommand* comm = new ComputeJobCommand(cm->job_name(), cm->job_id(),
      cm->read_set(), cm->write_set(), cm->before_set(), cm->after_set(),
      cm->params());
  pending_compute_jobs_.push_back(comm);
}

void SimpleScheduler::ProcessDefineDataCommand(DefineDataCommand* cm) {
  std::vector<job_id_t> j;
  id_maker_.GetNewJobID(&j, 1);
  ID<job_id_t> id(j[0]);
  IDSet<job_id_t> before, after;
  SchedulerCommand* comm = new CreateDataCommand(id, cm->data_name(),
      cm->logical_data_id(), cm->logical_data_id(), before, after);
  job_data_map_[id.elem()] = cm->logical_data_id().elem();
  create_data_[cm->logical_data_id().elem()] = false;
  std::cout << "Sending command: " << comm->toStringWTags() << std::endl;
  server_->SendCommand(*(server_->workers()->begin()), comm);
  delete comm;

  data_manager_->AddLogicalObject(cm->logical_data_id().elem(),
      cm->data_name(),
      cm->partition_id().elem());
}

void SimpleScheduler::ProcessJobDoneCommand(JobDoneCommand* cm) {
  SchedulerWorkerList::iterator iter = server_->workers()->begin();
  for (; iter != server_->workers()->end(); iter++) {
    std::cout << "Sending command: " << cm->toStringWTags() << std::endl;
    server_->SendCommand(*iter, cm);
  }
  if (job_data_map_.count(cm->job_id().elem()) != 0) {
    create_data_[job_data_map_[cm->job_id().elem()]] = true;
  }
}



