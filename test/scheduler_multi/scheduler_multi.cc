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
  * Simple Nimbus scheduler that is supposed to run the application over
  * multiple workers. This scheduler is not still smart enough to make the
  * dynamic decisions. It assigned each part of the simulation to a worker.
  * This test is intended to check the data exchange functionality and how the
  * system correctly enforces the correct program frlow.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "./scheduler_multi.h"

SimpleScheduler::SimpleScheduler(unsigned int p)
: Scheduler(p) {
}

void SimpleScheduler::SchedulerCoreProcessor() {
  Log::dbg_printLine("Simple Scheduler Core Processor");

  while (true) {
    sleep(1);
    int ready_num = 0;
    SchedulerWorkerList::iterator iter;
    for (iter = server_->workers()->begin();
        iter != server_->workers()->end(); iter++) {
      if ((*iter)->handshake_done()) {
        ready_num++;
      } else {
        ID<worker_id_t> worker_id;
        worker_id.set_elem((*iter)->worker_id());
        std::string ip("you-know");
        ID<port_t> port(0);
        HandshakeCommand cm(worker_id, ip, port);
        std::cout << "Sending command: " << cm.toStringWTags() << std::endl;
        server_->SendCommand(*iter, &cm);
      }
    }
    if (ready_num >= WORKER_NUM)
      break;

    std::cout << ready_num << " worker(s) are registered, waiting for " <<
      WORKER_NUM - ready_num << " more worker(s) to join ..."  << std::endl;

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
  IDSet<data_id_t> read, write;
  IDSet<job_id_t> before, after;
  std::string params;
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
      IDSet<data_id_t>::IDSetIter it;
      IDSet<data_id_t> read = comm->read_set();
      for (it = read.begin(); it != read.end(); it++) {
        if (create_data_.count(*it) == 0) {
          data_created = false;
          break;
        } else if (!create_data_[*it].second) {
          data_created = false;
          break;
        }
      }
      IDSet<data_id_t> write = comm->write_set();
      for (it = write.begin(); it != write.end(); it++) {
        if (create_data_.count(*it) == 0) {
          data_created = false;
          break;
        } else if (!create_data_[*it].second) {
          data_created = false;
          break;
        }
      }
      if (data_created) {
        SchedulerWorker *worker;
        if ((read.size() == 0) && (write.size() == 0)) {
          worker = *(server_->workers()->begin());
        } else if (read.size() != 0) {
          worker_id_t worker_id = create_data_[*(read.begin())].first;
          SchedulerWorkerList::iterator it = server_->workers()->begin();
          for (; it != server_->workers()->end(); it++) {
            if ((*it)->worker_id() == worker_id) {
              worker = *it;
              break;
            }
          }
        } else {
          worker_id_t worker_id = create_data_[*(write.begin())].first;
          SchedulerWorkerList::iterator it = server_->workers()->begin();
          for (; it != server_->workers()->end(); it++) {
            if ((*it)->worker_id() == worker_id) {
              worker = *it;
              break;
            }
          }
        }

        std::cout << "Sending command [to worker " << worker->worker_id()
          << "]: " << comm->toStringWTags() << std::endl;
        server_->SendCommand(worker, comm);
        pending_compute_jobs_.erase(iter++);
        delete comm;
      } else {
        ++iter;
      }
    }

    for (iter = pending_copy_jobs_.begin();
        iter != pending_copy_jobs_.end();) {
      SpawnCopyJobCommand* comm = reinterpret_cast<SpawnCopyJobCommand*>(*iter);
      bool data_created = true;
      data_id_t from_id = comm->from_id().elem();
      if (create_data_.count(from_id) == 0) {
        data_created = false;
      } else if (!create_data_[from_id].second) {
        data_created = false;
      }

      data_id_t to_id = comm->to_id().elem();
      if (create_data_.count(to_id) == 0) {
        data_created = false;
      } else if (!create_data_[to_id].second) {
        data_created = false;
      }
      if (data_created) {
        worker_id_t sender_id = create_data_[from_id].first;
        worker_id_t receiver_id = create_data_[to_id].first;

        if (sender_id == receiver_id) {
          SchedulerWorker* worker;
          SchedulerWorkerList::iterator it = server_->workers()->begin();
          for (; it != server_->workers()->end(); it++) {
            if ((*it)->worker_id() == sender_id) {
              worker = *it;
              break;
            }
          }
          LocalCopyCommand* cm = new LocalCopyCommand(comm->job_id(),
              comm->from_id(), comm->to_id(),
              comm->before_set(), comm->after_set());
          std::cout << "Sending command [to worker " << worker->worker_id()
            << "]: " << cm->toStringWTags() << std::endl;
          server_->SendCommand(worker, cm);
          delete cm;
        } else {
          SchedulerWorker* worker_sender;
          SchedulerWorkerList::iterator it = server_->workers()->begin();
          for (; it != server_->workers()->end(); it++) {
            if ((*it)->worker_id() == sender_id) {
              worker_sender = *it;
              break;
            }
          }
          SchedulerWorker* worker_receiver;
          it = server_->workers()->begin();
          for (; it != server_->workers()->end(); it++) {
            if ((*it)->worker_id() == receiver_id) {
              worker_receiver = *it;
              break;
            }
          }



          std::cout << "OMID********************" <<
            worker_receiver->worker_id() << std::endl;

          std::vector<job_id_t> j;
          id_maker_.GetNewJobID(&j, 1);
          ID<job_id_t> id(j[0]);
          RemoteCopyCommand* cm_s = new RemoteCopyCommand(id,
              comm->from_id(), comm->to_id(),
              ID<worker_id_t>(worker_receiver->worker_id()),
              worker_receiver->ip(), ID<port_t>(worker_receiver->port()),
              comm->before_set(), comm->after_set());
          std::cout << "Sending command [to worker " << worker_sender->worker_id()
            << "]: " << cm_s->toStringWTags() << std::endl;
          server_->SendCommand(worker_sender, cm_s);
          delete cm_s;

          RemoteCopyCommand* cm_r = new RemoteCopyCommand(comm->job_id(),
              comm->from_id(), comm->to_id(),
              ID<worker_id_t>(worker_receiver->worker_id()),
              worker_receiver->ip(), ID<port_t>(worker_receiver->port()),
              comm->before_set(), comm->after_set());
          std::cout << "Sending command [to worker " << worker_receiver->worker_id()
            << "]: " << cm_r->toStringWTags() << std::endl;
          server_->SendCommand(worker_receiver, cm_r);
          delete cm_r;
        }
        pending_copy_jobs_.erase(iter++);
        delete comm;
      } else {
        ++iter;
      }
    }
  }
}

void SimpleScheduler::ProcessSpawnComputeJobCommand(SpawnComputeJobCommand* cm) {
  SchedulerCommand* comm = new ComputeJobCommand(cm->job_name(), cm->job_id(),
      cm->read_set(), cm->write_set(), cm->before_set(), cm->after_set(),
      cm->params());
  pending_compute_jobs_.push_back(comm);
}

void SimpleScheduler::ProcessSpawnCopyJobCommand(SpawnCopyJobCommand* cm) {
  SchedulerCommand* comm = new SpawnCopyJobCommand(cm->job_id(),
      cm->from_id(), cm->to_id(), cm->before_set(), cm->after_set(),
      cm->params());
  pending_copy_jobs_.push_back(comm);
}

void SimpleScheduler::ProcessDefineDataCommand(DefineDataCommand* cm) {
  std::vector<job_id_t> j;
  id_maker_.GetNewJobID(&j, 1);
  ID<job_id_t> id(j[0]);
  IDSet<job_id_t> before, after;
  SchedulerCommand* comm = new CreateDataCommand(id, cm->data_name(),
      cm->data_id(), before, after);
  job_data_map_[id.elem()] = cm->data_id().elem();

  SchedulerWorker* worker = NULL;
  // TODO(omidm): do somthing smarter!!
  worker_id_t worker_id = (worker_id_t)cm->partition_id().elem();
  SchedulerWorkerList::iterator it = server_->workers()->begin();
  for (; it != server_->workers()->end(); it++) {
    if ((*it)->worker_id() == worker_id) {
      worker = *it;
      break;
    }
  }

  create_data_[cm->data_id().elem()] = std::make_pair(worker_id, false);
  std::cout << "Sending command [to worker " << worker->worker_id()
    << "]: " << comm->toStringWTags() << std::endl;
  server_->SendCommand(worker, comm);
  delete comm;
}

void SimpleScheduler::ProcessJobDoneCommand(JobDoneCommand* cm) {
  SchedulerWorkerList::iterator iter = server_->workers()->begin();
  for (; iter != server_->workers()->end(); iter++) {
    std::cout << "Sending command: " << cm->toStringWTags() << std::endl;
    server_->SendCommand(*iter, cm);
  }
  if (job_data_map_.count(cm->job_id().elem()) != 0) {
    create_data_[job_data_map_[cm->job_id().elem()]].second = true;
  }
}



