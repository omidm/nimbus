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
  * A Nimbus scheduler's abstraction of a worker.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "src/scheduler/scheduler_worker.h"
#include "src/shared/dbg.h"

namespace nimbus {

SchedulerWorker::SchedulerWorker(worker_id_t id,
                                 SchedulerServerConnection* conn,
                                 ApplicationGroup* app) {
  worker_id_ = id;
  connection_ = conn;
  application_ = app;
  is_alive_ = true;
  handshake_done_ = false;
}

SchedulerWorker::~SchedulerWorker() {
  delete connection_;
}

worker_id_t SchedulerWorker::worker_id() const {
  return worker_id_;
}

std::string SchedulerWorker::ip() const {
  return ip_;
}

void SchedulerWorker::set_ip(std::string ip) {
  ip_ = ip;
}

port_t SchedulerWorker::port() const {
  return port_;
}

void SchedulerWorker::set_port(port_t port) {
  port_ = port;
}

SchedulerServerConnection* SchedulerWorker::connection() {
  return connection_;
}

ApplicationGroup* SchedulerWorker::application() {
  return application_;
}

bool SchedulerWorker::is_alive() const {
  return is_alive_;
}

bool SchedulerWorker::handshake_done() const {
  return handshake_done_;
}

void SchedulerWorker::set_handshake_done(bool flag) {
  handshake_done_ = flag;
}

void SchedulerWorker::MarkDead() {
  is_alive_ = false;
}

}  // namespace nimbus
