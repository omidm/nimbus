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
  * Scheduler abstraction of a worker.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_SCHEDULER_WORKER_H_
#define NIMBUS_SCHEDULER_SCHEDULER_WORKER_H_

#include <list>
#include "lib/nimbus.h"
#include "lib/scheduler_server_connection.h"

namespace nimbus {

class Application;

class SchedulerWorker {
 public:
  SchedulerWorker(worker_id_t id,
                  SchedulerServerConnection* conn,
                  Application* app);
  virtual ~SchedulerWorker();

  virtual worker_id_t worker_id();
  virtual SchedulerServerConnection* connection();
  virtual Application* application();
  virtual bool is_alive();
  virtual void MarkDead();
  virtual char* read_buffer();
  virtual uint32_t existing_bytes();
  virtual void set_existing_bytes(uint32_t bytes);
  virtual uint32_t read_buffer_length();

 private:
  worker_id_t worker_id_;
  SchedulerServerConnection* connection_;
  Application* application_;
  bool is_alive_;
  char* read_buffer_;
  // How many bytes in the read buffer are valid before
  // a read.
  uint32_t existing_bytes_;
};

typedef std::list<SchedulerWorker*> SchedulerWorkerList;

}  // namespace nimbus

#endif  // NIMBUS_SCHEDULER_SCHEDULER_WORKER_H_
