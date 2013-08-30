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
 * Global declaration of Nimbus-wide types.
 * Author: Philip Levis <pal@cs.stanford.edu>
 */

#ifndef NIMBUS_SHARED_NIMBUS_TYPES_H_
#define NIMBUS_SHARED_NIMBUS_TYPES_H_

#include <inttypes.h>
#include "shared/address_book.h"

namespace nimbus {
  typedef uint32_t port_t;
  typedef uint32_t worker_id_t;
  typedef uint32_t app_id_t;
  typedef uint64_t data_id_t;
  typedef uint64_t job_id_t;
  typedef uint64_t command_id_t;
  typedef uint64_t partition_t;
  enum {
    WORKER_ID_NONE = 0,
    WORKER_ID_SCHEDULER = 1
  };

  enum JobType {
    JOB_COMP,
    JOB_SYNC
  };

  enum SchedulerCommandType {
    COMMAND_SPAWN_JOB,
    COMMAND_DEFINE_DATA,
    COMMAND_HANDSHAKE,
    COMMAND_JOB_DONE,
    COMMAND_COMPUTE_JOB,
    COMMAND_CREATE_DATA
  };

}  // namespace nimbus

#endif  // NIMBUS_SHARED_NIMBUS_TYPES_H_
