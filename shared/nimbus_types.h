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

#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif  // __STDC_FORMAT_MACROS

#include <inttypes.h>
#include <string>

#define NIMBUS_LEVELDB_PRIVATE_KEY "~/cloud/src/nimbus/scripts/nimbus-rsa-key-pair"

#define TCP_NODELAY_OPTION true

#define NIMBUS_TERMINATE_SUCCESS (exit_status_t)(0)
#define NIMBUS_TERMINATE_FAILURE (exit_status_t)(-1)

#define NIMBUS_INIT_DATA_VERSION (data_version_t)(1)
#define NIMBUS_UNDEFINED_DATA_VERSION (data_version_t)(0)
#define NIMBUS_INIT_JOB_DEPTH (job_depth_t)(0)

#define NIMBUS_SCHEDULER_ID (worker_id_t)(0)
#define NIMBUS_KERNEL_JOB_ID (job_id_t)(0)

#define JOB_ID_BATCH (job_id_t)(10000000000)
#define LOGICAL_DATA_ID_BATCH (logical_data_id_t)(10000000000)
#define PHYSICAL_DATA_ID_BATCH (physical_data_id_t)(10000000000)

#define STATIC_BINDING_RECORD (size_t)(0)

#define NIMBUS_KERNEL_JOB_NAME "kernel"
#define NIMBUS_MAIN_JOB_NAME "main"
#define NIMBUS_LOCAL_COPY_JOB_NAME "localcopy"
#define NIMBUS_REMOTE_COPY_SEND_JOB_NAME "remotecopysend"
#define NIMBUS_REMOTE_COPY_RECEIVE_JOB_NAME "remotecopyreceive"
#define NIMBUS_CREATE_DATA_JOB_NAME "createdata"
#define NIMBUS_SAVE_DATA_JOB_NAME "savedata"
#define NIMBUS_LOAD_DATA_JOB_NAME "loaddata"
#define NIMBUS_COMPLEX_JOB_NAME "complex"

#define NIMBUS_RECEIVER_KNOWN_IP "receiver_known_ip"

#define NIMBUS_FAULT_TOLERANCE_ACTIVE false
#define DEFAULT_CHECKPOINT_CREATION_RATE 30
#define NIMBUS_INIT_CHECKPOINT_ID 0

namespace nimbus {
  typedef uint32_t port_t;
  typedef uint32_t worker_id_t;
  typedef uint32_t app_id_t;
  typedef uint64_t physical_data_id_t;
  typedef uint64_t logical_data_id_t;
  typedef uint64_t job_id_t;
  typedef uint64_t checkpoint_id_t;
  typedef uint64_t command_id_t;
  typedef uint64_t partition_id_t;
  typedef uint64_t param_id_t;
  typedef uint64_t data_version_t;
  typedef uint64_t version_table_id_t;
  typedef uint64_t job_depth_t;
  typedef uint64_t counter_t;

  typedef int32_t exit_status_t;

  typedef uint32_t switch_id_t;  // Used in cluster map for network switches

  typedef int64_t int_dimension_t;
  typedef double  float_dimension_t;

  typedef uint64_t app_data_version_t;

  enum {
    WORKER_ID_NONE = 0,
    WORKER_ID_SCHEDULER = 1
  };

  enum JobType {
    JOB_COMP   = 1,
    JOB_COPY   = 2,
    JOB_CREATE = 3,
    JOB_SCHED  = 4,
    JOB_FUTURE = 5,
    JOB_SAVE   = 6,
    JOB_LOAD   = 7,
    JOB_CMPX   = 8,
    JOB_SHDW   = 9,
    JOB_TMPL   = 10
  };


}  // namespace nimbus

#endif  // NIMBUS_SHARED_NIMBUS_TYPES_H_
