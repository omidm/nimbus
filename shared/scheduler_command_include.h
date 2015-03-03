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

#ifndef NIMBUS_SHARED_SCHEDULER_COMMAND_INCLUDE_H_
#define NIMBUS_SHARED_SCHEDULER_COMMAND_INCLUDE_H_

#include "shared/scheduler_command.h"
#include "shared/handshake_command.h"
#include "shared/add_compute_job_command.h"
#include "shared/add_copy_job_command.h"
#include "shared/spawn_compute_job_command.h"
#include "shared/spawn_copy_job_command.h"
#include "shared/spawn_job_graph_command.h"
#include "shared/compute_job_command.h"
#include "shared/create_data_command.h"
#include "shared/remote_copy_send_command.h"
#include "shared/remote_copy_receive_command.h"
#include "shared/local_copy_command.h"
#include "shared/job_done_command.h"
#include "shared/define_data_command.h"
#include "shared/define_partition_command.h"
#include "shared/ldo_add_command.h"
#include "shared/ldo_remove_command.h"
#include "shared/partition_add_command.h"
#include "shared/partition_remove_command.h"
#include "shared/terminate_command.h"
#include "shared/profile_command.h"
#include "shared/start_template_command.h"
#include "shared/end_template_command.h"
#include "shared/defined_template_command.h"
#include "shared/spawn_template_command.h"
#include "shared/save_data_command.h"
#include "shared/load_data_command.h"
#include "shared/save_data_job_done_command.h"
#include "shared/prepare_rewind_command.h"
#include "shared/worker_down_command.h"
#include "shared/start_command_template_command.h"
// #include "shared/end_command_template_command.h"
// #include "shared/spawn_command_template_command.h"


#endif  // NIMBUS_SHARED_SCHEDULER_COMMAND_INCLUDE_H_
