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

#ifndef NIMBUS_SRC_SHARED_SCHEDULER_COMMAND_INCLUDE_H_
#define NIMBUS_SRC_SHARED_SCHEDULER_COMMAND_INCLUDE_H_

#include "src/shared/scheduler_command.h"
#include "src/shared/handshake_command.h"
#include "src/shared/add_compute_job_command.h"
#include "src/shared/add_copy_job_command.h"
#include "src/shared/spawn_compute_job_command.h"
#include "src/shared/spawn_copy_job_command.h"
#include "src/shared/spawn_job_graph_command.h"
#include "src/shared/compute_job_command.h"
#include "src/shared/create_data_command.h"
#include "src/shared/remote_copy_send_command.h"
#include "src/shared/remote_copy_receive_command.h"
#include "src/shared/mega_rcr_command.h"
#include "src/shared/mega_job_done_command.h"
#include "src/shared/local_copy_command.h"
#include "src/shared/job_done_command.h"
#include "src/shared/define_data_command.h"
#include "src/shared/define_partition_command.h"
#include "src/shared/ldo_add_command.h"
#include "src/shared/ldo_remove_command.h"
#include "src/shared/partition_add_command.h"
#include "src/shared/partition_remove_command.h"
#include "src/shared/terminate_command.h"
#include "src/shared/profile_command.h"
#include "src/shared/start_template_command.h"
#include "src/shared/end_template_command.h"
#include "src/shared/defined_template_command.h"
#include "src/shared/spawn_template_command.h"
#include "src/shared/save_data_command.h"
#include "src/shared/load_data_command.h"
#include "src/shared/save_data_job_done_command.h"
#include "src/shared/prepare_rewind_command.h"
#include "src/shared/worker_down_command.h"
#include "src/shared/start_command_template_command.h"
#include "src/shared/end_command_template_command.h"
#include "src/shared/spawn_command_template_command.h"
#include "src/shared/request_stat_command.h"
#include "src/shared/respond_stat_command.h"
#include "src/shared/print_stat_command.h"


#endif  // NIMBUS_SRC_SHARED_SCHEDULER_COMMAND_INCLUDE_H_
