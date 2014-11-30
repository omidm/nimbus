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
  * Scheduler Checkpoint Manager. This module serves the job manager by
  * creating and keeping track of the checkpoints in the system.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_CHECKPOINT_MANAGER_H_
#define NIMBUS_SCHEDULER_CHECKPOINT_MANAGER_H_

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/thread.hpp>
#include <utility>
#include <list>
#include <string>
#include "shared/nimbus_types.h"
#include "shared/dbg.h"
#include "scheduler/job_entry.h"
#include "scheduler/checkpoint_entry.h"

namespace nimbus {

class CheckpointEntry;

class CheckpointManager {
  public:
    typedef boost::unordered_map<checkpoint_id_t, CheckpointEntry*> Index;

    CheckpointManager();
    virtual ~CheckpointManager();

    bool CreatNewCheckpoint(checkpoint_id_t *checkpoint_id);

    bool AddJobToCheckpoint(checkpoint_id_t checkpoint_id,
                            const JobEntry *job);

    bool CompleteJobForCheckpoint(checkpoint_id_t checkpoint_id,
                                  const JobEntry *job);

    bool AddSaveDataJobToCheckpoint(checkpoint_id_t checkpoint_id,
                                    job_id_t job_id,
                                    logical_data_id_t ldid,
                                    data_version_t version);

    bool NotifySaveDataJobDoneForCheckpoint(checkpoint_id_t checkpoint_id,
                                            job_id_t job_id,
                                            std::string handle);

    bool GetCheckpointToRewind(checkpoint_id_t *checkpoint_id);

    bool GetJobListFromCheckpoint(checkpoint_id_t checkpoint_id,
                                  JobEntryList *list);

    bool GetHandleToLoadData(checkpoint_id_t checkpoint_id,
                             logical_data_id_t ldid,
                             data_version_t version,
                             std::string *handle);

    bool RemoveObsoleteCheckpoints(std::list<checkpoint_id_t> *removed_list);

  private:
    Log log_;
    Index index_;
    checkpoint_id_t checkpoint_id_;
};

}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_CHECKPOINT_MANAGER_H_
