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

#include "scheduler/checkpoint_manager.h"

using namespace nimbus; // NOLINT



CheckpointManager::CheckpointManager() {
  checkpoint_id_ = NIMBUS_INIT_CHECKPOINT_ID;
}

CheckpointManager::~CheckpointManager() {
}

bool CheckpointManager::CreateNewCheckpoint(checkpoint_id_t *checkpoint_id) {
  ++checkpoint_id_;
  index_[checkpoint_id_] = new CheckpointEntry(checkpoint_id_);
  *checkpoint_id = checkpoint_id_;
  return true;
}

bool CheckpointManager::AddJobToCheckpoint(checkpoint_id_t checkpoint_id,
                                           const JobEntry *job) {
  Index::iterator iter = index_.find(checkpoint_id);
  if (iter == index_.end()) {
    dbg(DBG_ERROR, "ERROR: checkpoint with id %lu does not exist!\n", checkpoint_id);
    exit(-1);
    return false;
  } else {
    iter->second->AddJob(job);
    return true;
  }
}

bool CheckpointManager::CompleteJobForCheckpoint(checkpoint_id_t checkpoint_id,
                                                 const JobEntry *job) {
  Index::iterator iter = index_.find(checkpoint_id);
  if (iter == index_.end()) {
    dbg(DBG_ERROR, "ERROR: checkpoint with id %lu does not exist!\n", checkpoint_id);
    exit(-1);
    return false;
  } else {
    iter->second->CompleteJob(job);
    return true;
  }
}

bool CheckpointManager::AddSaveDataJobToCheckpoint(checkpoint_id_t checkpoint_id,
                                                   job_id_t job_id,
                                                   logical_data_id_t ldid,
                                                   data_version_t version) {
  Index::iterator iter = index_.find(checkpoint_id);
  if (iter == index_.end()) {
    dbg(DBG_ERROR, "ERROR: checkpoint with id %lu does not exist!\n", checkpoint_id);
    exit(-1);
    return false;
  } else {
    iter->second->AddSaveDataJob(job_id, ldid, version);
    return true;
  }
}

bool CheckpointManager::NotifySaveDataJobDoneForCheckpoint(checkpoint_id_t checkpoint_id,
                                                           job_id_t job_id,
                                                           std::string handle) {
  Index::iterator iter = index_.find(checkpoint_id);
  if (iter == index_.end()) {
    dbg(DBG_ERROR, "ERROR: checkpoint with id %lu does not exist!\n", checkpoint_id);
    exit(-1);
    return false;
  } else {
    iter->second->NotifySaveDataJobDone(job_id, handle);
    return true;
  }
}

bool CheckpointManager::GetCheckpointToRewind(checkpoint_id_t *checkpoint_id) {
  return false;
}

bool CheckpointManager::GetJobListFromCheckpoint(checkpoint_id_t checkpoint_id,
                                                 JobEntryList *list) {
  Index::iterator iter = index_.find(checkpoint_id);
  if (iter == index_.end()) {
    dbg(DBG_ERROR, "ERROR: checkpoint with id %lu does not exist!\n", checkpoint_id);
    exit(-1);
    return false;
  } else {
    iter->second->GetJobList(list);
    return true;
  }
}

bool CheckpointManager::GetHandleToLoadData(checkpoint_id_t checkpoint_id,
                                            logical_data_id_t ldid,
                                            data_version_t version,
                                            std::string *handle) {
  Index::iterator iter = index_.find(checkpoint_id);
  if (iter == index_.end()) {
    dbg(DBG_ERROR, "ERROR: checkpoint with id %lu does not exist!\n", checkpoint_id);
    exit(-1);
    return false;
  } else {
    iter->second->GetHandleToLoadData(ldid, version, handle);
    return true;
  }
}

bool CheckpointManager::RemoveObsoleteCheckpoints(std::list<checkpoint_id_t> *removed_list) {
  return false;
}

