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
  * A simple implementation.
  * Author: Hang Qu <quhang@stanford.edu>
  */

#include "test/experimental_quh/job_index.h"

namespace nimbus {
JobIndex::JobIndex() {
  _job_table.clear();
}

JobIndex::~JobIndex() {
  while (!_job_table.empty()) {
    delete _job_table.back();
    _job_table.pop_back();
  }
}

int JobIndex::QueryJobEntry(
    const VariableRegionSet& read_set,
    const VariableRegionSet& write_set,
    JobIdSet* result) {
  int number = 0;
  VariableRegionSet temp_read_set(read_set);
  VariableRegionSet temp_write_set(write_set);
  result->clear();
  bool test_read_after_write, test_write_after_read, test_write_after_write;
  for (JobEntryTable::reverse_iterator index = _job_table.rbegin();
       index != _job_table.rend();
       index++) {
    test_read_after_write =
        temp_read_set.IntersectsAndDelete((*index)->write_set);
    test_write_after_write =
        temp_write_set.IntersectsAndDelete((*index)->write_set);
    test_write_after_read =
        temp_write_set.IntersectsTest((*index)->read_set);
    if (test_read_after_write || test_write_after_write ||
        test_write_after_read) {
      result->insert((*index)->job_id);
      number++;
    }
  }
  return number;
}

bool JobIndex::AddJobEntry(
    const job_id_t job_id,
    const VariableRegionSet& read_set,
    const VariableRegionSet& write_set) {
  JobEntry* temp = new JobEntry;
  temp->job_id = job_id;
  temp->read_set.CopyFrom(read_set);
  temp->write_set.CopyFrom(write_set);
  _job_table.push_back(temp);
  return true;
}

// TODO(quhang) Never run.
bool JobIndex::DeleteJobEntry(const job_id_t job_id) {
  for (JobEntryTable::iterator index = _job_table.begin();
       index != _job_table.end();
       index++)
    if ((*index)->job_id == job_id) {
      delete (*index);
      _job_table.erase(index);
      return true;
  }
  return false;
}

// TODO(quhang) Never run.
int JobIndex::AllJobEntries(JobIdSet* result) {
  int number = 0;
  for (JobEntryTable::reverse_iterator index = _job_table.rbegin();
       index != _job_table.rend();
       index++) {
    result->insert((*index)->job_id);
    number++;
  }
  return number;
}
}  // namespace nimbus
