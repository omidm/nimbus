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
  * StragglerMap is the data structure that keeps track of the load imbalance
  * among the workers. LoadBalancer uses this data structure to record the
  * information about the blamed workers that blocked a faster worker, and then
  * can query the data structure to fine the most imbalance workers
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#include "scheduler/straggler_map.h"

namespace nimbus {

StragglerMap::StragglerMap() {
}

StragglerMap::~StragglerMap() {
  MapIter iter = map_.begin();
  for (; iter != map_.end(); ++iter) {
    delete iter->second;
  }
}

void StragglerMap::AddRecord(const worker_id_t& suffered,
                             const worker_id_t& blamed) {
  MapIter iter = map_.find(suffered);
  if (iter == map_.end()) {
    Table *table = new Table();
    table->operator[](blamed) = 1;
    map_[suffered] = table;
  } else {
    Table *table = iter->second;
    TableIter it = table->find(blamed);
    if (it == table->end()) {
      table->operator[](blamed) = 1;
    } else {
      it->second++;
    }
  }
}

void StragglerMap::ClearRecords() {
  MapIter iter = map_.begin();
  for (; iter != map_.end(); ++iter) {
    delete iter->second;
  }
  map_.clear();
}

bool StragglerMap::GetMostImbalanceWorkers(worker_id_t *fast,
                                           worker_id_t *slow) {
  if (map_.size() == 0) {
    dbg(DBG_ERROR, "ERROR: StragglerMap: map is emty, cannot search for most imbalanced workers.");
    return false;
  }

  int diff;
  worker_id_t f_id, s_id;
  MapIter iter = map_.begin();
  Table *table = iter->second;
  assert(table->size());
  TableIter it = table->begin();
  f_id = iter->first;
  s_id = it->first;
  diff = it->second - LookUp(s_id, f_id);

  for (; iter != map_.end(); ++iter) {
    Table *table = iter->second;
    assert(table->size());
    TableIter it = table->begin();
    int diff_temp = it->second - LookUp(it->first, iter->first);
    if (diff_temp > diff) {
      f_id = iter->first;
      s_id = it->first;
      diff = diff_temp;
    }
  }

  *fast = f_id;
  *slow = s_id;
  return true;
}

size_t StragglerMap::LookUp(const worker_id_t& suffered,
                            const worker_id_t& blamed) {
  MapIter iter = map_.find(suffered);
  if (iter == map_.end()) {
    return 0;
  } else {
    Table *table = iter->second;
    TableIter it = table->find(blamed);
    if (it == table->end()) {
      return 0;
    } else {
      return it->second;
    }
  }
}

}  // namespace nimbus
