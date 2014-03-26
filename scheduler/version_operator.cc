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
  * Scheduler Version Operator. It is the main class that performs operations
  * over version tables including merging and making roots out of normal nodes.
  * This module provides caching to expedite the merging operations over known
  * previously calculated merges.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/version_operator.h"

#define CACHE_SEED_ 123
#define VERSION_OPERATOR_CACHE_SIZE 10

using namespace nimbus; // NOLINT

VersionOperator::VersionOperator() {
  max_cache_size_ = VERSION_OPERATOR_CACHE_SIZE;
  assert(max_cache_size_ >= 1);
  cache_seed_ = CACHE_SEED_;
}

VersionOperator::~VersionOperator() {
}

bool VersionOperator::MergeVersionTables(
    std::vector<boost::shared_ptr<VersionTable> > tables,
    boost::shared_ptr<VersionTable> *result) {
  size_t count = tables.size();

  if (count == 0) {
    return false;
  }

  if (count == 1) {
    *result = tables[0];
    return true;
  }

  std::vector<boost::shared_ptr<VersionTable> > reduced;
  if ((count % 2) == 1) {
    reduced.push_back(tables[count - 1]);
    --count;
  }

  for (size_t i = 0; i < count; ++i) {
    boost::shared_ptr<VersionTable> merged;
    if (MergeTwoVersionTables(tables[i], tables[i+1], &merged)) {
      reduced.push_back(merged);
    } else {
      return false;
    }
  }

  return MergeVersionTables(reduced, result);
}

bool VersionOperator::MergeTwoVersionTables(
    boost::shared_ptr<VersionTable> t1,
    boost::shared_ptr<VersionTable> t2,
    boost::shared_ptr<VersionTable> *result) {
  std::set<version_table_id_t> ids;
  ids.insert(t1->id());
  ids.insert(t2->id());
  if (LookUpCache(ids, result)) {
    dbg(DBG_SCHED, "Version Operator: hit cache.\n");
    return true;
  } else {
    dbg(DBG_SCHED, "Version Operator: missed cache.\n");
    boost::shared_ptr<VersionTable> merged(new VersionTable());
    merged->set_id(GetNewVersionTableId());
    if (t1->root_raw()->id() == t2->root_raw()->id()) {
      dbg(DBG_SCHED, "Version Operator: roots are the same for merge.\n");
      merged->set_root(t1->root());
      VersionTable::Map content = t1->content();
      VersionTable::MapConstIter it;
      for (it = t2->content_p()->begin(); it != t2->content_p()->end(); ++it) {
        VersionTable::MapIter it_p = content.find(it->first);
        if (it_p == content.end()) {
          content[it->first] = it->second;
        } else if (it->second > it_p->second) {
          it_p->second = it->second;
        }
      }
      merged->set_content(content);
    } else {
      dbg(DBG_SCHED, "Version Operator: roots are different for merge.\n");
      // TODO(omidm): implement!
      exit(-1);
    }

    CacheMergeResult(ids, merged);
    *result = merged;
    return true;
  }
}

bool VersionOperator::MakeRootTable(
    boost::shared_ptr<VersionTable> table,
    boost::shared_ptr<VersionTable> *result) {
  // TODO(omidm): implement!
  return false;
}

bool VersionOperator::MakeVersionTableOut(
    boost::shared_ptr<VersionTable> tabel_in,
    IDSet<logical_data_id_t> write_set,
    boost::shared_ptr<VersionTable> *teble_out) {
  // TODO(omidm): implement!
  return false;
}

bool VersionOperator::LookUpCache(
    std::set<version_table_id_t> ids,
    boost::shared_ptr<VersionTable>* result) {
  Cache::iterator iter = cache_.find(ids);
  if (iter != cache_.end()) {
    *result = iter->second;
    return true;
  }
  return false;
}

bool VersionOperator::CacheMergeResult(
    std::set<version_table_id_t> ids,
    boost::shared_ptr<VersionTable> merged) {
  if (cache_.size() < max_cache_size_) {
    cache_[ids] = merged;
  } else {
    Cache::iterator iter = cache_.begin();
    std::advance(iter, rand_r(&cache_seed_) % max_cache_size_);
    cache_.erase(iter);
    cache_[ids] = merged;
  }
  return true;
}

version_table_id_t VersionOperator::GetNewVersionTableId() {
  static version_table_id_t id = 0;
  return ++id;
}

