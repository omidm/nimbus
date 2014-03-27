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
    std::vector<boost::shared_ptr<const VersionTable> > tables,
    boost::shared_ptr<VersionTable> *result) {
  size_t count = tables.size();

  if (count == 0) {
    return false;
  }

  if (count == 1) {
    *result = boost::shared_ptr<VersionTable>(new VersionTable(GetNewVersionTableId()));
    (*result)->set_root(tables[0]->root());
    (*result)->set_content(tables[0]->content());
    return true;
  }

  std::vector<boost::shared_ptr<const VersionTable> > reduced;
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
    boost::shared_ptr<const VersionTable> t_1,
    boost::shared_ptr<const VersionTable> t_2,
    boost::shared_ptr<VersionTable> *result) {
  std::set<version_table_id_t> ids;
  ids.insert(t_1->id());
  ids.insert(t_2->id());
  if (LookUpCache(ids, result)) {
    dbg(DBG_SCHED, "Version Operator: hit cache.\n");
    return true;
  } else {
    dbg(DBG_SCHED, "Version Operator: missed cache.\n");
    boost::shared_ptr<VersionTable> merged(new VersionTable(GetNewVersionTableId()));
    boost::shared_ptr<const VersionTable::Map> root;
    VersionTable::Map content;
    VersionTable::MapIter it;
    VersionTable::MapConstIter it_p;

    if (t_1->root() == t_2->root()) {
      dbg(DBG_SCHED, "Version Operator: roots are the same for merge.\n");
      root = t_1->root();
      content = t_1->content();
      for (it_p = t_2->content_p()->begin(); it_p != t_2->content_p()->end(); ++it_p) {
        it = content.find(it_p->first);
        if (it == content.end()) {
          content[it_p->first] = it_p->second;
        } else if (it_p->second > it->second) {
          it->second = it_p->second;
        }
      }
    } else {
      dbg(DBG_SCHED, "Version Operator: roots are different for merge.\n");
      boost::shared_ptr<const VersionTable> t_r;
      boost::shared_ptr<VersionTable::Map> content_p;
      if (CompareRootDominance(t_1->root(), t_2->root())) {
        t_r = t_1;
        FlatenVersionTable(t_2, &content_p);
      } else {
        t_r = t_2;
        FlatenVersionTable(t_1, &content_p);
      }

      root = t_r->root();
      content = t_r->content();

      for (it_p = content_p->begin(); it_p != content_p->end(); ++it_p) {
        data_version_t version;
        if (t_r->query_entry(it_p->first, &version)) {
          if (it_p->second > version) {
            content[it_p->first] = it_p->second;
          }
        } else {
          content[it_p->first] = it_p->second;
        }
      }
    }

    merged->set_content(content);
    merged->set_root(root);
    CacheMergedResult(ids, merged);
    *result = merged;
    return true;
  }
}

bool VersionOperator::MakeVersionTableOut(
    boost::shared_ptr<const VersionTable> table_in,
    const IDSet<logical_data_id_t>& write_set,
    boost::shared_ptr<VersionTable> *table_out) {
  boost::shared_ptr<VersionTable> result(new VersionTable(GetNewVersionTableId()));
  result->set_root(table_in->root());
  result->set_content(table_in->content());
  IDSet<logical_data_id_t>::ConstIter iter;
  for (iter = write_set.begin(); iter != write_set.end(); ++iter) {
    data_version_t version;
    if (result->query_entry(*iter, &version)) {
      result->set_entry(*iter, ++version);
    } else {
      dbg(DBG_ERROR, "ERROR: Version Operator, writing to unknown logical id in the context.\n");
      exit(-1);
      return false;
    }
  }
  *table_out = result;
  return true;
}

bool VersionOperator::RecomputeRootForVersionTable(
    std::vector<boost::shared_ptr<VersionTable> > tables) {
  // TODO(omidm): implement!
  FlushCache();
  return false;
}

bool VersionOperator::CompareRootDominance(
    boost::shared_ptr<const VersionTable::Map> r1,
    boost::shared_ptr<const VersionTable::Map> r2) {
  VersionTable::MapConstIter iter_1;
  VersionTable::MapConstIter iter_2;

  int count_1 = 0;
  for (iter_1 = r1->begin(); iter_1 != r1->end(); ++iter_1) {
    ++count_1;
    iter_2 = r2->find(iter_1->first);
    if (iter_2 != r2->end()) {
      if (iter_2->second > iter_1->second) {
        --count_1;
      }
    }
  }

  int count_2 = 0;
  for (iter_2 = r2->begin(); iter_2 != r2->end(); ++iter_2) {
    ++count_2;
    iter_1 = r1->find(iter_2->first);
    if (iter_1 != r1->end()) {
      if (iter_1->second > iter_2->second) {
        --count_2;
      }
    }
  }

  return (count_1 >= count_2);
}

bool FlatenVersionTable(
    boost::shared_ptr<const VersionTable> table,
    boost::shared_ptr<VersionTable::Map> *result) {
  VersionTable::MapConstIter iter;
  boost::shared_ptr<VersionTable::Map> content(new VersionTable::Map());

  for (iter = table->root()->begin(); iter != table->root()->end(); ++iter) {
    content->operator[](iter->first) = iter->second;
  }
  for (iter = table->content_p()->begin(); iter != table->content_p()->end(); ++iter) {
    content->operator[](iter->first) = iter->second;
  }

  *result = content;
  return true;
}

bool VersionOperator::LookUpCache(
    const std::set<version_table_id_t>& ids,
    boost::shared_ptr<VersionTable>* result) {
  Cache::iterator iter = cache_.find(ids);
  if (iter != cache_.end()) {
    *result = iter->second;
    return true;
  }
  return false;
}

bool VersionOperator::CacheMergedResult(
    const std::set<version_table_id_t>& ids,
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

void VersionOperator::FlushCache() {
  cache_.clear();
}

version_table_id_t VersionOperator::GetNewVersionTableId() {
  static version_table_id_t id = 0;
  return ++id;
}

