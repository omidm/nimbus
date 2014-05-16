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
  * Meta before set information used to keep track of immediate and indirect
  * before set of a job.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "scheduler/meta_before_set.h"

namespace nimbus {

MetaBeforeSet::MetaBeforeSet() {
  is_root_ = false;
}

MetaBeforeSet::MetaBeforeSet(const MetaBeforeSet& other) {
  table_ = other.table_;
  positive_query_ = other.positive_query_;
  negative_query_ = other.negative_query_;
  is_root_ = other.is_root_;
}

MetaBeforeSet::~MetaBeforeSet() {
}

MetaBeforeSet::Table MetaBeforeSet::table() const {
  return table_;
}

const MetaBeforeSet::Table* MetaBeforeSet::table_p() const {
  return &table_;
}

MetaBeforeSet::Table* MetaBeforeSet::table_p() {
  return &table_;
}

MetaBeforeSet::Cache MetaBeforeSet::positive_query() const {
  return positive_query_;
}

const MetaBeforeSet::Cache* MetaBeforeSet::positive_query_p() const {
  return &positive_query_;
}

MetaBeforeSet::Cache* MetaBeforeSet::positive_query_p() {
  return &positive_query_;
}

MetaBeforeSet::Cache MetaBeforeSet::negative_query() const {
  return negative_query_;
}

const MetaBeforeSet::Cache* MetaBeforeSet::negative_query_p() const {
  return &negative_query_;
}

MetaBeforeSet::Cache* MetaBeforeSet::negative_query_p() {
  return &negative_query_;
}

bool MetaBeforeSet::is_root() const {
  return is_root_;
}

job_depth_t MetaBeforeSet::job_depth() {
  return job_depth_;
}

void MetaBeforeSet::set_table(const Table& table) {
  table_ = table;
}

void MetaBeforeSet::set_positive_query(const Cache& positive_query) {
  positive_query_ = positive_query;
}

void MetaBeforeSet::set_negative_query(const Cache& negative_query) {
  negative_query_ = negative_query;
}

void MetaBeforeSet::set_is_root(bool flag) {
  is_root_ = flag;
}

void MetaBeforeSet::set_job_depth(job_depth_t job_depth) {
  job_depth_ = job_depth;
}

MetaBeforeSet& MetaBeforeSet::operator= (const MetaBeforeSet& right) {
  table_ = right.table_;
  positive_query_ = right.positive_query_;
  negative_query_ = right.negative_query_;
  is_root_ = right.is_root_;
  return *this;
}

bool MetaBeforeSet::LookUpBeforeSetChain(job_id_t job_id, job_depth_t job_depth) {
  if (is_root_ || (job_depth_ <= job_depth)) {
    return false;
  }

  if (table_.count(job_id) > 0) {
    return true;
  }

  if (positive_query_.count(job_id) > 0) {
    return true;
  }

  if (negative_query_.count(job_id) > 0) {
    return false;
  }

  Table::iterator iter;
  for (iter = table_.begin(); iter != table_.end(); ++iter) {
    if (iter->second->LookUpBeforeSetChain(job_id, job_depth)) {
      positive_query_.insert(job_id);
      return true;
    }
  }

  negative_query_.insert(job_id);
  return false;
}

void MetaBeforeSet::InvalidateNegativeQueryCache() {
  negative_query_.clear();
}


}  // namespace nimbus
