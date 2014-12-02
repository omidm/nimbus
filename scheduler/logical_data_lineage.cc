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
  * This class keeps the lineage information about how logical data evolves as
  * it is written by jobs.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */


#include "scheduler/logical_data_lineage.h"

namespace nimbus {

LogicalDataLineage::LogicalDataLineage() {
  last_version_ = NIMBUS_UNDEFINED_DATA_VERSION;
}

LogicalDataLineage::LogicalDataLineage(
    const logical_data_id_t& ldid) :
  ldid_(ldid) {
  last_version_ = NIMBUS_UNDEFINED_DATA_VERSION;
}

LogicalDataLineage::LogicalDataLineage(
    const logical_data_id_t& ldid,
    const Chain& chain,
    const Index& parents_index) :
  ldid_(ldid), chain_(chain),
  parents_index_(parents_index) {
  last_version_ = NIMBUS_UNDEFINED_DATA_VERSION;
}

LogicalDataLineage::LogicalDataLineage(
    const LogicalDataLineage& other) {
  ldid_ = other.ldid_;
  chain_ = other.chain_;
  parents_index_ = other.parents_index_;
  last_version_ = other.last_version_;
}

void LogicalDataLineage::Reinitialize() {
  chain_.clear();
  parents_index_.clear();
  last_version_ = NIMBUS_UNDEFINED_DATA_VERSION;
}


LogicalDataLineage::~LogicalDataLineage() {
}

logical_data_id_t LogicalDataLineage::ldid() const {
  return ldid_;
}

LogicalDataLineage::Chain LogicalDataLineage::chain() const {
  return chain_;
}

const LogicalDataLineage::Chain* LogicalDataLineage::chain_p() const {
  return &chain_;
}

LogicalDataLineage::Chain* LogicalDataLineage::chain_p() {
  return &chain_;
}

LogicalDataLineage::Index LogicalDataLineage::parents_index() const {
  return parents_index_;
}

const LogicalDataLineage::Index* LogicalDataLineage::parents_index_p() const {
  return &parents_index_;
}

LogicalDataLineage::Index* LogicalDataLineage::parents_index_p() {
  return &parents_index_;
}

void LogicalDataLineage::set_ldid(const logical_data_id_t& ldid) {
  ldid_ = ldid;
}

void LogicalDataLineage::set_chain(const Chain& chain) {
  chain_ = chain;
}

void LogicalDataLineage::set_parents_index(const Index& parents_index) {
  parents_index_ = parents_index;
}

LogicalDataLineage& LogicalDataLineage::operator= (
    const LogicalDataLineage& right) {
  ldid_ = right.ldid_;
  chain_ = right.chain_;
  parents_index_ = right.parents_index_;
  last_version_ = right.last_version_;
  return *this;
}

bool LogicalDataLineage::AppendLdlEntry(
    const job_id_t& job_id,
    const data_version_t& version,
    const job_depth_t& job_depth,
    const bool& sterile) {
  assert((last_version_ + 1) == version);
  last_version_ = version;

  chain_.push_back(LdlEntry(job_id, version, job_depth, sterile));

//  if (!sterile) {
//    parents_index_.push_back(--(chain_.end()));
//  }

  return true;
}

bool LogicalDataLineage::InsertCheckpointLdlEntry(
    const job_id_t& job_id,
    const data_version_t& version,
    const job_depth_t& job_depth) {
  if (version > last_version_) {
    last_version_ = version;
  }

  Chain::reverse_iterator it = chain_.rbegin();
  for (; it != chain_.rend(); ++it) {
    if (it->version() <= version) {
      break;
    }
  }

  chain_.insert(it.base(), LdlEntry(job_id, version, job_depth, false));

//  Index::reverse_iterator iit = parents_index_.rbegin();
//  for (; iit != parents_index_.rend(); ++iit) {
//    if ((*iit)->version() <= version) {
//      break;
//    }
//  }
//
//  parents_index_.insert(iit.base(), --it.base());
//
//  // TODO(omidm) : remove this check!
//  for (iit = parents_index_.rbegin(); iit != parents_index_.rend(); ++iit) {
//    assert(!(*iit)->sterile());
//  }

  return true;
}

bool LogicalDataLineage::CleanChain(const IDSet<job_id_t>& snap_shot) {
  if (snap_shot.size() == 0) {
    return true;
  }

  bool found = false;
  Chain::iterator it = chain_.begin();
  for (; it != chain_.end(); ++it) {
    if (snap_shot.contains(it->job_id())) {
      found = true;
      break;
    }
  }

  assert(found);
  chain_.erase(chain_.begin(), it);

//  Index::reverse_iterator iit = parents_index_.rbegin();
//  for (; iit != parents_index_.rend(); ++iit) {
//    temp.remove((*iit)->job_id());
//    if (temp.size() == 0) {
//      break;
//    }
//  }
//
//  assert(temp.size() == 0);
//
//  parents_index_.erase(parents_index_.begin(), --(iit.base()));
//  chain_.erase(chain_.begin(), *iit);

  return true;
}

data_version_t LogicalDataLineage::last_version() {
  return last_version_;
}

bool LogicalDataLineage::LookUpVersion(
    boost::shared_ptr<MetaBeforeSet> mbs,
    data_version_t *version) {
  Chain::const_reverse_iterator iter;
  for (iter = chain_.rbegin(); iter != chain_.rend(); ++iter) {
    if (mbs->LookUpBeforeSetChain(iter->job_id(), iter->job_depth())) {
      *version = iter->version();
      return true;
    }
  }
  return false;
}

}  // namespace nimbus
