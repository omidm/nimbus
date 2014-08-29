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


#include "scheduler/ldl_map.h"

namespace nimbus {

LdlMap::LdlMap() {
}

LdlMap::LdlMap(const Table& table)
  : table_(table) {
}

LdlMap::LdlMap(const LdlMap& other) {
  table_ = other.table_;
}

LdlMap::~LdlMap() {
}

LdlMap::Table LdlMap::table() const {
  return table_;
}

const LdlMap::Table* LdlMap::table_p() const {
  return &table_;
}

LdlMap::Table* LdlMap::table_p() {
  return &table_;
}

void LdlMap::set_table(const Table& table) {
  table_ = table;
}

LdlMap& LdlMap::operator= (
    const LdlMap& right) {
  table_ = right.table_;
  return *this;
}

bool LdlMap::DefineData(
    const logical_data_id_t ldid,
    const job_id_t& job_id,
    const job_depth_t& job_depth,
    const bool& sterile) {
  if (table_.count(ldid) == 0) {
    LogicalDataLineage* ldl = new LogicalDataLineage(ldid);
    ldl->AppendLdlEntry(job_id, NIMBUS_INIT_DATA_VERSION, job_depth, sterile);
    table_[ldid] = ldl;
    return true;
  } else {
    dbg(DBG_ERROR, "ERROR: defining logical data id %lu, which already exist.\n", ldid);
    exit(-1);
    return false;
  }
}

bool LdlMap::AppendLdlEntry(
    const logical_data_id_t ldid,
    const job_id_t& job_id,
    const data_version_t& version,
    const job_depth_t& job_depth,
    const bool& sterile) {
  Table::iterator iter;
  iter = table_.find(ldid);
  if (iter != table_.end()) {
    return iter->second->AppendLdlEntry(job_id, version, job_depth, sterile);
  } else {
    dbg(DBG_ERROR, "ERROR: ldid %lu does not exist in the ldl_map.\n");
    exit(-1);
    return false;
  }
}

bool LdlMap::InsertParentLdlEntry(
    const logical_data_id_t ldid,
    const job_id_t& job_id,
    const data_version_t& version,
    const job_depth_t& job_depth,
    const bool& sterile) {
  Table::iterator iter;
  iter = table_.find(ldid);
  if (iter != table_.end()) {
    return iter->second->InsertParentLdlEntry(job_id, version, job_depth);
  } else {
    dbg(DBG_ERROR, "ERROR: ldid %lu does not exist in the ldl_map");
    exit(-1);
    return false;
  }
}

bool LdlMap::CleanTable(
    const IDSet<job_id_t>& live_parents) {
  Table::iterator iter;
  for (iter = table_.begin(); iter != table_.end(); ++iter) {
    iter->second->CleanChain(live_parents);
  }
  return true;
}

bool LdlMap::LookUpVersion(
    logical_data_id_t ldid,
    boost::shared_ptr<MetaBeforeSet> mbs,
    data_version_t *version) {
  return table_[ldid]->LookUpVersion(mbs, version);
}

}  // namespace nimbus
