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

#ifndef NIMBUS_SCHEDULER_VERSION_OPERATOR_H_
#define NIMBUS_SCHEDULER_VERSION_OPERATOR_H_

#include <boost/shared_ptr.hpp>
#include <map>
#include <set>
#include <vector>
#include "shared/nimbus_types.h"
#include "scheduler/version_table.h"
#include "scheduler/job_entry.h"
#include "shared/dbg.h"

namespace nimbus {

class VersionOperator {
  public:
    typedef std::map<std::set<version_table_id_t>, boost::shared_ptr<VersionTable> > Cache;
    typedef std::map<version_table_id_t, size_t> RankTable;

    VersionOperator();
    virtual ~VersionOperator();

    bool MergeVersionTables(
        std::vector<boost::shared_ptr<VersionTable> > tables,
        boost::shared_ptr<VersionTable> *result);

    bool MergeTwoVersionTables(
        boost::shared_ptr<VersionTable> first,
        boost::shared_ptr<VersionTable> second,
        boost::shared_ptr<VersionTable> *result);

    bool MakeRootTable(
        boost::shared_ptr<VersionTable> table,
        boost::shared_ptr<VersionTable> *result);

    bool MakeVersionTableOut(
        boost::shared_ptr<VersionTable> tabel_in,
        IDSet<logical_data_id_t> write_set,
        boost::shared_ptr<VersionTable> *teble_out);

    bool LookUpCache(std::set<version_table_id_t> ids,
        boost::shared_ptr<VersionTable>* result);

    bool CacheMergeResult(
        std::set<version_table_id_t> ids,
        boost::shared_ptr<VersionTable> merged);

    version_table_id_t GetNewVersionTableId();

  private:
    Cache cache_;
    size_t cache_size_;
    RankTable ranks_;
};


}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_VERSION_OPERATOR_H_
