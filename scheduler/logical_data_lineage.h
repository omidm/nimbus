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

#ifndef NIMBUS_SCHEDULER_LOGICAL_DATA_LINEAGE_H_
#define NIMBUS_SCHEDULER_LOGICAL_DATA_LINEAGE_H_

#include <boost/unordered_map.hpp>
#include <list>
#include <utility>
#include "shared/nimbus_types.h"
#include "shared/idset.h"
#include "shared/dbg.h"
#include "scheduler/ldl_entry.h"
#include "scheduler/meta_before_set.h"

namespace nimbus {

  class LogicalDataLineage {
  public:
    typedef std::list<LdlEntry> Chain;
    typedef std::list<Chain::iterator> Index;

    LogicalDataLineage();
    explicit LogicalDataLineage(
        const logical_data_id_t& ldid);
    LogicalDataLineage(
        const logical_data_id_t& ldid,
        const Chain& chain,
        const Index& parents_index);

    LogicalDataLineage(
        const LogicalDataLineage& other);

    virtual ~LogicalDataLineage();

    logical_data_id_t ldid() const;
    Chain chain() const;
    const Chain* chain_p() const;
    Chain* chain_p();
    Index parents_index() const;
    const Index* parents_index_p() const;
    Index* parents_index_p();

    void set_ldid(const logical_data_id_t& ldid);
    void set_chain(const Chain& chain);
    void set_parents_index(const Index& parents_index);

    LogicalDataLineage& operator= (
        const LogicalDataLineage& right);

    bool AppendLdlEntry(
        const job_id_t& job_id,
        const data_version_t& version,
        const job_depth_t& job_depth,
        const bool& sterile);

    bool InsertCheckpointLdlEntry(
        const job_id_t& job_id,
        const data_version_t& version,
        const job_depth_t& job_depth);

    bool CleanChain(const IDSet<job_id_t>& snap_shot);


    bool LookUpVersion(
        boost::shared_ptr<MetaBeforeSet> mbs,
        data_version_t *version);

    data_version_t last_version();

  private:
    logical_data_id_t ldid_;
    Chain chain_;
    Index parents_index_;
    data_version_t last_version_;
  };

}  // namespace nimbus

#endif  // NIMBUS_SCHEDULER_LOGICAL_DATA_LINEAGE_H_
