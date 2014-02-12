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
  * Scheduler Version Entry. It holds the meta data for each version of the
  * data containing logical data id, version number, and a pointer to a job
  * related to the version, and the type of relation meaning whether job needs
  * the version as input or dumps it as output.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_VERSION_ENTRY_H_
#define NIMBUS_SCHEDULER_VERSION_ENTRY_H_

#include <list>
#include "shared/nimbus_types.h"
#include "shared/dbg.h"
#include "scheduler/job_entry.h"

namespace nimbus {

class VersionEntry {
  public:
    enum Relation {IN, OUT};

    VersionEntry();
    VersionEntry(logical_data_id_t logical_id, data_version_t version,
    JobEntry* job_entry, Relation relation);
    VersionEntry(const VersionEntry& ve);

    virtual ~VersionEntry();

    logical_data_id_t logical_id() const;
    data_version_t version() const;
    JobEntry* job_entry() const;
    Relation relation() const;

    VersionEntry& operator= (const VersionEntry& right);

  private:
    logical_data_id_t logical_id_;
    data_version_t version_;
    JobEntry* job_entry_;
    Relation relation_;
};

typedef std::list<VersionEntry> VersionEntryList;

}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_VERSION_ENTRY_H_
