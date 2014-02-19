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
  * Job entry in the job table of the job manager. Each entry holds the
  * meta data of the compute and copy jobs received by the scheduler.   
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SCHEDULER_JOB_ENTRY_H_
#define NIMBUS_SCHEDULER_JOB_ENTRY_H_

#include <vector>
#include <string>
#include <set>
#include <list>
#include <utility>
#include <map>
#include "worker/data.h"
#include "shared/idset.h"
#include "shared/parameter.h"
#include "shared/nimbus_types.h"

namespace nimbus {

class JobEntry;
typedef std::map<job_id_t, JobEntry*> JobEntryMap;
typedef std::map<job_id_t, JobEntry*> JobEntryTable;
typedef std::list<JobEntry*> JobEntryList;
typedef std::vector<Data*> DataArray;

class JobEntry {
  public:
    typedef std::map<logical_data_id_t, data_version_t> VersionTable;
    typedef std::map<logical_data_id_t, physical_data_id_t> PhysicalTable;
    typedef VersionTable::iterator VTIter;

    JobEntry();
    JobEntry(const JobType& job_type,
        const std::string& job_name,
        const job_id_t& job_id,
        const IDSet<logical_data_id_t>& read_set,
        const IDSet<logical_data_id_t>& write_set,
        const IDSet<job_id_t>& before_set,
        const IDSet<job_id_t>& after_set,
        const job_id_t& parent_job_id,
        const Parameter& params,
        const bool& is_parent);

    JobEntry(const JobType& job_type,
        const std::string& job_name,
        const job_id_t& job_id,
        const job_id_t& parent_job_id,
        const bool& is_parent,
        const bool& versioned,
        const bool& assigned);

    virtual ~JobEntry();

    JobType job_type();
    std::string job_name();
    job_id_t job_id();
    IDSet<logical_data_id_t> read_set();
    IDSet<logical_data_id_t> write_set();
    IDSet<logical_data_id_t> union_set();
    IDSet<job_id_t> before_set();
    IDSet<job_id_t> after_set();
    job_id_t parent_job_id();
    Parameter params();
    VersionTable version_table_in();
    VersionTable version_table_out();
    PhysicalTable physical_table();
    bool is_parent();
    bool versioned();
    bool assigned();
    bool done();

    void set_job_name(std::string job_name);
    void set_job_id(job_id_t job_id);
    void set_before_set(IDSet<job_id_t> before_set);
    void set_after_set(IDSet<job_id_t> after_set);
    void set_version_table_in(VersionTable version_table);
    void set_version_table_out(VersionTable version_table);
    void set_physical_table(PhysicalTable physical_table);
    void set_is_parent(bool flag);
    void set_versioned(bool flag);
    void set_assigned(bool flag);
    void set_done(bool flag);

    bool GetPhysicalReadSet(IDSet<physical_data_id_t>* set);
    bool GetPhysicalWriteSet(IDSet<physical_data_id_t>* set);

  private:
    JobType job_type_;
    std::string job_name_;
    job_id_t job_id_;
    IDSet<logical_data_id_t> read_set_;
    IDSet<logical_data_id_t> write_set_;
    IDSet<logical_data_id_t> union_set_;
    IDSet<job_id_t> before_set_;
    IDSet<job_id_t> after_set_;
    job_id_t parent_job_id_;
    Parameter params_;
    VersionTable version_table_in_;
    VersionTable version_table_out_;
    PhysicalTable physical_table_;
    bool is_parent_;
    bool versioned_;
    bool assigned_;
    bool done_;
};

typedef std::map<job_id_t, JobEntry*> JobEntryTable;

class RemoteCopySendJobEntry : public JobEntry {
  public:
    explicit RemoteCopySendJobEntry() {}
    ~RemoteCopySendJobEntry() {}

    ID<job_id_t> receive_job_id() {return receive_job_id_;}
    ID<worker_id_t> to_worker_id() {return to_worker_id_;}
    std::string to_ip() {return to_ip_;}
    ID<port_t> to_port() {return to_port_;}

    void set_receive_job_id(ID<job_id_t> receive_job_id) {}
    void set_to_worker_id(ID<worker_id_t> worker_id) {}
    void set_to_ip(std::string ip) {}
    void set_to_port(ID<port_t> port) {}

  private:
    ID<job_id_t> receive_job_id_;
    ID<worker_id_t> to_worker_id_;
    std::string to_ip_;
    ID<port_t> to_port_;
};

class RemoteCopyReceiveJobEntry : public JobEntry {
  public:
    RemoteCopyReceiveJobEntry() {}
    ~RemoteCopyReceiveJobEntry() {}
};

class LocalCopyJobEntry : public JobEntry {
  public:
    LocalCopyJobEntry() {}
    ~LocalCopyJobEntry() {}
};

class CreateDataJobEntry : public JobEntry {
  public:
    CreateDataJobEntry() {}
    ~CreateDataJobEntry() {}
};

}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_JOB_ENTRY_H_


