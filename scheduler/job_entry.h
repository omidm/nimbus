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
#include "scheduler/version_map.h"
#include "scheduler/ancestor_chain.h"
#include "scheduler/data_manager.h"
#include "scheduler/meta_before_set.h"
#include "scheduler/logical_data_lineage.h"

namespace nimbus {

class JobEntry;
typedef std::map<job_id_t, JobEntry*> JobEntryMap;
typedef std::map<job_id_t, JobEntry*> JobEntryTable;
typedef std::list<JobEntry*> JobEntryList;
typedef std::vector<Data*> DataArray;

class JobEntry {
  public:
    typedef std::map<logical_data_id_t, physical_data_id_t> PhysicalTable;

    JobEntry();

    explicit JobEntry(const job_id_t& job_id);

    virtual ~JobEntry();

    JobType job_type() const;
    std::string job_name() const;
    job_id_t job_id() const;
    IDSet<logical_data_id_t> read_set() const;
    IDSet<logical_data_id_t> write_set() const;
    IDSet<logical_data_id_t> union_set() const;
    IDSet<job_id_t> before_set() const;
    IDSet<job_id_t> after_set() const;
    job_id_t parent_job_id() const;
    Parameter params() const;
    boost::shared_ptr<VersionMap> vmap_read() const;
    boost::shared_ptr<VersionMap> vmap_write() const;
    boost::shared_ptr<VersionMap> vmap_partial() const;
    boost::shared_ptr<MetaBeforeSet> meta_before_set() const;
    job_depth_t job_depth() const;
    PhysicalTable physical_table() const;
    IDSet<job_id_t> jobs_passed_versions() const;
    IDSet<job_id_t> need_set() const;
    worker_id_t assigned_worker() const;
    bool sterile() const;
    bool partial_versioned() const;
    bool versioned() const;
    bool assigned() const;
    bool done() const;
    bool future() const;

    const IDSet<logical_data_id_t>* read_set_p() const;
    const IDSet<logical_data_id_t>* write_set_p() const;
    const IDSet<logical_data_id_t>* union_set_p() const;
    const IDSet<job_id_t>* before_set_p() const;
    IDSet<job_id_t>* before_set_p();

    void set_job_type(JobType job_type);
    void set_job_name(std::string job_name);
    void set_job_id(job_id_t job_id);
    void set_read_set(IDSet<logical_data_id_t> read_set);
    void set_write_set(IDSet<logical_data_id_t> write_set);
    void set_before_set(IDSet<job_id_t> before_set, bool update_dependencies = false);
    void set_after_set(IDSet<job_id_t> after_set);
    void set_parent_job_id(job_id_t parent_job_id, bool update_dependencies = false);
    void set_params(Parameter params);
    void set_vmap_read(boost::shared_ptr<VersionMap> vmap_read);
    void set_vmap_write(boost::shared_ptr<VersionMap> vmap_write);
    void set_vmap_partial(boost::shared_ptr<VersionMap> vmap_partial);
    void set_meta_before_set(boost::shared_ptr<MetaBeforeSet> meta_before_set);
    void set_job_depth(job_depth_t job_depth);
    void set_physical_table(PhysicalTable physical_table);
    void set_jobs_passed_versions(IDSet<job_id_t> jobs);
    void add_job_passed_versions(job_id_t job_id);
    void set_assigned_worker(worker_id_t worker_id);
    void set_sterile(bool flag);
    void set_partial_versioned(bool flag);
    void set_versioned(bool flag);
    void set_assigned(bool flag);
    void set_done(bool flag);
    void set_future(bool flag);

    void set_physical_table_entry(logical_data_id_t l_id, physical_data_id_t p_id);

    bool GetPhysicalReadSet(IDSet<physical_data_id_t>* set);
    bool GetPhysicalWriteSet(IDSet<physical_data_id_t>* set);

    bool IsReadyToAssign();
    void remove_assignment_dependency(job_id_t job_id);

    bool IsReadyForCompleteVersioning();
    void remove_versioning_dependency(job_id_t job_id);

    bool GetRegion(DataManager *data_manager, GeometricRegion *region);

  protected:
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
    boost::shared_ptr<VersionMap> vmap_read_;
    boost::shared_ptr<VersionMap> vmap_write_;
    boost::shared_ptr<VersionMap> vmap_partial_;
    boost::shared_ptr<MetaBeforeSet> meta_before_set_;
    job_depth_t job_depth_;
    PhysicalTable physical_table_;
    IDSet<job_id_t> jobs_passed_versions_;
    IDSet<job_id_t> assignment_dependencies_;
    IDSet<job_id_t> versioning_dependencies_;
    GeometricRegion region_;
    worker_id_t assigned_worker_;
    bool sterile_;
    bool partial_versioned_;
    bool versioned_;
    bool assigned_;
    bool done_;
    bool future_;
    bool region_valid_;

  private:
    void Initialize();
};

typedef std::map<job_id_t, JobEntry*> JobEntryTable;


class ComputeJobEntry : public JobEntry {
  public:
    ComputeJobEntry(
        const std::string& job_name,
        const job_id_t& job_id,
        const IDSet<logical_data_id_t>& read_set,
        const IDSet<logical_data_id_t>& write_set,
        const IDSet<job_id_t>& before_set,
        const IDSet<job_id_t>& after_set,
        const job_id_t& parent_job_id,
        const Parameter& params,
        const bool& sterile);
    ~ComputeJobEntry();
};


class KernelJobEntry : public JobEntry {
  public:
    KernelJobEntry();
    ~KernelJobEntry();
};


class MainJobEntry : public JobEntry {
  public:
    explicit MainJobEntry(const job_id_t& job_id);
    ~MainJobEntry();
};


class FutureJobEntry : public JobEntry {
  public:
    explicit FutureJobEntry(const job_id_t& job_id);
    ~FutureJobEntry();
};


class LocalCopyJobEntry : public JobEntry {
  public:
    explicit LocalCopyJobEntry(const job_id_t& job_id);
    ~LocalCopyJobEntry();
};


class CreateDataJobEntry : public JobEntry {
  public:
    explicit CreateDataJobEntry(const job_id_t& job_id);
    ~CreateDataJobEntry();
};


class RemoteCopyReceiveJobEntry : public JobEntry {
  public:
    explicit RemoteCopyReceiveJobEntry(const job_id_t& job_id);
    ~RemoteCopyReceiveJobEntry();
};


class RemoteCopySendJobEntry : public JobEntry {
  public:
    explicit RemoteCopySendJobEntry(const job_id_t& job_id);
    ~RemoteCopySendJobEntry();

    ID<job_id_t> receive_job_id();
    ID<worker_id_t> to_worker_id();
    std::string to_ip();
    ID<port_t> to_port();

    void set_receive_job_id(ID<job_id_t> receive_job_id);
    void set_to_worker_id(ID<worker_id_t> worker_id);
    void set_to_ip(std::string ip);
    void set_to_port(ID<port_t> port);

  private:
    ID<job_id_t> receive_job_id_;
    ID<worker_id_t> to_worker_id_;
    std::string to_ip_;
    ID<port_t> to_port_;
};

}  // namespace nimbus
#endif  // NIMBUS_SCHEDULER_JOB_ENTRY_H_


