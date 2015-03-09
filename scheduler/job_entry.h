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

#include <boost/unordered_map.hpp>
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
#include "shared/geometric_region.h"
#include "scheduler/version_map.h"
#include "scheduler/data_manager.h"
#include "scheduler/binding_template.h"
#include "scheduler/scheduler_worker.h"
#include "scheduler/meta_before_set.h"
#include "scheduler/logical_data_lineage.h"

namespace nimbus {

class JobEntry;
typedef boost::unordered_map<job_id_t, JobEntry*> JobEntryMap;
typedef boost::unordered_map<job_id_t, JobEntry*> JobEntryTable;
typedef std::list<JobEntry*> JobEntryList;
typedef std::vector<Data*> DataArray;

class TemplateJobEntry;

class JobEntry {
  public:
    typedef boost::unordered_map<logical_data_id_t, physical_data_id_t> PhysicalTable;

    JobEntry();

    explicit JobEntry(const job_id_t& job_id);

    virtual ~JobEntry();

    virtual JobType job_type() const;
    virtual std::string job_name() const;
    virtual std::string parent_job_name() const;
    virtual std::string grand_parent_job_name() const;
    virtual job_id_t job_id() const;
    virtual IDSet<logical_data_id_t> read_set() const;
    virtual IDSet<logical_data_id_t> write_set() const;
    virtual IDSet<logical_data_id_t> union_set() const;
    virtual IDSet<job_id_t> before_set() const;
    virtual IDSet<job_id_t> after_set() const;
    virtual job_id_t parent_job_id() const;
    virtual job_id_t future_job_id() const;
    virtual Parameter params() const;
    virtual boost::shared_ptr<VersionMap> vmap_read() const;
    virtual boost::shared_ptr<VersionMap> vmap_write() const;
    virtual boost::shared_ptr<VersionMap> vmap_partial() const;
    virtual boost::shared_ptr<MetaBeforeSet> meta_before_set() const;
    virtual job_depth_t job_depth() const;
    virtual PhysicalTable physical_table() const;
    virtual IDSet<job_id_t> jobs_passed_versions() const;
    virtual IDSet<job_id_t> need_set() const;
    virtual worker_id_t assigned_worker_id() const;
    virtual SchedulerWorker* assigned_worker() const;
    virtual bool sterile() const;
    virtual bool memoize() const;
    virtual bool memoize_binding() const;
    virtual bool to_finalize_binding_template() const;
    virtual TemplateJobEntry* template_job() const;
    virtual BindingTemplate* binding_template() const;
    virtual GeometricRegion region() const;
    virtual bool partial_versioned() const;
    virtual bool versioned() const;
    virtual bool versioned_for_pattern() const;
    virtual bool versioned_entire_context() const;
    virtual bool assigned() const;
    virtual bool done() const;
    virtual bool future() const;
    virtual checkpoint_id_t checkpoint_id() const;

    virtual const IDSet<logical_data_id_t>* read_set_p() const;
    virtual const IDSet<logical_data_id_t>* write_set_p() const;
    virtual const IDSet<logical_data_id_t>* union_set_p() const;
    virtual const IDSet<job_id_t>* before_set_p() const;
    virtual IDSet<job_id_t>* before_set_p();

    virtual void set_job_type(JobType job_type);
    virtual void set_job_name(std::string job_name);
    virtual void set_parent_job_name(std::string parent_job_name);
    virtual void set_grand_parent_job_name(std::string grand_parent_job_name);
    virtual void set_job_id(job_id_t job_id);
    virtual void set_read_set(IDSet<logical_data_id_t> read_set);
    virtual void set_write_set(IDSet<logical_data_id_t> write_set);
    virtual void set_before_set(IDSet<job_id_t> before_set, bool update_dependencies = false);
    virtual void set_after_set(IDSet<job_id_t> after_set);
    virtual void set_parent_job_id(job_id_t parent_job_id, bool update_dependencies = false);
    virtual void set_future_job_id(job_id_t future_job_id);
    virtual void set_params(Parameter params);
    virtual void set_vmap_read(boost::shared_ptr<VersionMap> vmap_read);
    virtual void set_vmap_write(boost::shared_ptr<VersionMap> vmap_write);
    virtual void set_vmap_partial(boost::shared_ptr<VersionMap> vmap_partial);
    virtual void set_meta_before_set(boost::shared_ptr<MetaBeforeSet> meta_before_set);
    virtual void set_job_depth(job_depth_t job_depth);
    virtual void set_physical_table(PhysicalTable physical_table);
    virtual void set_jobs_passed_versions(IDSet<job_id_t> jobs);
    virtual void add_job_passed_versions(job_id_t job_id);
    virtual void set_assigned_worker_id(worker_id_t assigned_worker_id);
    virtual void set_assigned_worker(SchedulerWorker* assigned_worker);
    virtual void set_sterile(bool flag);
    virtual void set_memoize(bool flag);
    virtual void set_memoize_binding(bool flag);
    virtual void set_to_finalize_binding_template(bool flag);
    virtual void set_template_job(TemplateJobEntry* template_job);
    virtual void set_binding_template(BindingTemplate* binding_template);
    virtual void set_region(GeometricRegion region);
    virtual void set_partial_versioned(bool flag);
    virtual void set_versioned(bool flag);
    virtual void set_versioned_for_pattern(bool flag);
    virtual void set_versioned_entire_context(bool flag);
    virtual void set_assigned(bool flag);
    virtual void set_done(bool flag);
    virtual void set_future(bool flag);
    virtual void set_checkpoint_id(checkpoint_id_t checkpoint_id);

    virtual void set_physical_table_entry(logical_data_id_t l_id, physical_data_id_t p_id);

    virtual bool GetPhysicalReadSet(IDSet<physical_data_id_t>* set);
    virtual bool GetPhysicalWriteSet(IDSet<physical_data_id_t>* set);

    virtual bool IsReadyToAssign();
    virtual void remove_assignment_dependency(job_id_t job_id);

    virtual bool IsReadyForCompleteVersioning();
    virtual void remove_versioning_dependency(job_id_t job_id);

    virtual void MarkJobAsCompletelyResolved();

    virtual bool LookUpMetaBeforeSet(JobEntry* job);

    virtual bool GetRegion(GeometricRegion *region);
    virtual bool GetReadSetRegion(DataManager *data_manager, GeometricRegion *region);
    virtual bool GetWriteSetRegion(DataManager *data_manager, GeometricRegion *region);
    virtual bool GetUnionSetRegion(DataManager *data_manager, GeometricRegion *region);

  protected:
    JobType job_type_;
    std::string job_name_;
    std::string parent_job_name_;
    std::string grand_parent_job_name_;
    job_id_t job_id_;
    IDSet<logical_data_id_t> read_set_;
    IDSet<logical_data_id_t> write_set_;
    IDSet<logical_data_id_t> union_set_;
    IDSet<job_id_t> before_set_;
    IDSet<job_id_t> after_set_;
    job_id_t parent_job_id_;
    job_id_t future_job_id_;
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
    GeometricRegion read_region_;
    GeometricRegion write_region_;
    GeometricRegion union_region_;
    worker_id_t assigned_worker_id_;
    SchedulerWorker *assigned_worker_;
    bool sterile_;
    bool memoize_;
    bool memoize_binding_;
    bool to_finalize_binding_template_;
    TemplateJobEntry* template_job_;
    BindingTemplate* binding_template_;
    GeometricRegion region_;
    bool partial_versioned_;
    bool versioned_;
    bool versioned_for_pattern_;
    bool versioned_entire_context_;
    bool assigned_;
    bool done_;
    bool future_;
    bool read_region_valid_;
    bool write_region_valid_;
    bool union_region_valid_;
    checkpoint_id_t checkpoint_id_;

  private:
    void Initialize();
};


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
        const job_id_t& future_job_id,
        const bool& sterile,
        const GeometricRegion& region,
        const Parameter& params);
    ~ComputeJobEntry();
};

class SaveDataJobEntry : public JobEntry {
  public:
    SaveDataJobEntry(
        const job_id_t& job_id,
        const IDSet<logical_data_id_t>& read_set);
    ~SaveDataJobEntry();
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


