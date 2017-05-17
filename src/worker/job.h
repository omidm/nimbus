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
  * A Nimbus job. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SRC_WORKER_JOB_H_
#define NIMBUS_SRC_WORKER_JOB_H_

#include <vector>
#include <string>
#include <set>
#include <map>
#include <list>
#include "src/shared/helpers.h"
#include "src/shared/nimbus_types.h"
#include "src/shared/id.h"
#include "src/shared/idset.h"
#include "src/shared/geometric_region.h"
#include "src/shared/serialized_data.h"
#include "src/shared/worker_data_exchanger.h"
#include "src/shared/distributed_db.h"
#include "src/worker/app_data_manager.h"
#include "src/worker/data.h"
#include "src/worker/ldo_index_cache.h"
#include "src/worker/dependency_query.h"
#include "src/worker/worker_ldo_map.h"

namespace nimbus {

class ExecutionTemplate;
class Application;
class Job;
class WorkerThread;
class StaticConfigManager;
typedef std::list<Job*> JobList;
typedef std::map<job_id_t, Job*> JobMap;
typedef std::map<std::string, Job*> JobTable;

class Job {
  public:
    Job();
    virtual ~Job();

    // TODO(omidm) should remove this later. left it now so the tests
    // that use it still pass.
    Job(Application* app, JobType type);

    virtual void Execute(Parameter params, const DataArray& da) {}
    virtual Job* Clone();
    virtual void Sleep() {}
    virtual void Cancel() {}

    enum SpawnState {
      INIT,
      NORMAL,
      START_KNOWN_TEMPLATE,
      START_UNKNOWN_TEMPLATE,
      END_TEMPLATE
    };

    bool SpawnComputeJob(const std::string& name,
                         const job_id_t& id,
                         const IDSet<logical_data_id_t>& read,
                         const IDSet<logical_data_id_t>& write,
                         const IDSet<job_id_t>& before,
                         const IDSet<job_id_t>& after,
                         const Parameter& params,
                         const bool& sterile = false,
                         const GeometricRegion& region = GeometricRegion(),
                         const job_id_t& future_job_id = 0);

    bool SpawnComputeJob(const std::string& name,
                         const job_id_t& id,
                         const IDSet<logical_data_id_t>& read,
                         const IDSet<logical_data_id_t>& write,
                         const IDSet<logical_data_id_t>& scratch,
                         const IDSet<logical_data_id_t>& reduce,
                         const IDSet<job_id_t>& before,
                         const IDSet<job_id_t>& after,
                         const Parameter& params,
                         const bool& sterile = false,
                         const GeometricRegion& region = GeometricRegion(),
                         const job_id_t& future_job_id = 0);

    bool SpawnComputeJob(const std::string& name,
                         const job_id_t& id,
                         const IDSet<logical_data_id_t>& read,
                         const IDSet<logical_data_id_t>& write,
                         const IDSet<logical_data_id_t>& scratch,
                         const IDSet<logical_data_id_t>& reduce,
                         const IDSet<job_id_t>& before,
                         const IDSet<job_id_t>& after,
                         const Parameter& params,
                         const std::string& combiner,
                         const bool& sterile = false,
                         const GeometricRegion& region = GeometricRegion(),
                         const job_id_t& future_job_id = 0);

    bool SpawnCopyJob(const job_id_t& id,
                      const logical_data_id_t& from_logical_id,
                      const logical_data_id_t& to_logical_id,
                      const IDSet<job_id_t>& before,
                      const IDSet<job_id_t>& after,
                      const Parameter& params);

    bool DefineData(const std::string& name,
                    const logical_data_id_t& logical_data_id,
                    const partition_id_t& partition_id,
                    const IDSet<partition_id_t>& neighbor_partition);

    bool DefinePartition(const ID<partition_id_t>& partition_id,
                         const GeometricRegion& r);

    bool DefineJobGraph(const std::string& job_graph_name);

    bool TerminateApplication(const exit_status_t& exit_status_id = NIMBUS_TERMINATE_SUCCESS);

    bool GetNewJobID(std::vector<job_id_t>* result, size_t req_num);

    bool GetNewLogicalDataID(std::vector<logical_data_id_t>* result, size_t req_num);

    bool GetPartition(partition_id_t id, GeometricRegion* r) const;

    const LogicalDataObject* GetLogicalObject(logical_data_id_t id) const;

    int GetCoveredLogicalObjects(CLdoVector* result,
                                 const std::string& variable,
                                 const GeometricRegion* r);

    int GetAdjacentLogicalObjects(CLdoVector* result,
                                  const std::string& variable,
                                  const GeometricRegion* r);

    int AddIntersectingLdoIds(const std::string& variable,
                              const nimbus::GeometricRegion& region,
                              IDSet<logical_data_id_t>* result);

    int GetIntersectingLogicalObjects(CLdoVector* result,
                                      const std::string& variable,
                                      const GeometricRegion* r);

    void LoadLdoIdsInSet(IDSet<logical_data_id_t>* set,
                         const nimbus::GeometricRegion& region,
                         ...);

    bool StageJobAndLoadBeforeSet(IDSet<job_id_t> *before_set,
                                  const std::string& name,
                                  const job_id_t& id,
                                  const IDSet<logical_data_id_t>& read,
                                  const IDSet<logical_data_id_t>& write,
                                  const bool barrier = false);

    bool StageJobAndLoadBeforeSet(IDSet<job_id_t> *before_set,
                                  const std::string& name,
                                  const job_id_t& id,
                                  const IDSet<logical_data_id_t>& read,
                                  const IDSet<logical_data_id_t>& write,
                                  const IDSet<logical_data_id_t>& scratch,
                                  const IDSet<logical_data_id_t>& reduce,
                                  const bool barrier = false);

    bool MarkEndOfStage();

    void StartTemplate(const std::string& template_name);

    void EndTemplate(const std::string& template_name);

    bool IsTemplateDefined(const std::string& template_name);

    std::string name() const;
    ID<job_id_t> id() const;
    IDSet<physical_data_id_t> read_set() const;
    IDSet<physical_data_id_t> write_set() const;
    IDSet<physical_data_id_t> scratch_set() const;
    IDSet<physical_data_id_t> reduce_set() const;
    const IDSet<physical_data_id_t>& get_read_set() const;
    const IDSet<physical_data_id_t>& get_write_set() const;
    const IDSet<physical_data_id_t>& get_scratch_set() const;
    const IDSet<physical_data_id_t>& get_reduce_set() const;
    IDSet<job_id_t> before_set() const;
    const IDSet<job_id_t>* before_set_p() const;
    IDSet<job_id_t> after_set() const;
    ID<job_id_t> future_job_id() const;
    job_id_t shadow_job_id() const;
    ExecutionTemplate* execution_template() const;

    Parameter parameters() const;
    Application* application() const;
    bool  sterile() const;
    GeometricRegion region() const;
    double run_time() const;
    double wait_time() const;
    size_t max_alloc() const {
      return max_alloc_;
    }

    void set_name(const std::string& name);
    void set_id(const ID<job_id_t>& id);
    void set_read_set(const IDSet<physical_data_id_t>& read_set);
    void set_write_set(const IDSet<physical_data_id_t>& write_set);
    void set_scratch_set(const IDSet<physical_data_id_t>& scratch_set);
    void set_reduce_set(const IDSet<physical_data_id_t>& reduce_set);
    void set_before_set(const IDSet<job_id_t>& before_set);
    void set_after_set(const IDSet<job_id_t>& after_set);
    void set_parameters(const Parameter& parameters);
    void set_application(Application* app);
    void set_sterile(const bool& sterile);
    void set_region(const GeometricRegion& region);
    void set_future_job_id(const ID<job_id_t>& future_job_id);
    void set_shadow_job_id(const job_id_t& id);
    void set_execution_template(ExecutionTemplate *execution_template);
    void clear_template_variables();
    void set_run_time(const double& run_time);
    void set_wait_time(const double& wait_time);
    void set_max_alloc(const size_t& max_alloc) {
      max_alloc_ = max_alloc;
    }
    // TODO(quhang) should add accesssors.
    DataArray data_array;

    void refresh_dependency_query();

    AppDataManager* GetAppDataManager() const;
    StaticConfigManager* GetStaticConfigManager() const;

    void set_core_quota(int core_quota) {
      core_quota_ = core_quota;
    }
    int core_quota() {
      return core_quota_;
    }
    void set_use_threading(bool use_threading) {
      use_threading_ = use_threading;
    }
    bool use_threading() {
      return use_threading_;
    }
    virtual bool SupportMultiThread() const {
      return false;
    }
    WorkerThread* worker_thread() const {
      return worker_thread_;
    }
    void set_worker_thread(WorkerThread* worker_thread) {
      worker_thread_ = worker_thread;
    }

  private:
    std::map<std::string, LdoIndexCache> query_cache_;
    int core_quota_;
    bool use_threading_;
    std::string name_;
    ID<job_id_t> id_;
    IDSet<physical_data_id_t> read_set_;
    IDSet<physical_data_id_t> write_set_;
    IDSet<physical_data_id_t> scratch_set_;
    IDSet<physical_data_id_t> reduce_set_;
    IDSet<job_id_t> before_set_;
    IDSet<job_id_t> after_set_;
    ID<job_id_t> future_job_id_;
    Parameter parameters_;
    Application* application_;
    bool sterile_;
    GeometricRegion region_;
    bool app_is_set_;
    double run_time_;
    double wait_time_;
    size_t max_alloc_;
    WorkerThread* worker_thread_;

    job_id_t shadow_job_id_;
    ExecutionTemplate *execution_template_;

    SpawnState spawn_state_;
    std::string template_name_;
    std::vector<job_id_t> template_inner_job_ids_;
    std::vector<job_id_t> template_outer_job_ids_;
    std::vector<Parameter> template_parameters_;

    DependencyQuery *dependency_query_;

  protected:
    // TODO(omidm) should remove it later; left them now so the tests
    // that use it still pass.
    JobType type_;
};

class RemoteCopySendJob : public Job {
  public:
    explicit RemoteCopySendJob(WorkerDataExchanger* de, Application *app);
    ~RemoteCopySendJob();

    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job* Clone();
    virtual void Sleep() {}
    virtual void Cancel() {}

    ID<job_id_t> receive_job_id();
    ID<job_id_t> mega_rcr_job_id();
    ID<worker_id_t> to_worker_id();
    std::string to_ip();
    ID<port_t> to_port();

    void set_receive_job_id(ID<job_id_t> receive_job_id);
    void set_mega_rcr_job_id(ID<job_id_t> mega_rcr_job_id);
    void set_to_worker_id(ID<worker_id_t> worker_id);
    void set_to_ip(std::string ip);
    void set_to_port(ID<port_t> port);
    void set_template_generation_id(template_id_t id);

  private:
    ID<job_id_t> receive_job_id_;
    ID<job_id_t> mega_rcr_job_id_;
    ID<worker_id_t> to_worker_id_;
    std::string to_ip_;
    ID<port_t> to_port_;
    WorkerDataExchanger* data_exchanger_;
    template_id_t template_generation_id_;
};

class RemoteCopyReceiveJob : public Job {
  public:
    explicit RemoteCopyReceiveJob(Application *app);
    ~RemoteCopyReceiveJob();

    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job* Clone();
    virtual void Sleep() {}
    virtual void Cancel() {}

    void set_serialized_data(SerializedData* ser_data);
    void set_data_version(data_version_t version);

  private:
    SerializedData * serialized_data_;
    data_version_t data_version_;
};

class MegaRCRJob : public Job {
  public:
    explicit MegaRCRJob(Application *app);
    MegaRCRJob(Application *app,
               const std::vector<job_id_t>& receive_job_ids,
               const std::vector<physical_data_id_t>& to_phy_ids);
    ~MegaRCRJob();

    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job* Clone() {assert(false); return NULL;}
    virtual void Sleep() {}
    virtual void Cancel() {}

    const std::vector<job_id_t>* receive_job_ids_p();
    const std::vector<physical_data_id_t>* to_phy_ids_p();

    void set_receive_job_ids(const std::vector<job_id_t>& receive_job_ids);
    void set_to_phy_ids(const std::vector<physical_data_id_t>& to_phy_ids);
    void set_serialized_data(job_id_t job_id, SerializedData* ser_data);
    void set_serialized_data_map(const std::map<job_id_t, SerializedData*>& map);
    void clear_serialized_data_map();

    bool AllDataReceived();

  private:
    std::vector<job_id_t> receive_job_ids_;
    std::vector<physical_data_id_t> to_phy_ids_;
    std::map<job_id_t, SerializedData*> serialized_data_map_;
};

class SaveDataJob : public Job {
  public:
    explicit SaveDataJob(DistributedDB *ddb, Application *app);
    ~SaveDataJob();

    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job* Clone();
    virtual void Sleep() {}
    virtual void Cancel() {}

    std::string handle();
    checkpoint_id_t checkpoint_id();

    void set_checkpoint_id(checkpoint_id_t checkpoint_id);

  private:
    std::string handle_;
    checkpoint_id_t checkpoint_id_;
    DistributedDB *ddb_;
};

class LoadDataJob : public Job {
  public:
    explicit LoadDataJob(DistributedDB *ddb, Application *app);
    ~LoadDataJob();

    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job* Clone();
    virtual void Sleep() {}
    virtual void Cancel() {}

    void set_handle(std::string handle);
    void set_data_version(data_version_t version);

  private:
    std::string handle_;
    data_version_t data_version_;
    DistributedDB *ddb_;
};

class LocalCopyJob : public Job {
  public:
    explicit LocalCopyJob(Application *app);
    ~LocalCopyJob();

    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job* Clone();
    virtual void Sleep() {}
    virtual void Cancel() {}

    static void PrintTimeProfile();

  private:
    static double copy_ghost_time_;
    static double copy_central_time_;
    static int64_t copy_ghost_count_;
    static int64_t copy_central_count_;
};

class CreateDataJob : public Job {
  public:
    CreateDataJob();
    ~CreateDataJob();

    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job* Clone();
    virtual void Sleep() {}
    virtual void Cancel() {}

  private:
};




}  // namespace nimbus
#endif  // NIMBUS_SRC_WORKER_JOB_H_


