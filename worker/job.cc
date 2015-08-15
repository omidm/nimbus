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

#include <time.h>
#include "shared/fast_log.hh"
#include "worker/job.h"
#include "worker/application.h"
#include "worker/app_data_manager.h"

using namespace nimbus; // NOLINT

Job::Job() {
  app_is_set_ = false;
  sterile_ = true;
  use_threading_ = false;
  core_quota_ = 1;
  run_time_ = 0;
  wait_time_ = 0;
  max_alloc_ = 0;
  worker_thread_ = NULL;
  spawn_state_ = INIT;
  dependency_query_ = new DependencyQuery();
}

Job::~Job() {
  delete dependency_query_;
}

// TODO(omidm) should remove this later. left it now so the tests
// that use it still pass.
Job::Job(Application* app, JobType type) {
  application_ = app;
  type_ = type;
  run_time_ = 0;
  wait_time_ = 0;
  max_alloc_ = 0;
}

Job* Job::Clone() {
  std::cout << "cloning the base job\n";
  Job* j = new Job();
  return j;
}

bool Job::SpawnComputeJob(const std::string& name,
                          const job_id_t& id,
                          const IDSet<logical_data_id_t>& read,
                          const IDSet<logical_data_id_t>& write,
                          const IDSet<job_id_t>& before,
                          const IDSet<job_id_t>& after,
                          const Parameter& params,
                          const bool& sterile,
                          const GeometricRegion& region,
                          const job_id_t& future_job_id) {
  if (sterile_) {
    dbg(DBG_ERROR, "ERROR: the job is sterile, it cannot spawn jobs.\n");
    exit(-1);
    return false;
  }

  if (!app_is_set_) {
    dbg(DBG_ERROR, "ERROR: SpawnComputeJob, application has not been set.\n");
    exit(-1);
    return false;
  }

  switch (spawn_state_) {
    case START_KNOWN_TEMPLATE:
      template_inner_job_ids_.push_back(id);
      template_parameters_.push_back(params);
      return true;
      break;
    case END_TEMPLATE:
      dbg(DBG_ERROR, "ERROR: currently we do not support both normal jobs and templates in same non-sterile job!\n"); // NOLINT
      exit(-1);
      return false;
      break;
    case INIT:
      spawn_state_ = NORMAL;
    case NORMAL:
    case START_UNKNOWN_TEMPLATE:
      application_->SpawnComputeJob(name,
                                    id,
                                    read,
                                    write,
                                    before,
                                    after,
                                    id_.elem(),
                                    future_job_id,
                                    sterile,
                                    region,
                                    params);
      return true;
      break;
  }

  return true;
}

bool Job::SpawnCopyJob(const job_id_t& id,
    const logical_data_id_t& from_logical_id,
    const logical_data_id_t& to_logical_id,
    const IDSet<job_id_t>& before,
    const IDSet<job_id_t>& after,
    const Parameter& params) {
  if (sterile_) {
    dbg(DBG_ERROR, "ERROR: the job is sterile, it cannot spawn jobs.\n");
    return false;
  }
  if (app_is_set_) {
    application_->SpawnCopyJob(id, from_logical_id, to_logical_id, before,
                               after, id_.elem());
    return true;
  } else {
    std::cout << "ERROR: SpawnCopyJob, application has not been set." <<
      std::endl;
    return false;
  }
}

bool Job::DefineData(const std::string& name,
                     const logical_data_id_t& logical_data_id,
                     const partition_id_t& partition_id,
                     const IDSet<partition_id_t>& neighbor_partition) {
  if (app_is_set_) {
    query_cache_.clear();
    application_->DefineData(name, logical_data_id, partition_id,
                             neighbor_partition, id_.elem());
    return true;
  } else {
    std::cout << "ERROR: DefineData, application has not been set." <<
      std::endl;
    return false;
  }
}

bool Job::DefinePartition(const ID<partition_id_t>& partition_id,
                          const GeometricRegion& r) {
    if (app_is_set_) {
        application_->DefinePartition(partition_id, r);
        return true;
    } else {
        std::cout << "ERROR: DefinePartition, application has not been set." <<
            std::endl;
        return false;
    }
}

bool Job::DefineJobGraph(const std::string& job_graph_name) {
  if (app_is_set_) {
    application_->DefineJobGraph(job_graph_name);
    return true;
  } else {
    std::cout << "ERROR: DefineJobGraph, application has not been set." << std::endl;
    return false;
  }
}

bool Job::TerminateApplication(const exit_status_t& exit_status) {
    if (app_is_set_) {
        application_->TerminateApplication(exit_status);
        return true;
    } else {
        std::cout << "ERROR: TerminateApplication, application has not been set." <<
            std::endl;
        return false;
    }
}

bool Job::GetNewJobID(std::vector<job_id_t>* result, size_t req_num) {
  if (app_is_set_) {
    return application_->GetNewJobID(result, req_num);
  } else {
    std::cout << "ERROR: GetNewJobID, application has not been set." <<
      std::endl;
    return false;
  }
}

bool Job::GetNewLogicalDataID(std::vector<logical_data_id_t>* result, size_t req_num) {
  if (app_is_set_) {
    return application_->GetNewLogicalDataID(result, req_num);
  } else {
    std::cout << "ERROR: GetNewDataID, application has not been set." <<
      std::endl;
    return false;
  }
}

bool Job::GetPartition(partition_id_t id, GeometricRegion* r) const {
  if (app_is_set_) {
      return application_->GetPartition(id, r);
  } else {
      std::cout << "Error: GetLogicalObject, application has not been set." << std::endl;
      exit(-1);
  }
}

const LogicalDataObject* Job::GetLogicalObject(logical_data_id_t id) const {
  if (app_is_set_) {
      return application_->GetLogicalObject(id);
  } else {
      std::cout << "Error: GetLogicalObject, application has not been set." << std::endl;
      exit(-1);
  }
}

int Job::GetCoveredLogicalObjects(CLdoVector* result,
    const std::string& variable,
    const GeometricRegion* r) {
  if (app_is_set_) {
      return application_->GetCoveredLogicalObjects(result, variable, r);
  } else {
      std::cout << "Error: GetCoveredLogicalObjects, application has not been set." << std::endl;
      return -1;
  }
}

int Job::GetAdjacentLogicalObjects(CLdoVector* result,
    const std::string& variable,
    const GeometricRegion* r) {
  if (app_is_set_) {
      return application_->GetAdjacentLogicalObjects(result, variable, r);
  } else {
      std::cout << "Error: GetAdjacentLogicalObjects, application has not been set." << std::endl;
      return -1;
  }
}

int Job::AddIntersectingLdoIds(const std::string& variable,
                               const nimbus::GeometricRegion& region,
                               IDSet<logical_data_id_t>* result) {
  if (app_is_set_) {
    LdoIndexCache* cache = NULL;
    if (query_cache_.find(variable) == query_cache_.end()) {
      cache = &(query_cache_[variable]);
      cache->Initialize(application_, variable);
    } else {
      cache = &(query_cache_[variable]);
    }
    cache->GetResult(region, result);
    return result->size();
  } else {
    std::cout << "Error: GetAdjacentLogicalObjects, application has not been set." << std::endl;
    return -1;
  }
}

int Job::GetIntersectingLogicalObjects(CLdoVector* result,
    const std::string& variable,
    const GeometricRegion* r) {
  if (app_is_set_) {
      return application_->GetIntersectingLogicalObjects(result, variable, r);
  } else {
      std::cout << "Error: GetAdjacentLogicalObjects, application has not been set." << std::endl;
      return -1;
  }
}



void Job::LoadLdoIdsInSet(IDSet<logical_data_id_t>* set,
                          const nimbus::GeometricRegion& region,
                          ...) {
  switch (spawn_state_) {
    case START_KNOWN_TEMPLATE:
      // Neutralize the call - omidm
      break;
    case INIT:
    case NORMAL:
    case END_TEMPLATE:
    case START_UNKNOWN_TEMPLATE:
      CLdoVector result;
      va_list vl;
      va_start(vl, region);
      char* arg = va_arg(vl, char*);
      while (arg != NULL) {
        // AddIntersectingLdoIds(arg, region, set);
        GetIntersectingLogicalObjects(&result, arg, &region);
        for (size_t i = 0; i < result.size(); ++i) {
          set->insert(result[i]->id());
        }
        arg = va_arg(vl, char*);
      }
      va_end(vl);
      break;
  }
}

bool Job::StageJobAndLoadBeforeSet(IDSet<job_id_t> *before_set,
                                   const std::string& name,
                                   const job_id_t& id,
                                   const IDSet<logical_data_id_t>& read,
                                   const IDSet<logical_data_id_t>& write,
                                   const bool barrier) {
  switch (spawn_state_) {
    case START_KNOWN_TEMPLATE:
      // Neutralize the call - omidm
      break;
    case END_TEMPLATE:
      dbg(DBG_ERROR, "ERROR: currently dependency quey is not valid if template is already ended in a non-sterile job!\n"); // NOLINT
      exit(-1);
      return false;
      break;
    case INIT:
    case NORMAL:
    case START_UNKNOWN_TEMPLATE:
      return dependency_query_->StageJobAndLoadBeforeSet(before_set,
                                                         name,
                                                         id,
                                                         read,
                                                         write,
                                                         barrier);
      break;
  }

  return true;
}

bool Job::MarkEndOfStage() {
  switch (spawn_state_) {
    case START_KNOWN_TEMPLATE:
      // Neutralize the call - omidm
      break;
    case END_TEMPLATE:
      dbg(DBG_ERROR, "ERROR: currently dependency quey is not valid if template is already ended in a non-sterile job!\n"); // NOLINT
      exit(-1);
      return false;
      break;
    case INIT:
    case NORMAL:
    case START_UNKNOWN_TEMPLATE:
      return dependency_query_->MarkEndOfStage();
      break;
  }

  return true;
}

void Job::StartTemplate(const std::string& template_name) {
  if (sterile_) {
    dbg(DBG_ERROR, "ERROR: the job is sterile, it cannot start a template.\n");
    exit(-1);
  }

  if (!app_is_set_) {
      dbg(DBG_ERROR, "ERROR: StartTEmplate, application has not been set.\n");
      exit(-1);
  }

  switch (spawn_state_) {
    case START_KNOWN_TEMPLATE:
    case START_UNKNOWN_TEMPLATE:
      dbg(DBG_ERROR, "ERROR: Cannot start a template with in another template!\n");
      exit(-1);
      break;
    case END_TEMPLATE:
      dbg(DBG_ERROR, "ERROR: currently we do not support spawning two templates in same non-sterile job!\n"); // NOLINT
      exit(-1);
      break;
    case NORMAL:
      dbg(DBG_ERROR, "ERROR: currently we do not support both normal jobs and templates in same non-sterile job!\n"); // NOLINT
      exit(-1);
      break;
    case INIT:
      template_name_ = template_name;
      if (IsTemplateDefined(template_name)) {
        spawn_state_ = START_KNOWN_TEMPLATE;
      } else {
        spawn_state_ = START_UNKNOWN_TEMPLATE;
        application_->StartTemplate(template_name, id_.elem());
      }
      break;
  }
}

void Job::EndTemplate(const std::string& template_name) {
  if (sterile_) {
    dbg(DBG_ERROR, "ERROR: the job is sterile, it cannot end a template.\n");
    exit(-1);
  }

  if (!app_is_set_) {
      dbg(DBG_ERROR, "ERROR: EndTemplate, application has not been set.\n");
      exit(-1);
  }

  if (template_name_ != template_name) {
      dbg(DBG_ERROR, "ERROR: template name in end does not match the name in start mark!\n");
      exit(-1);
  }

  switch (spawn_state_) {
    case INIT:
    case NORMAL:
    case END_TEMPLATE:
      dbg(DBG_ERROR, "ERROR: Unbalanced end and start template marks!\n");
      exit(-1);
      break;
    case START_UNKNOWN_TEMPLATE:
      spawn_state_ = END_TEMPLATE;
      application_->EndTemplate(template_name, id_.elem());
      break;
    case START_KNOWN_TEMPLATE:
      spawn_state_ = END_TEMPLATE;
      application_->SpawnTemplate(template_name,
                                  template_inner_job_ids_,
                                  template_outer_job_ids_,
                                  template_parameters_,
                                  id_.elem());
      break;
  }
}

bool Job::IsTemplateDefined(const std::string& template_name) {
  if (!app_is_set_) {
      dbg(DBG_ERROR, "ERROR: IsTemplateDefined, application has not been set.\n");
      exit(-1);
      return false;
  }

  return application_->IsTemplateDefined(template_name);
}

std::string Job::name() const {
  return name_;
}

ID<job_id_t> Job::id() const {
  return id_;
}

IDSet<physical_data_id_t> Job::read_set() const {
  return read_set_;
}

IDSet<physical_data_id_t> Job::write_set() const {
  return write_set_;
}

const IDSet<physical_data_id_t>& Job::get_read_set() const {
  return read_set_;
}

const IDSet<physical_data_id_t>& Job::get_write_set() const {
  return write_set_;
}

IDSet<job_id_t> Job::before_set() const {
  return before_set_;
}

const IDSet<job_id_t>* Job::before_set_p() const {
  return &before_set_;
}

IDSet<job_id_t> Job::after_set() const {
  return after_set_;
}

Parameter Job::parameters() const {
  return parameters_;
}

Application* Job::application() const {
  return application_;
}

bool Job::sterile() const {
  return sterile_;
}

GeometricRegion Job::region() const {
  return region_;
}

ID<job_id_t> Job::future_job_id() const {
  return future_job_id_;
}

double Job::run_time() const {
  return run_time_;
}

double Job::wait_time() const {
  return wait_time_;
}

void Job::set_name(std::string name) {
  name_ = name;
}

void Job::set_id(ID<job_id_t> id) {
  id_ = id;
}

void Job::set_read_set(IDSet<physical_data_id_t> read_set) {
  read_set_ = read_set;
}

void Job::set_write_set(IDSet<physical_data_id_t> write_set) {
  write_set_ = write_set;
}

void Job::set_before_set(IDSet<job_id_t> before_set) {
  before_set_ = before_set;
}

void Job::set_after_set(IDSet<job_id_t> after_set) {
  after_set_ = after_set;
}

void Job::set_parameters(Parameter parameters) {
  parameters_ = parameters;
}

void Job::set_application(Application* app) {
  application_ = app;
  app_is_set_ = true;
}

void Job::set_sterile(bool sterile) {
  sterile_ = sterile;
}

void Job::set_region(GeometricRegion region) {
  region_ = region;
}

void Job::set_future_job_id(ID<job_id_t> future_job_id) {
  future_job_id_ = future_job_id;
}

void Job::set_run_time(double run_time) {
  run_time_ = run_time;
}

void Job::set_wait_time(double wait_time) {
  wait_time_ = wait_time;
}

AppDataManager* Job::GetAppDataManager() const {
  return application_->app_data_manager();
}

StaticConfigManager* Job::GetStaticConfigManager() const {
  return application_->static_config_manager();
}

RemoteCopySendJob::RemoteCopySendJob(WorkerDataExchanger* da, Application *app) {
  data_exchanger_ = da;
  set_application(app);
}

RemoteCopySendJob::~RemoteCopySendJob() {
}

// TODO(quhang) data exchanger is thread-safe?
void RemoteCopySendJob::Execute(Parameter params, const DataArray& da) {
  AppDataManager *am = GetAppDataManager();
  am->SyncData(da[0]);
  SerializedData ser_data;
  da[0]->Serialize(&ser_data);
  data_exchanger_->SendSerializedData(receive_job_id().elem(),
                                      mega_rcr_job_id().elem(),
                                      to_worker_id_.elem(),
                                      ser_data,
                                      da[0]->version());
  // delete ser_data.data_ptr(); // Not needed with shared pointer.
}

Job* RemoteCopySendJob::Clone() {
  return new RemoteCopySendJob(data_exchanger_, application());
}

ID<job_id_t> RemoteCopySendJob::receive_job_id() {
  return receive_job_id_;
}

ID<job_id_t> RemoteCopySendJob::mega_rcr_job_id() {
  return mega_rcr_job_id_;
}

ID<worker_id_t> RemoteCopySendJob::to_worker_id() {
  return to_worker_id_;
}

std::string RemoteCopySendJob::to_ip() {
  return to_ip_;
}

ID<port_t> RemoteCopySendJob::to_port() {
  return to_port_;
}


void RemoteCopySendJob::set_receive_job_id(ID<job_id_t> receive_job_id) {
  receive_job_id_ = receive_job_id;
}

void RemoteCopySendJob::set_mega_rcr_job_id(ID<job_id_t> mega_rcr_job_id) {
  mega_rcr_job_id_ = mega_rcr_job_id;
}

void RemoteCopySendJob::set_to_worker_id(ID<worker_id_t> worker_id) {
  to_worker_id_ = worker_id;
}

void RemoteCopySendJob::set_to_ip(std::string ip) {
  to_ip_ = ip;
}

void RemoteCopySendJob::set_to_port(ID<port_t> port) {
  to_port_ = port;
}


RemoteCopyReceiveJob::RemoteCopyReceiveJob(Application *app) {
  set_application(app);
}

RemoteCopyReceiveJob::~RemoteCopyReceiveJob() {
}

void RemoteCopyReceiveJob::Execute(Parameter params, const DataArray& da) {
  AppDataManager *am = GetAppDataManager();
  timer::StartTimer(timer::kInvalidateMappings);
  am->InvalidateMappings(da[0]);
  timer::StopTimer(timer::kInvalidateMappings);
  Data * data_copy = NULL;
  da[0]->DeSerialize(*serialized_data_, &data_copy);
  timer::StartTimer(timer::kRCRCopy);
  da[0]->Copy(data_copy);
  da[0]->set_version(data_version_);

  data_copy->Destroy();
  timer::StopTimer(timer::kRCRCopy);
  // delete serialized_data_->data_ptr(); // Not needed with shared pointer.
  delete serialized_data_;
}

Job* RemoteCopyReceiveJob::Clone() {
  return new RemoteCopyReceiveJob(application());
}

void RemoteCopyReceiveJob::set_serialized_data(SerializedData* ser_data) {
  serialized_data_ = ser_data;
}

void RemoteCopyReceiveJob::set_data_version(data_version_t version) {
  data_version_ = version;
}

MegaRCRJob::MegaRCRJob(Application *app) {
  set_application(app);
}

MegaRCRJob::MegaRCRJob(Application *app,
                    const std::vector<job_id_t>& receive_job_ids,
                    const std::vector<physical_data_id_t>& to_phy_ids)
  : receive_job_ids_(receive_job_ids),
    to_phy_ids_(to_phy_ids) {
  assert(receive_job_ids_.size() == to_phy_ids_.size());
  set_application(app);
}

MegaRCRJob::~MegaRCRJob() {
}

void MegaRCRJob::Execute(Parameter params, const DataArray& da) {
  assert(AllDataReceived());
  AppDataManager *am = GetAppDataManager();
  size_t idx = 0;
  std::vector<job_id_t>::iterator iter = receive_job_ids_.begin();
  for (; iter != receive_job_ids_.end(); ++iter, ++idx) {
    timer::StartTimer(timer::kInvalidateMappings);
    am->InvalidateMappings(da[idx]);
    timer::StopTimer(timer::kInvalidateMappings);
    Data * data_copy = NULL;
    std::map<job_id_t, SerializedData*>::iterator it =
      serialized_data_map_.find(*iter);
    assert(it != serialized_data_map_.end());
    SerializedData *ser_data = it->second;
    da[idx]->DeSerialize(*ser_data, &data_copy);
    da[idx]->Copy(data_copy);

    data_copy->Destroy();
    // delete ser_data->data_ptr(); // Not needed with shared pointer.
    delete ser_data;
  }
}

const std::vector<job_id_t>* MegaRCRJob::receive_job_ids_p() {
  return &receive_job_ids_;
}

const std::vector<physical_data_id_t>* MegaRCRJob::to_phy_ids_p() {
  return &to_phy_ids_;
}

void MegaRCRJob::set_receive_job_ids(const std::vector<job_id_t>& receive_job_ids) {
  receive_job_ids_ = receive_job_ids;
}

void MegaRCRJob::set_to_phy_ids(const std::vector<physical_data_id_t>& to_phy_ids) {
  to_phy_ids_ = to_phy_ids;
}

void MegaRCRJob::set_serialized_data(job_id_t job_id, SerializedData* ser_data) {
  serialized_data_map_[job_id] = ser_data;
}

void MegaRCRJob::set_serialized_data_map(const std::map<job_id_t, SerializedData*>& map) {
  serialized_data_map_ = map;
}

void MegaRCRJob::clear_serialized_data_map() {
  serialized_data_map_.clear();
}

bool MegaRCRJob::AllDataReceived() {
  return serialized_data_map_.size() == receive_job_ids_.size();
}

SaveDataJob::SaveDataJob(DistributedDB *ddb, Application *app) {
  ddb_ = ddb;
  set_application(app);
}

SaveDataJob::~SaveDataJob() {
}

void SaveDataJob::Execute(Parameter params, const DataArray& da) {
  AppDataManager *am = GetAppDataManager();
  am->SyncData(da[0]);

  SerializedData ser_data;
  da[0]->Serialize(&ser_data);

  std::string key = int2string(id().elem());
  std::string value(ser_data.data_ptr().get(), ser_data.size());
  if (!ddb_->Put(key, value, checkpoint_id_, &handle_)) {
    dbg(DBG_ERROR, "ERROR: could not save the data in ddb!\n");
    exit(-1);
  }

  // delete ser_data.data_ptr(); // Not needed with shared pointer.
}

Job* SaveDataJob::Clone() {
  return new SaveDataJob(ddb_, application());
}

std::string SaveDataJob::handle() {
  return handle_;
}

checkpoint_id_t SaveDataJob::checkpoint_id() {
  return checkpoint_id_;
}

void SaveDataJob::set_checkpoint_id(checkpoint_id_t checkpoint_id) {
  checkpoint_id_ = checkpoint_id;
}

LoadDataJob::LoadDataJob(DistributedDB *ddb, Application *app) {
  ddb_ = ddb;
  set_application(app);
}

LoadDataJob::~LoadDataJob() {
}

void LoadDataJob::Execute(Parameter params, const DataArray& da) {
  AppDataManager *am = GetAppDataManager();
  am->InvalidateMappings(da[0]);

  std::string value;
  if (!ddb_->Get(handle_, &value)) {
    dbg(DBG_ERROR, "ERROR: could not load the data from ddb!\n");
    exit(-1);
  }
  SerializedData serialized_data_(value);

  Data * data_copy = NULL;
  da[0]->DeSerialize(serialized_data_, &data_copy);
  da[0]->Copy(data_copy);
  da[0]->set_version(data_version_);

  data_copy->Destroy();
  // delete serialized_data_->data_ptr(); // Not needed with shared pointer.
  // delete serialized_data; // Not needed, goes out of context.
}

Job* LoadDataJob::Clone() {
  return new LoadDataJob(ddb_, application());
}

void LoadDataJob::set_handle(std::string handle) {
  handle_ = handle;
}

void LoadDataJob::set_data_version(data_version_t version) {
  data_version_ = version;
}

LocalCopyJob::LocalCopyJob(Application *app) {
  set_application(app);
}

LocalCopyJob::~LocalCopyJob() {
}

Job* LocalCopyJob::Clone() {
  return new LocalCopyJob(application());
}

void LocalCopyJob::Execute(Parameter params, const DataArray& da) {
  struct timespec start_time;
  struct timespec t;
  clock_gettime(CLOCK_REALTIME, &start_time);
  AppDataManager *am = GetAppDataManager();
  am->SyncData(da[0]);
  am->InvalidateMappings(da[1]);
  da[1]->Copy(da[0]);
  da[1]->set_version(da[0]->version());
  clock_gettime(CLOCK_REALTIME, &t);
  GeometricRegion region = da[1]->region();
  // TODO(quhang): Is this needed?
  // const int_dimension_t kMargin = 10;
  // if (region.dx() < kMargin || region.dy() < kMargin || region.dz() < kMargin) {
  //   copy_ghost_count_++;
  //   copy_ghost_time_ += difftime(t.tv_sec, start_time.tv_sec)
  //       + .000000001 * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));
  // } else {
  //   copy_central_count_++;
  //   copy_central_time_ += difftime(t.tv_sec, start_time.tv_sec)
  //       + .000000001 * (static_cast<double>(t.tv_nsec - start_time.tv_nsec));
  //   // TODO(quhang): temporary use. Should use dbg instead.
  //   // printf("[PROFILE] Central Copy %s, %s\n", da[1]->name().c_str(),
  //   //        region.ToNetworkData().c_str());
  // }
}

void LocalCopyJob::PrintTimeProfile() {
  // TODO(quhang): temporary use. Should use dbg instead.
  // printf("[PROFILE] copy of ghost: %lld, %f seconds.\n"
  //        "[PROFILE] copy of central: %lld, %f seconds.\n",
  //        copy_ghost_count_, copy_ghost_time_,
  //        copy_central_count_, copy_central_time_);
}

double LocalCopyJob::copy_ghost_time_ = 0;
double LocalCopyJob::copy_central_time_ = 0;
int64_t LocalCopyJob::copy_ghost_count_ = 0;
int64_t LocalCopyJob::copy_central_count_ = 0;

CreateDataJob::CreateDataJob() {
}

CreateDataJob::~CreateDataJob() {
}

Job* CreateDataJob::Clone() {
  return new CreateDataJob();
}

void CreateDataJob::Execute(Parameter params, const DataArray& da) {
  da[0]->Create();
  da[0]->set_version(NIMBUS_INIT_DATA_VERSION);
}


