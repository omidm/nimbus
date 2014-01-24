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

#ifndef NIMBUS_WORKER_JOB_H_
#define NIMBUS_WORKER_JOB_H_

#include <vector>
#include <string>
#include <set>
#include <map>
#include <list>
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"
#include "shared/id.h"
#include "shared/idset.h"
#include "shared/serialized_data.h"
#include "shared/worker_data_exchanger.h"
#include "worker/data.h"
#include "worker/worker_ldo_map.h"

namespace nimbus {

class Application;
class Job;
typedef std::list<Job*> JobList;
typedef std::map<job_id_t, Job*> JobMap;
typedef std::map<std::string, Job*> JobTable;
typedef std::vector<Data*> DataArray;

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

    bool SpawnComputeJob(const std::string& name,
        const job_id_t& id,
        const IDSet<logical_data_id_t>& read,
        const IDSet<logical_data_id_t>& write,
        const IDSet<job_id_t>& before,
        const IDSet<job_id_t>& after,
        const Parameter& params);

    bool SpawnCopyJob(const job_id_t& id,
        const logical_data_id_t& from_logical_id,
        const logical_data_id_t& to_logical_id,
        const IDSet<job_id_t>& before,
        const IDSet<job_id_t>& after,
        const Parameter& params);

    bool DefineData(const std::string& name,
        const logical_data_id_t& logical_data_id,
        const partition_id_t& partition_id,
        const IDSet<partition_id_t>& neighbor_partition,
        const Parameter& params);

    bool DefinePartition(const ID<partition_id_t>& partition_id,
         const GeometricRegion& r,
         const Parameter& params);

    bool TerminateApplication(const exit_status_t& exit_status_id = NIMBUS_TERMINATE_SUCCESS);

    bool GetNewJobID(std::vector<job_id_t>* result, size_t req_num);
    bool GetNewLogicalDataID(std::vector<logical_data_id_t>* result, size_t req_num);

    GeometricRegion GetPartition(partition_id_t id) const;

    const LogicalDataObject* GetLogicalObject(logical_data_id_t id) const;
    int GetCoveredLogicalObjects(CLdoVector* result,
         const std::string& variable,
         const GeometricRegion* r);
    int GetAdjacentLogicalObjects(CLdoVector* result,
         const std::string& variable,
         const GeometricRegion* r);

    std::string name();
    ID<job_id_t> id();
    IDSet<physical_data_id_t> read_set();
    IDSet<physical_data_id_t> write_set();
    IDSet<job_id_t> before_set();
    IDSet<job_id_t> after_set();
    Parameter parameters();
    Application* application();

    void set_name(std::string name);
    void set_id(ID<job_id_t> id);
    void set_read_set(IDSet<physical_data_id_t> read_set);
    void set_write_set(IDSet<physical_data_id_t> write_set);
    void set_before_set(IDSet<job_id_t> before_set);
    void set_after_set(IDSet<job_id_t> after_set);
    void set_parameters(Parameter parameters);
    void set_application(Application* app);

  private:
    std::string name_;
    ID<job_id_t> id_;
    IDSet<physical_data_id_t> read_set_;
    IDSet<physical_data_id_t> write_set_;
    IDSet<job_id_t> before_set_;
    IDSet<job_id_t> after_set_;
    Parameter parameters_;
    Application* application_;
    bool app_is_set_;

  protected:
    // TODO(omidm) should remove it later; left them now so the tests
    // that use it still pass.
    JobType type_;
};

class RemoteCopySendJob : public Job {
  public:
    explicit RemoteCopySendJob(WorkerDataExchanger* de);
    ~RemoteCopySendJob();

    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job* Clone();
    virtual void Sleep() {}
    virtual void Cancel() {}

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
    WorkerDataExchanger* data_exchanger_;
};

class RemoteCopyReceiveJob : public Job {
  public:
    RemoteCopyReceiveJob();
    ~RemoteCopyReceiveJob();

    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job* Clone();
    virtual void Sleep() {}
    virtual void Cancel() {}

    void set_serialized_data(SerializedData* ser_data);

  private:
    SerializedData * serialized_data_;
};

class LocalCopyJob : public Job {
  public:
    LocalCopyJob();
    ~LocalCopyJob();

    virtual void Execute(Parameter params, const DataArray& da);
    virtual Job* Clone();
    virtual void Sleep() {}
    virtual void Cancel() {}

  private:
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
#endif  // NIMBUS_WORKER_JOB_H_


