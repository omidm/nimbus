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
#include "shared/nimbus_types.h"
#include "worker/data.h"
#include "shared/id.h"
#include "shared/idset.h"
#include "shared/serialized_data.h"
#include "shared/worker_data_exchanger.h"

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

    virtual void Execute(std::string params, const DataArray& da) {}
    virtual Job* Clone();
    virtual void Sleep() {}
    virtual void Cancel() {}

    ID<job_id_t> id();
    IDSet<data_id_t> read_set();
    IDSet<data_id_t> write_set();
    IDSet<job_id_t> before_set();
    IDSet<job_id_t> after_set();
    std::string parameters();

    void set_id(ID<job_id_t> id);
    void set_read_set(IDSet<data_id_t> read_set);
    void set_write_set(IDSet<data_id_t> write_set);
    void set_before_set(IDSet<job_id_t> before_set);
    void set_after_set(IDSet<job_id_t> after_set);
    void set_parapeters(std::string parameters);


  protected:
    ID<job_id_t> id_;
    IDSet<data_id_t> read_set_;
    IDSet<data_id_t> write_set_;
    IDSet<job_id_t> before_set_;
    IDSet<job_id_t> after_set_;
    std::string parameters_;

    // TODO(omidm) should remove bothe of them, later; left them now so the tests
    // that use them still pass.
    Application* application_;
    JobType type_;
};

class RemoteCopySendJob : public Job {
  public:
    explicit RemoteCopySendJob(WorkerDataExchanger* de);
    ~RemoteCopySendJob();

    virtual void Execute(std::string params, const DataArray& da);
    virtual Job* Clone();
    virtual void Sleep() {}
    virtual void Cancel() {}

    ID<worker_id_t> to_worker_id();
    std::string to_ip();
    ID<port_t> to_port();

    void set_to_worker_id(ID<worker_id_t> worker_id);
    void set_to_ip(std::string ip);
    void set_to_port(ID<port_t> port);

  private:
    ID<worker_id_t> to_worker_id_;
    std::string to_ip_;
    ID<port_t> to_port_;
    WorkerDataExchanger* data_exchanger_;
};

class RemoteCopyReceiveJob : public Job {
  public:
    RemoteCopyReceiveJob();
    ~RemoteCopyReceiveJob();

    virtual void Execute(std::string params, const DataArray& da);
    virtual Job* Clone();
    virtual void Sleep() {}
    virtual void Cancel() {}

    void set_serialized_data(SerializedData* ser_data);

  private:
    SerializedData * serialized_data_;
};





}  // namespace nimbus
#endif  // NIMBUS_WORKER_JOB_H_


