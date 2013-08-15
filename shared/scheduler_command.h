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
  * Object representation of a scheduler command. Used by workers to
  * send commands to server and server to send commands down to workers. The
  * super class SchedulerCommand is inherited by its children implemented here.
  * Each child represents a specific command exchanged between scheduler and
  * worker.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SHARED_SCHEDULER_COMMAND_H_
#define NIMBUS_SHARED_SCHEDULER_COMMAND_H_


#include <sstream> // NOLINT
#include <string>
#include <vector>
#include <list>
#include <map>
#include <set>
#include "worker/job.h"
#include "shared/parser.h"
#include "shared/idset.h"
#include "shared/nimbus_types.h"

namespace nimbus {

class CommandParameter {
 public:
  CommandParameter();
  explicit CommandParameter(std::string parameter);
  CommandParameter(std::string name, std::string value, const IDSet<job_id_t>& set);
  virtual ~CommandParameter();

  virtual std::string toString();
  virtual std::string name();
  virtual std::string value();
  virtual IDSet<job_id_t>* identifier_set();

 private:
  std::string name_;
  std::string value_;
  IDSet<job_id_t> identifier_set_;
};

typedef std::map<std::string, CommandParameter> CommandParameterList;

class SchedulerCommand {
 public:
  SchedulerCommand();
  explicit SchedulerCommand(std::string command);
  SchedulerCommand(std::string name,
                   const CommandParameterList& parameters);
  virtual ~SchedulerCommand();

  virtual void addParameter(CommandParameter parameter);
  virtual std::string toString();
  virtual std::string name();
  virtual CommandParameterList* parameters();

  virtual worker_id_t worker_id();
  virtual void set_worker_id(worker_id_t id);

 protected:
  std::string name_;

 private:
  CommandParameterList parameters_;
  worker_id_t worker_id_;
};

typedef std::vector<SchedulerCommand*> SchedulerCommandVector;
typedef std::list<SchedulerCommand*> SchedulerCommandList;


class SpawnJobCommand : public SchedulerCommand {
  public:
    SpawnJobCommand();
    SpawnJobCommand(std::string job_name,
        const IDSet<job_id_t>& job_id,
        const IDSet<data_id_t>& read, const IDSet<data_id_t>& write,
        const IDSet<job_id_t>& before, const IDSet<job_id_t>& after,
        JobType job_type, std::string params);
    ~SpawnJobCommand();

    virtual std::string toString();
    std::string job_name();
    JobType job_type();
    IDSet<job_id_t> job_id();
    IDSet<data_id_t> read_set();
    IDSet<data_id_t> write_set();
    IDSet<job_id_t> before_set();
    IDSet<job_id_t> after_set();
    std::string params();

  private:
    std::string job_name_;
    IDSet<job_id_t> job_id_;
    IDSet<data_id_t> read_set_;
    IDSet<data_id_t> write_set_;
    IDSet<job_id_t> before_set_;
    IDSet<job_id_t> after_set_;
    JobType job_type_;
    std::string params_;
};

class JobDoneCommand : public SchedulerCommand {
  public:
    JobDoneCommand();
    JobDoneCommand(job_id_t job_id, std::string params);
    ~JobDoneCommand();

    virtual std::string toString();
    job_id_t job_id();
    std::string params();

  private:
    job_id_t job_id_;
    std::string params_;
};

class KillJobCommand : public SchedulerCommand {
  public:
    KillJobCommand();
    KillJobCommand(job_id_t job_id, std::string params);
    ~KillJobCommand();

    virtual std::string toString();
    job_id_t job_id();
    std::string params();

  private:
    job_id_t job_id_;
    std::string params_;
};

class PauseJobCommand : public SchedulerCommand {
  public:
    PauseJobCommand();
    PauseJobCommand(job_id_t job_id, std::string params);
    ~PauseJobCommand();

    virtual std::string toString();
    job_id_t job_id();
    std::string params();

  private:
    job_id_t job_id_;
    std::string params_;
};

class DefineDataCommand : public SchedulerCommand {
  public:
    DefineDataCommand();
    DefineDataCommand(std::string data_name, data_id_t data_id,
        IDSet<data_id_t>* neighbors, std::string params);
    ~DefineDataCommand();

    virtual std::string toString();
    std::string data_name();
    data_id_t data_id();
    IDSet<data_id_t>* neighbor_set();
    std::string params();

  private:
    std::string data_name_;
    data_id_t data_id_;
    IDSet<data_id_t> neighbor_set_;
    std::string params_;
};

class DataCreatedCommand : public SchedulerCommand {
  public:
    DataCreatedCommand();
    DataCreatedCommand(data_id_t data_id, std::string params);
    ~DataCreatedCommand();

    virtual std::string toString();
    data_id_t data_id();
    std::string params();

  private:
    data_id_t data_id_;
    std::string params_;
};

class CopyDataCommand : public SchedulerCommand {
  public:
    CopyDataCommand();
    CopyDataCommand(data_id_t data_id_from, data_id_t data_id_to,
        std::string params);
    ~CopyDataCommand();

    virtual std::string toString();
    data_id_t data_id_from();
    data_id_t data_id_to();
    std::string params();

  private:
    data_id_t data_id_from_;
    data_id_t data_id_to_;
    std::string params_;
};

class MigrateDataCommand : public SchedulerCommand {
  public:
    MigrateDataCommand();
    MigrateDataCommand(data_id_t data_id_from, data_id_t data_id_to,
        std::string params);
    ~MigrateDataCommand();

    virtual std::string toString();
    data_id_t data_id_from();
    data_id_t data_id_to();
    std::string params();

  private:
    data_id_t data_id_from_;
    data_id_t data_id_to_;
    // TODO(omidm): add the information of the remote worker holding the data.
    std::string params_;
};


}  // namespace nimbus

#endif  // NIMBUS_SHARED_SCHEDULER_COMMAND_H_