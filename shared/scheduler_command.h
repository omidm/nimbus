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
#include "shared/id.h"
#include "shared/idset.h"
#include "shared/escaper.h"
#include "shared/parameter.h"
#include "shared/geometric_region.h"
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

  enum Type {
    BASE,
    SPAWN_JOB,
    SPAWN_COMPUTE_JOB,
    SPAWN_COPY_JOB,
    DEFINE_DATA,
    HANDSHAKE,
    JOB_DONE,
    COMPUTE_JOB,
    CREATE_DATA,
    REMOTE_COPY_SEND,
    REMOTE_COPY_RECEIVE,
    LOCAL_COPY,
    DEFINE_PARTITION
  };

  typedef std::set<Type> TypeSet;

  virtual void addParameter(CommandParameter parameter);
  virtual std::string toString();
  virtual std::string toStringWTags();
  virtual std::string name();
  virtual Type type();
  virtual CommandParameterList* parameters();

  static std::string GetNameFromType(Type type);
  static bool GenerateSchedulerCommandChild(const std::string& input,
      TypeSet* command_set,
      SchedulerCommand*& ptr_generated_command);


  // virtual worker_id_t worker_id();
  // virtual void set_worker_id(worker_id_t id);

 protected:
  Type type_;
  std::string name_;
  static const std::string BASE_NAME;
  static const std::string SPAWN_JOB_NAME;
  static const std::string SPAWN_COMPUTE_JOB_NAME;
  static const std::string SPAWN_COPY_JOB_NAME;
  static const std::string DEFINE_DATA_NAME;
  static const std::string HANDSHAKE_NAME;
  static const std::string JOB_DONE_NAME;
  static const std::string COMPUTE_JOB_NAME;
  static const std::string CREATE_DATA_NAME;
  static const std::string REMOTE_COPY_SEND_NAME;
  static const std::string REMOTE_COPY_RECEIVE_NAME;
  static const std::string LOCAL_COPY_NAME;
  static const std::string DEFINE_PARTITION_NAME;

 private:
  CommandParameterList parameters_;
  // worker_id_t worker_id_;
};

// ************************************************

typedef std::vector<SchedulerCommand*> SchedulerCommandVector;
typedef std::list<SchedulerCommand*> SchedulerCommandList;

class HandshakeCommand : public SchedulerCommand {
  public:
    HandshakeCommand();
    explicit HandshakeCommand(const ID<worker_id_t>& worker_id,
        const std::string& ip, const ID<port_t>& port);
    ~HandshakeCommand();

    virtual std::string toString();
    virtual std::string toStringWTags();
    ID<worker_id_t> worker_id();
    std::string ip();
    ID<port_t> port();

  private:
    ID<worker_id_t> worker_id_;
    std::string ip_;
    ID<port_t> port_;
};


class SpawnJobCommand : public SchedulerCommand {
  public:
    SpawnJobCommand();
    SpawnJobCommand(const std::string& job_name,
        const IDSet<job_id_t>& job_id,
        const IDSet<data_id_t>& read, const IDSet<data_id_t>& write,
        const IDSet<job_id_t>& before, const IDSet<job_id_t>& after,
        const JobType& job_type, const Parameter& params);
    ~SpawnJobCommand();

    virtual std::string toString();
    virtual std::string toStringWTags();
    std::string job_name();
    JobType job_type();
    IDSet<job_id_t> job_id();
    IDSet<data_id_t> read_set();
    IDSet<data_id_t> write_set();
    IDSet<job_id_t> before_set();
    IDSet<job_id_t> after_set();
    Parameter params();

  private:
    std::string job_name_;
    IDSet<job_id_t> job_id_;
    IDSet<data_id_t> read_set_;
    IDSet<data_id_t> write_set_;
    IDSet<job_id_t> before_set_;
    IDSet<job_id_t> after_set_;
    JobType job_type_;
    Parameter params_;
};


class SpawnComputeJobCommand : public SchedulerCommand {
  public:
    SpawnComputeJobCommand();
    SpawnComputeJobCommand(const std::string& job_name,
        const ID<job_id_t>& job_id,
        const IDSet<data_id_t>& read, const IDSet<data_id_t>& write,
        const IDSet<job_id_t>& before, const IDSet<job_id_t>& after,
        const Parameter& params);
    ~SpawnComputeJobCommand();

    virtual std::string toString();
    virtual std::string toStringWTags();
    std::string job_name();
    ID<job_id_t> job_id();
    IDSet<data_id_t> read_set();
    IDSet<data_id_t> write_set();
    IDSet<job_id_t> before_set();
    IDSet<job_id_t> after_set();
    Parameter params();

  private:
    std::string job_name_;
    ID<job_id_t> job_id_;
    IDSet<data_id_t> read_set_;
    IDSet<data_id_t> write_set_;
    IDSet<job_id_t> before_set_;
    IDSet<job_id_t> after_set_;
    Parameter params_;
};

class SpawnCopyJobCommand : public SchedulerCommand {
  public:
    SpawnCopyJobCommand();
    SpawnCopyJobCommand(const ID<job_id_t>& job_id,
        const ID<data_id_t>& from_id, const ID<data_id_t>& to_id,
        const IDSet<job_id_t>& before, const IDSet<job_id_t>& after,
        const Parameter& params);
    ~SpawnCopyJobCommand();

    virtual std::string toString();
    virtual std::string toStringWTags();
    ID<job_id_t> job_id();
    ID<data_id_t> from_id();
    ID<data_id_t> to_id();
    IDSet<job_id_t> before_set();
    IDSet<job_id_t> after_set();
    Parameter params();

  private:
    ID<job_id_t> job_id_;
    ID<data_id_t> from_id_;
    ID<data_id_t> to_id_;
    IDSet<job_id_t> before_set_;
    IDSet<job_id_t> after_set_;
    Parameter params_;
};

class ComputeJobCommand : public SchedulerCommand {
  public:
    ComputeJobCommand();
    ComputeJobCommand(const std::string& job_name,
        const ID<job_id_t>& job_id,
        const IDSet<data_id_t>& read, const IDSet<data_id_t>& write,
        const IDSet<job_id_t>& before, const IDSet<job_id_t>& after,
        const Parameter& params);
    ~ComputeJobCommand();

    virtual std::string toString();
    virtual std::string toStringWTags();
    std::string job_name();
    ID<job_id_t> job_id();
    IDSet<data_id_t> read_set();
    IDSet<data_id_t> write_set();
    IDSet<job_id_t> before_set();
    IDSet<job_id_t> after_set();
    Parameter params();

  private:
    std::string job_name_;
    ID<job_id_t> job_id_;
    IDSet<data_id_t> read_set_;
    IDSet<data_id_t> write_set_;
    IDSet<job_id_t> before_set_;
    IDSet<job_id_t> after_set_;
    Parameter params_;
};



class CreateDataCommand : public SchedulerCommand {
  public:
    CreateDataCommand();
    CreateDataCommand(const ID<job_id_t>& job_id,
        const std::string& data_name, const ID<data_id_t>& data_id,
        const IDSet<job_id_t>& before, const IDSet<job_id_t>& after);
    ~CreateDataCommand();

    virtual std::string toString();
    virtual std::string toStringWTags();
    ID<job_id_t> job_id();
    std::string data_name();
    ID<data_id_t> data_id();
    IDSet<job_id_t> before_set();
    IDSet<job_id_t> after_set();

  private:
    ID<job_id_t> job_id_;
    std::string data_name_;
    ID<data_id_t> data_id_;
    IDSet<job_id_t> before_set_;
    IDSet<job_id_t> after_set_;
};


class RemoteCopySendCommand : public SchedulerCommand {
  public:
    RemoteCopySendCommand();
    RemoteCopySendCommand(const ID<job_id_t>& job_id,
        const ID<job_id_t>& receive_job_id,
        const ID<data_id_t>& from_data_id,
        const ID<worker_id_t>& to_worker_id,
        const std::string to_ip, const ID<port_t>& to_port,
        const IDSet<job_id_t>& before, const IDSet<job_id_t>& after);
    ~RemoteCopySendCommand();

    virtual std::string toString();
    virtual std::string toStringWTags();
    ID<job_id_t> job_id();
    ID<job_id_t> receive_job_id();
    ID<data_id_t> from_data_id();
    ID<worker_id_t> to_worker_id();
    std::string to_ip();
    ID<port_t> to_port();
    IDSet<job_id_t> before_set();
    IDSet<job_id_t> after_set();

  private:
    ID<job_id_t> job_id_;
    ID<job_id_t> receive_job_id_;
    ID<data_id_t> from_data_id_;
    ID<worker_id_t> to_worker_id_;
    std::string to_ip_;
    ID<port_t> to_port_;
    IDSet<job_id_t> before_set_;
    IDSet<job_id_t> after_set_;
};

class RemoteCopyReceiveCommand : public SchedulerCommand {
  public:
    RemoteCopyReceiveCommand();
    RemoteCopyReceiveCommand(const ID<job_id_t>& job_id,
        const ID<data_id_t>& to_data_id,
        const IDSet<job_id_t>& before, const IDSet<job_id_t>& after);
    ~RemoteCopyReceiveCommand();

    virtual std::string toString();
    virtual std::string toStringWTags();
    ID<job_id_t> job_id();
    ID<data_id_t> to_data_id();
    IDSet<job_id_t> before_set();
    IDSet<job_id_t> after_set();

  private:
    ID<job_id_t> job_id_;
    ID<data_id_t> to_data_id_;
    IDSet<job_id_t> before_set_;
    IDSet<job_id_t> after_set_;
};





class LocalCopyCommand : public SchedulerCommand {
  public:
    LocalCopyCommand();
    LocalCopyCommand(const ID<job_id_t>& jog_id,
        const ID<data_id_t>& from_data_id,
        const ID<data_id_t>& to_data_id,
        const IDSet<job_id_t>& before, const IDSet<job_id_t>& after);
    ~LocalCopyCommand();

    virtual std::string toString();
    virtual std::string toStringWTags();
    ID<job_id_t> job_id();
    ID<data_id_t> from_data_id();
    ID<data_id_t> to_data_id();
    IDSet<job_id_t> before_set();
    IDSet<job_id_t> after_set();

  private:
    ID<job_id_t> job_id_;
    ID<data_id_t> from_data_id_;
    ID<data_id_t> to_data_id_;
    IDSet<job_id_t> before_set_;
    IDSet<job_id_t> after_set_;
};




class JobDoneCommand : public SchedulerCommand {
  public:
    JobDoneCommand();
    JobDoneCommand(const ID<job_id_t>& job_id,
        const IDSet<job_id_t>& after_set,
        const Parameter& params);
    ~JobDoneCommand();

    virtual std::string toString();
    virtual std::string toStringWTags();
    ID<job_id_t> job_id();
    IDSet<job_id_t> after_set();
    Parameter params();


  private:
    ID<job_id_t> job_id_;
    IDSet<job_id_t> after_set_;
    Parameter params_;
};

class KillJobCommand : public SchedulerCommand {
  public:
    KillJobCommand();
    KillJobCommand(job_id_t job_id, Parameter params);
    ~KillJobCommand();

    virtual std::string toString();
    virtual std::string toStringWTags();
    job_id_t job_id();
    Parameter params();

  private:
    job_id_t job_id_;
    Parameter params_;
};

class PauseJobCommand : public SchedulerCommand {
  public:
    PauseJobCommand();
    PauseJobCommand(job_id_t job_id, Parameter params);
    ~PauseJobCommand();

    virtual std::string toString();
    virtual std::string toStringWTags();
    job_id_t job_id();
    Parameter params();

  private:
    job_id_t job_id_;
    Parameter params_;
};

class DefineDataCommand : public SchedulerCommand {
  public:
    DefineDataCommand();
    DefineDataCommand(const std::string& data_name,
                      const ID<data_id_t>& data_id,
                      const ID<partition_id_t>& partition_id,
                      const IDSet<partition_id_t>& neighbor_partition,
                      const Parameter& params);
    ~DefineDataCommand();

    virtual std::string toString();
    virtual std::string toStringWTags();
    std::string data_name();
    ID<data_id_t> data_id();
    ID<partition_id_t> partition_id();
    IDSet<partition_id_t> neighbor_partitions();
    Parameter params();

  private:
    std::string data_name_;
    ID<data_id_t> data_id_;
    ID<partition_id_t> partition_id_;
    IDSet<partition_id_t> neighbor_partitions_;
    Parameter params_;
};

class DefinePartitionCommand : public SchedulerCommand {
 public:
  DefinePartitionCommand(const ID<partition_id_t>& partition_id,
                         const GeometricRegion& r,
                         const Parameter& params);
  ~DefinePartitionCommand();

  virtual std::string toString();
  virtual std::string toStringWTags();

  ID<partition_id_t> partition_id();
  const GeometricRegion* region();
  Parameter params();

 private:
  ID<partition_id_t> id_;
  GeometricRegion region_;
  Parameter params_;
};

/*

class DataCreatedCommand : public SchedulerCommand {
  public:
    DataCreatedCommand();
    DataCreatedCommand(data_id_t data_id, std::string params);
    ~DataCreatedCommand();

    virtual std::string toString();
    virtual std::string toStringWTags();
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
    virtual std::string toStringWTags();
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
    virtual std::string toStringWTags();
    data_id_t data_id_from();
    data_id_t data_id_to();
    std::string params();

  private:
    data_id_t data_id_from_;
    data_id_t data_id_to_;
    // TODO(omidm): add the information of the remote worker holding the data.
    std::string params_;
};

*/

}  // namespace nimbus

#endif  // NIMBUS_SHARED_SCHEDULER_COMMAND_H_
