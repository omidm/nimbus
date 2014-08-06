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
  * super class SchedulerCommand is inherited by its children implemented in
  * other files. Each child represents a specific command exchanged between
  * scheduler and worker.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SHARED_SCHEDULER_COMMAND_H_
#define NIMBUS_SHARED_SCHEDULER_COMMAND_H_


#include <boost/tokenizer.hpp>
#include <sstream> // NOLINT
#include <string>
#include <vector>
#include <list>
#include <map>
#include <set>
#include "shared/id.h"
#include "shared/idset.h"
#include "shared/dbg.h"
#include "shared/escaper.h"
#include "shared/parameter.h"
#include "shared/geometric_region.h"
#include "shared/nimbus_types.h"
#include "shared/protobuf_compiled/commands.pb.h"

namespace nimbus {

class SchedulerCommand {
 public:
  SchedulerCommand();
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
    DEFINE_PARTITION,
    LDO_ADD,
    LDO_REMOVE,
    PARTITION_ADD,
    PARTITION_REMOVE,
    TERMINATE,
    PROFILE
  };

  typedef std::set<Type> TypeSet;
  typedef std::map<uint16_t, SchedulerCommand*> PrototypeTable;

  virtual SchedulerCommand* Clone();

  // Parsing is a bit tricky for historical reasons.
  // The first Parse method takes a string, and it assumes
  // that this string is only the *subtype* of command. So,
  // if this is a SubmitCopyJobCommand, the string should be
  // the result of SubmitCopyPBuf.ToString.
  // THe second Parse method takes a SchedulerPBuf, and it
  // assumes that this is a complete scheduler command, including
  // the type information. That is, expects it to a SubmitCopyPBuf
  // inside a SchedulerPBuf.
  virtual bool Parse(const std::string& param_segment);
  virtual bool Parse(const SchedulerPBuf& buf);
  virtual std::string toString();
  virtual std::string toStringWTags();
  virtual std::string name();
  virtual Type type();

  static std::string GetNameFromType(Type type);
  static bool GenerateSchedulerCommandChild(const std::string& input,
      PrototypeTable* command_table,
      SchedulerCommand*& generated_command);


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
  static const std::string LDO_ADD_NAME;
  static const std::string LDO_REMOVE_NAME;
  static const std::string PARTITION_ADD_NAME;
  static const std::string PARTITION_REMOVE_NAME;
  static const std::string TERMINATE_NAME;
  static const std::string PROFILE_NAME;

 private:
};


typedef std::vector<SchedulerCommand*> SchedulerCommandVector;
typedef std::list<SchedulerCommand*> SchedulerCommandList;


/*

class KillJobCommand : public SchedulerCommand {
  public:
    KillJobCommand();
    KillJobCommand(job_id_t job_id, Parameter params);
    ~KillJobCommand();

    virtual SchedulerCommand* Clone();
    virtual bool Parse(const std::string& param_segment);
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

    virtual SchedulerCommand* Clone();
    virtual bool Parse(const std::string& param_segment);
    virtual std::string toString();
    virtual std::string toStringWTags();
    job_id_t job_id();
    Parameter params();

  private:
    job_id_t job_id_;
    Parameter params_;
};





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
