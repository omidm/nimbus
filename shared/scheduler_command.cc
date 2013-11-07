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

#include "shared/scheduler_command.h"
#include "shared/dbg.h"

using namespace nimbus;  // NOLINT

SchedulerCommand::SchedulerCommand() {
  name_ = BASE_NAME;
  type_ = BASE;
}

// TODO(omidm): remove this.
SchedulerCommand::SchedulerCommand(std::string n,
    const CommandParameterList& p) {
  name_ = n;
  CommandParameterList::const_iterator  iter = p.begin();
  for (; iter != p.end(); ++iter)
    addParameter(iter->second);
}

// TODO(omidm): remove this.
#include "shared/parser.h"
SchedulerCommand::SchedulerCommand(std::string command) {
  std::vector<std::string> string_params;
  parseCommandFromString(command, name_, string_params);
  std::vector<std::string>::iterator  iter = string_params.begin();
  for (; iter != string_params.end(); ++iter) {
    dbg(DBG_USR1, "Adding parameter %s\n", (*iter).c_str());
    addParameter(CommandParameter(*iter));
  }
}

SchedulerCommand::~SchedulerCommand() {}

// TODO(omidm): remove this.
void SchedulerCommand::addParameter(CommandParameter cm) {
  parameters_[cm.name()] = cm;
}

std::string SchedulerCommand::toString() {
  std::string rval = name_;
  CommandParameterList::iterator iter = parameters_.begin();
  for (; iter != parameters_.end(); ++iter) {
    rval += " ";
    rval += (iter->second).toString();
  }
  return rval;
}

std::string SchedulerCommand::toStringWTags() {
  return toString();
}

std::string SchedulerCommand::name() {
  return name_;
}

SchedulerCommand::Type SchedulerCommand::type() {
  return type_;
}

CommandParameterList* SchedulerCommand::parameters() {
  return &parameters_;
}

// worker_id_t SchedulerCommand::worker_id() {
//   return worker_id_;
// }

// void SchedulerCommand::set_worker_id(worker_id_t id) {
//   worker_id_ = id;
// }

CommandParameter::CommandParameter() {
  name_ = "empty-field";
}

CommandParameter::CommandParameter(std::string n,
    std::string v, const IDSet<job_id_t>& s) {
  name_ = n;
  value_ = v;
  identifier_set_ = s;
}

CommandParameter::CommandParameter(std::string parameter) {
  std::string string_set;
  parseParameterFromString(parameter, name_, value_, string_set);
  if (isSet(string_set))
    identifier_set_ = IDSet<job_id_t>(string_set);
}

CommandParameter::~CommandParameter() {}

std::string CommandParameter::toString() {
  std::string rval = name_;
  rval += ":";
  if (value_ == "")
    rval += identifier_set_.toString();
  else
    rval += value_;

  return rval;
}

std::string CommandParameter::name() {
  return name_;
}

std::string CommandParameter::value() {
  return value_;
}

IDSet<job_id_t>* CommandParameter::identifier_set() {
  return &identifier_set_;
}

// ************************************************

const std::string SchedulerCommand::BASE_NAME = "nop";
const std::string SchedulerCommand::SPAWN_JOB_NAME = "spawnjob";
const std::string SchedulerCommand::SPAWN_COMPUTE_JOB_NAME = "spawncomputejob";
const std::string SchedulerCommand::SPAWN_COPY_JOB_NAME = "spawncopyjob";
const std::string SchedulerCommand::DEFINE_DATA_NAME = "definedata";
const std::string SchedulerCommand::HANDSHAKE_NAME = "handshake";
const std::string SchedulerCommand::JOB_DONE_NAME = "jobdone";
const std::string SchedulerCommand::COMPUTE_JOB_NAME = "computejob";
const std::string SchedulerCommand::CREATE_DATA_NAME = "createdata";
const std::string SchedulerCommand::REMOTE_COPY_SEND_NAME = "remotecopysend";
const std::string SchedulerCommand::REMOTE_COPY_RECEIVE_NAME = "remotecopyreceive";
const std::string SchedulerCommand::LOCAL_COPY_NAME = "localcopy";
const std::string SchedulerCommand::DEFINE_PARTITION_NAME = "definepartition";

std::string SchedulerCommand::GetNameFromType(SchedulerCommand::Type type) {
  std::string str;
  switch (type) {
    case BASE:
      str = BASE_NAME;
      break;
    case SPAWN_JOB:
      str = SPAWN_JOB_NAME;
      break;
    case SPAWN_COMPUTE_JOB:
      str = SPAWN_COMPUTE_JOB_NAME;
      break;
    case SPAWN_COPY_JOB:
      str = SPAWN_COPY_JOB_NAME;
      break;
    case DEFINE_DATA:
      str = DEFINE_DATA_NAME;
      break;
    case HANDSHAKE:
      str = HANDSHAKE_NAME;
      break;
    case JOB_DONE:
      str = JOB_DONE_NAME;
      break;
    case COMPUTE_JOB:
      str = COMPUTE_JOB_NAME;
      break;
    case CREATE_DATA:
      str = CREATE_DATA_NAME;
      break;
    case REMOTE_COPY_SEND:
      str = REMOTE_COPY_SEND_NAME;
      break;
    case REMOTE_COPY_RECEIVE:
      str = REMOTE_COPY_RECEIVE_NAME;
      break;
    case LOCAL_COPY:
      str = LOCAL_COPY_NAME;
      break;
    case DEFINE_PARTITION:
      str = DEFINE_PARTITION_NAME;
      break;
  }
  return str;
}

bool SchedulerCommand::GenerateSchedulerCommandChild(const std::string& input,
    SchedulerCommand::TypeSet* command_set,
    SchedulerCommand*& generated) {
  std::string name, param_segment;
  SchedulerCommand::Type type;
  if (!ParseSchedulerCommand(input, command_set, name, param_segment, type)) {
    std::cout << "ERROR: Could not detect valid scheduler command." << std::endl;
    return false;
  } else {
    if (type == SPAWN_JOB) {
      std::string job_name;
      IDSet<job_id_t> job_id;
      IDSet<data_id_t> read;
      IDSet<data_id_t> write;
      IDSet<job_id_t> before;
      IDSet<job_id_t> after;
      JobType job_type;
      Parameter params;

      bool cond = ParseSpawnJobCommand(param_segment, job_name, job_id,
          read, write, before, after, job_type, params);
      if (!cond) {
        std::cout << "ERROR: Could not detect valid spawnjob." << std::endl;
        return false;
      } else {
        generated = new SpawnJobCommand(job_name, job_id,
            read, write, before, after, job_type, params);
      }
    } else if (type == SPAWN_COMPUTE_JOB) {
        std::string job_name;
        ID<job_id_t> job_id;
        IDSet<data_id_t> read;
        IDSet<data_id_t> write;
        IDSet<job_id_t> before;
        IDSet<job_id_t> after;
        Parameter params;

        bool cond = ParseSpawnComputeJobCommand(param_segment, job_name, job_id,
            read, write, before, after, params);
        if (!cond) {
          std::cout << "ERROR: Could not detect valid spawncomputejob." << std::endl;
          return false;
        } else {
          generated = new SpawnComputeJobCommand(job_name, job_id,
              read, write, before, after, params);
        }
    } else if (type == SPAWN_COPY_JOB) {
      ID<job_id_t> job_id;
      ID<data_id_t> from_id;
      ID<data_id_t> to_id;
      IDSet<job_id_t> before;
      IDSet<job_id_t> after;
      Parameter params;

      bool cond = ParseSpawnCopyJobCommand(param_segment, job_id,
          from_id, to_id, before, after, params);
      if (!cond) {
        std::cout << "ERROR: Could not detect valid spawncopyjob." << std::endl;
        return false;
      } else {
        generated = new SpawnCopyJobCommand(job_id,
            from_id, to_id, before, after, params);
      }
    } else if (type == COMPUTE_JOB) {
      std::string job_name;
      ID<job_id_t> job_id;
      IDSet<data_id_t> read;
      IDSet<data_id_t> write;
      IDSet<job_id_t> before;
      IDSet<job_id_t> after;
      Parameter params;

      bool cond = ParseComputeJobCommand(param_segment, job_name, job_id,
          read, write, before, after, params);
      if (!cond) {
        std::cout << "ERROR: Could not detect valid computejob." << std::endl;
        return false;
      } else {
        generated = new ComputeJobCommand(job_name, job_id,
            read, write, before, after, params);
      }
    } else if (type == CREATE_DATA) {
      std::string data_name;
      ID<job_id_t> job_id;
      ID<data_id_t> data_id;
      IDSet<job_id_t> before;
      IDSet<job_id_t> after;

      bool cond = ParseCreateDataCommand(param_segment, job_id,
          data_name, data_id, before, after);
      if (!cond) {
        std::cout << "ERROR: Could not detect valid createdata." << std::endl;
        return false;
      } else {
        generated = new CreateDataCommand(job_id,
            data_name, data_id, before, after);
      }
    } else if (type == REMOTE_COPY_SEND) {
      ID<job_id_t> job_id;
      ID<job_id_t> receive_job_id;
      ID<data_id_t> from_data_id;
      ID<worker_id_t> to_worker_id;
      std::string to_ip;
      ID<port_t> to_port;
      IDSet<job_id_t> before;
      IDSet<job_id_t> after;

      bool cond = ParseRemoteCopySendCommand(param_segment, job_id, receive_job_id,
          from_data_id, to_worker_id, to_ip, to_port, before, after);
      if (!cond) {
        std::cout << "ERROR: Could not detect valid remotecopysend." << std::endl;
        return false;
      } else {
        generated = new RemoteCopySendCommand(job_id, receive_job_id,
            from_data_id, to_worker_id, to_ip, to_port, before, after);
      }
    } else if (type == REMOTE_COPY_RECEIVE) {
      ID<job_id_t> job_id;
      ID<data_id_t> to_data_id;
      IDSet<job_id_t> before;
      IDSet<job_id_t> after;

      bool cond = ParseRemoteCopyReceiveCommand(param_segment, job_id,
          to_data_id, before, after);
      if (!cond) {
        std::cout << "ERROR: Could not detect valid remotecopyreceive." << std::endl;
        return false;
      } else {
        generated = new RemoteCopyReceiveCommand(job_id,
            to_data_id, before, after);
      }
    } else if (type == LOCAL_COPY) {
      ID<job_id_t> job_id;
      ID<data_id_t> from_data_id;
      ID<data_id_t> to_data_id;
      IDSet<job_id_t> before;
      IDSet<job_id_t> after;

      bool cond = ParseLocalCopyCommand(param_segment, job_id, from_data_id,
          to_data_id, before, after);
      if (!cond) {
        std::cout << "ERROR: Could not detect valid localcopy." << std::endl;
        return false;
      } else {
        generated = new LocalCopyCommand(job_id, from_data_id,
            to_data_id, before, after);
      }
    } else if (type == DEFINE_DATA) {
      std::string data_name;
      ID<data_id_t> data_id;
      ID<partition_id_t> partition_id;
      IDSet<partition_id_t> neighbor_partitions;
      Parameter params;

      bool cond = ParseDefineDataCommand(param_segment, data_name,
          data_id, partition_id, neighbor_partitions, params);
      if (!cond) {
        std::cout << "ERROR: Could not detect valid definedata." << std::endl;
        return false;
      } else {
        generated = new DefineDataCommand(data_name, data_id,
            partition_id, neighbor_partitions, params);
      }
    } else if (type == HANDSHAKE) {
      ID<worker_id_t> worker_id;
      std::string ip;
      ID<port_t> port;

      bool cond = ParseHandshakeCommand(param_segment, worker_id, ip, port);

      if (!cond) {
        std::cout << "ERROR: Could not detect valid handshake." << std::endl;
        return false;
      } else {
        generated = new HandshakeCommand(worker_id, ip, port);
      }
    } else if (type == JOB_DONE) {
      ID<job_id_t> job_id;
      IDSet<job_id_t> after_set;
      Parameter params;

      bool cond = ParseJobDoneCommand(param_segment, job_id, after_set, params);

      if (!cond) {
        std::cout << "ERROR: Could not detect valid jobdone." << std::endl;
        return false;
      } else {
        generated = new JobDoneCommand(job_id, after_set, params);
      }
    } else {
      std::cout << "ERROR: Unknown command." << std::endl;
      return false;
    }
    return true;
  }
}


HandshakeCommand::HandshakeCommand() {
  name_ = HANDSHAKE_NAME;
  type_ = HANDSHAKE;
}

HandshakeCommand::HandshakeCommand(const ID<worker_id_t>& worker_id,
    const std::string& ip, const ID<port_t>& port)
: worker_id_(worker_id), ip_(ip), port_(port) {
  name_ = HANDSHAKE_NAME;
  type_ = HANDSHAKE;
}

HandshakeCommand::~HandshakeCommand() {
}

std::string HandshakeCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (worker_id_.toString() + " ");
  str += (ip_ + " ");
  str += port_.toString();
  return str;
}

std::string HandshakeCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("worker-id:" + worker_id_.toString() + " ");
  str += ("ip:" + ip_ + " ");
  str += ("port:" + port_.toString());
  return str;
}


ID<worker_id_t> HandshakeCommand::worker_id() {
  return worker_id_;
}

std::string HandshakeCommand::ip() {
  return ip_;
}

ID<port_t> HandshakeCommand::port() {
  return port_;
}




SpawnJobCommand::SpawnJobCommand() {
  name_ = SPAWN_JOB_NAME;
  type_ = SPAWN_JOB;
}

SpawnJobCommand::SpawnJobCommand(const std::string& job_name,
    const IDSet<job_id_t>& job_id,
    const IDSet<data_id_t>& read, const IDSet<data_id_t>& write,
    const IDSet<job_id_t>& before, const IDSet<job_id_t>& after,
    const JobType& job_type, const Parameter& params)
: job_name_(job_name), job_id_(job_id),
  read_set_(read), write_set_(write),
  before_set_(before), after_set_(after),
  job_type_(job_type), params_(params) {
  name_ = SPAWN_JOB_NAME;
  type_ = SPAWN_JOB;
}

SpawnJobCommand::~SpawnJobCommand() {
}

std::string SpawnJobCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (job_name_ + " ");
  str += (job_id_.toString() + " ");
  str += (read_set_.toString() + " ");
  str += (write_set_.toString() + " ");
  str += (before_set_.toString() + " ");
  str += (after_set_.toString() + " ");
  if (job_type_ == JOB_COMP)
    str += "COMP ";
  else
    str += "SYNC ";
  str += params_.toString();

  return str;
}

std::string SpawnJobCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("name:" + job_name_ + " ");
  str += ("id:" + job_id_.toString() + " ");
  str += ("read:" + read_set_.toString() + " ");
  str += ("write:" + write_set_.toString() + " ");
  str += ("before:" + before_set_.toString() + " ");
  str += ("after:" + after_set_.toString() + " ");
  if (job_type_ == JOB_COMP)
    str += "type:COMP ";
  else
    str += "type:SYNC ";
  str += ("params:" + params_.toString());
  return str;
}

std::string SpawnJobCommand::job_name() {
  return job_name_;
}

JobType SpawnJobCommand::job_type() {
  return job_type_;
}

IDSet<job_id_t> SpawnJobCommand::job_id() {
  return job_id_;
}

IDSet<data_id_t> SpawnJobCommand::read_set() {
  return read_set_;
}

IDSet<data_id_t> SpawnJobCommand::write_set() {
  return write_set_;
}

IDSet<job_id_t> SpawnJobCommand::after_set() {
  return after_set_;
}
IDSet<job_id_t> SpawnJobCommand::before_set() {
  return before_set_;
}

Parameter SpawnJobCommand::params() {
  return params_;
}


SpawnComputeJobCommand::SpawnComputeJobCommand() {
  name_ = SPAWN_COMPUTE_JOB_NAME;
  type_ = SPAWN_COMPUTE_JOB;
}

SpawnComputeJobCommand::SpawnComputeJobCommand(const std::string& job_name,
    const ID<job_id_t>& job_id,
    const IDSet<data_id_t>& read, const IDSet<data_id_t>& write,
    const IDSet<job_id_t>& before, const IDSet<job_id_t>& after,
    const Parameter& params)
: job_name_(job_name), job_id_(job_id),
  read_set_(read), write_set_(write),
  before_set_(before), after_set_(after),
  params_(params) {
  name_ = SPAWN_COMPUTE_JOB_NAME;
  type_ = SPAWN_COMPUTE_JOB;
}

SpawnComputeJobCommand::~SpawnComputeJobCommand() {
}

std::string SpawnComputeJobCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (job_name_ + " ");
  str += (job_id_.toString() + " ");
  str += (read_set_.toString() + " ");
  str += (write_set_.toString() + " ");
  str += (before_set_.toString() + " ");
  str += (after_set_.toString() + " ");
  str += params_.toString();

  return str;
}

std::string SpawnComputeJobCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("name:" + job_name_ + " ");
  str += ("id:" + job_id_.toString() + " ");
  str += ("read:" + read_set_.toString() + " ");
  str += ("write:" + write_set_.toString() + " ");
  str += ("before:" + before_set_.toString() + " ");
  str += ("after:" + after_set_.toString() + " ");
  str += ("params:" + params_.toString());
  return str;
}

std::string SpawnComputeJobCommand::job_name() {
  return job_name_;
}

ID<job_id_t> SpawnComputeJobCommand::job_id() {
  return job_id_;
}

IDSet<data_id_t> SpawnComputeJobCommand::read_set() {
  return read_set_;
}

IDSet<data_id_t> SpawnComputeJobCommand::write_set() {
  return write_set_;
}

IDSet<job_id_t> SpawnComputeJobCommand::after_set() {
  return after_set_;
}
IDSet<job_id_t> SpawnComputeJobCommand::before_set() {
  return before_set_;
}

Parameter SpawnComputeJobCommand::params() {
  return params_;
}


SpawnCopyJobCommand::SpawnCopyJobCommand() {
  name_ = SPAWN_COPY_JOB_NAME;
  type_ = SPAWN_COPY_JOB;
}

SpawnCopyJobCommand::SpawnCopyJobCommand(const ID<job_id_t>& job_id,
    const ID<data_id_t>& from_id, const ID<data_id_t>& to_id,
    const IDSet<job_id_t>& before, const IDSet<job_id_t>& after,
    const Parameter& params)
: job_id_(job_id),
  from_id_(from_id), to_id_(to_id),
  before_set_(before), after_set_(after),
  params_(params) {
  name_ = SPAWN_COPY_JOB_NAME;
  type_ = SPAWN_COPY_JOB;
}

SpawnCopyJobCommand::~SpawnCopyJobCommand() {
}

std::string SpawnCopyJobCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (job_id_.toString() + " ");
  str += (from_id_.toString() + " ");
  str += (to_id_.toString() + " ");
  str += (before_set_.toString() + " ");
  str += (after_set_.toString() + " ");
  str += params_.toString();

  return str;
}

std::string SpawnCopyJobCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("id:" + job_id_.toString() + " ");
  str += ("from:" + from_id_.toString() + " ");
  str += ("to:" + to_id_.toString() + " ");
  str += ("before:" + before_set_.toString() + " ");
  str += ("after:" + after_set_.toString() + " ");
  str += ("params:" + params_.toString());
  return str;
}

ID<job_id_t> SpawnCopyJobCommand::job_id() {
  return job_id_;
}

ID<data_id_t> SpawnCopyJobCommand::from_id() {
  return from_id_;
}

ID<data_id_t> SpawnCopyJobCommand::to_id() {
  return to_id_;
}

IDSet<job_id_t> SpawnCopyJobCommand::after_set() {
  return after_set_;
}
IDSet<job_id_t> SpawnCopyJobCommand::before_set() {
  return before_set_;
}

Parameter SpawnCopyJobCommand::params() {
  return params_;
}

ComputeJobCommand::ComputeJobCommand() {
  name_ = COMPUTE_JOB_NAME;
  type_ = COMPUTE_JOB;
}

ComputeJobCommand::ComputeJobCommand(const std::string& job_name,
    const ID<job_id_t>& job_id,
    const IDSet<data_id_t>& read, const IDSet<data_id_t>& write,
    const IDSet<job_id_t>& before, const IDSet<job_id_t>& after,
    const Parameter& params)
: job_name_(job_name), job_id_(job_id),
  read_set_(read), write_set_(write),
  before_set_(before), after_set_(after),
  params_(params) {
  name_ = COMPUTE_JOB_NAME;
  type_ = COMPUTE_JOB;
}

ComputeJobCommand::~ComputeJobCommand() {
}

std::string ComputeJobCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (job_name_ + " ");
  str += (job_id_.toString() + " ");
  str += (read_set_.toString() + " ");
  str += (write_set_.toString() + " ");
  str += (before_set_.toString() + " ");
  str += (after_set_.toString() + " ");
  str += params_.toString();

  return str;
}

std::string ComputeJobCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("name:" + job_name_ + " ");
  str += ("id:" + job_id_.toString() + " ");
  str += ("read:" + read_set_.toString() + " ");
  str += ("write:" + write_set_.toString() + " ");
  str += ("before:" + before_set_.toString() + " ");
  str += ("after:" + after_set_.toString() + " ");
  str += ("params:" + params_.toString());
  return str;
}

std::string ComputeJobCommand::job_name() {
  return job_name_;
}

ID<job_id_t> ComputeJobCommand::job_id() {
  return job_id_;
}

IDSet<data_id_t> ComputeJobCommand::read_set() {
  return read_set_;
}

IDSet<data_id_t> ComputeJobCommand::write_set() {
  return write_set_;
}

IDSet<job_id_t> ComputeJobCommand::after_set() {
  return after_set_;
}
IDSet<job_id_t> ComputeJobCommand::before_set() {
  return before_set_;
}

Parameter ComputeJobCommand::params() {
  return params_;
}




CreateDataCommand::CreateDataCommand() {
  name_ = CREATE_DATA_NAME;
  type_ = CREATE_DATA;
}

CreateDataCommand::CreateDataCommand(const ID<job_id_t>& job_id,
    const std::string& data_name, const ID<data_id_t>& data_id,
    const IDSet<job_id_t>& before, const IDSet<job_id_t>& after)
: job_id_(job_id),
  data_name_(data_name), data_id_(data_id),
  before_set_(before), after_set_(after) {
  name_ = CREATE_DATA_NAME;
  type_ = CREATE_DATA;
}

CreateDataCommand::~CreateDataCommand() {
}

std::string CreateDataCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (job_id_.toString() + " ");
  str += (data_name_ + " ");
  str += (data_id_.toString() + " ");
  str += (before_set_.toString() + " ");
  str += after_set_.toString();
  return str;
}

std::string CreateDataCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("job-id:" + job_id_.toString() + " ");
  str += ("name:" + data_name_ + " ");
  str += ("data-id:" + data_id_.toString() + " ");
  str += ("before:" + before_set_.toString() + " ");
  str += ("after:" + after_set_.toString());
  return str;
}

ID<job_id_t> CreateDataCommand::job_id() {
  return job_id_;
}

std::string CreateDataCommand::data_name() {
  return data_name_;
}

ID<data_id_t> CreateDataCommand::data_id() {
  return data_id_;
}

IDSet<job_id_t> CreateDataCommand::after_set() {
  return after_set_;
}
IDSet<job_id_t> CreateDataCommand::before_set() {
  return before_set_;
}


RemoteCopySendCommand::RemoteCopySendCommand() {
  name_ = REMOTE_COPY_SEND_NAME;
  type_ = REMOTE_COPY_SEND;
}

RemoteCopySendCommand::RemoteCopySendCommand(const ID<job_id_t>& job_id,
    const ID<job_id_t>& receive_job_id,
    const ID<data_id_t>& from_data_id,
    const ID<worker_id_t>& to_worker_id,
    const std::string to_ip, const ID<port_t>& to_port,
    const IDSet<job_id_t>& before, const IDSet<job_id_t>& after)
: job_id_(job_id),
  receive_job_id_(receive_job_id),
  from_data_id_(from_data_id),
  to_worker_id_(to_worker_id),
  to_ip_(to_ip), to_port_(to_port),
  before_set_(before), after_set_(after) {
  name_ = REMOTE_COPY_SEND_NAME;
  type_ = REMOTE_COPY_SEND;
}

RemoteCopySendCommand::~RemoteCopySendCommand() {
}

std::string RemoteCopySendCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (job_id_.toString() + " ");
  str += (receive_job_id_.toString() + " ");
  str += (from_data_id_.toString() + " ");
  str += (to_worker_id_.toString() + " ");
  str += (to_ip_ + " ");
  str += (to_port_.toString() + " ");
  str += (before_set_.toString() + " ");
  str += after_set_.toString();
  return str;
}

std::string RemoteCopySendCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("job-id:" + job_id_.toString() + " ");
  str += ("receive-job-id:" + receive_job_id_.toString() + " ");
  str += ("from-data-id:" + from_data_id_.toString() + " ");
  str += ("to-worker-id:" + to_worker_id_.toString() + " ");
  str += ("to-ip:" + to_ip_ + " ");
  str += ("to-port:" + to_port_.toString() + " ");
  str += ("before:" + before_set_.toString() + " ");
  str += ("after:" + after_set_.toString());
  return str;
}

ID<job_id_t> RemoteCopySendCommand::job_id() {
  return job_id_;
}

ID<job_id_t> RemoteCopySendCommand::receive_job_id() {
  return receive_job_id_;
}

ID<data_id_t> RemoteCopySendCommand::from_data_id() {
  return from_data_id_;
}

ID<worker_id_t> RemoteCopySendCommand::to_worker_id() {
  return to_worker_id_;
}

std::string RemoteCopySendCommand::to_ip() {
  return to_ip_;
}

ID<port_t> RemoteCopySendCommand::to_port() {
  return to_port_;
}

IDSet<job_id_t> RemoteCopySendCommand::after_set() {
  return after_set_;
}

IDSet<job_id_t> RemoteCopySendCommand::before_set() {
  return before_set_;
}



RemoteCopyReceiveCommand::RemoteCopyReceiveCommand() {
  name_ = REMOTE_COPY_RECEIVE_NAME;
  type_ = REMOTE_COPY_RECEIVE;
}

RemoteCopyReceiveCommand::RemoteCopyReceiveCommand(const ID<job_id_t>& job_id,
    const ID<data_id_t>& to_data_id,
    const IDSet<job_id_t>& before, const IDSet<job_id_t>& after)
: job_id_(job_id),
  to_data_id_(to_data_id),
  before_set_(before), after_set_(after) {
  name_ = REMOTE_COPY_RECEIVE_NAME;
  type_ = REMOTE_COPY_RECEIVE;
}

RemoteCopyReceiveCommand::~RemoteCopyReceiveCommand() {
}

std::string RemoteCopyReceiveCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (job_id_.toString() + " ");
  str += (to_data_id_.toString() + " ");
  str += (before_set_.toString() + " ");
  str += after_set_.toString();
  return str;
}

std::string RemoteCopyReceiveCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("job-id:" + job_id_.toString() + " ");
  str += ("to-data-id:" + to_data_id_.toString() + " ");
  str += ("before:" + before_set_.toString() + " ");
  str += ("after:" + after_set_.toString());
  return str;
}

ID<job_id_t> RemoteCopyReceiveCommand::job_id() {
  return job_id_;
}

ID<data_id_t> RemoteCopyReceiveCommand::to_data_id() {
  return to_data_id_;
}

IDSet<job_id_t> RemoteCopyReceiveCommand::after_set() {
  return after_set_;
}

IDSet<job_id_t> RemoteCopyReceiveCommand::before_set() {
  return before_set_;
}


LocalCopyCommand::LocalCopyCommand() {
  name_ = LOCAL_COPY_NAME;
  type_ = LOCAL_COPY;
}

LocalCopyCommand::LocalCopyCommand(const ID<job_id_t>& job_id,
    const ID<data_id_t>& from_data_id,
    const ID<data_id_t>& to_data_id,
    const IDSet<job_id_t>& before, const IDSet<job_id_t>& after)
: job_id_(job_id),
  from_data_id_(from_data_id),
  to_data_id_(to_data_id),
  before_set_(before), after_set_(after) {
  name_ = LOCAL_COPY_NAME;
  type_ = LOCAL_COPY;
}

LocalCopyCommand::~LocalCopyCommand() {
}

std::string LocalCopyCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (job_id_.toString() + " ");
  str += (from_data_id_.toString() + " ");
  str += (to_data_id_.toString() + " ");
  str += (before_set_.toString() + " ");
  str += after_set_.toString();
  return str;
}

std::string LocalCopyCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("job-id:" + job_id_.toString() + " ");
  str += ("from-data-id:" + from_data_id_.toString() + " ");
  str += ("to-data-id:" + to_data_id_.toString() + " ");
  str += ("before:" + before_set_.toString() + " ");
  str += ("after:" + after_set_.toString());
  return str;
}

ID<job_id_t> LocalCopyCommand::job_id() {
  return job_id_;
}

ID<data_id_t> LocalCopyCommand::from_data_id() {
  return from_data_id_;
}

ID<data_id_t> LocalCopyCommand::to_data_id() {
  return to_data_id_;
}

IDSet<job_id_t> LocalCopyCommand::after_set() {
  return after_set_;
}

IDSet<job_id_t> LocalCopyCommand::before_set() {
  return before_set_;
}






DefineDataCommand::DefineDataCommand() {
  name_ = DEFINE_DATA_NAME;
  type_ = DEFINE_DATA;
}

DefineDataCommand::DefineDataCommand(const std::string& data_name,
    const ID<data_id_t>& data_id,
    const ID<partition_id_t>& partition_id,
    const IDSet<partition_id_t>& neighbor_partitions,
    const Parameter& params)
: data_name_(data_name), data_id_(data_id),
  partition_id_(partition_id),
  neighbor_partitions_(neighbor_partitions),
  params_(params) {
  name_ = DEFINE_DATA_NAME;
  type_ = DEFINE_DATA;
}

DefineDataCommand::~DefineDataCommand() {
}

std::string DefineDataCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (data_name_ + " ");
  str += (data_id_.toString() + " ");
  str += (partition_id_.toString() + " ");
  str += (neighbor_partitions_.toString() + " ");
  str += params_.toString();

  return str;
}

std::string DefineDataCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("name:" + data_name_ + " ");
  str += ("id:" + data_id_.toString() + " ");
  str += ("partition-id:" + partition_id_.toString() + " ");
  str += ("neighbor-partitions:" + neighbor_partitions_.toString() + " ");
  str += ("params:" + params_.toString());
  return str;
}

std::string DefineDataCommand::data_name() {
  return data_name_;
}


ID<data_id_t> DefineDataCommand::data_id() {
  return data_id_;
}

ID<partition_id_t> DefineDataCommand::partition_id() {
  return partition_id_;
}

IDSet<partition_id_t> DefineDataCommand::neighbor_partitions() {
  return neighbor_partitions_;
}

Parameter DefineDataCommand::params() {
  return params_;
}




JobDoneCommand::JobDoneCommand() {
  name_ = JOB_DONE_NAME;
  type_ = JOB_DONE;
}

JobDoneCommand::JobDoneCommand(const ID<job_id_t>& job_id,
    const IDSet<job_id_t>& after_set,
    const Parameter& params)
: job_id_(job_id), after_set_(after_set), params_(params) {
  name_ = JOB_DONE_NAME;
  type_ = JOB_DONE;
}

JobDoneCommand::~JobDoneCommand() {
}

std::string JobDoneCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (job_id_.toString() + " ");
  str += (after_set_.toString() + " ");
  str += params_.toString();

  return str;
}

std::string JobDoneCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("id:" + job_id_.toString() + " ");
  str += ("after:" + after_set_.toString() + " ");
  str += ("params:" + params_.toString());
  return str;
}

ID<job_id_t> JobDoneCommand::job_id() {
  return job_id_;
}

IDSet<job_id_t> JobDoneCommand::after_set() {
  return after_set_;
}

Parameter JobDoneCommand::params() {
  return params_;
}

DefinePartitionCommand::DefinePartitionCommand(const ID<partition_id_t>& part,
                                               const GeometricRegion& r,
                                               const Parameter& params):
  id_(part), region_(r), params_(params) {
  name_ = DEFINE_PARTITION_NAME;
  type_ = DEFINE_PARTITION;
}
DefinePartitionCommand::~DefinePartitionCommand() {}

std::string DefinePartitionCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (id_.toString() + " ");
  str += (region_.toString() + " ");
  str += params_.toString();
  return str;
}
std::string DefinePartitionCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("id:" + id_.toString() + " ");
  str += ("region:" + region_.toString() + " ");
  str += ("params:" + params_.toString());
  return str;
}

ID<partition_id_t> DefinePartitionCommand::partition_id() {
  return id_;
}

const GeometricRegion* DefinePartitionCommand::region() {
  return &region_;
}

Parameter DefinePartitionCommand::params() {
  return params_;
}
