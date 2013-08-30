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

using namespace nimbus; // NOLINT

SchedulerCommand::SchedulerCommand() {
  name_ = "no-op";
}

SchedulerCommand::SchedulerCommand(std::string n,
    const CommandParameterList& p) {
  name_ = n;
  CommandParameterList::const_iterator  iter = p.begin();
  for (; iter != p.end(); ++iter)
    addParameter(iter->second);
}

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

bool SchedulerCommand::GenerateSchedulerCommandChild(const std::string& input,
    CommandSet* command_set,
    SchedulerCommand*& generated) {
  std::string name, param_segment;
  SchedulerCommandType type;
  if (!ParseSchedulerCommand(input, command_set, name, param_segment, type)) {
    std::cout << "ERROR: Could not detect valid scheduler command." << std::endl;
    return false;
  } else {
    if (type == COMMAND_SPAWN_JOB) {
      std::string job_name;
      IDSet<job_id_t> job_id;
      IDSet<data_id_t> read;
      IDSet<data_id_t> write;
      IDSet<job_id_t> before;
      IDSet<job_id_t> after;
      JobType job_type;
      std::string params;

      bool cond = ParseSpawnJobCommand(param_segment, job_name, job_id,
          read, write, before, after, job_type, params);
      if (!cond) {
        std::cout << "ERROR: Could not detect valid spawnjob." << std::endl;
        return false;
      } else {
        generated = new SpawnJobCommand(job_name, job_id,
            read, write, before, after, job_type, params);
      }
    } else if (type == COMMAND_DEFINE_DATA) {
      std::string data_name;
      IDSet<data_id_t> data_id;
      IDSet<partition_t> partition_id;
      IDSet<partition_t> neighbor_partitions;
      std::string params;

      bool cond = ParseDefineDataCommand(param_segment, data_name,
          data_id, partition_id, neighbor_partitions, params);
      if (!cond) {
        std::cout << "ERROR: Could not detect valid definedata." << std::endl;
        return false;
      } else {
        generated = new DefineDataCommand(data_name, data_id,
            partition_id, neighbor_partitions, params);
      }
    } else if (type == COMMAND_HANDSHAKE) {
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
    } else if (type == COMMAND_JOBDONE) {
      ID<job_id_t> job_id;
      std::string params;

      bool cond = ParseJobDoneCommand(param_segment, job_id, params);

      if (!cond) {
        std::cout << "ERROR: Could not detect valid jobdone." << std::endl;
        return false;
      } else {
        generated = new JobDoneCommand(job_id, params);
      }
    } else {
      std::cout << "ERROR: Unknown command." << std::endl;
      return false;
    }
    return true;
  }
}


HandshakeCommand::HandshakeCommand() {
  name_ = "handshake";
}

HandshakeCommand::HandshakeCommand(const ID<worker_id_t>& worker_id,
    const std::string& ip, const ID<port_t>& port)
: worker_id_(worker_id), ip_(ip), port_(port) {
    name_ = "handshake";
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
  str += ("worker-id:" + worker_id_.toString());
  str += ("ip:" + ip_ + " ");
  str += ("port" + port_.toString());
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
  name_ = "spawnjob";
}

SpawnJobCommand::SpawnJobCommand(const std::string& job_name,
    const IDSet<job_id_t>& job_id,
    const IDSet<data_id_t>& read, const IDSet<data_id_t>& write,
    const IDSet<job_id_t>& before, const IDSet<job_id_t>& after,
    const JobType& job_type, const std::string& params)
: job_name_(job_name), job_id_(job_id),
  read_set_(read), write_set_(write),
  before_set_(before), after_set_(after),
  job_type_(job_type), params_(params) {
    name_ = "spawnjob";
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
  str += params_;
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
  str += ("params:" + params_);
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

std::string SpawnJobCommand::params() {
  return params_;
}




DefineDataCommand::DefineDataCommand() {
  name_ = "definedata";
}

DefineDataCommand::DefineDataCommand(const std::string& data_name,
    const IDSet<data_id_t>& data_id,
    const IDSet<partition_t>& partition_id,
    const IDSet<partition_t>& neighbor_partitions,
    const std::string& params)
: data_name_(data_name), data_id_(data_id),
  partition_id_(partition_id),
  neighbor_partitions_(neighbor_partitions),
  params_(params) {
    name_ = "definedata";
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
  str += params_;
  return str;
}

std::string DefineDataCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("name:" + data_name_ + " ");
  str += ("id:" + data_id_.toString() + " ");
  str += ("partition-id:" + data_id_.toString() + " ");
  str += ("neighbor-partitions:" + neighbor_partitions_.toString() + " ");
  str += ("params:" + params_);
  return str;
}

std::string DefineDataCommand::data_name() {
  return data_name_;
}


IDSet<data_id_t> DefineDataCommand::data_id() {
  return data_id_;
}

IDSet<partition_t> DefineDataCommand::partition_id() {
  return data_id_;
}

IDSet<partition_t> DefineDataCommand::neighbor_partitions() {
  return neighbor_partitions_;
}

std::string DefineDataCommand::params() {
  return params_;
}




JobDoneCommand::JobDoneCommand() {
  name_ = "jobdone";
}

JobDoneCommand::JobDoneCommand(const ID<job_id_t>& job_id,
    const std::string& params)
: job_id_(job_id), params_(params) {
  name_ = "jobdone";
}

JobDoneCommand::~JobDoneCommand() {
}

std::string JobDoneCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (job_id_.toString() + " ");
  str += params_;
  return str;
}

std::string JobDoneCommand::toStringWTags() {
  std::string str;
  str += ("name:" + name_ + " ");
  str += ("id:" + job_id_.toString() + " ");
  str += ("params:" + params_);
  return str;
}
ID<job_id_t> JobDoneCommand::job_id() {
  return job_id_;
}

std::string JobDoneCommand::params() {
  return params_;
}


