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

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;


SchedulerCommand::SchedulerCommand() {
  name_ = BASE_NAME;
  type_ = BASE;
}

SchedulerCommand::~SchedulerCommand() {}

std::string SchedulerCommand::toString() {
  std::string rval = name_;
  rval += " PARAMETER_SEGMENT";
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

// worker_id_t SchedulerCommand::worker_id() {
//   return worker_id_;
// }

// void SchedulerCommand::set_worker_id(worker_id_t id) {
//   worker_id_ = id;
// }


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
const std::string SchedulerCommand::LDO_ADD_NAME = "ldoadd";
const std::string SchedulerCommand::LDO_REMOVE_NAME = "ldoremove";

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
    case LDO_ADD:
      str = LDO_ADD_NAME;
      break;
    case LDO_REMOVE:
      str = LDO_REMOVE_NAME;
      break;
  }
  return str;
}

SchedulerCommand* SchedulerCommand::Clone() {
  std::cout << "WARNING: Base Scheduler Command Cloned." << std::endl;
  return new SchedulerCommand();
}

bool SchedulerCommand::Parse(const std::string& param_segment) {
  std::cout << "WARNING: Base Scheduler Command Parsed." << std::endl;
  return false;
}

bool SchedulerCommand::ParseCommandType(const std::string& input,
    SchedulerCommand::PrototypeTable* command_table,
    SchedulerCommand*& generated_command,
    std::string& param_segment) {
  char_separator<char> separator(" \n\t\r");
  tokenizer<char_separator<char> > tokens(input, separator);
  tokenizer<char_separator<char> >::iterator iter = tokens.begin();
  if (iter == tokens.end()) {
    std::cout << "ERROR: Command is empty." << std::endl;
    return false;
  }
  std::string name = *iter;
  bool name_is_valid = false;
  SchedulerCommand::PrototypeTable::iterator itr = command_table->begin();
  for (; itr != command_table->end(); itr++) {
    if (name == (*itr)->name()) {
      name_is_valid = true;
      generated_command = (*itr)->Clone();
      break;
    }
  }
  if (!name_is_valid) {
    std::cout << "ERROR: Command name is unknown: " << name << std::endl;
    return false;
  }

  param_segment = input.substr(name.length());
  return true;
}

bool SchedulerCommand::GenerateSchedulerCommandChild(const std::string& input,
    SchedulerCommand::PrototypeTable* command_table,
    SchedulerCommand*& generated_command) {
  std::string param_segment;
  if (!ParseCommandType(input, command_table, generated_command, param_segment)) {
    std::cout << "ERROR: Could not detect valid scheduler command type." << std::endl;
    return false;
  } else {
    if (!generated_command->Parse(param_segment)) {
      std::cout << "ERROR: Could not parse valid " << generated_command->name()
        << "command" << std::endl;
      delete generated_command;
      return false;
    }

    return true;
  }
}
