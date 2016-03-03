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

#include "src/shared/scheduler_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;


SchedulerCommand::SchedulerCommand() {
  name_ = BASE_NAME;
  type_ = BASE;
}

SchedulerCommand::~SchedulerCommand() {}

std::string SchedulerCommand::ToNetworkData() {
  std::string rval = name_;
  rval += " PARAMETER_SEGMENT";
  return rval;
}

std::string SchedulerCommand::ToString() {
  return ToNetworkData();
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
const std::string SchedulerCommand::ADD_COMPUTE_NAME = "addcomputejob";
const std::string SchedulerCommand::ADD_COPY_NAME = "addcopyjob";
const std::string SchedulerCommand::SPAWN_JOB_GRAPH_NAME = "spawnjobgraph";
const std::string SchedulerCommand::SPAWN_COMPUTE_NAME = "spawncomputejob";
const std::string SchedulerCommand::SPAWN_COPY_NAME = "spawncopyjob";
const std::string SchedulerCommand::DEFINE_DATA_NAME = "definedata";
const std::string SchedulerCommand::HANDSHAKE_NAME = "handshake";
const std::string SchedulerCommand::JOB_DONE_NAME = "jobdone";
const std::string SchedulerCommand::EXECUTE_COMPUTE_NAME = "computejob";
const std::string SchedulerCommand::CREATE_DATA_NAME = "createdata";
const std::string SchedulerCommand::REMOTE_SEND_NAME = "remotecopysend";
const std::string SchedulerCommand::REMOTE_RECEIVE_NAME = "remotecopyreceive";
const std::string SchedulerCommand::LOCAL_COPY_NAME = "localcopy";
const std::string SchedulerCommand::DEFINE_PARTITION_NAME = "definepartition";
const std::string SchedulerCommand::LDO_ADD_NAME = "ldoadd";
const std::string SchedulerCommand::LDO_REMOVE_NAME = "ldoremove";
const std::string SchedulerCommand::PARTITION_ADD_NAME = "partitionadd";
const std::string SchedulerCommand::PARTITION_REMOVE_NAME = "partitionremove";
const std::string SchedulerCommand::TERMINATE_NAME = "terminate";
const std::string SchedulerCommand::PROFILE_NAME = "profile";
const std::string SchedulerCommand::START_TEMPLATE_NAME = "starttemplate";
const std::string SchedulerCommand::END_TEMPLATE_NAME = "endtemplate";
const std::string SchedulerCommand::DEFINED_TEMPLATE_NAME = "definedtemplate";
const std::string SchedulerCommand::SPAWN_TEMPLATE_NAME = "spawntemplate";
const std::string SchedulerCommand::SAVE_DATA_NAME = "savedata";
const std::string SchedulerCommand::LOAD_DATA_NAME = "loaddata";
const std::string SchedulerCommand::SAVE_DATA_JOB_DONE_NAME = "savedatajobdone";
const std::string SchedulerCommand::PREPARE_REWIND_NAME = "preparerewind";
const std::string SchedulerCommand::WORKER_DOWN_NAME = "workerdown";
const std::string SchedulerCommand::START_COMMAND_TEMPLATE_NAME = "startcommandtemplate";
const std::string SchedulerCommand::END_COMMAND_TEMPLATE_NAME = "endcommandtemplate";
const std::string SchedulerCommand::SPAWN_COMMAND_TEMPLATE_NAME = "spawncommandtemplate";
const std::string SchedulerCommand::REQUEST_STAT_NAME = "requeststat";
const std::string SchedulerCommand::RESPOND_STAT_NAME = "respondstat";
const std::string SchedulerCommand::PRINT_STAT_NAME = "printstat";
const std::string SchedulerCommand::MEGA_RCR_NAME = "megarcr";
const std::string SchedulerCommand::MEGA_JOB_DONE_NAME = "megajobdone";
const std::string SchedulerCommand::EXECUTE_COMBINE_NAME = "combinejob";

std::string SchedulerCommand::GetNameFromType(SchedulerCommand::Type type) {
  std::string str;
  switch (type) {
    case BASE:
      str = BASE_NAME;
      break;
    case ADD_COMPUTE:
      str = ADD_COMPUTE_NAME;
      break;
    case ADD_COPY:
      str = ADD_COPY_NAME;
      break;
    case SPAWN_JOB_GRAPH:
      str = SPAWN_JOB_GRAPH_NAME;
      break;
    case SPAWN_COMPUTE:
      str = SPAWN_COMPUTE_NAME;
      break;
    case SPAWN_COPY:
      str = SPAWN_COPY_NAME;
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
    case EXECUTE_COMPUTE:
      str = EXECUTE_COMPUTE_NAME;
      break;
    case CREATE_DATA:
      str = CREATE_DATA_NAME;
      break;
    case REMOTE_SEND:
      str = REMOTE_SEND_NAME;
      break;
    case REMOTE_RECEIVE:
      str = REMOTE_RECEIVE_NAME;
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
    case PARTITION_ADD:
      str = LDO_ADD_NAME;
      break;
    case PARTITION_REMOVE:
      str = LDO_REMOVE_NAME;
      break;
    case TERMINATE:
      str = TERMINATE_NAME;
      break;
    case PROFILE:
      str = PROFILE_NAME;
      break;
    case START_TEMPLATE:
      str = START_TEMPLATE_NAME;
      break;
    case END_TEMPLATE:
      str = END_TEMPLATE_NAME;
      break;
    case DEFINED_TEMPLATE:
      str = DEFINED_TEMPLATE_NAME;
      break;
    case SPAWN_TEMPLATE:
      str = SPAWN_TEMPLATE_NAME;
      break;
    case SAVE_DATA:
      str = SAVE_DATA_NAME;
      break;
    case LOAD_DATA:
      str = LOAD_DATA_NAME;
      break;
    case SAVE_DATA_JOB_DONE:
      str = SAVE_DATA_JOB_DONE_NAME;
      break;
    case PREPARE_REWIND:
      str = PREPARE_REWIND_NAME;
      break;
    case WORKER_DOWN:
      str = WORKER_DOWN_NAME;
      break;
    case START_COMMAND_TEMPLATE:
      str = START_COMMAND_TEMPLATE_NAME;
      break;
    case END_COMMAND_TEMPLATE:
      str = END_COMMAND_TEMPLATE_NAME;
      break;
    case SPAWN_COMMAND_TEMPLATE:
      str = SPAWN_COMMAND_TEMPLATE_NAME;
      break;
    case REQUEST_STAT:
      str = REQUEST_STAT_NAME;
      break;
    case RESPOND_STAT:
      str = RESPOND_STAT_NAME;
      break;
    case PRINT_STAT:
      str = PRINT_STAT_NAME;
      break;
    case MEGA_RCR:
      str = MEGA_RCR_NAME;
      break;
    case MEGA_JOB_DONE:
      str = MEGA_JOB_DONE_NAME;
      break;
    case EXECUTE_COMBINE:
      str = EXECUTE_COMBINE_NAME;
      break;
    default:
      std::cout << "Type did not found\n";
      exit(-1);
  }
  return str;
}

SchedulerCommand* SchedulerCommand::Clone() {
  std::cout << "WARNING: Base Scheduler Command Cloned." << std::endl;
  return new SchedulerCommand();
}

bool SchedulerCommand::Parse(const std::string& param_segment) {
  dbg(DBG_ERROR, "ERROR: Base Scheduler Command Parsed: should be parsing a subclass.\n");
  dbg(DBG_ERROR, "       A Parse(const std::string) hasn't been redefined.\n");  return false;
}

bool SchedulerCommand::Parse(const SchedulerPBuf& param_segment) {
  dbg(DBG_ERROR, "ERROR: Base Scheduler Command Parsed: should be parsing a subclass.\n");
  dbg(DBG_ERROR, "       A Parse(const SchedulerPBuf&) hasn't been redefined.\n");
  return false;
}

bool SchedulerCommand::GenerateSchedulerCommandChild(const std::string& input,
                                                     SchedulerCommand::PrototypeTable* command_table, // NOLINT
                                                     SchedulerCommand*& generated_command) {
  SchedulerPBuf pBuf;
  bool result = pBuf.ParseFromString(input);
  if (!result) {
    dbg(DBG_ERROR, "ERROR: could not parse command.\n");
    return false;
  }

  SchedulerCommand::PrototypeTable::iterator it = command_table->find((uint16_t)pBuf.type());
  if (it == command_table->end()) {
    dbg(DBG_ERROR, "ERROR: Could not find command type %d in command table.\n", (uint16_t)pBuf.type()); // NOLINT
    exit(-1);
  }

  SchedulerCommand* prototype = command_table->at((uint16_t)pBuf.type());
  generated_command = prototype->Clone();
  if (!generated_command->Parse(pBuf)) {
    dbg(DBG_ERROR, "ERROR: Could not parse valid %s command.\n", generated_command->name().c_str());
    delete generated_command;
    return false;
  }
  return true;
}
