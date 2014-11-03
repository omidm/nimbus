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
  * Object representation of a network command between a controller
  * and workers (either direction). The super class SchedulerCommand
  * is inherited by its children implemented in other files. Each
  * child represents a specific command. Commands are encoded for the
  * network using Google Protocol Buffers.
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
    ADD_COMPUTE      = SchedulerPBuf_Type_ADD_COMPUTE,
    ADD_COPY         = SchedulerPBuf_Type_ADD_COPY,
    SPAWN_JOB_GRAPH  = SchedulerPBuf_Type_SPAWN_JOB_GRAPH,
    SPAWN_COMPUTE    = SchedulerPBuf_Type_SPAWN_COMPUTE,
    SPAWN_COPY       = SchedulerPBuf_Type_SPAWN_COPY,
    DEFINE_DATA      = SchedulerPBuf_Type_DEFINE_DATA,
    HANDSHAKE        = SchedulerPBuf_Type_HANDSHAKE,
    JOB_DONE         = SchedulerPBuf_Type_JOB_DONE,
    EXECUTE_COMPUTE  = SchedulerPBuf_Type_EXECUTE_COMPUTE,
    CREATE_DATA      = SchedulerPBuf_Type_CREATE_DATA,
    REMOTE_SEND      = SchedulerPBuf_Type_REMOTE_SEND,
    REMOTE_RECEIVE   = SchedulerPBuf_Type_REMOTE_RECEIVE,
    LOCAL_COPY       = SchedulerPBuf_Type_LOCAL_COPY,
    DEFINE_PARTITION = SchedulerPBuf_Type_DEFINE_PARTITION,
    LDO_ADD          = SchedulerPBuf_Type_LDO_ADD,
    LDO_REMOVE       = SchedulerPBuf_Type_LDO_REMOVE,
    PARTITION_ADD    = SchedulerPBuf_Type_PARTITION_ADD,
    PARTITION_REMOVE = SchedulerPBuf_Type_PARTITION_REMOVE,
    TERMINATE        = SchedulerPBuf_Type_TERMINATE,
    PROFILE          = SchedulerPBuf_Type_PROFILE,
    START_TEMPLATE   = SchedulerPBuf_Type_START_TEMPLATE,
    END_TEMPLATE     = SchedulerPBuf_Type_END_TEMPLATE,
    DEFINED_TEMPLATE = SchedulerPBuf_Type_DEFINED_TEMPLATE
  };

  typedef std::set<Type> TypeSet;
  typedef std::map<uint16_t, SchedulerCommand*> PrototypeTable;
  typedef uint32_t length_field_t;

  virtual SchedulerCommand* Clone();

  // Parsing is a bit tricky for historical reasons.
  //
  // The first Parse method takes a string, and it assumes
  // that this string is only the *subtype* of command. So,
  // if this is a SubmitCopyJobCommand, the string should be
  // the result of SubmitCopyPBuf.ToString.
  //
  // The second Parse method takes a SchedulerPBuf, and it
  // assumes that this is a complete scheduler command, including
  // the type information. That is, expects it to a SubmitCopyPBuf
  // inside a SchedulerPBuf.
  //
  // This should be cleaned up and we should get rid of the first
  // Parse eventuall.
  virtual bool Parse(const std::string& param_segment);
  virtual bool Parse(const SchedulerPBuf& buf);

  // ToNetworkData and ToString are very different.
  // ToNetworkData returns an encoded series of bytes for an on-wire
  // format. ToString returns a human-readable ASCII string.
  virtual std::string ToNetworkData();
  virtual std::string ToString();
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
  static const std::string ADD_COMPUTE_NAME;
  static const std::string ADD_COPY_NAME;
  static const std::string SPAWN_JOB_GRAPH_NAME;
  static const std::string SPAWN_COMPUTE_NAME;
  static const std::string SPAWN_COPY_NAME;
  static const std::string DEFINE_DATA_NAME;
  static const std::string HANDSHAKE_NAME;
  static const std::string JOB_DONE_NAME;
  static const std::string EXECUTE_COMPUTE_NAME;
  static const std::string CREATE_DATA_NAME;
  static const std::string REMOTE_SEND_NAME;
  static const std::string REMOTE_RECEIVE_NAME;
  static const std::string LOCAL_COPY_NAME;
  static const std::string DEFINE_PARTITION_NAME;
  static const std::string LDO_ADD_NAME;
  static const std::string LDO_REMOVE_NAME;
  static const std::string PARTITION_ADD_NAME;
  static const std::string PARTITION_REMOVE_NAME;
  static const std::string TERMINATE_NAME;
  static const std::string PROFILE_NAME;
  static const std::string START_TEMPLATE_NAME;
  static const std::string END_TEMPLATE_NAME;
  static const std::string DEFINED_TEMPLATE_NAME;

 private:
};


typedef std::vector<SchedulerCommand*> SchedulerCommandVector;
typedef std::list<SchedulerCommand*> SchedulerCommandList;

}  // namespace nimbus

#endif  // NIMBUS_SHARED_SCHEDULER_COMMAND_H_
