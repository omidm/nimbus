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
  * Parser for Nimbus scheduler protocol. 
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#ifndef NIMBUS_SHARED_PARSER_H_
#define NIMBUS_SHARED_PARSER_H_

#include <stdint.h>
#include <boost/tokenizer.hpp>
#include <boost/thread.hpp>
#include <iostream> // NOLINT
#include <fstream> // NOLINT
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <utility>
#include "shared/id.h"
#include "shared/idset.h"
#include "shared/nimbus_types.h"

namespace nimbus {

typedef std::set<std::string> CmSet;
typedef std::set<std::pair<std::string, SchedulerCommandType> > CommandSet;

void parseCommand(const std::string& string,
                  const CmSet& commandSet,
                  std::string& command,
                  std::vector<int>& arguments);

int parseCommandFile(const std::string& fname,
                     CmSet& cs);

/** Returns true if there was a valid command in the string,
    false if no valid command. */
bool parseCommandFromString(const std::string& input,
    std::string& command,
    std::vector<std::string>& parameters);

void parseParameterFromString(const std::string& input, std::string& tag,
    std::string& args, std::string& string_set);

void parseIDSetFromString(const std::string& input, std::set<uint64_t>& set);

void parseIDSetFromString(const std::string& input, std::set<uint32_t>& set);

bool isSet(const std::string& tag);

int countOccurence(std::string str, std::string substr);

// ********************************************************

bool ParseSchedulerCommand(const std::string& input,
    CommandSet* command_set,
    std::string& name, std::string& param_segment,
    SchedulerCommandType& command_type);

bool ParseSpawnJobCommand(const std::string& input,
    std::string& job_name,
    IDSet<job_id_t>& job_id,
    IDSet<data_id_t>& read, IDSet<data_id_t>& write,
    IDSet<job_id_t>& before, IDSet<job_id_t>& after,
    JobType& job_type, std::string& params);

bool ParseSpawnComputeJobCommand(const std::string& input,
    std::string& job_name,
    ID<job_id_t>& job_id,
    IDSet<data_id_t>& read, IDSet<data_id_t>& write,
    IDSet<job_id_t>& before, IDSet<job_id_t>& after,
    std::string& params);

bool ParseComputeJobCommand(const std::string& input,
    std::string& job_name,
    ID<job_id_t>& job_id,
    IDSet<data_id_t>& read, IDSet<data_id_t>& write,
    IDSet<job_id_t>& before, IDSet<job_id_t>& after,
    std::string& params);

bool ParseCreateDataCommand(const std::string& input,
    ID<job_id_t>& job_id,
    std::string& data_name, ID<data_id_t>& data_id,
    IDSet<job_id_t>& before, IDSet<job_id_t>& after);

bool ParseRemoteCopyCommand(const std::string& input,
    ID<job_id_t>& job_id,
    ID<data_id_t>& from_data_id,
    ID<data_id_t>& to_data_id,
    ID<worker_id_t>& to_worker_id,
    std::string& to_ip,
    ID<port_t>& to_port,
    IDSet<job_id_t>& before, IDSet<job_id_t>& after);

bool ParseLocalCopyCommand(const std::string& input,
    ID<job_id_t>& job_id,
    ID<data_id_t>& from_data_id,
    ID<data_id_t>& to_data_id,
    IDSet<job_id_t>& before, IDSet<job_id_t>& after);

bool ParseDefineDataCommand(const std::string& input,
    std::string& data_name,
    IDSet<data_id_t>& data_id,
    IDSet<partition_t>& partition_id,
    IDSet<partition_t>& neighbor_partitions,
    std::string& params);

bool ParseHandshakeCommand(const std::string& input,
    ID<worker_id_t>& worker_id,
    std::string& ip, ID<port_t>& port);

bool ParseJobDoneCommand(const std::string& input,
    ID<job_id_t>& job_id,
    IDSet<job_id_t>& after_set,
    std::string& params);

bool ParseIDSet(const std::string& input, std::set<uint64_t>& set);

bool ParseIDSet(const std::string& input, std::set<uint32_t>& set);

bool ParseID(const std::string& input, uint64_t& elem);

bool ParseID(const std::string& input, uint32_t& elem);

bool ParseWorkerDataHeader(const std::string& input,
    job_id_t& job_id, size_t& data_length);









}  // namespace nimbus
#endif  // NIMBUS_SHARED_PARSER_H_
