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
  * Define data command to define a logical region from the application point
  * of view.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/define_data_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

DefineDataCommand::DefineDataCommand() {
  name_ = DEFINE_DATA_NAME;
  type_ = DEFINE_DATA;
}

DefineDataCommand::DefineDataCommand(const std::string& data_name,
    const ID<logical_data_id_t>& logical_data_id,
    const ID<partition_id_t>& partition_id,
    const IDSet<partition_id_t>& neighbor_partitions,
    const Parameter& params)
: data_name_(data_name), logical_data_id_(logical_data_id),
  partition_id_(partition_id),
  neighbor_partitions_(neighbor_partitions),
  params_(params) {
  name_ = DEFINE_DATA_NAME;
  type_ = DEFINE_DATA;
}

DefineDataCommand::~DefineDataCommand() {
}

SchedulerCommand* DefineDataCommand::Clone() {
  return new DefineDataCommand();
}

bool DefineDataCommand::Parse(const std::string& params) {
  int num = 5;

  char_separator<char> separator(" \n\t\r");
  tokenizer<char_separator<char> > tokens(params, separator);
  tokenizer<char_separator<char> >::iterator iter = tokens.begin();
  for (int i = 0; i < num; i++) {
    if (iter == tokens.end()) {
      std::cout << "ERROR: DefineDataCommand has only " << i <<
        " parameters (expected " << num << ")." << std::endl;
      return false;
    }
    iter++;
  }
  if (iter != tokens.end()) {
    std::cout << "ERROR: DefineDataCommand has more than "<<
      num << " parameters." << std::endl;
    return false;
  }

  iter = tokens.begin();
  data_name_ = *iter;

  iter++;
  if (!logical_data_id_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid logical data id." << std::endl;
    return false;
  }

  iter++;
  if (!partition_id_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid partiiton id." << std::endl;
    return false;
  }

  iter++;
  if (!neighbor_partitions_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid partition neighbor set." << std::endl;
    return false;
  }

  iter++;
  if (!params_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid parameter." << std::endl;
    return false;
  }

  return true;
}

std::string DefineDataCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (data_name_ + " ");
  str += (logical_data_id_.toString() + " ");
  str += (partition_id_.toString() + " ");
  str += (neighbor_partitions_.toString() + " ");
  str += params_.toString();

  return str;
}

std::string DefineDataCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("name:" + data_name_ + " ");
  str += ("logical-id:" + logical_data_id_.toString() + " ");
  str += ("partition-id:" + partition_id_.toString() + " ");
  str += ("neighbor-partitions:" + neighbor_partitions_.toString() + " ");
  str += ("params:" + params_.toString());
  return str;
}

std::string DefineDataCommand::data_name() {
  return data_name_;
}


ID<logical_data_id_t> DefineDataCommand::logical_data_id() {
  return logical_data_id_;
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


