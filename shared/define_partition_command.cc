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
  * Define partition command to define a geometrical region from application
  * point of view.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#include "shared/define_partition_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

DefinePartitionCommand::DefinePartitionCommand() {
  name_ = DEFINE_PARTITION_NAME;
  type_ = DEFINE_PARTITION;
}

DefinePartitionCommand::DefinePartitionCommand(const ID<partition_id_t>& part,
                                               const GeometricRegion& r,
                                               const Parameter& params):
  id_(part), region_(r), params_(params) {
  name_ = DEFINE_PARTITION_NAME;
  type_ = DEFINE_PARTITION;
}

DefinePartitionCommand::~DefinePartitionCommand() {}

SchedulerCommand* DefinePartitionCommand::Clone() {
  return new DefinePartitionCommand();
}

bool DefinePartitionCommand::Parse(const std::string& params) {
  std::cout << "WARNING: Has not been built yet!" << std::endl;
  return false;
}

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
