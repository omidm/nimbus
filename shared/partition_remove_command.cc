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
  * add partition command to add partition definition to the worker from
  * scheduler.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/partition_remove_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

PartitionRemoveCommand::PartitionRemoveCommand() {
  name_ = PARTITION_REMOVE_NAME;
  type_ = PARTITION_REMOVE;
}

PartitionRemoveCommand::PartitionRemoveCommand(const ID<partition_id_t>& part,
                                               const GeometricRegion& r):
  id_(part), region_(r) {
  name_ = PARTITION_REMOVE_NAME;
  type_ = PARTITION_REMOVE;
}

PartitionRemoveCommand::~PartitionRemoveCommand() {}

SchedulerCommand* PartitionRemoveCommand::Clone() {
  return new PartitionRemoveCommand();
}

bool PartitionRemoveCommand::Parse(const std::string& params) {
  int num = 2;

  char_separator<char> separator(" \n\t\r");
  tokenizer<char_separator<char> > tokens(params, separator);
  tokenizer<char_separator<char> >::iterator iter = tokens.begin();
  for (int i = 0; i < num; i++) {
    if (iter == tokens.end()) {
      std::cout << "ERROR: PartitionRemoveCommand has only " << i <<
        " parameters (expected " << num << ")." << std::endl;
      return false;
    }
    iter++;
  }
  if (iter != tokens.end()) {
    std::cout << "ERROR: PartitionRemoveCommand has more than "<<
      num << " parameters." << std::endl;
    return false;
  }

  iter = tokens.begin();
  if (!id_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid partition id." << std::endl;
    return false;
  }

  iter++;
  if (!region_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid region." << std::endl;
    return false;
  }

  return true;
}

std::string PartitionRemoveCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (id_.toString() + " ");
  str += region_.toString();
  return str;
}
std::string PartitionRemoveCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("id:" + id_.toString() + " ");
  str += ("region:" + region_.toString());
  return str;
}

ID<partition_id_t> PartitionRemoveCommand::id() {
  return id_;
}

const GeometricRegion* PartitionRemoveCommand::region() {
  return &region_;
}

