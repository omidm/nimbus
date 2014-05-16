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
  * Profile command with memory usage statistics.  
  *
  * Author: Andrew Lim <alim16@stanford.edu>
  */

#include "shared/profile_command.h"
#include "boost/lexical_cast.hpp"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

ProfileCommand::ProfileCommand() {
  name_ = PROFILE_NAME;
  type_ = PROFILE;
}

ProfileCommand::ProfileCommand(const ID<worker_id_t>& worker_id,
    const uint64_t total_virtual,
    const uint64_t used_virtual,
    const uint64_t proc_virtual,
    const uint64_t total_physical,
    const uint64_t used_physical,
    const uint64_t proc_physical)
: worker_id_(worker_id), total_virtual_(total_virtual),
  used_virtual_(used_virtual), proc_virtual_(proc_virtual),
  total_physical_(total_physical), used_physical_(used_physical),
  proc_physical_(proc_physical) {
  name_ = PROFILE_NAME;
  type_ = PROFILE;
}

ProfileCommand::~ProfileCommand() {
}

SchedulerCommand* ProfileCommand::Clone() {
  return new ProfileCommand();
}

bool ProfileCommand::Parse(const std::string& params) {
  int num = 7;

  char_separator<char> separator(" \n\t\r");
  tokenizer<char_separator<char> > tokens(params, separator);
  tokenizer<char_separator<char> >::iterator iter = tokens.begin();
  for (int i = 0; i < num; i++) {
    if (iter == tokens.end()) {
      std::cout << "ERROR: ProfileCommand has only " << i <<
        " parameters (expected " << num << ")." << std::endl;
      return false;
    }
    iter++;
  }
  if (iter != tokens.end()) {
    std::cout << "ERROR: ProfileCommand has more than "<<
      num << " parameters." << std::endl;
    return false;
  }

  iter = tokens.begin();
  if (!worker_id_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid worker id." << std::endl;
    return false;
  }

  iter++;
  total_virtual_ = boost::lexical_cast<double>(*iter);

  iter++;
  used_virtual_ = boost::lexical_cast<double>(*iter);

  iter++;
  proc_virtual_ = boost::lexical_cast<double>(*iter);

  iter++;
  total_physical_ = boost::lexical_cast<double>(*iter);

  iter++;
  used_physical_ = boost::lexical_cast<double>(*iter);

  iter++;
  proc_physical_ = boost::lexical_cast<double>(*iter);

  return true;
}

std::string ProfileCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += (worker_id_.toString() + " ");
  str += (boost::lexical_cast<std::string>(total_virtual_) + " ");
  str += (boost::lexical_cast<std::string>(used_virtual_) + " ");
  str += (boost::lexical_cast<std::string>(proc_virtual_) + " ");
  str += (boost::lexical_cast<std::string>(total_physical_) + " ");
  str += (boost::lexical_cast<std::string>(used_physical_) + " ");
  str += boost::lexical_cast<std::string>(proc_physical_);
  return str;
}

std::string ProfileCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("worker_id: " + worker_id_.toString() + " ");
  str += ("total_virtual: " + boost::lexical_cast<std::string>(total_virtual_) + " ");
  str += ("used_virtual: " + boost::lexical_cast<std::string>(used_virtual_) + " ");
  str += ("proc_virtual: " + boost::lexical_cast<std::string>(proc_virtual_) + " ");
  str += ("total_physical: " + boost::lexical_cast<std::string>(total_physical_) + " ");
  str += ("used_physical: " + boost::lexical_cast<std::string>(used_physical_) + " ");
  str += ("proc_physical: " + boost::lexical_cast<std::string>(proc_physical_));
  return str;
}

ID<worker_id_t> ProfileCommand::worker_id() {
  return worker_id_;
}

uint64_t ProfileCommand::total_virtual() {
  return total_virtual_;
}

uint64_t ProfileCommand::used_virtual() {
  return used_virtual_;
}

uint64_t ProfileCommand::proc_virtual() {
  return proc_virtual_;
}

uint64_t ProfileCommand::total_physical() {
  return total_physical_;
}

uint64_t ProfileCommand::used_physical() {
  return used_physical_;
}

uint64_t ProfileCommand::proc_physical() {
  return proc_physical_;
}




