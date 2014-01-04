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
  * Terminate command used to signal the scheduler that the application is
  * complete and there  would be no more spawned jobs.
  * Also used by the scheduler to terminate the workers.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "shared/terminate_command.h"

using namespace nimbus;  // NOLINT
using boost::tokenizer;
using boost::char_separator;

TerminateCommand::TerminateCommand() {
  name_ = TERMINATE_NAME;
  type_ = TERMINATE;
}

TerminateCommand::TerminateCommand(const ID<exit_status_t>& exit_status)
: exit_status_(exit_status) {
  name_ = TERMINATE_NAME;
  type_ = TERMINATE;
}

TerminateCommand::~TerminateCommand() {
}

SchedulerCommand* TerminateCommand::Clone() {
  return new TerminateCommand();
}

bool TerminateCommand::Parse(const std::string& params) {
  int num = 1;

  char_separator<char> separator(" \n\t\r");
  tokenizer<char_separator<char> > tokens(params, separator);
  tokenizer<char_separator<char> >::iterator iter = tokens.begin();
  for (int i = 0; i < num; i++) {
    if (iter == tokens.end()) {
      std::cout << "ERROR: TerminateCommand has only " << i <<
        " parameters (expected " << num << ")." << std::endl;
      return false;
    }
    iter++;
  }
  if (iter != tokens.end()) {
    std::cout << "ERROR: TerminateCommand has more than "<<
      num << " parameters." << std::endl;
    return false;
  }

  iter = tokens.begin();
  if (!exit_status_.Parse(*iter)) {
    std::cout << "ERROR: Could not detect valid exit status." << std::endl;
    return false;
  }

  return true;
}

std::string TerminateCommand::toString() {
  std::string str;
  str += (name_ + " ");
  str += exit_status_.toString();
  return str;
}

std::string TerminateCommand::toStringWTags() {
  std::string str;
  str += (name_ + " ");
  str += ("exit-status:" + exit_status_.toString());
  return str;
}


ID<exit_status_t> TerminateCommand::exit_status() {
  return exit_status_;
}

