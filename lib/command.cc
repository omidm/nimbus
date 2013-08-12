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
  * Object representation of a command. Used by workers to send commands to
  * server and server to send commands down to workers. There is a base class
  * Command to define a command in general. Children of this superclass
  * implement each of the commands used in Nimbus. 
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  */

#include "lib/command.h"

using namespace nimbus; // NOLINT

Command::Command() {
  name_ = "no-op";
}

Command::Command(std::string n,
    const CommandParameterList& p) {
  name_ = n;
  CommandParameterList::const_iterator  iter = p.begin();
  for (; iter != p.end(); ++iter)
    addParameter(iter->second);
}

Command::Command(std::string command) {
  std::vector<std::string> string_params;
  parseCommandFromString(command, name_, string_params);
  std::vector<std::string>::iterator  iter = string_params.begin();
  for (; iter != string_params.end(); ++iter)
    addParameter(CommandParameter(*iter));
}

Command::~Command() {}

void Command::addParameter(CommandParameter cm) {
  parameters_[cm.name()] = cm;
}

std::string Command::toString() {
  std::string rval = name_;
  CommandParameterList::iterator iter = parameters_.begin();
  for (; iter != parameters_.end(); ++iter) {
    rval += " ";
    rval += (iter->second).toString();
  }
  return rval;
}

std::string Command::name() {
  return name_;
}

CommandParameterList* Command::parameters() {
  return &parameters_;
}

worker_id_t Command::worker_id() {
  return worker_id_;
}

void Command::set_worker_id(worker_id_t id) {
  worker_id_ = id;
}

CommandParameter::CommandParameter() {
  name_ = "empty-field";
}

CommandParameter::CommandParameter(std::string n,
    std::string v, const IDSet& s) {
  name_ = n;
  value_ = v;
  identifier_set_ = s;
}

CommandParameter::CommandParameter(std::string parameter) {
  std::string string_set;
  parseParameterFromString(parameter, name_, value_, string_set);
  if (isSet(string_set))
    identifier_set_ = IDSet(string_set);
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

IDSet* CommandParameter::identifier_set() {
  return &identifier_set_;
}






