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
  * A Nimbus scheduler command.
  *
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#include "lib/scheduler_command.h"

using namespace nimbus; // NOLINT

SchedulerCommand::SchedulerCommand() {
  name_ = "no-op";
}

SchedulerCommand::SchedulerCommand(std::string n,
    const CommandParameterList& p) {
  name_ = n;
  CommandParameterList::const_iterator  iter = p.begin();
  for (; iter != p.end(); ++iter)
    addParameter(iter->second);
}

SchedulerCommand::SchedulerCommand(std::string command) {
  std::vector<std::string> string_params;
  parseCommandFromString(command, name_, string_params);
  std::vector<std::string>::iterator  iter = string_params.begin();
  for (; iter != string_params.end(); ++iter)
    addParameter(CommandParameter(*iter));
}

SchedulerCommand::~SchedulerCommand() {}

void SchedulerCommand::addParameter(CommandParameter cm) {
  parameters_[cm.name()] = cm;
}

std::string SchedulerCommand::toString() {
  std::string rval = name_;
  CommandParameterList::iterator iter = parameters_.begin();
  for (; iter != parameters_.end(); ++iter) {
    rval += " ";
    rval += (iter->second).toString();
  }
  return rval;
}

std::string SchedulerCommand::name() {
  return name_;
}

CommandParameterList* SchedulerCommand::parameters() {
  return &parameters_;
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






IDSet::IDSet() {}

IDSet::IDSet(std::string s) {
  parseIDSetFromString(s, identifiers_);
}

IDSet::~IDSet() {}

std::string IDSet::toString() {
  bool empty = true;
  std::string rval = "{";
  IDSetIter iter =  identifiers_.begin();
  for (; iter !=  identifiers_.end(); ++iter) {
    empty = false;
    std::ostringstream ss;
    ss << *iter;
    rval += ss.str();
    rval += ",";
  }
  if (empty)
    rval += "}";
  else
    rval[rval.length() - 1] = '}';
  return rval;
}

void IDSet::insert(int n) {
  identifiers_.insert(n);
}

void IDSet::clear() {
  identifiers_.clear();
}

int IDSet::size() {
  return identifiers_.size();
}


IDSet::IDSetIter IDSet::begin() {
  return identifiers_.begin();
}

IDSet::IDSetIter IDSet::end() {
  return identifiers_.end();
}


