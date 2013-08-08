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
#include "lib/parser.h"

SchedulerCommand::SchedulerCommand() {
  name = "no-op";
}

SchedulerCommand::SchedulerCommand(std::string n,
    const CommandParameterList& p) {
  name = n;
  CommandParameterList::const_iterator  iter = p.begin();
  for (; iter != p.end(); ++iter)
    addParameter(iter->second);
}

SchedulerCommand::SchedulerCommand(std::string command) {
  std::vector<std::string> string_params;
  parseCommandFromString(command, name, string_params);
  std::vector<std::string>::iterator  iter = string_params.begin();
  for (; iter != string_params.end(); ++iter)
    addParameter(CommandParameter(*iter));
}

SchedulerCommand::~SchedulerCommand() {}

void SchedulerCommand::addParameter(CommandParameter cm) {
  parameters[cm.getTag()] = cm;
}

std::string SchedulerCommand::toString() {
  std::string rval = name;
  CommandParameterList::iterator iter = parameters.begin();
  for (; iter != parameters.end(); ++iter) {
    rval += " ";
    rval += (iter->second).toString();
  }
  return rval;
}

std::string SchedulerCommand::getName() {
  return name;
}

CommandParameterList SchedulerCommand::getParameters() {
  return parameters;
}





CommandParameter::CommandParameter() {
  tag = "empty-field";
}

CommandParameter::CommandParameter(std::string t, std::string a,
    const IDSet& s) {
  tag = t;
  arg = a;
  set = s;
}

CommandParameter::CommandParameter(std::string parameter) {
  std::string string_set;
  parseParameterFromString(parameter, tag, arg, string_set);
  if (isSet(string_set))
    set = IDSet(string_set);
}

CommandParameter::~CommandParameter() {}

std::string CommandParameter::toString() {
  std::string rval = tag;
  rval += ":";
  if (arg == "")
    rval += set.toString();
  else
    rval += arg;

  return rval;
}

std::string CommandParameter::getTag() {
  return tag;
}

std::string CommandParameter::getArg() {
  return arg;
}

IDSet CommandParameter::getIDSet() {
  return set;
}






IDSet::IDSet() {}

IDSet::IDSet(std::string s) {
  ParseIDSetFromString(s, set);
}

IDSet::~IDSet() {}

std::string IDSet::toString() {
  bool empty = true;
  std::string rval = "{";
  std::set<int>::iterator iter = set.begin();
  for (; iter != set.end(); ++iter) {
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
  set.insert(n);
}

void IDSet::clear() {
  set.clear();
}

int IDSet::size() {
  return set.size();
}


