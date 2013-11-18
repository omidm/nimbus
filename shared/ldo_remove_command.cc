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

/***********************************************************************
 * AUTHOR: Philip Levis <pal> *   FILE: .//ldo_remove_command.cc
 *   DATE: Fri Nov 15 12:25:27 2013
 *  DESCR:
 ***********************************************************************/
#include "shared/ldo_remove_command.h"
#include "shared/escaper.h"

namespace nimbus {

/**
 * \fn nimbus::LdoRemoveCommand::LdoRemoveCommand()
 * \brief Brief description.
 * \return
*/
  nimbus::LdoRemoveCommand::LdoRemoveCommand() {
  name_ = LDO_REMOVE_NAME;
  type_ = LDO_REMOVE;
  region_ = NULL;
  object_ = NULL;
}


/**
 * \fn nimbus::LdoRemoveCommand::LdoRemoveCommand()
 * \brief Brief description.
 * \return
*/
nimbus::LdoRemoveCommand::LdoRemoveCommand(LogicalDataObject* obj) {
  name_ = LDO_REMOVE_NAME;
  type_ = LDO_REMOVE;
  region_ = new GeometricRegion(*obj->region());
  object_ = new LogicalDataObject(obj->id(), obj->variable(), region_);
}

nimbus::LdoRemoveCommand::~LdoRemoveCommand() {
  delete object_;
}

/**
 * \fn SchedulerCommand * nimbus::LdoRemoveCommand::Clone()
 * \brief Brief description.
 * \return
*/
SchedulerCommand * nimbus::LdoRemoveCommand::Clone() {
  return new LdoRemoveCommand();
}


/**
 * \fn bool nimbus::LdoRemoveCommand::Parse(const std::string &param_segment)
 * \brief Brief description.
 * \param param_segment
 * \return
*/
bool nimbus::LdoRemoveCommand::Parse(const std::string &param_segment) {
  std::string strCopy = param_segment;
  object_ = new LogicalDataObject();
  UnescapeString(&strCopy);
  object_->Parse(strCopy);
  return true;
}


/**
 * \fn std::string nimbus::LdoRemoveCommand::toString()
 * \brief Brief description.
 * \return
*/
std::string nimbus::LdoRemoveCommand::toString() {
  std::string str;
  std::string payload;
  str += (name_ + " ");
  object_->SerializeToString(&payload);
  str += payload;
  EscapeString(&str);
  return str;
}


/**
 * \fn std::string nimbus::LdoRemoveCommand::toStringWTags()
 * \brief Brief description.
 * \return
*/
std::string nimbus::LdoRemoveCommand::toStringWTags() {
  std::string str;
  std::string payload;
  str += (name_ + " ");
  object_->SerializeToString(&payload);
  str += ("object:" + payload);
  EscapeString(&str);
  return str;
}


/**
 * \fn LogicalDataObject * nimbus::LdoRemoveCommand::object()
 * \brief Brief description.
 * \return
*/
LogicalDataObject * nimbus::LdoRemoveCommand::object() {
  return object_;
}


}  // namespace nimbus
