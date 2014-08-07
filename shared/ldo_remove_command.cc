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

using namespace nimbus; // NOLINT

/**
 * \fn LdoRemoveCommand::LdoRemoveCommand()
 * \brief Brief description.
 * \return
*/
LdoRemoveCommand::LdoRemoveCommand() {
  name_ = LDO_REMOVE_NAME;
  type_ = LDO_REMOVE;
  region_ = NULL;
  object_ = NULL;
}


/**
 * \fn LdoRemoveCommand::LdoRemoveCommand()
 * \brief Brief description.
 * \return
*/
LdoRemoveCommand::LdoRemoveCommand(LogicalDataObject* obj) {
  name_ = LDO_REMOVE_NAME;
  type_ = LDO_REMOVE;
  region_ = new GeometricRegion(*obj->region());
  object_ = new LogicalDataObject(obj->id(), obj->variable(), region_);
}

LdoRemoveCommand::~LdoRemoveCommand() {
  delete object_;
}

/**
 * \fn SchedulerCommand * LdoRemoveCommand::Clone()
 * \brief Brief description.
 * \return
*/
SchedulerCommand * LdoRemoveCommand::Clone() {
  return new LdoRemoveCommand();
}


/**
 * \fn bool LdoRemoveCommand::Parse(const std::string &param_segment)
 * \brief Brief description.
 * \param param_segment
 * \return
*/
bool LdoRemoveCommand::Parse(const std::string &data) {
  LdoRemovePBuf buf;
  bool result = buf.ParseFromString(data);

  if (!result) {
    dbg(DBG_ERROR, "ERROR: Failed to parse LdoRemoveCommand from string.\n");
    return false;
  } else {
    ReadFromProtobuf(buf);
    return true;
  }
}


bool LdoRemoveCommand::Parse(const SchedulerPBuf& buf) {
  if (!buf.has_ldo_remove()) {
    dbg(DBG_ERROR, "ERROR: Failed to parse LdoRemoveCommand from SchedulerPBuf.\n");
    return false;
  } else {
    return ReadFromProtobuf(buf.ldo_remove());
  }
}

  /**
 * \fn std::string LdoRemoveCommand::toString()
 * \brief Brief description.
 * \return
*/
std::string LdoRemoveCommand::toString() {
  std::string result;

  // First we construct a general scheduler buffer, then
  // add the spawn compute field to it, then serialize.
  SchedulerPBuf buf;
  buf.set_type(SchedulerPBuf_Type_LDO_REMOVE);
  LdoRemovePBuf* cmd = buf.mutable_ldo_remove();
  WriteToProtobuf(cmd);

  buf.SerializeToString(&result);

  return result;
}


/**
 * \fn std::string LdoRemoveCommand::toStringWTags()
 * \brief Brief description.
 * \return
*/
std::string LdoRemoveCommand::toStringWTags() {
  std::ostringstream sstream;
  sstream << (name_ + ",");
  sstream << "logical_id:" << object()->id() << ",";
  sstream << "variable:" << object()->variable() << ",";
  sstream << object()->region()->toString();
  return sstream.str();
}


/**
 * \fn LogicalDataObject * LdoRemoveCommand::object()
 * \brief Brief description.
 * \return
*/
LogicalDataObject * LdoRemoveCommand::object() {
  return object_;
}

bool LdoRemoveCommand::ReadFromProtobuf(const LdoRemovePBuf& buf) {
  if (object_ != NULL) {
    delete object_;
    object_ = NULL;
    region_ = NULL;
  }
  region_ = new GeometricRegion();
  region_->FillInValues(&buf.ldo().region());
  object_ = new LogicalDataObject(buf.ldo().data_id(),
                                  buf.ldo().variable(),
                                  region_);
  return true;
}

bool LdoRemoveCommand::WriteToProtobuf(LdoRemovePBuf* buf) {
  if (object_ == NULL) {
    return false;
  }
  object_->FillInMessage(buf->mutable_ldo());
  return true;
}
